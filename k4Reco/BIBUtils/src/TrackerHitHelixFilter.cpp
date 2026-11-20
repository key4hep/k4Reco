/*
 * Copyright (c) 2020-2024 Key4hep-Project.
 *
 * This file is part of Key4hep.
 * See https://key4hep.github.io/key4hep-doc/ for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#include "TrackerHitHelixFilter.h"

#include "LinkUtils.h"

#include <edm4hep/utils/Helix.h>
#include <edm4hep/utils/vector_utils.h>

#include <podio/LinkNavigator.h>
#include <podio/ObjectID.h>

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <unordered_set>

TrackerHitHelixFilter::TrackerHitHelixFilter(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(name, svcLoc,
                       {KeyValue("MCParticleCollection", "MCParticle"),
                        KeyValue("TrackerHitInputCollections", "VBTrackerHits"),
                        KeyValue("TrackerHitInputRelations", "VBTrackerHitsRelations")},
                       {KeyValue("TrackerHitOutputCollections", "VBTrackerHitsConed"),
                        KeyValue("TrackerSimHitOutputCollections", "VertexBarrelCollectionConed"),
                        KeyValue("TrackerHitOutputRelations", "VBTrackerHitsRelationsConed")}) {}

StatusCode TrackerHitHelixFilter::initialize() {
  m_geoSvc = serviceLocator()->service(m_geoSvcName);
  if (!m_geoSvc) {
    error() << "Unable to retrieve the GeoSvc" << endmsg;
    return StatusCode::FAILURE;
  }

  // Get the value of the magnetic field along z at the origin, in Tesla.
  double bfield[3] = {0., 0., 0.};
  m_geoSvc->getDetector()->field().magneticField({0., 0., 0.}, bfield);
  m_bField = bfield[2] / dd4hep::tesla;
  info() << "Magnetic field at the origin: Bz = " << m_bField << " T" << endmsg;

  return StatusCode::SUCCESS;
}

std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::SimTrackerHitCollection,
           edm4hep::TrackerHitSimTrackerHitLinkCollection>
TrackerHitHelixFilter::operator()(const edm4hep::MCParticleCollection& mcParticles,
                                  const edm4hep::TrackerHitPlaneCollection& trackerHits,
                                  const edm4hep::TrackerHitSimTrackerHitLinkCollection& trackerHitLinks) const {
  edm4hep::TrackerHitPlaneCollection outHits;
  outHits.setSubsetCollection();
  edm4hep::SimTrackerHitCollection outSimHits;
  outSimHits.setSubsetCollection();
  edm4hep::TrackerHitSimTrackerHitLinkCollection outLinks;

  const podio::LinkNavigator linkNavigator(trackerHitLinks);
  // A sim hit can be linked to more than one accepted hit, but is written out once.
  std::unordered_set<podio::ObjectID> writtenSimHits;

  // Track which hits have already been accepted so a hit close to two particles
  // is not written out twice.
  std::unordered_set<podio::ObjectID> acceptedHits;

  for (const auto& part : mcParticles) {
    const int genStatus = part.getGeneratorStatus();
    if (std::find(m_coneAroundStatus.begin(), m_coneAroundStatus.end(), genStatus) == m_coneAroundStatus.end()) {
      continue;
    }
    // Neutral particles are not relevant for tracking.
    if (part.getCharge() == 0.f) {
      continue;
    }

    const auto& mom = part.getMomentum();
    const auto& vtx = part.getVertex();
    const double pt = std::hypot(mom.x, mom.y);
    if (pt == 0. || !edm4hep::utils::finite(mom) || !edm4hep::utils::finite(vtx) || !std::isfinite(part.getCharge())) {
      continue;
    }
    const double partP = std::hypot(pt, mom.z);
    const edm4hep::utils::Helix helix(vtx, mom, part.getCharge(), m_bField, vtx);

    // A curling track may never reach the tracker outer cylinder.
    const auto intersectionLength = helix.lengthToRadius(m_trackerOuterRadius);

    for (const auto& hit : trackerHits) {
      const auto key = hit.getObjectID();
      if (acceptedHits.count(key)) {
        continue;
      }

      const auto& pos = hit.getPosition();

      // Skip hits in the opposite hemisphere w.r.t. the particle direction.
      if (!edm4hep::utils::finite(pos) || (pos.x * mom.x + pos.y * mom.y + pos.z * mom.z) < 0.) {
        continue;
      }

      const auto hitDistance = helix.getDistanceToPoint({pos.x, pos.y, pos.z});
      double transverseLength = hitDistance.path;
      // For a planar circle, z cannot distinguish the turns. Keep the first forward projection.
      if (transverseLength < 0. && helix.getTanLambda() == 0. && helix.getOmega() != 0.) {
        transverseLength += 2. * std::numbers::pi / std::abs(helix.getOmega());
      }
      // Reject hits before production or after the first exit from the tracker.
      if (transverseLength < 0. || (intersectionLength && transverseLength > *intersectionLength)) {
        continue;
      }

      const double timeAtPCA = transverseLength / pt;
      const double pathLength = partP * timeAtPCA;
      const double hitAngle = std::atan2(hitDistance.threeDimensional, pathLength);

      if (m_fillHistos) {
        ++m_histograms[hDistXY][hitDistance.transverse];
        ++m_histograms[hDistZ][hitDistance.longitudinal];
        ++m_histograms[hDist3D][hitDistance.threeDimensional];
        ++m_histograms[hAngle][hitAngle];
        ++m_histograms[hTime][timeAtPCA];
        ++m_histograms[hPathLength][pathLength];
      }

      bool save = false;
      if (m_deltaRCut > 0. && hitAngle < m_deltaRCut) {
        save = true;
      }
      if (m_dist3DCut > 0. && hitDistance.threeDimensional < m_dist3DCut) {
        save = true;
      }
      if (!save) {
        continue;
      }

      const auto linkedSimHits = linkNavigator.getLinkedTo(hit);
      if (linkedSimHits.empty()) {
        continue;
      }

      outHits.push_back(hit);
      k4reco::bibutils::copyLinks(hit, linkedSimHits, outLinks, outSimHits, writtenSimHits);
      acceptedHits.insert(key);
    }
  }

  debug() << acceptedHits.size() << " / " << trackerHits.size() << " tracker hits kept inside the cone" << endmsg;

  return std::make_tuple(std::move(outHits), std::move(outSimHits), std::move(outLinks));
}

DECLARE_COMPONENT(TrackerHitHelixFilter)
