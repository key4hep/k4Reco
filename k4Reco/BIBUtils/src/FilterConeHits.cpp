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
#include "FilterConeHits.h"

#include "BIBUtilsHelpers.h"
#include "TrackHelix.h"

#include <edm4hep/SimTrackerHit.h>

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"

#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <unordered_set>

using k4reco::bibutils::objectKey;
using k4reco::bibutils::TrackHelix;

FilterConeHits::FilterConeHits(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(name, svcLoc,
                       {KeyValues("MCParticleCollection", {"MCParticle"}),
                        KeyValues("TrackerHitInputCollections", {"VBTrackerHits"}),
                        KeyValues("TrackerHitInputRelations", {"VBTrackerHitsRelations"})},
                       {KeyValues("TrackerHitOutputCollections", {"VBTrackerHitsConed"}),
                        KeyValues("TrackerSimHitOutputCollections", {"VertexBarrelCollectionConed"}),
                        KeyValues("TrackerHitOutputRelations", {"VBTrackerHitsRelationsConed"})}) {}

StatusCode FilterConeHits::initialize() {
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

  m_histograms[hDistXY].reset(new Gaudi::Accumulators::StaticRootHistogram<1>{
      this, "distXY", "hit-to-helix XY distance;d_{XY} [mm]", {1000, 0., 1000.}});
  m_histograms[hDistZ].reset(new Gaudi::Accumulators::StaticRootHistogram<1>{
      this, "distZ", "hit-to-helix Z distance;d_{Z} [mm]", {1000, 0., 1000.}});
  m_histograms[hDist3D].reset(new Gaudi::Accumulators::StaticRootHistogram<1>{
      this, "dist3D", "hit-to-helix 3D distance;d_{3D} [mm]", {1000, 0., 1000.}});
  m_histograms[hAngle].reset(new Gaudi::Accumulators::StaticRootHistogram<1>{
      this, "angle", "angle between hit and particle;angle [rad]", {1000, 0., 1.}});
  m_histograms[hTime].reset(new Gaudi::Accumulators::StaticRootHistogram<1>{
      this, "time", "time at the point of closest approach;T [mm/GeV]", {1000, 0., 2000.}});
  m_histograms[hPathLength].reset(new Gaudi::Accumulators::StaticRootHistogram<1>{
      this, "pathLength", "pathlength at the point of closest approach;L [mm]", {1000, 0., 12000.}});

  return StatusCode::SUCCESS;
}

std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::SimTrackerHitCollection,
           edm4hep::TrackerHitSimTrackerHitLinkCollection>
FilterConeHits::operator()(const edm4hep::MCParticleCollection& mcParticles,
                           const edm4hep::TrackerHitPlaneCollection& trackerHits,
                           const edm4hep::TrackerHitSimTrackerHitLinkCollection& trackerHitLinks) const {
  edm4hep::TrackerHitPlaneCollection outHits;
  outHits.setSubsetCollection();
  edm4hep::SimTrackerHitCollection outSimHits;
  outSimHits.setSubsetCollection();
  edm4hep::TrackerHitSimTrackerHitLinkCollection outLinks;

  // Map each reconstructed hit to its simulated hit through the input links.
  std::unordered_map<std::uint64_t, edm4hep::SimTrackerHit> hitToSim;
  hitToSim.reserve(trackerHitLinks.size());
  for (const auto& link : trackerHitLinks) {
    hitToSim.emplace(objectKey(link.getFrom().getObjectID()), link.getTo());
  }

  // Track which hits have already been accepted so a hit close to two particles
  // is not written out twice.
  std::unordered_set<std::uint64_t> acceptedHits;

  for (const auto& part : mcParticles) {
    const int genStatus = part.getGeneratorStatus();
    if (std::find(m_coneAroundStatus.begin(), m_coneAroundStatus.end(), genStatus) == m_coneAroundStatus.end()) {
      continue;
    }

    const auto& mom = part.getMomentum();
    const auto& vtx = part.getVertex();
    double momentum[3] = {mom.x, mom.y, mom.z};
    double vertex[3] = {vtx.x, vtx.y, vtx.z};
    const double partP = std::sqrt(mom.x * mom.x + mom.y * mom.y + mom.z * mom.z);

    TrackHelix helix;
    helix.initialize(vertex, momentum, static_cast<double>(part.getCharge()), m_bField);

    // Intersection of the helix with the tracker outer cylinder. If the particle
    // spirals and never reaches it, getPointOnCircle returns -1e20.
    const double intersectionTime = helix.getPointOnCircle(m_trackerOuterRadius, vertex);

    for (const auto& hit : trackerHits) {
      const std::uint64_t key = objectKey(hit.getObjectID());
      if (acceptedHits.count(key)) {
        continue;
      }

      const auto& pos = hit.getPosition();

      // Skip hits in the opposite hemisphere w.r.t. the particle direction.
      if ((pos.x * mom.x + pos.y * mom.y + pos.z * mom.z) < 0.) {
        continue;
      }

      double hitPosition[3] = {pos.x, pos.y, pos.z};
      double hitDistance[3] = {0., 0., 0.};
      const double timeAtPCA = helix.getDistanceToPoint(hitPosition, hitDistance);

      // Exclude the opposite side of the helix w.r.t. the production vertex.
      if (timeAtPCA < 0.) {
        continue;
      }
      // Avoid the helix of central trajectories re-entering the tracker.
      if (timeAtPCA > intersectionTime && intersectionTime != -1.e20) {
        continue;
      }

      const double pathLength = partP * timeAtPCA;
      const double hitAngle = std::atan2(hitDistance[2], pathLength);

      if (m_fillHistos) {
        ++(*m_histograms[hDistXY])[hitDistance[0]];
        ++(*m_histograms[hDistZ])[hitDistance[1]];
        ++(*m_histograms[hDist3D])[hitDistance[2]];
        ++(*m_histograms[hAngle])[hitAngle];
        ++(*m_histograms[hTime])[timeAtPCA];
        ++(*m_histograms[hPathLength])[pathLength];
      }

      bool save = false;
      if (m_deltaRCut > 0. && hitAngle < m_deltaRCut) {
        save = true;
      }
      if (m_dist3DCut > 0. && hitDistance[2] < m_dist3DCut) {
        save = true;
      }
      if (!save) {
        continue;
      }

      const auto simIt = hitToSim.find(key);
      if (simIt == hitToSim.end()) {
        continue;
      }

      outHits.push_back(hit);
      outSimHits.push_back(simIt->second);
      auto link = outLinks.create();
      link.setFrom(hit);
      link.setTo(simIt->second);
      link.setWeight(1.0);
      acceptedHits.insert(key);
    }
  }

  debug() << acceptedHits.size() << " / " << trackerHits.size() << " tracker hits kept inside the cone" << endmsg;

  return std::make_tuple(std::move(outHits), std::move(outSimHits), std::move(outLinks));
}

DECLARE_COMPONENT(FilterConeHits)
