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
#include "SplitCollectionByPolarAngle.h"

#include "BIBUtilsHelpers.h"

#include <edm4hep/SimTrackerHit.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <unordered_map>

using k4reco::bibutils::objectKey;

SplitCollectionByPolarAngle::SplitCollectionByPolarAngle(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(name, svcLoc,
                       {KeyValues("TrackerHitInputCollections", {"VBTrackerHits"}),
                        KeyValues("TrackerHitInputRelations", {"VBTrackerHitsRelations"})},
                       {KeyValues("TrackerHitOutputCollections", {"VBTrackerHitsSplit"}),
                        KeyValues("TrackerSimHitOutputCollections", {"VertexBarrelCollectionSplit"}),
                        KeyValues("TrackerHitOutputRelations", {"VBTrackerHitsRelationsSplit"})}) {}

StatusCode SplitCollectionByPolarAngle::initialize() {
  m_histograms[hTheta].reset(new Gaudi::Accumulators::StaticRootHistogram<1>{
      this, "theta", "polar angle of the hit;#theta [rad]", {1000, 0., M_PI}});

  return StatusCode::SUCCESS;
}

std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::SimTrackerHitCollection,
           edm4hep::TrackerHitSimTrackerHitLinkCollection>
SplitCollectionByPolarAngle::operator()(const edm4hep::TrackerHitPlaneCollection& trackerHits,
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

  std::size_t nKept = 0;

  for (const auto& hit : trackerHits) {
    const auto& pos = hit.getPosition();

    // Polar angle theta = acos(z / r) of the hit, in radians.
    const double posMod = std::sqrt(pos.x * pos.x + pos.y * pos.y + pos.z * pos.z);
    if (posMod == 0.) {
      continue;
    }
    const double hitTheta = std::acos(pos.z / posMod);
    const double hitThetaDeg = hitTheta * 180. / M_PI;

    if (hitThetaDeg < m_thetaMin || hitThetaDeg > m_thetaMax) {
      continue;
    }

    if (m_fillHistos) {
      ++(*m_histograms[hTheta])[hitTheta];
    }

    const auto simIt = hitToSim.find(objectKey(hit.getObjectID()));
    if (simIt == hitToSim.end()) {
      continue;
    }

    outHits.push_back(hit);
    outSimHits.push_back(simIt->second);
    auto link = outLinks.create();
    link.setFrom(hit);
    link.setTo(simIt->second);
    link.setWeight(1.0);
    ++nKept;
  }

  debug() << nKept << " / " << trackerHits.size() << " tracker hits kept within the polar-angle window" << endmsg;

  return std::make_tuple(std::move(outHits), std::move(outSimHits), std::move(outLinks));
}

DECLARE_COMPONENT(SplitCollectionByPolarAngle)
