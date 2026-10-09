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

#include "LinkUtils.h"

#include <edm4hep/utils/vector_utils.h>

#include <podio/ObjectID.h>

#include <cmath>
#include <cstddef>
#include <unordered_set>

SplitCollectionByPolarAngle::SplitCollectionByPolarAngle(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(name, svcLoc,
                       {KeyValue("TrackerHitInputCollections", "VBTrackerHits"),
                        KeyValue("TrackerHitInputRelations", "VBTrackerHitsRelations")},
                       {KeyValue("TrackerHitOutputCollections", "VBTrackerHitsSplit"),
                        KeyValue("TrackerSimHitOutputCollections", "VertexBarrelCollectionSplit"),
                        KeyValue("TrackerHitOutputRelations", "VBTrackerHitsRelationsSplit")}) {}

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

  const auto linksByHit = k4reco::bibutils::linksByFrom(trackerHitLinks);
  // A sim hit can be linked to more than one accepted hit, but is written out once.
  std::unordered_set<podio::ObjectID> writtenSimHits;

  std::size_t nKept = 0;

  for (const auto& hit : trackerHits) {
    // Polar angle of the hit, in radians.
    const double hitTheta = edm4hep::utils::anglePolar(hit.getPosition());
    const double hitThetaDeg = hitTheta * 180. / M_PI;

    if (hitThetaDeg < m_thetaMin || hitThetaDeg > m_thetaMax) {
      continue;
    }

    if (m_fillHistos) {
      ++(*m_histograms[hTheta])[hitTheta];
    }

    const auto linksIt = linksByHit.find(hit.getObjectID());
    if (linksIt == linksByHit.end()) {
      continue;
    }

    outHits.push_back(hit);
    k4reco::bibutils::copyLinks(linksIt->second, outLinks, outSimHits, writtenSimHits);
    ++nKept;
  }

  debug() << nKept << " / " << trackerHits.size() << " tracker hits kept within the polar-angle window" << endmsg;

  return std::make_tuple(std::move(outHits), std::move(outSimHits), std::move(outLinks));
}

DECLARE_COMPONENT(SplitCollectionByPolarAngle)
