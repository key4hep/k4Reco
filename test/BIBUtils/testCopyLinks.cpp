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

#include "LinkUtils.h"

#include <edm4hep/SimTrackerHitCollection.h>
#include <edm4hep/TrackerHitPlaneCollection.h>
#include <edm4hep/TrackerHitSimTrackerHitLinkCollection.h>
#include <podio/LinkNavigator.h>

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <unordered_set>

int main() {
  edm4hep::TrackerHitPlaneCollection recoHits;
  const edm4hep::TrackerHit from = recoHits.create();

  // The first hits in separate collections have the same ID until the collections are put in a frame.
  edm4hep::SimTrackerHitCollection barrelHits, endcapHits;
  const edm4hep::SimTrackerHit barrelHit = barrelHits.create();
  const edm4hep::SimTrackerHit endcapHit = endcapHits.create();
  if (barrelHit.getObjectID() != endcapHit.getObjectID()) {
    std::cerr << "Test input does not contain an ObjectID collision\n";
    return 1;
  }

  edm4hep::TrackerHitSimTrackerHitLinkCollection inputLinks;
  const auto addLink = [&](const edm4hep::SimTrackerHit& to, float weight) {
    auto link = inputLinks.create();
    link.setFrom(from);
    link.setTo(to);
    link.setWeight(weight);
  };
  addLink(barrelHit, 0.25f);
  addLink(endcapHit, 0.5f);
  addLink(barrelHit, 0.75f);

  const podio::LinkNavigator navigator(inputLinks);
  const auto linkedTo = navigator.getLinkedTo(from);
  edm4hep::TrackerHitSimTrackerHitLinkCollection outLinks;
  edm4hep::SimTrackerHitCollection outHits;
  outHits.setSubsetCollection();
  std::unordered_set<edm4hep::SimTrackerHit> seen;
  k4reco::bibutils::copyLinks(from, linkedTo, outLinks, outHits, seen);

  if (outHits.size() != 2 || std::count(outHits.begin(), outHits.end(), barrelHit) != 1 ||
      std::count(outHits.begin(), outHits.end(), endcapHit) != 1) {
    std::cerr << "Expected each distinct target in the subset exactly once\n";
    return 1;
  }
  if (outLinks.size() != inputLinks.size()) {
    std::cerr << "Expected all links to be copied, including repeated targets\n";
    return 1;
  }
  for (std::size_t i = 0; i < linkedTo.size(); ++i) {
    const auto& [to, weight] = linkedTo[i];
    const auto link = outLinks[i];
    if (link.getFrom() != from || link.getTo() != to || link.getWeight() != weight) {
      std::cerr << "Copied link endpoints or weight differ from the input\n";
      return 1;
    }
  }
}
