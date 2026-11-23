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

#include <catch2/catch_test_macros.hpp>
#include <edm4hep/SimTrackerHitCollection.h>
#include <edm4hep/TrackerHitPlaneCollection.h>
#include <edm4hep/TrackerHitSimTrackerHitLinkCollection.h>
#include <podio/LinkNavigator.h>

#include <algorithm>
#include <cstddef>
#include <unordered_set>

TEST_CASE("copyLinks keeps distinct targets with colliding ObjectIDs", "[BIBUtils][copyLinks]") {
  edm4hep::TrackerHitPlaneCollection recoHits;
  const edm4hep::TrackerHit from = recoHits.create();

  // The first hits in separate collections have the same ID until the collections are put in a frame.
  edm4hep::SimTrackerHitCollection barrelHits, endcapHits;
  const edm4hep::SimTrackerHit barrelHit = barrelHits.create();
  const edm4hep::SimTrackerHit endcapHit = endcapHits.create();
  REQUIRE(barrelHit.getObjectID() == endcapHit.getObjectID());

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

  CHECK(outHits.size() == 2);
  CHECK(std::count(outHits.begin(), outHits.end(), barrelHit) == 1);
  CHECK(std::count(outHits.begin(), outHits.end(), endcapHit) == 1);

  REQUIRE(outLinks.size() == inputLinks.size());
  for (std::size_t i = 0; i < linkedTo.size(); ++i) {
    CAPTURE(i);
    const auto& [to, weight] = linkedTo[i];
    const auto link = outLinks[i];
    CHECK(link.getFrom() == from);
    CHECK(link.getTo() == to);
    CHECK(link.getWeight() == weight);
  }
}
