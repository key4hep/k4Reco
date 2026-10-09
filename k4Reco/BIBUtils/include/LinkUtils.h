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
#ifndef K4RECO_BIBUTILS_LINKUTILS_H
#define K4RECO_BIBUTILS_LINKUTILS_H 1

#include <podio/ObjectID.h>

#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace k4reco::bibutils {

/// Groups the links of a reco -> sim link collection by their reco ("from") object,
/// keeping all of them, since a reco hit can be linked to several sim hits.
template <typename LinkCollT>
std::unordered_map<podio::ObjectID, std::vector<typename LinkCollT::value_type>> linksByFrom(const LinkCollT& links) {
  std::unordered_map<podio::ObjectID, std::vector<typename LinkCollT::value_type>> byFrom;
  byFrom.reserve(links.size());
  for (const auto& link : links) {
    byFrom[link.getFrom().getObjectID()].push_back(link);
  }
  return byFrom;
}

/// Copies the links of an accepted object into outLinks, keeping their targets and weights.
template <typename LinkCollT>
void copyLinks(const std::vector<typename LinkCollT::value_type>& links, LinkCollT& outLinks) {
  for (const auto& link : links) {
    outLinks.push_back(link.clone());
  }
}

/// As above, and also appends the target of each link to the subset collection outTo,
/// once per target: seenTo records the targets already written.
template <typename LinkCollT, typename ToCollT>
void copyLinks(const std::vector<typename LinkCollT::value_type>& links, LinkCollT& outLinks, ToCollT& outTo,
               std::unordered_set<podio::ObjectID>& seenTo) {
  for (const auto& link : links) {
    outLinks.push_back(link.clone());
    const auto to = link.getTo();
    if (seenTo.insert(to.getObjectID()).second) {
      outTo.push_back(to);
    }
  }
}

} // namespace k4reco::bibutils

#endif
