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

#include <unordered_set>

namespace k4reco::bibutils {

/// Links an accepted object to each of its targets, as returned by podio::LinkNavigator::getLinkedTo,
/// keeping their weights.
template <typename LinkCollT>
void copyLinks(const typename LinkCollT::from_type& from, const auto& linkedTo, LinkCollT& outLinks) {
  for (const auto& [to, weight] : linkedTo) {
    auto link = outLinks.create();
    link.setFrom(from);
    link.setTo(to);
    link.setWeight(weight);
  }
}

/// As above, and also appends each target to the subset collection outTo,
/// once per target: seenTo records the targets already written.
template <typename LinkCollT, typename ToCollT>
void copyLinks(const typename LinkCollT::from_type& from, const auto& linkedTo, LinkCollT& outLinks, ToCollT& outTo,
               std::unordered_set<podio::ObjectID>& seenTo) {
  copyLinks(from, linkedTo, outLinks);
  for (const auto& [to, weight] : linkedTo) {
    if (seenTo.insert(to.getObjectID()).second) {
      outTo.push_back(to);
    }
  }
}

} // namespace k4reco::bibutils

#endif
