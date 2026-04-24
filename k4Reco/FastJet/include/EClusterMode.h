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
#ifndef K4RECO_FASTJET_ECLUSTERMODE_H
#define K4RECO_FASTJET_ECLUSTERMODE_H

#include <iostream>

// The enum, name and value of the enum for the Cluster Mode
namespace k4Reco::FastJet {
enum EClusterMode {
  NONE = 0,
  FJ_exclusive_yCut = 1,     // exclusive clustering mode implemented in FastJet
  FJ_exclusive_nJets = 2,    // exclusive clustering mode implemented in FastJet
  FJ_inclusive = 4,          // inclusive "-"
  OWN_inclusiveIteration = 8 // use FJ inclusive Clustering, but iterate until we have the desired number of jets
};
std::ostream& operator<<(std::ostream& out, EClusterMode& m) {
  switch (m) {
  case OWN_inclusiveIteration:
    out << "InclusiveIterativeNJets";
    break;
  case FJ_inclusive:
    out << "Inclusive";
    break;
  case FJ_exclusive_nJets:
    out << "ExclusiveNJets";
    break;
  case FJ_exclusive_yCut:
    out << "ExclusiveYCut";
    break;
  default:
    out << "unknown";
    break;
  }
  return out;
}

} // namespace k4Reco::FastJet
#endif
