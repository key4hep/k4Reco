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
#include "lcfiplus.h"

#include <algorithm>

#include "TMatrixDSym.h"
#include "TVectorD.h"

namespace lcfiplus {

Jet::Jet(const Track* track) : TLorentzVector(*track), _id(-1) { _tracks.push_back(track); }

Jet::Jet(const Neutral* neutral) : TLorentzVector(*neutral), _id(-1) { _neutrals.push_back(neutral); }

void Jet::add(const Jet& jet) {
  *this += jet;
  _tracks.insert(_tracks.end(), jet.getTracks().begin(), jet.getTracks().end());
  _neutrals.insert(_neutrals.end(), jet.getNeutrals().begin(), jet.getNeutrals().end());
  _vertices.insert(_vertices.end(), jet.getVertices().begin(), jet.getVertices().end());
}

Jet::Jet(const Jet& from, bool extractVertex) : TLorentzVector(from), _id(-1) {
  _tracks = from.getTracks();
  _neutrals = from.getNeutrals();
  if (extractVertex) {
    for (const auto* vertex : from.getVertices())
      _tracks.insert(_tracks.end(), vertex->getTracks().begin(), vertex->getTracks().end());
  } else {
    _vertices = from.getVertices();
  }
  _params = from.params();
}

double Jet::sphericity() const {
  TMatrixDSym sphericityMatrix(3);
  const TVector3 jetBoost = BoostVector();
  for (const auto* track : getAllTracks()) {
    TLorentzVector fourMomentum(*track);
    fourMomentum.Boost(-jetBoost);
    sphericityMatrix(0, 0) += fourMomentum.X() * fourMomentum.X();
    sphericityMatrix(0, 1) += fourMomentum.X() * fourMomentum.Y();
    sphericityMatrix(0, 2) += fourMomentum.X() * fourMomentum.Z();
    sphericityMatrix(1, 1) += fourMomentum.Y() * fourMomentum.Y();
    sphericityMatrix(1, 2) += fourMomentum.Y() * fourMomentum.Z();
    sphericityMatrix(2, 2) += fourMomentum.Z() * fourMomentum.Z();
  }
  const double normalisation = sphericityMatrix(0, 0) + sphericityMatrix(1, 1) + sphericityMatrix(2, 2);
  if (normalisation == 0.)
    return 0.;
  sphericityMatrix *= 1.5 / normalisation;
  TVectorD eigenvalues;
  sphericityMatrix.EigenVectors(eigenvalues);
  return eigenvalues(1) + eigenvalues(2);
}

std::vector<const Track*> Jet::getAllTracks(bool withoutV0) const {
  std::vector<const Track*> tracks = _tracks;
  for (const auto* vertex : _vertices) {
    // Identifying V0s requires the event-store primary vertex, which is not
    // available in the EDM4hep-native data model. The native vertex-refiner
    // performs its V0 selection explicitly before this method is used.
    if (withoutV0)
      continue;
    tracks.insert(tracks.end(), vertex->getTracks().begin(), vertex->getTracks().end());
  }
  return tracks;
}

std::vector<const Vertex*> Jet::getVerticesForFT() const {
  std::vector<const Vertex*> vertices;
  for (const auto* vertex : _vertices) {
    if (vertex->getTracks().size() > 1)
      vertices.push_back(vertex);
  }
  std::sort(vertices.begin(), vertices.end(),
            [](const Vertex* lhs, const Vertex* rhs) { return lhs->getPos().Mag() < rhs->getPos().Mag(); });
  return vertices;
}

} // namespace lcfiplus
