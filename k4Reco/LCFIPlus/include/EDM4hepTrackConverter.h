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
#ifndef K4RECO_EDM4HEPTRACKCONVERTER_H
#define K4RECO_EDM4HEPTRACKCONVERTER_H

#include "lcfiplus.h"

#include <edm4hep/ReconstructedParticleCollection.h>
#include <edm4hep/TrackState.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <utility>
#include <vector>

namespace lcfiplus::edm4hep_utils {

struct ConvertedTracks {
  std::vector<std::pair<Track, edm4hep::ReconstructedParticle>> tracks;
  std::vector<const Track*> pointers;
  unsigned int skippedNoIPState{0};
};

inline std::vector<edm4hep::ReconstructedParticle>
sortedPFOsByEnergy(const edm4hep::ReconstructedParticleCollection& pfos) {
  std::vector<edm4hep::ReconstructedParticle> sorted{pfos.begin(), pfos.end()};
  // LCIOStorer presents every LCFIPlus algorithm with PFOs sorted in
  // descending energy order. Several vertex candidate tie-breaks depend on
  // that order, so preserve it in the EDM4hep bridge.
  std::sort(sorted.begin(), sorted.end(),
            [](const auto& first, const auto& second) { return first.getEnergy() > second.getEnergy(); });
  return sorted;
}

inline ConvertedTracks convertTracks(const edm4hep::ReconstructedParticleCollection& pfos) {
  ConvertedTracks converted;
  const auto sortedPfos = sortedPFOsByEnergy(pfos);
  std::size_t numberOfTracks = 0;
  for (const auto& pfo : sortedPfos) {
    if (pfo.getCharge() != 0.0f)
      numberOfTracks += pfo.tracks_size();
  }
  converted.tracks.reserve(numberOfTracks);
  converted.pointers.reserve(numberOfTracks);

  for (const auto& pfo : sortedPfos) {
    if (pfo.getCharge() == 0.0f)
      continue;

    // Match LCIOStorer: Pandora can attach parent tracks to a PFO, so retain
    // only the helix whose momentum best agrees with the PFO momentum.
    edm4hep::Track selectedTrack{};
    bool hasSelectedTrack = false;
    edm4hep::TrackState selectedState{};
    double bestMomentumDifference = 1.e300;
    const auto momentum = pfo.getMomentum();
    const double pfoMomentum = std::hypot(momentum[0], momentum[1], momentum[2]);
    const double eB = std::abs(lcfiplus::Globals::Instance()->getBField()) * 2.99792458e-4;
    for (const auto& candidate : pfo.getTracks()) {
      edm4hep::TrackState state{};
      bool hasIPState = false;
      for (const auto& candidateState : candidate.getTrackStates()) {
        if (candidateState.location == edm4hep::TrackState::AtIP) {
          state = candidateState;
          hasIPState = true;
          break;
        }
      }
      if (!hasIPState)
        continue;
      const double omega = state.omega;
      const double candidateMomentum =
          omega == 0. ? 1.e300 : eB * std::sqrt(1. + state.tanLambda * state.tanLambda) / std::abs(omega);
      const double difference = std::abs(candidateMomentum - pfoMomentum);
      if (difference < bestMomentumDifference) {
        bestMomentumDifference = difference;
        // Relation ranges yield handles by value. Keep our own handle instead
        // of a pointer to the loop variable, which becomes dangling once the
        // range iteration advances.
        selectedTrack = candidate;
        hasSelectedTrack = true;
        selectedState = state;
      }
    }
    if (!hasSelectedTrack) {
      ++converted.skippedNoIPState;
      continue;
    }
    {
      const auto& edmTrack = selectedTrack;
      const auto& state = selectedState;

      Track track;
      track.SetPxPyPzE(pfo.getMomentum()[0], pfo.getMomentum()[1], pfo.getMomentum()[2], pfo.getEnergy());
      track.setCharge(pfo.getCharge());
      track.setPDG(pfo.getPDG());
      // k4EDM4hep2LcioConv derives LCIO's radiusOfInnermostHit from the
      // AtFirstHit reference point. Reproduce that bridge exactly so that
      // LCFIPlus track selections have the same input in both workflows.
      double innermostHitRadius = -1.0;
      for (const auto& trackState : edmTrack.getTrackStates()) {
        if (trackState.location == edm4hep::TrackState::AtFirstHit) {
          const auto referencePoint = trackState.referencePoint;
          innermostHitRadius = std::hypot(referencePoint.x, referencePoint.y);
          break;
        }
      }
      track.setRadiusOfInnermostHit(innermostHitRadius);
      track.setHelix(state.D0, state.Z0, state.phi, state.omega, state.tanLambda);
      track.setChi2(edmTrack.getChi2());
      track.setNdf(edmTrack.getNdf());

      std::array<double, tpar::covN> covariance{};
      // EDM4hep stores the first five helix parameters in the same lower-triangular
      // order as LCFIPlus: d0, phi, omega, z0, tan(lambda).
      for (unsigned int i = 0; i < covariance.size(); ++i)
        covariance[i] = state.covMatrix[i];
      track.setCovMatrix(covariance.data());

      std::array<int, tpar::hitN> hits{};
      // CLDConfig uses LCFIPlus' CLICdet hit ordering. EDM4hep preserves the
      // LCIO subdetector-hit vector, whose even entries hold fitted-hit counts.
      // Keep the same mapping as LCIOStorer::storeTracks(..., ordering=2).
      const auto subdetectorHits = edmTrack.getSubdetectorHitNumbers();
      if (subdetectorHits.size() >= 11) {
        hits[tpar::VTX] = subdetectorHits[0];
        hits[tpar::FTD] = subdetectorHits[2];
        hits[tpar::TPC] = subdetectorHits[4] + subdetectorHits[6] + subdetectorHits[8] + subdetectorHits[10];
      } else {
        // The generic EDM4hep fallback remains useful for inputs not produced
        // by the Marlin LCIO converter.
        hits[tpar::VTX] = static_cast<int>(edmTrack.trackerHits_size());
      }
      track.setTrackHits(hits.data());
      std::array<double, tpar::caloN> caloEnergy{};
      track.setCaloEdep(caloEnergy.data());

      converted.tracks.emplace_back(std::move(track), pfo);
      converted.pointers.push_back(&converted.tracks.back().first);
    }
  }

  return converted;
}

} // namespace lcfiplus::edm4hep_utils

#endif // K4RECO_EDM4HEPTRACKCONVERTER_H
