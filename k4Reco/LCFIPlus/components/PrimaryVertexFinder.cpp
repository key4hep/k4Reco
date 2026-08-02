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
#include "PrimaryVertexFinder.h"

#include "EDM4hepTrackConverter.h"
#include "VertexFinderTearDown.h"

#include <edm4hep/ReconstructedParticleCollection.h>

PrimaryVertexFinder::PrimaryVertexFinder(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(name, svcLoc,
                       {
                           KeyValues("PFOCollection", {"PFOsFromJets"}),
                       },
                       {
                           KeyValues("VertexNames", {"PrimaryVertices"}),
                       }) {}

// Copied and pasted from k4EDM4hep2LcioConv
static int find_ndf(double chi2, double prob) {
  if (chi2 < 0 || prob < 0 || prob > 1) {
    throw std::invalid_argument("Invalid input for find_ndf. Either chi2 < 0, prob < 0 or prob > 1 in a LCIO Vertex.");
  }
  int lower = 0;
  // Initial guess for the upper bound. If it's not enough, it will be increased
  int upper = 100;
  while (TMath::Prob(chi2, upper) < prob) {
    lower = upper;
    upper *= 2;
  }
  while (lower < upper - 1) {
    int mid = (lower + upper) / 2;
    if (TMath::Prob(chi2, mid) < prob) {
      lower = mid;
    } else {
      upper = mid;
    }
  }
  if (std::abs(TMath::Prob(chi2, lower) - prob) < std::abs(TMath::Prob(chi2, upper) - prob)) {
    return lower;
  } else {
    return upper;
  }
}

StatusCode PrimaryVertexFinder::initialize() {
  _priVtxCfg.maxD0 = m_trackMaxD0;
  _priVtxCfg.maxZ0 = m_trackMaxZ0;
  _priVtxCfg.minVtxPlusFtdHits = m_trackMinVtxPlusFtdHits;
  _priVtxCfg.minTpcHits = m_trackMinTpcHits;
  _priVtxCfg.minTpcHitsMinPt = m_trackMinTpcHitsMinPt;
  _priVtxCfg.minFtdHits = m_trackMinFtdHits;
  _priVtxCfg.minVtxHits = m_trackMinVxdHits;
  _chi2th = m_chi2th;
  _beamspotConstraint = m_beamspotConstraint;
  _beamspotSmearing = m_beamspotSmearing;
  auto* globals = lcfiplus::Globals::Instance();
  globals->setBField(m_magneticField);
  globals->setBeamSizeX(m_beamSizeX);
  globals->setBeamSizeY(m_beamSizeY);
  globals->setBeamSizeZ(m_beamSizeZ);
  return StatusCode::SUCCESS;
}

std::tuple<edm4hep::VertexCollection>
PrimaryVertexFinder::operator()(const edm4hep::ReconstructedParticleCollection& pfos) const {

  auto converted = lcfiplus::edm4hep_utils::convertTracks(pfos);
  if (converted.skippedNoIPState != 0)
    warning() << "Skipped " << converted.skippedNoIPState << " track(s) without an AtIP track state." << endmsg;

  const auto passedTracks = lcfiplus::TrackSelector()(converted.pointers, _priVtxCfg);
  debug() << "PrimaryVertexFinder track selection: " << passedTracks.size() << "/" << converted.tracks.size()
          << " accepted." << endmsg;

  const lcfiplus::Vertex* vtx = findPrimaryVertex(passedTracks, _chi2th, _beamspotConstraint, _beamspotSmearing);
  auto vtxColl = edm4hep::VertexCollection();
  if (!vtx)
    return std::make_tuple(std::move(vtxColl));

  for (const auto* vertexTrack : vtx->getTracks()) {
    const auto it = std::ranges::find_if(converted.tracks,
                                         [vertexTrack](const auto& track) { return &track.first == vertexTrack; });
    if (it != converted.tracks.end())
      debug() << "Primary track PFO=" << it->second.getObjectID().index << " chi2=" << vtx->getChi2Track(vertexTrack)
              << endmsg;
  }

  auto vtxObj = vtxColl.create();
  vtxObj.setPrimary(true);
  vtxObj.setChi2(vtx->getChi2());
  vtxObj.setNdf(find_ndf(vtx->getChi2(), vtx->getProb()));
  vtxObj.setPosition(
      {static_cast<float>(vtx->getX()), static_cast<float>(vtx->getY()), static_cast<float>(vtx->getZ())});
  vtxObj.setCovMatrix(
      {vtx->getCov()[0], vtx->getCov()[1], vtx->getCov()[2], vtx->getCov()[3], vtx->getCov()[4], vtx->getCov()[5]});

  for (const auto* vertexTrack : vtx->getTracks()) {
    const auto it = std::ranges::find_if(converted.tracks,
                                         [vertexTrack](const auto& track) { return &track.first == vertexTrack; });
    if (it != converted.tracks.end())
      vtxObj.addToParticles(it->second);
  }

  return std::make_tuple(std::move(vtxColl));
}

DECLARE_COMPONENT(PrimaryVertexFinder)
