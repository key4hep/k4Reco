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
#include "BuildUpVertex.h"

#include "EDM4hepTrackConverter.h"
#include "VertexFinderSuehara.h"
#include "VertexFitterSimple.h"

#include <edm4hep/ReconstructedParticleCollection.h>

#include <ranges>
#include <string>

namespace {
int find_ndf(const double chi2, const double prob) {
  if (chi2 < 0 || prob < 0 || prob > 1)
    return 0;
  int lower = 0;
  int upper = 100;
  while (TMath::Prob(chi2, upper) < prob) {
    lower = upper;
    upper *= 2;
  }
  while (lower < upper - 1) {
    const int middle = (lower + upper) / 2;
    if (TMath::Prob(chi2, middle) < prob)
      lower = middle;
    else
      upper = middle;
  }
  return std::abs(TMath::Prob(chi2, lower) - prob) < std::abs(TMath::Prob(chi2, upper) - prob) ? lower : upper;
}

void copy_vertices(const std::vector<lcfiplus::Vertex*>& source,
                   const lcfiplus::edm4hep_utils::ConvertedTracks& converted, edm4hep::VertexCollection& destination) {
  for (const auto* vertex : source) {
    auto output = destination.create();
    output.setChi2(vertex->getChi2());
    output.setNdf(find_ndf(vertex->getChi2(), vertex->getProb()));
    output.setPosition(
        {static_cast<float>(vertex->getX()), static_cast<float>(vertex->getY()), static_cast<float>(vertex->getZ())});
    output.setCovMatrix({vertex->getCov()[0], vertex->getCov()[1], vertex->getCov()[2], vertex->getCov()[3],
                         vertex->getCov()[4], vertex->getCov()[5]});
    for (const auto* vertexTrack : vertex->getTracks()) {
      const auto it = std::ranges::find_if(converted.tracks,
                                           [vertexTrack](const auto& track) { return &track.first == vertexTrack; });
      if (it != converted.tracks.end())
        output.addToParticles(it->second);
    }
  }
}
} // namespace

BuildUpVertex::BuildUpVertex(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(name, svcLoc,
                       {
                           KeyValues("PFOCollection", {"PFOsFromJets"}),
                           KeyValues("PrimaryVertexCollectionName", {"PrimaryVertices"}),
                       },
                       {
                           KeyValues("BuildUpVertexCollectionName", {"BuildUpVertices"}),
                           KeyValues("V0VertexCollectionName", {"BuildUpVertices_V0"}),
                       }) {}

StatusCode BuildUpVertex::initialize() {
  _secVtxCfg.maxD0 = m_trackMaxD0;
  _secVtxCfg.maxZ0 = m_trackMaxZ0;
  _secVtxCfg.minPt = m_trackMinPt;
  _secVtxCfg.maxInnermostHitRadius = 1e10;
  _secVtxCfg.maxD0Err = m_trackMaxD0Err;
  _secVtxCfg.maxZ0Err = m_trackMaxZ0Err;
  _secVtxCfg.minTpcHits = m_trackMinTpcHits;
  _secVtxCfg.minTpcHitsMinPt = m_trackMinTpcHitsMinPt;
  _secVtxCfg.minFtdHits = m_trackMinFtdHits;
  _secVtxCfg.minVtxHits = m_trackMinVxdHits;
  _secVtxCfg.minVtxPlusFtdHits = m_trackMinVxdFtdHits;
  _chi2thpri = m_primaryChi2Threshold;
  _chi2thsec = m_secondaryChi2Threshold;
  _massth = m_massThreshold;
  _posth = m_minDistFromIP;
  _chi2orderinglimit = m_maxChi2ForDistOrder;
  _doassoc = m_assocIPTracks;
  _minimumdist = m_assocIPTracksMinDist;
  _chi2ratio = m_assocIPTracksChi2RatioSecToPri;
  _v0sel = m_useV0Selection;
  _avf = m_useAVF;
  _temperature = m_aVFTemperature;
  _beamspotConstraint = m_beamspotConstraint;
  _beamspotSmearing = m_beamspotSmearing;
  return StatusCode::SUCCESS;
}

std::tuple<edm4hep::VertexCollection, edm4hep::VertexCollection>
BuildUpVertex::operator()(const edm4hep::ReconstructedParticleCollection& pfos,
                          const edm4hep::VertexCollection& primaryVertices) const {
  auto secondaryVertices = edm4hep::VertexCollection();
  auto v0Vertices = edm4hep::VertexCollection();
  if (primaryVertices.empty()) {
    warning() << "No primary vertex supplied; skipping secondary vertex finding." << endmsg;
    return std::make_tuple(std::move(secondaryVertices), std::move(v0Vertices));
  }

  auto converted = lcfiplus::edm4hep_utils::convertTracks(pfos);
  if (converted.skippedNoIPState != 0)
    warning() << "Skipped " << converted.skippedNoIPState << " track(s) without an AtIP track state." << endmsg;

  const auto& primaryVertex = primaryVertices[0];
  double covariance[6]{};
  for (int i = 0; i < 6; ++i)
    covariance[i] = primaryVertex.getCovMatrix()[i];
  auto* primary = new lcfiplus::Vertex(
      primaryVertex.getChi2(), TMath::Prob(primaryVertex.getChi2(), primaryVertex.getNdf()),
      primaryVertex.getPosition()[0], primaryVertex.getPosition()[1], primaryVertex.getPosition()[2], covariance, true);
  // The legacy processor copies the primary vertex, including its tracks.
  // They are subsequently considered by associateIPTracks. EDM4hep stores the
  // associated PFOs, so recover their corresponding transient LCFI tracks.
  for (const auto& primaryParticle : primaryVertex.getParticles()) {
    const auto primaryIndex = primaryParticle.getObjectID().index;
    const auto track = std::ranges::find_if(converted.tracks, [primaryIndex](const auto& candidate) {
      return candidate.second.getObjectID().index == primaryIndex;
    });
    if (track != converted.tracks.end()) {
      // VertexFitterSimple records this value in Vertex::_chi2Tracks. It is
      // not representable in edm4hep::Vertex, but associateIPTracks relies on
      // it to decide whether a primary track belongs to a secondary vertex.
      const auto chi2 =
          lcfiplus::VertexFitterSimple_V{}.getChi2(primary, &track->first, 0, lcfiplus::PointBase::PRIVTX);
      primary->add(&track->first, chi2);
    }
  }
  const auto selectedTracks = lcfiplus::TrackSelector()(converted.pointers, _secVtxCfg, primary);

  lcfiplus::VertexFinderSuehara::VertexFinderSueharaConfig cfg;
  cfg.chi2th = _chi2thsec;
  cfg.massth = _massth;
  cfg.v0selVertex.minpos = _posth;
  cfg.chi2orderinglimit = _chi2orderinglimit;
  if (!_v0sel) {
    cfg.v0selTrack.setNoV0Cut();
    cfg.v0selVertex.setNoV0Cut();
  }
  cfg.minimumdistIP = _minimumdist;
  cfg.chi2ratioIP = _chi2ratio;
  cfg.avf = _avf;
  cfg.temperature = _temperature;
  cfg.beamspotConstraint = _beamspotConstraint;
  cfg.beamspotSmearing = _beamspotSmearing;

  std::vector<lcfiplus::Vertex*> lcfiplusSecondaryVertices;
  std::vector<lcfiplus::Vertex*> lcfiplusV0Vertices;
  lcfiplus::VertexFinderSuehara::buildUp(selectedTracks, lcfiplusSecondaryVertices, lcfiplusV0Vertices, _chi2thpri, cfg,
                                         &primary);
  if (_doassoc) {
    if (_avf)
      lcfiplus::VertexFinderSuehara::associateIPTracksAVF(lcfiplusSecondaryVertices, primary, cfg);
    else
      lcfiplus::VertexFinderSuehara::associateIPTracks(lcfiplusSecondaryVertices, primary, cfg);
  }

  copy_vertices(lcfiplusSecondaryVertices, converted, secondaryVertices);
  copy_vertices(lcfiplusV0Vertices, converted, v0Vertices);
  for (auto* vertex : lcfiplusSecondaryVertices)
    delete vertex;
  for (auto* vertex : lcfiplusV0Vertices)
    delete vertex;
  delete primary;

  return std::make_tuple(std::move(secondaryVertices), std::move(v0Vertices));
}

DECLARE_COMPONENT(BuildUpVertex)
