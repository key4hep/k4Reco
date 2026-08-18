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
#include "JetClustering.h"

#include "EDM4hepTrackConverter.h"
#include "JetFinder.h"
#include "VertexFitterSimple.h"

#include <TMath.h>

#include <algorithm>
#include <array>
#include <memory>
#include <ranges>
#include <utility>
#include <vector>

namespace {
struct ConvertedNeutrals {
  std::vector<std::pair<lcfiplus::Neutral, edm4hep::ReconstructedParticle>> particles;
  std::vector<const lcfiplus::Neutral*> pointers;
};

ConvertedNeutrals convert_neutrals(const edm4hep::ReconstructedParticleCollection& pfos) {
  ConvertedNeutrals converted;
  const auto sortedPfos = lcfiplus::edm4hep_utils::sortedPFOsByEnergy(pfos);
  const auto numberOfNeutrals =
      std::count_if(sortedPfos.begin(), sortedPfos.end(), [](const auto& pfo) { return pfo.getCharge() == 0.f; });
  // Pointers handed to JetFinder refer to particle storage. Reserve before
  // constructing it so vector growth cannot invalidate them.
  converted.particles.reserve(numberOfNeutrals);
  converted.pointers.reserve(numberOfNeutrals);
  for (const auto& pfo : sortedPfos) {
    if (pfo.getCharge() != 0.f)
      continue;
    lcfiplus::Neutral neutral;
    neutral.SetPxPyPzE(pfo.getMomentum()[0], pfo.getMomentum()[1], pfo.getMomentum()[2], pfo.getEnergy());
    neutral.setPDG(pfo.getPDG());
    std::array<double, lcfiplus::tpar::caloN> deposits{};
    neutral.setCaloEdep(deposits.data());
    converted.particles.emplace_back(std::move(neutral), pfo);
    converted.pointers.push_back(&converted.particles.back().first);
  }
  return converted;
}

std::vector<std::unique_ptr<lcfiplus::Vertex>>
convert_vertices(const edm4hep::VertexCollection& input, const lcfiplus::edm4hep_utils::ConvertedTracks& tracks) {
  std::vector<std::unique_ptr<lcfiplus::Vertex>> converted;
  converted.reserve(input.size());
  for (const auto& edmVertex : input) {
    double covariance[6]{};
    for (int i = 0; i < 6; ++i)
      covariance[i] = edmVertex.getCovMatrix()[i];
    auto vertex = std::make_unique<lcfiplus::Vertex>(
        edmVertex.getChi2(), TMath::Prob(edmVertex.getChi2(), edmVertex.getNdf()), edmVertex.getPosition()[0],
        edmVertex.getPosition()[1], edmVertex.getPosition()[2], covariance, false);
    for (const auto& particle : edmVertex.getParticles()) {
      const auto found =
          std::ranges::find_if(tracks.tracks, [&particle](const auto& track) { return track.second == particle; });
      if (found != tracks.tracks.end())
        vertex->add(&found->first);
    }
    converted.push_back(std::move(vertex));
  }
  return converted;
}

void copy_primary_tracks(const edm4hep::Vertex& input, const lcfiplus::edm4hep_utils::ConvertedTracks& tracks,
                         lcfiplus::Vertex& output) {
  for (const auto& particle : input.getParticles()) {
    const auto found =
        std::ranges::find_if(tracks.tracks, [&particle](const auto& track) { return track.second == particle; });
    if (found != tracks.tracks.end()) {
      const auto chi2 =
          lcfiplus::VertexFitterSimple_V{}.getChi2(&output, &found->first, 0, lcfiplus::PointBase::PRIVTX);
      output.add(&found->first, chi2);
    }
  }
}
} // namespace

JetClustering::JetClustering(const std::string& name, ISvcLocator* svcLoc)
    : Transformer(name, svcLoc,
                  {KeyValues("PFOCollection", {"PFOsFromJets"}),
                   KeyValues("PrimaryVertexCollectionName", {"PrimaryVertices"}),
                   KeyValues("InputVertexCollectionName", {"BuildUpVertices"})},
                  {KeyValues("OutputJetCollectionName", {"VertexJets"})}) {}

StatusCode JetClustering::initialize() {
  // Match the initialization of CLDConfig's separate
  // JetClusteringAndRefiner LcfiplusProcessor. LCFIPlus keeps these values in
  // a process-wide singleton, so this must happen before event processing,
  // after PrimaryVertexFinder has initialized the first processor's globals.
  auto* globals = lcfiplus::Globals::Instance();
  globals->setBField(m_magneticField);
  globals->setBeamSizeX(m_beamSizeX);
  globals->setBeamSizeY(m_beamSizeY);
  globals->setBeamSizeZ(m_beamSizeZ);
  return StatusCode::SUCCESS;
}

edm4hep::ReconstructedParticleCollection
JetClustering::operator()(const edm4hep::ReconstructedParticleCollection& pfos,
                          const edm4hep::VertexCollection& primaryVertices,
                          const edm4hep::VertexCollection& secondaryVertices) const {
  edm4hep::ReconstructedParticleCollection output;
  auto tracks = lcfiplus::edm4hep_utils::convertTracks(pfos);
  auto neutrals = convert_neutrals(pfos);
  auto vertices = convert_vertices(secondaryVertices, tracks);
  std::vector<const lcfiplus::Vertex*> vertexPointers;
  vertexPointers.reserve(vertices.size());
  for (const auto& vertex : vertices)
    vertexPointers.push_back(vertex.get());

  std::unique_ptr<lcfiplus::Vertex> primary;
  if (!primaryVertices.empty()) {
    const auto& input = primaryVertices[0];
    double covariance[6]{};
    for (int i = 0; i < 6; ++i)
      covariance[i] = input.getCovMatrix()[i];
    primary = std::make_unique<lcfiplus::Vertex>(input.getChi2(), TMath::Prob(input.getChi2(), input.getNdf()),
                                                 input.getPosition()[0], input.getPosition()[1], input.getPosition()[2],
                                                 covariance, true);
    copy_primary_tracks(input, tracks, *primary);
  }

  lcfiplus::JetConfig config;
  config.nJet = m_nJetsRequested;
  config.Ycut = m_yCut;
  config.algo = m_jetAlgorithm;
  config.useBeamJets = m_useBeamJets;
  config.rParameter = m_rParameter;
  config.alphaParameter = m_alphaParameter;
  config.betaParameter = m_betaParameter;
  config.gammaParameter = m_gammaParameter;
  config.YaddVV = m_yAddedForJetVertexVertex;
  config.YaddVL = m_yAddedForJetVertexLepton;
  config.YaddLL = m_yAddedForJetLeptonLepton;
  config.useMuonID = m_useMuonID;
  config.muonIDExternal = m_muonIDExternal;
  config.muonIDMinEnergy = m_muonIDMinimumEnergy;
  config.muonIDMinD0Sig = m_muonIDMinimumD0Significance;
  config.muonIDMinZ0Sig = m_muonIDMinimumZ0Significance;
  config.muonIDMaxDist = m_muonIDMaximum3DImpactParameter;
  config.muonIDMinProb = m_muonIDMinimumProbability;

  lcfiplus::JetFinder finder(config, primary.get());
  const lcfiplus::TrackVec& trackPointers = tracks.pointers;
  const lcfiplus::NeutralVec& neutralPointers = neutrals.pointers;
  const lcfiplus::VertexVec& secondaryPointers = vertexPointers;
  auto jets = finder.run(trackPointers, neutralPointers, secondaryPointers);
  for (const auto* jet : jets) {
    auto edmJet = output.create();
    edmJet.setMomentum({static_cast<float>(jet->Px()), static_cast<float>(jet->Py()), static_cast<float>(jet->Pz())});
    edmJet.setEnergy(jet->E());
    edmJet.setMass(jet->M());
    std::vector<edm4hep::ReconstructedParticle> constituents;
    const auto addConstituent = [&edmJet, &constituents](const edm4hep::ReconstructedParticle& particle) {
      if (std::ranges::find(constituents, particle) == constituents.end()) {
        edmJet.addToParticles(particle);
        constituents.push_back(particle);
      }
    };
    for (const auto* track : jet->getAllTracks()) {
      const auto found =
          std::ranges::find_if(tracks.tracks, [track](const auto& entry) { return &entry.first == track; });
      if (found != tracks.tracks.end())
        addConstituent(found->second);
    }
    for (const auto* neutral : jet->getNeutrals()) {
      const auto found =
          std::ranges::find_if(neutrals.particles, [neutral](const auto& entry) { return &entry.first == neutral; });
      if (found != neutrals.particles.end())
        addConstituent(found->second);
    }
    delete jet;
  }
  return output;
}

DECLARE_COMPONENT(JetClustering)
