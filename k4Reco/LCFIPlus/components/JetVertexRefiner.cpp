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
#include "JetVertexRefiner.h"

#include "EDM4hepTrackConverter.h"
#include "JetFinder.h"
#include "TrackSelector.h"
#include "VertexFinderSuehara.h"
#include "VertexFinderTearDown.h"
#include "VertexFitterSimple.h"
#include "algoEtc.h"

#include <TMath.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <ranges>

namespace {
using Tracks = lcfiplus::edm4hep_utils::ConvertedTracks;

int find_ndf(const double chi2, const double probability) {
  if (chi2 < 0 || probability < 0 || probability > 1)
    return 0;
  int lower = 0;
  int upper = 100;
  while (TMath::Prob(chi2, upper) < probability) {
    lower = upper;
    upper *= 2;
  }
  while (lower < upper - 1) {
    const int middle = (lower + upper) / 2;
    if (TMath::Prob(chi2, middle) < probability)
      lower = middle;
    else
      upper = middle;
  }
  return std::abs(TMath::Prob(chi2, lower) - probability) < std::abs(TMath::Prob(chi2, upper) - probability) ? lower
                                                                                                             : upper;
}

std::vector<std::unique_ptr<lcfiplus::Vertex>> convertVertices(const edm4hep::VertexCollection& input,
                                                               const Tracks& tracks) {
  std::vector<std::unique_ptr<lcfiplus::Vertex>> vertices;
  for (const auto vertex : input) {
    std::vector<const lcfiplus::Track*> associatedTracks;
    for (const auto& particle : vertex.getParticles())
      for (const auto& [track, pfo] : tracks.tracks)
        if (pfo == particle)
          associatedTracks.push_back(&track);

    double covariance[6]{};
    for (int index = 0; index < 6; ++index)
      covariance[index] = vertex.getCovMatrix()[index];
    auto converted = std::make_unique<lcfiplus::Vertex>(
        vertex.getChi2(), TMath::Prob(vertex.getChi2(), vertex.getNdf()), vertex.getPosition()[0],
        vertex.getPosition()[1], vertex.getPosition()[2], covariance, false);
    for (const auto* track : associatedTracks)
      converted->add(track);
    vertices.push_back(std::move(converted));
  }
  return vertices;
}

void copyPrimaryTracks(const edm4hep::Vertex& input, const Tracks& tracks, lcfiplus::Vertex& output) {
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

bool hasSameParticles(const lcfiplus::Vertex& candidate, const edm4hep::Vertex& input, const Tracks& tracks) {
  const auto particles = input.getParticles();
  if (candidate.getTracks().size() != particles.size())
    return false;
  for (const auto* candidateTrack : candidate.getTracks()) {
    const auto track = std::ranges::find_if(
        tracks.tracks, [candidateTrack](const auto& converted) { return &converted.first == candidateTrack; });
    if (track == tracks.tracks.end() || std::ranges::find(particles, track->second) == particles.end())
      return false;
  }
  return true;
}

std::vector<std::unique_ptr<lcfiplus::Vertex>> reproduceBuildUpVertices(const lcfiplus::Vertex& geometricPrimary,
                                                                        const Tracks& tracks) {
  auto* primary =
      new lcfiplus::Vertex(geometricPrimary.getChi2(), geometricPrimary.getProb(), geometricPrimary.getX(),
                           geometricPrimary.getY(), geometricPrimary.getZ(), geometricPrimary.getCov(), true);
  for (const auto* track : geometricPrimary.getTracks())
    primary->add(track, geometricPrimary.getChi2Track(track));

  lcfiplus::TrackSelectorConfig trackConfiguration;
  trackConfiguration.maxD0 = 10.;
  trackConfiguration.maxZ0 = 20.;
  trackConfiguration.minPt = .1;
  trackConfiguration.maxInnermostHitRadius = 1e10;
  trackConfiguration.maxD0Err = .1;
  trackConfiguration.maxZ0Err = .1;
  trackConfiguration.minTpcHits = 1;
  trackConfiguration.minTpcHitsMinPt = 999999;
  trackConfiguration.minFtdHits = 1;
  trackConfiguration.minVtxHits = 1;
  trackConfiguration.minVtxPlusFtdHits = 1;
  const auto selectedTracks = lcfiplus::TrackSelector{}(tracks.pointers, trackConfiguration, primary);

  lcfiplus::VertexFinderSuehara::VertexFinderSueharaConfig configuration;
  configuration.chi2th = 9.;
  configuration.massth = 10.;
  configuration.v0selVertex.minpos = .3;
  configuration.chi2orderinglimit = 1.;
  configuration.minimumdistIP = 0.;
  configuration.chi2ratioIP = 2.;
  configuration.avf = false;
  configuration.temperature = 5.;
  configuration.beamspotConstraint = true;
  configuration.beamspotSmearing = false;

  std::vector<lcfiplus::Vertex*> vertices;
  std::vector<lcfiplus::Vertex*> v0Vertices;
  lcfiplus::VertexFinderSuehara::buildUp(selectedTracks, vertices, v0Vertices, 25., configuration, &primary);
  lcfiplus::VertexFinderSuehara::associateIPTracks(vertices, primary, configuration);

  std::vector<std::unique_ptr<lcfiplus::Vertex>> anchors;
  anchors.reserve(vertices.size());
  for (auto* vertex : vertices)
    anchors.emplace_back(vertex);
  for (auto* vertex : v0Vertices)
    delete vertex;
  delete primary;
  return anchors;
}

std::unique_ptr<lcfiplus::Vertex> reproducePrimaryVertex(const Tracks& tracks) {
  lcfiplus::TrackSelectorConfig configuration;
  configuration.maxD0 = 20.;
  configuration.maxZ0 = 20.;
  configuration.minVtxPlusFtdHits = 1;
  configuration.minTpcHits = 999999;
  configuration.minTpcHitsMinPt = 999999;
  configuration.minFtdHits = 999999;
  configuration.minVtxHits = 999999;

  const auto selectedTracks = lcfiplus::TrackSelector{}(tracks.pointers, configuration);
  auto primary = std::unique_ptr<lcfiplus::Vertex>{lcfiplus::findPrimaryVertex(selectedTracks, 25., true, false)};
  return primary;
}
} // namespace

JetVertexRefiner::JetVertexRefiner(const std::string& name, ISvcLocator* serviceLocator)
    : MultiTransformer(name, serviceLocator,
                       {KeyValues("PFOCollection", {"PFOsFromJets"}),
                        KeyValues("InputJetCollectionName", {"VertexJets"}),
                        KeyValues("PrimaryVertexCollectionName", {"PrimaryVertices"}),
                        KeyValues("InputVertexCollectionName", {"BuildUpVertices"}),
                        KeyValues("V0VertexCollectionName", {"BuildUpVertices_V0"})},
                       {KeyValues("OutputJetCollectionName", {"RefinedVertexJets"}),
                        KeyValues("OutputVertexCollectionName", {"RefinedVertices"})}) {}

std::tuple<edm4hep::ReconstructedParticleCollection, edm4hep::VertexCollection> JetVertexRefiner::operator()(
    const edm4hep::ReconstructedParticleCollection& pfos, const edm4hep::ReconstructedParticleCollection& inputJets,
    const edm4hep::VertexCollection& primaryVertices, const edm4hep::VertexCollection& inputVertices,
    const edm4hep::VertexCollection& v0Vertices) const {
  edm4hep::ReconstructedParticleCollection outputJets;
  edm4hep::VertexCollection outputVertices;
  const auto tracks = lcfiplus::edm4hep_utils::convertTracks(pfos);
  const auto v0 = convertVertices(v0Vertices, tracks);
  const auto geometricPrimary = reproducePrimaryVertex(tracks);
  const auto geometricAnchors = geometricPrimary ? reproduceBuildUpVertices(*geometricPrimary, tracks)
                                                 : std::vector<std::unique_ptr<lcfiplus::Vertex>>{};

  std::unique_ptr<lcfiplus::Vertex> primary;
  if (!primaryVertices.empty()) {
    const auto vertex = primaryVertices[0];
    double covariance[6]{};
    for (int index = 0; index < 6; ++index)
      covariance[index] = vertex.getCovMatrix()[index];
    primary = std::make_unique<lcfiplus::Vertex>(vertex.getChi2(), TMath::Prob(vertex.getChi2(), vertex.getNdf()),
                                                 vertex.getPosition()[0], vertex.getPosition()[1],
                                                 vertex.getPosition()[2], covariance, true);
    copyPrimaryTracks(vertex, tracks, *primary);
  }

  for (const auto inputJet : inputJets) {
    lcfiplus::Jet jet;
    for (const auto& particle : inputJet.getParticles())
      for (const auto& [track, pfo] : tracks.tracks)
        if (pfo == particle)
          jet.add(&track);

    // The refiner takes ownership of vertices it changes or produces.  Start with
    // event-local copies for every input jet, then release only the vertices that
    // were assigned to this jet before passing them to the legacy routine.
    auto secondary = convertVertices(inputVertices, tracks);
    std::vector<lcfiplus::Vertex*> secondaryPointers;
    for (const auto& vertex : secondary)
      secondaryPointers.push_back(vertex.get());

    lcfiplus::JetVec jets{&jet};
    std::vector<std::vector<lcfiplus::Vertex*>> verticesByJet;
    std::vector<std::vector<const lcfiplus::Track*>> residualTracks;
    lcfiplus::algoEtc::connectVerticesToJets(jets, secondaryPointers, verticesByJet, residualTracks, primary.get());

    std::vector<const lcfiplus::Vertex*> anchorPointers;
    anchorPointers.reserve(verticesByJet[0].size());
    for (const auto* vertex : verticesByJet[0]) {
      const auto found =
          std::ranges::find_if(secondary, [vertex](const auto& candidate) { return candidate.get() == vertex; });
      if (found == secondary.end()) {
        anchorPointers.push_back(vertex);
      } else {
        const auto index = static_cast<std::size_t>(found - secondary.begin());
        const auto anchor =
            std::ranges::find_if(geometricAnchors, [&inputVertices, &tracks, index](const auto& candidate) {
              return hasSameParticles(*candidate, inputVertices[index], tracks);
            });
        anchorPointers.push_back(anchor != geometricAnchors.end() ? anchor->get() : vertex);
      }
    }

    for (auto& vertex : secondary)
      if (std::find(verticesByJet[0].begin(), verticesByJet[0].end(), vertex.get()) != verticesByJet[0].end())
        vertex.release();

    std::vector<const lcfiplus::Vertex*> vertexPointers;
    std::vector<const lcfiplus::Vertex*> v0Pointers;
    for (const auto* vertex : verticesByJet[0])
      vertexPointers.push_back(vertex);
    for (const auto& vertex : v0)
      v0Pointers.push_back(vertex.get());
    lcfiplus::VertexVec vertices = vertexPointers;
    lcfiplus::VertexVec v0VerticesForJet = v0Pointers;
    lcfiplus::VertexVec anchors = anchorPointers;
    lcfiplus::TrackVec residual = residualTracks[0];
    lcfiplus::VertexFinderSuehara::VertexFinderSueharaConfig configuration;
    auto singleTrackVertices = lcfiplus::VertexFinderSuehara::makeSingleTrackVertices(
        anchors, residual, v0VerticesForJet, geometricPrimary ? geometricPrimary.get() : primary.get(), configuration);
    lcfiplus::VertexFinderSuehara::recombineVertices(verticesByJet[0], singleTrackVertices);

    auto outputJet = outputJets.create();
    outputJet.setMomentum(inputJet.getMomentum());
    outputJet.setEnergy(inputJet.getEnergy());
    outputJet.setMass(inputJet.getMass());
    for (const auto& particle : inputJet.getParticles())
      outputJet.addToParticles(particle);

    for (const auto* vertex : verticesByJet[0]) {
      auto outputVertex = outputVertices.create();
      outputVertex.setPosition(
          {static_cast<float>(vertex->getX()), static_cast<float>(vertex->getY()), static_cast<float>(vertex->getZ())});
      outputVertex.setChi2(vertex->getChi2());
      outputVertex.setNdf(find_ndf(vertex->getChi2(), vertex->getProb()));
      outputVertex.setCovMatrix({vertex->getCov()[0], vertex->getCov()[1], vertex->getCov()[2], vertex->getCov()[3],
                                 vertex->getCov()[4], vertex->getCov()[5]});
      for (const auto* track : vertex->getTracks())
        for (const auto& [convertedTrack, pfo] : tracks.tracks)
          if (&convertedTrack == track)
            outputVertex.addToParticles(pfo);
      delete vertex;
    }
  }
  return {std::move(outputJets), std::move(outputVertices)};
}

DECLARE_COMPONENT(JetVertexRefiner)
