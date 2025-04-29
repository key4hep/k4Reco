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
#ifndef K4RECO_CLONESANDSPLITTRACKSFINDER_H
#define K4RECO_CLONESANDSPLITTRACKSFINDER_H 1

#include <edm4hep/TrackCollection.h>

#include <k4FWCore/Transformer.h>

#include <Gaudi/Property.h>

struct ClonesAndSplitTracksFinder final
    : k4FWCore::Transformer<edm4hep::TrackCollection(const edm4hep::TrackCollection&)> {
  ClonesAndSplitTracksFinder(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;

  edm4hep::TrackCollection operator()(const edm4hep::TrackCollection&) const override;

  // Checks for overlapping hits
  size_t overlappingHits(const edm4hep::Track&, const edm4hep::Track&) const;

  // Picks up the best track between two clones (based on chi2 and length requirements)
  edm4hep::Track bestInClones(const edm4hep::Track&, const edm4hep::Track&, size_t) const;

  // Removes doubles (from clone treatments and track merging) and filters multiple connections (clones and mergeable
  // tracks treated differently)
  void filterClonesAndMergedTracks(std::multimap<size_t, std::pair<size_t, edm4hep::Track>>&,
                                   const edm4hep::TrackCollection&, edm4hep::TrackCollection&, bool) const;

  // Contains the whole clone skimming procedure (calls bestInClones and filterClonesAndMergedTracks(bool true))
  edm4hep::TrackCollection removeClones(const edm4hep::TrackCollection&) const;

  Gaudi::Property<bool> m_MSOn{this, "MultipleScatteringOn", true, "Use MultipleScattering in Fit"};
  Gaudi::Property<bool> m_ElossOn{this, "EnergyLossOn", true, "Use Energy Loss in Fit"};
  Gaudi::Property<bool> m_SmoothOn{this, "SmoothOn", false, "Smooth All Measurement Sites in Fit"};
  Gaudi::Property<bool> m_extrapolateForward{this, "extrapolateForward", true,
                                             "if true extrapolation in the forward direction (in-out), otherwise "
                                             "backward (out-in)"};
  Gaudi::Property<double> m_minPt{this, "minTrackPt", 1.0, "minimum track pt for merging (in GeV/c)"};
  Gaudi::Property<double> m_maxSignificanceTheta{this, "maxSignificanceTheta", 3.0,
                                                 "maximum significance separation in tanLambda"};
  Gaudi::Property<double> m_maxSignificancePhi{this, "maxSignificancePhi", 3.0,
                                               "maximum significance separation in phi"};
  Gaudi::Property<double> m_maxSignificancePt{this, "maxSignificancePt", 2.0, "maximum significance separation in pt"};
  Gaudi::Property<bool> m_mergeSplitTracks{this, "mergeSplitTracks", false,
                                           "if true, the merging of split tracks is performed"};
};

#endif

DECLARE_COMPONENT(ClonesAndSplitTracksFinder)
