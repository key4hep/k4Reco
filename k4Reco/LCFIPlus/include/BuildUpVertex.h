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
#ifndef K4RECO_BUILDUPOVERTEX_H
#define K4RECO_BUILDUPOVERTEX_H

#include "TrackSelector.h"

#include <edm4hep/ReconstructedParticleCollection.h>
#include <edm4hep/VertexCollection.h>

#include <Gaudi/Property.h>

#include <k4FWCore/Transformer.h>

struct BuildUpVertex final
    : k4FWCore::MultiTransformer<std::tuple<edm4hep::VertexCollection, edm4hep::VertexCollection>(
          const edm4hep::ReconstructedParticleCollection&, const edm4hep::VertexCollection&)> {
  BuildUpVertex(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;
  std::tuple<edm4hep::VertexCollection, edm4hep::VertexCollection>
  operator()(const edm4hep::ReconstructedParticleCollection& pfos,
             const edm4hep::VertexCollection& headers) const override;

  // parameters
  // vertex formation limits
  double _chi2thpri;
  double _chi2thsec;
  double _massth;
  double _posth;
  double _chi2orderinglimit;
  int _v0sel;

  // association parameters
  bool _doassoc;
  double _minimumdist;
  double _chi2ratio;

  // track cut parameters
  lcfiplus::TrackSelectorConfig _secVtxCfg; //!

  // AVF parameters
  bool _avf;           // flag AVF/chi2
  double _temperature; // AVF parameter

  // primary vertex finding parameters
  bool _beamspotConstraint;
  bool _beamspotSmearing;

  Gaudi::Property<double> m_trackMaxD0{this, "TrackMaxD0", 10., "Maximum D0 of the track"};
  Gaudi::Property<double> m_trackMaxZ0{this, "TrackMaxZ0", 20., "Maximum Z0 of the track"};
  Gaudi::Property<double> m_trackMinPt{this, "TrackMinPt", 0.1, "Minimum Pt of the track"};

  Gaudi::Property<double> m_trackMaxD0Err{this, "TrackMaxD0Err", 0.1, "Maximum D0 error of the track"};
  Gaudi::Property<double> m_trackMaxZ0Err{this, "TrackMaxZ0Err", 0.1, "Maximum Z0 error of the track"};

  Gaudi::Property<int> m_trackMinTpcHits{this, "TrackMinTpcHits", 999999, "Minimum TPC hits of the track"};
  Gaudi::Property<int> m_trackMinTpcHitsMinPt{this, "TrackMinTpcHitsMinPt", 999999,
                                              "Minimum TPC hits of the track at min Pt"};
  Gaudi::Property<int> m_trackMinFtdHits{this, "TrackMinFtdHits", 999999, "Minimum FTD hits of the track"};
  Gaudi::Property<int> m_trackMinVxdHits{this, "TrackMinVxdHits", 999999, "Minimum VXD hits of the track"};
  Gaudi::Property<int> m_trackMinVxdFtdHits{this, "TrackMinVxdFtdHits", 1, "Minimum VXD+FTD hits of the track"};

  Gaudi::Property<double> m_primaryChi2Threshold{this, "PrimaryChi2Threshold", 25.,
                                                 "Maximum Chi2 threshold for primary vertex"};
  Gaudi::Property<double> m_secondaryChi2Threshold{this, "SecondaryChi2Threshold", 9.,
                                                   "Maximum Chi2 threshold for secondary vertex"};
  Gaudi::Property<double> m_massThreshold{this, "MassThreshold", 10., "Maximum mass threshold for secondary vertex"};
  Gaudi::Property<double> m_minDistFromIP{this, "MinDistFromIP", 0.3, "Minimum distance from IP for secondary vertex"};
  Gaudi::Property<double> m_maxChi2ForDistOrder{this, "MaxChi2ForDistOrder", 1.0, "Maximum Chi2 for distance ordering"};

  Gaudi::Property<bool> m_assocIPTracks{this, "AssocIPTracks", true, "Flag for associating IP tracks"};
  Gaudi::Property<double> m_assocIPTracksMinDist{this, "AssocIPTracksMinDist", 0.,
                                                 "Minimum distance for associating IP tracks"};
  Gaudi::Property<double> m_assocIPTracksChi2RatioSecToPri{this, "AssocIPTracksChi2RatioSecToPri", 2.0,
                                                           "Chi2 ratio for associating IP tracks"};
  Gaudi::Property<bool> m_useV0Selection{this, "UseV0Selection", true, "Flag for using V0 selection"};

  Gaudi::Property<bool> m_useAVF{this, "UseAVF", true, "Flag for using AVF"};
  Gaudi::Property<double> m_aVFTemperature{this, "AVFTemperature", 5.0, "Temperature for AVF"};
  Gaudi::Property<bool> m_beamspotConstraint{this, "BeamspotConstraint", true, "Beamspot constraint for refitting"};
  Gaudi::Property<bool> m_beamspotSmearing{this, "BeamspotSmearing", true, "Beamspot smearing for refitting"};
};

#endif // K4RECO_BUILDUPOVERTEX_H
