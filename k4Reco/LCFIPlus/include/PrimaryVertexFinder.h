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
#ifndef K4RECO_PRIMARYVERTEXFINDER_H
#define K4RECO_PRIMARYVERTEXFINDER_H

#include "TrackSelector.h"

#include <Gaudi/Property.h>

#include <edm4hep/ReconstructedParticleCollection.h>
#include <edm4hep/VertexCollection.h>

#include <k4FWCore/Transformer.h>

struct PrimaryVertexFinder final : k4FWCore::MultiTransformer<std::tuple<edm4hep::VertexCollection>(
                                       const edm4hep::ReconstructedParticleCollection&)> {
  PrimaryVertexFinder(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;
  std::tuple<edm4hep::VertexCollection> operator()(const edm4hep::ReconstructedParticleCollection& pfos) const override;

private:
  // parameters
  double _chi2th{25.};
  bool _beamspotConstraint{true};
  bool _beamspotSmearing{true};

  // track cut parameters
  lcfiplus::TrackSelectorConfig _priVtxCfg;

  Gaudi::Property<double> m_trackMaxD0{this, "TrackMaxD0", 20., "Maximum D0 of the track"};
  Gaudi::Property<double> m_trackMaxZ0{this, "TrackMaxZ0", 20., "Maximum Z0 of the track"};
  Gaudi::Property<double> m_trackMinVtxPlusFtdHits{this, "TrackMinVtxFtdHits", 1., "Minimum Vtx+Ftd Hits of the track"};
  Gaudi::Property<double> m_trackMinTpcHits{this, "TrackMinTpcHits", 999999., "Minimum TPC Hits of the track"};
  Gaudi::Property<double> m_trackMinTpcHitsMinPt{this, "TrackMinTpcHitsMinPt", 999999.,
                                                 "Minimum TPC Hits Min Pt of the track"};
  Gaudi::Property<double> m_trackMinFtdHits{this, "TrackMinFtdHits", 999999., "Minimum FTD Hits of the track"};
  Gaudi::Property<double> m_trackMinVxdHits{this, "TrackMinVxdHits", 999999., "Minimum VXD Hits of the track"};

  Gaudi::Property<double> m_chi2th{this, "Chi2Threshold", 25., "Chi2 Threshold for the track"};
  Gaudi::Property<bool> m_beamspotConstraint{this, "BeamspotConstraint", true, "Beamspot Constraint"};
  Gaudi::Property<bool> m_beamspotSmearing{this, "BeamspotSmearing", true, "Beamspot Smearing"};
  Gaudi::Property<double> m_magneticField{this, "MagneticField", 2., "Magnetic field at the interaction point [T]"};
  Gaudi::Property<double> m_beamSizeX{this, "BeamSizeX", 38.2e-3, "Beam size along x [mm]"};
  Gaudi::Property<double> m_beamSizeY{this, "BeamSizeY", 68.e-6, "Beam size along y [mm]"};
  Gaudi::Property<double> m_beamSizeZ{this, "BeamSizeZ", 1.97, "Beam size along z [mm]"};
};

#endif // K4RECO_PRIMARYVERTEXFINDER_H
