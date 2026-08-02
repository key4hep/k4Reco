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
#ifndef K4RECO_JETCLUSTERING_H
#define K4RECO_JETCLUSTERING_H

#include <edm4hep/ReconstructedParticleCollection.h>
#include <edm4hep/VertexCollection.h>

#include <Gaudi/Property.h>

#include <k4FWCore/Transformer.h>

struct JetClustering final : k4FWCore::Transformer<edm4hep::ReconstructedParticleCollection(
                                 const edm4hep::ReconstructedParticleCollection&, const edm4hep::VertexCollection&,
                                 const edm4hep::VertexCollection&)> {
  JetClustering(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;
  edm4hep::ReconstructedParticleCollection
  operator()(const edm4hep::ReconstructedParticleCollection& pfos, const edm4hep::VertexCollection& primaryVertices,
             const edm4hep::VertexCollection& secondaryVertices) const override;

  Gaudi::Property<int> m_nJetsRequested{this, "NJetsRequested", 2, "Requested number of jets (zero enables YCut)"};
  Gaudi::Property<double> m_yCut{this, "YCut", 0., "Clustering Y cut when NJetsRequested is zero"};
  Gaudi::Property<std::string> m_jetAlgorithm{this, "JetAlgorithm", "ValenciaVertex", "LCFIPlus jet algorithm"};
  Gaudi::Property<bool> m_useBeamJets{this, "UseBeamJets", true, "Enable beam jets"};
  Gaudi::Property<double> m_rParameter{this, "RParameter", 1., "Jet-radius parameter"};
  Gaudi::Property<double> m_alphaParameter{this, "AlphaParameter", 1., "Durham alpha parameter"};
  Gaudi::Property<double> m_betaParameter{this, "BetaParameter", 1., "Valencia beta parameter"};
  Gaudi::Property<double> m_gammaParameter{this, "GammaParameter", 1., "Valencia gamma parameter"};
  Gaudi::Property<double> m_yAddedForJetVertexVertex{this, "YAddedForJetVertexVertex", 100.,
                                                     "Vertex-vertex merging penalty"};
  Gaudi::Property<double> m_yAddedForJetVertexLepton{this, "YAddedForJetVertexLepton", 100.,
                                                     "Vertex-lepton merging penalty"};
  Gaudi::Property<double> m_yAddedForJetLeptonLepton{this, "YAddedForJetLeptonLepton", 100.,
                                                     "Lepton-lepton merging penalty"};
  Gaudi::Property<bool> m_useMuonID{this, "UseMuonID", true, "Enable LCFIPlus secondary-muon identification"};
  Gaudi::Property<bool> m_muonIDExternal{this, "MuonIDExternal", false, "Use externally supplied muon identification"};
  Gaudi::Property<double> m_muonIDMinimumEnergy{this, "MuonIDMinimumEnergy", 0., "Minimum muon energy"};
  Gaudi::Property<double> m_muonIDMinimumD0Significance{this, "MuonIDMinimumD0Significance", 5.,
                                                        "Minimum muon d0 significance"};
  Gaudi::Property<double> m_muonIDMinimumZ0Significance{this, "MuonIDMinimumZ0Significance", 5.,
                                                        "Minimum muon z0 significance"};
  Gaudi::Property<double> m_muonIDMaximum3DImpactParameter{this, "MuonIDMaximum3DImpactParameter", 5.,
                                                           "Maximum muon impact parameter"};
  Gaudi::Property<double> m_muonIDMinimumProbability{this, "MuonIDMinimumProbability", .5,
                                                     "Minimum external muon probability"};

  // JetClustering is the first algorithm run by CLDConfig's second
  // LcfiplusProcessor. That processor resets the process-wide LCFIPlus
  // globals during initialization. Keep the legacy defaults here so the
  // native chain reproduces the same subsequent vertex-fit constraint.
  Gaudi::Property<double> m_magneticField{this, "MagneticField", 2., "Magnetic field at the interaction point [T]"};
  Gaudi::Property<double> m_beamSizeX{this, "BeamSizeX", 639.e-6, "Beam size along x [mm]"};
  Gaudi::Property<double> m_beamSizeY{this, "BeamSizeY", 5.7e-6, "Beam size along y [mm]"};
  Gaudi::Property<double> m_beamSizeZ{this, "BeamSizeZ", 9.13e-2, "Beam size along z [mm]"};
};

#endif
