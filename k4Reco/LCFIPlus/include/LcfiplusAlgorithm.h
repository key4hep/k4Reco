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
#ifndef K4RECO_LCFIPLUSPROCESSOR_H
#define K4RECO_LCFIPLUSPROCESSOR_H 1

#include <Gaudi/Property.h>

#include <k4FWCore/Transformer.h>

#include <edm4hep/MCParticleCollection.h>
#include <edm4hep/RecoMCParticleLinkCollection.h>
#include <edm4hep/ReconstructedParticleCollection.h>

#include <string>
#include <vector>

/**
 *  Marlin processor for LCFIPlus.
 *
 * @author Tomohiko Tanabe, ICEPP, The University of Tokyo
 * @author Taikan Suehara, ICEPP, The University of Tokyo
 * @version $Id$
 */

struct LcfiplusAlgorithm final
    : k4FWCore::MultiTransformer<std::tuple<edm4hep::ReconstructedParticleCollection>(
          const edm4hep::ReconstructedParticleCollection&, const edm4hep::MCParticleCollection&,
          const edm4hep::RecoMCParticleLinkCollection&)> {

  std::tuple<edm4hep::ReconstructedParticleCollection>
  operator()(const edm4hep::ReconstructedParticleCollection&, const edm4hep::MCParticleCollection&,
             const edm4hep::RecoMCParticleLinkCollection&) const override;

  LcfiplusAlgorithm(const std::string& name, ISvcLocator* svcLoc);
  StatusCode initialize() override;
  StatusCode finalize() override;

private:
  Gaudi::Property<std::vector<std::string>> _algonames{this, "Algorithms", {}, "LCFIPlus algorithms to run"};
  Gaudi::Property<int> _readSubdetectorEnergies{this, "ReadSubdetectorEnergies", 1, "Read subdetector energies (ILD)"};
  Gaudi::Property<int> _trackHitOrdering{this, "TrackHitOrdering", 0,
                                         "Track hit ordering: 0=ILD-LOI (default), 1=ILD-DBD, 2=CLICdet"};
  Gaudi::Property<int> _updateVertexRPDaughters{
      this, "UpdateVertexRPDaughters", 1,
      "Writing back obtained vertices to input RP collections (which must be writable)"};
  Gaudi::Property<int> _ignoreLackOfVertexRP{this, "IgnoreLackOfVertexRP", 0,
                                             "Keep running even if vertex RP collection is not present"};
  Gaudi::Property<int> _printPeriod{this, "PrintEventNumber", 0,
                                    "Event number printing period in std output: 0 = no printing"};
  Gaudi::Property<std::string> _pidAlgoName{this, "PIDAlgorithmName", "LikelihoodPID", "ParticleID Algorithm Name"};
  Gaudi::Property<float> _magneticField{this, "MagneticField", 0.0,
                                        "Manually set magnetic field, overriding the value from DD4hep [T]"};
  Gaudi::Property<float> _beamSizeX{this, "BeamSizeX", 639e-6, "Bunch size in the X direction [mm]"};
  Gaudi::Property<float> _beamSizeY{this, "BeamSizeY", 5.7e-6, "Bunch size in the Y direction [mm]"};
  Gaudi::Property<float> _beamSizeZ{this, "BeamSizeZ", 9.13e-2, "Bunch size in the Z direction [mm]"};

  // lciostorer singleton
  // static lcfiplus::LCIOStorer* _lcio;
  bool _lcioowner;

  int _useMcp;

  /** Input collection name.
   */
  std::string _pfoCollectionName;
  std::string _mcpCollectionName;
  std::string _mcpfoRelationName;

  // std::vector<lcfiplus::Algorithm*> _algos;
  // lcfiplus::Parameters* _param;

  int _nRun;
  int _nEvt;

  // collections to register
  std::vector<std::string> _vertexColNamesToWrite;
  std::vector<std::string> _jetColNamesToWrite;
  std::vector<int> _vertexColNamesToWriteFlags;
  std::vector<int> _jetColNamesToWriteFlags;

  bool _inInit;
};

#endif
