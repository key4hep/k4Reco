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
#ifndef K4RECO_GAUDILUMICALCLUSTERER_H
#define K4RECO_GAUDILUMICALCLUSTERER_H

#include "LumiCalClusterer.h"

#include <edm4hep/CalorimeterHitCollection.h>
#include <edm4hep/ClusterCollection.h>
#include <edm4hep/ReconstructedParticleCollection.h>
#include <edm4hep/SimCalorimeterHitCollection.h>

#include <Gaudi/Property.h>

#include <k4FWCore/Transformer.h>
#include <k4Interface/IGeoSvc.h>

#include <map>
#include <string>
#include <tuple>
#include <vector>

class LCCluster;
class ClusterClass;

typedef std::map<int, std::vector<int>> MapIntVInt;

struct GaudiLumiCalClusterer
    : k4FWCore::MultiTransformer<std::tuple<edm4hep::CalorimeterHitCollection, edm4hep::ClusterCollection,
                                            edm4hep::ReconstructedParticleCollection>(
          const edm4hep::SimCalorimeterHitCollection&)> {

public:
  GaudiLumiCalClusterer(const std::string& name, ISvcLocator* svcLoc);

  // LumiCalClusterer keeps an internal state, so it cannot be reentrant
  bool isReEntrant() const override { return false; }

  StatusCode initialize() override;

  std::tuple<edm4hep::CalorimeterHitCollection, edm4hep::ClusterCollection, edm4hep::ReconstructedParticleCollection>
  operator()(const edm4hep::SimCalorimeterHitCollection&) const override;

  double m_BeamCrossingAngle;

  mutable LumiCalClustererClass m_lumiCalClusterer;

  // void CreateClusters(std::map<int, MapIntPClusterClass>& clusterClassMapP, EVENT::LCEvent* evt);
  // std::tuple<ClusterImpl*, ReconstructedParticleImpl*> getLCIOObjects(LCCluster const& clusterInfo) const;
  // void writeRootInfo(EVENT::LCEvent* evt);

  // void storeMCParticleInfo(LCEvent* evt, int clusterInFlag);

  Gaudi::Property<int> MaxRecordNumber{this, "MaxRecordNumber", 10, "Number of event to work with"};
  Gaudi::Property<std::string> OutDirName{this, "OutDirName", "rootOut", "Name of output directory"};
  Gaudi::Property<std::string> OutRootFileName{
      this, "OutRootFileName", "",
      "Name of output ROOT file ( without suffix) NO DEFAULT. If no name provided root file will not be generated"};
  Gaudi::Property<int> NumEventsTree{this, "NumEventsTree", 500, "Number of events in memory resident ROOT tree."};
  Gaudi::Property<int> MemoryResidentTree{this, "MemoryResidentTree", 0, "Place for ROOT tree memory(1) or disk(0)"};

  Gaudi::Property<double> m_logWeigthConstant{this, "LogWeigthConstant", 6.,
                                              "Sets minimum for logarithmic energy weights"};
  Gaudi::Property<double> m_zLayerPhiOffset{
      this, "ZLayerPhiOffset", 3.75, "Relative offset of LCal z-layers [deg] default is half of the phi sector size"};
  Gaudi::Property<double> m_EnergyCalibConst{
      this, "EnergyCalibConst", 0.0105, "Calibration const E_dep = EnergyCalibConst*E_primary ( default for LCal ILD)"};
  Gaudi::Property<double> m_rMoliere{this, "MoliereRadius", 16.,
                                     "Moliere radius, controls clusters separation distance [mm]"};
  Gaudi::Property<double> m_minClusterEngy{this, "MinClusterEngy", 2.,
                                           "Sets minimum energy deposit for cluster to be accepted [GeV]"};
  Gaudi::Property<std::string> m_WeightingMethod{this, "WeightingMethod", "LogMethod",
                                                 "Hit positions weighting method (LogMthod, EnergyMethod )"};
  Gaudi::Property<double> m_minHitEnergy{this, "MinHitEnergy", 5.e-6, "Hit energy cut [Mev]"};
  Gaudi::Property<int> m_ClusterMinNumHits{this, "ClusterMinNumHits", 15, "Minimal number of hits in cluster"};
  Gaudi::Property<double> m_ElementsPercentInShowerPeakLayer{this, "ElementsPercentInShowerPeakLayer", 0.03,
                                                             "BP: Not sure what it is"};
  Gaudi::Property<double> m_MiddleEnergyHitBoundFrac{this, "MiddleEnergyHitBoundFrac", 0.01,
                                                     "BP: see explanation in LumiCalClusterer.cpp"};
  Gaudi::Property<int> m_NumOfNearNeighbor{this, "NumOfNearNeighbor", 6, "Number of neighbor hits to consider"};
  Gaudi::Property<bool> m_cutOnFiducialVolume{this, "CutOnFiducuialVolume", false,
                                              "Whether to cut clusters outside of the fiducial volume or not"};

  Gaudi::Property<std::string> m_geoSvcName{this, "GeoSvcName", "GeoSvc", "The name of the GeoSvc instance"};
  SmartIF<IGeoSvc> m_geoSvc;
};

DECLARE_COMPONENT(GaudiLumiCalClusterer)

#endif
