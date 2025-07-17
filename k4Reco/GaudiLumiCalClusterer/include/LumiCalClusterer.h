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
#ifndef K4RECO_LUMICALCLUSTERER_H
#define K4RECO_LUMICALCLUSTERER_H 1

// clustering options
#define _CLUSTER_MIDDLE_RANGE_ENGY_HITS 1
#define _MOLIERE_RADIUS_CORRECTIONS 1
#define _CLUSTER_MIXING_ENERGY_CORRECTIONS 1

// verbosity
#define _GENERAL_CLUSTERER_DEBUG 0
#define _CLUSTER_BUILD_DEBUG 0
#define _VIRTUALCLUSTER_BUILD_DEBUG 0
#define _MOL_RAD_CORRECT_DEBUG 0

#include "Global.h"
#include "GlobalMethodsClass.h"
#include "LCCluster.h"

#include <DDSegmentation/BitFieldCoder.h>

#include <Gaudi/Algorithm.h>

#include <edm4hep/CalorimeterHitCollection.h>
#include <edm4hep/SimCalorimeterHitCollection.h>

#include <memory>
#include <string>

enum RETVAL {
  NOK = 0,
  OK = 1,
};

class LumiCalClustererClass {

public:
  // Constructor
  LumiCalClustererClass(const Gaudi::Algorithm* alg);
  LumiCalClustererClass(LumiCalClustererClass const& rhs) = delete;
  LumiCalClustererClass& operator=(LumiCalClustererClass const& rhs) = delete;
  LumiCalClustererClass(LumiCalClustererClass&& rhs) = delete;
  LumiCalClustererClass& operator=(LumiCalClustererClass&& rhs) = default;
  ~LumiCalClustererClass() = default;

  // initialization routine - Called at the begining of the job.
  void init(GlobalMethodsClass const& gmc);

  /// set the cutOnFiducialVolume flag
  void setCutOnFiducialVolume(bool cutFlag) { m_cutOnFiducialVolume = cutFlag; }

  void createDecoder(const std::string& decoderString);

  // main actions in each event -Called for every event - the working horse.
  std::pair<RETVAL, edm4hep::CalorimeterHitCollection> processEvent(edm4hep::SimCalorimeterHitCollection& col);

  MapIntMapIntVInt m_superClusterIdToCellId;
  MapIntMapIntVDouble m_superClusterIdToCellEngy;
  MapIntMapIntLCCluster m_superClusterIdClusterInfo;

protected:
  // Processor Parameters
  int m_clusterMinNumHits;
  double m_hitMinEnergy;

  // global variables
  int m_numEventsPerTree, m_resetRootTrees;
  size_t m_maxLayerToAnalyse;
  double m_zFirstLayer, m_zLayerThickness, m_zLayerPhiOffset, m_rMin, m_rMax, m_rCellLength, m_phiCellLength;
  double m_beamCrossingAngle;
  double m_elementsPercentInShowerPeakLayer;
  double m_logWeightConst;
  int m_nNearNeighbor;
  int m_cellRMax, m_cellPhiMax;
  double m_middleEnergyHitBoundFrac;
  GlobalMethodsClass::WeightingMethod_t m_methodCM;
  double m_moliereRadius;
  double m_thetaContainmentBounds[2];
  double m_minSeparationDistance, m_minClusterEngyGeV;

  MapIntDouble m_totEngyArm;
  MapIntInt m_numHitsInArm;
  //  VInt _armsToCluster;

  std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder> m_mydecoder{};

  GlobalMethodsClass m_gmc;

  bool m_cutOnFiducialVolume = false;

  const Gaudi::Algorithm* m_alg;

  std::pair<int, edm4hep::CalorimeterHitCollection> getCalHits(const edm4hep::SimCalorimeterHitCollection& col,
                                                               MapIntMapIntVCalHit& calHits);

  edm4hep::CalorimeterHitCollection createCaloHitCollection(const edm4hep::SimCalorimeterHitCollection& input) const;

  int buildClusters(const MapIntVCalHit& calHits, MapIntCalHit& calHitsCellIdGlobal, MapIntVInt& superClusterIdToCellId,
                    MapIntVDouble& superClusterIdToCellEngy, MapIntLCCluster& superClusterCM, const int detectorArm);

  int initialClusterBuild(MapIntCalHit const& calHitsCellId, MapIntInt& cellIdToClusterId,
                          MapIntVInt& clusterIdToCellId, MapIntLCCluster& clusterCM, VInt const& controlVar);

  int initialLowEngyClusterBuild(MapIntCalHit const& calHitsSmallEngyCellId, MapIntCalHit& calHitsCellId,
                                 MapIntInt& cellIdToClusterId, MapIntVInt& clusterIdToCellId,
                                 MapIntLCCluster& clusterCM);

  int virtualCMClusterBuild(MapIntCalHit const& calHitsCellId, MapIntInt& cellIdToClusterId,
                            MapIntVInt& clusterIdToCellId, MapIntLCCluster& clusterCM,
                            MapIntVirtualCluster const& virtualClusterCM);

  int virtualCMPeakLayersFix(MapIntCalHit const& calHitsCellId, MapIntInt& cellIdToClusterId,
                             MapIntVInt& clusterIdToCellId, MapIntLCCluster& clusterCM,
                             MapIntVirtualCluster virtualClusterCM);

  int buildSuperClusters(MapIntCalHit& calHitsCellIdGlobal, VMapIntCalHit const& calHitsCellId,
                         VMapIntVInt const& clusterIdToCellId, VMapIntLCCluster const& clusterCM,
                         VMapIntVirtualCluster const& virtualClusterCM, MapIntInt& cellIdToSuperClusterId,
                         MapIntVInt& superClusterIdToCellId, MapIntLCCluster& superClusterCM);

  int engyInMoliereCorrections(MapIntCalHit const& calHitsCellIdGlobal, MapIntVCalHit const& calHits,
                               VMapIntCalHit const& calHitsCellIdLayer, VMapIntVInt& clusterIdToCellId,
                               VMapIntLCCluster& clusterCM, VMapIntInt& cellIdToClusterId,
                               MapIntInt& cellIdToSuperClusterId, MapIntVInt& superClusterIdToCellId,
                               MapIntLCCluster& superClusterCM, double middleEnergyHitBound, int detectorArm);

  void energyCorrections(MapIntVInt& superClusterIdToCellId, MapIntVDouble& superClusterIdToCellEngy,
                         MapIntLCCluster& superClusterCM, MapIntCalHit const& calHitsCellIdGlobal);

  void clusterMerger(MapIntVDouble& clusterIdToCellEngy, MapIntVInt& clusterIdToCellId, MapIntLCCluster& clusterCM,
                     MapIntCalHit& calHitsCellIdGlobal);

  void fiducialVolumeCuts(MapIntVInt& superClusterIdToCellId, MapIntVDouble& superClusterIdToCellEngy,
                          MapIntLCCluster& superClusterCM);

  int getNeighborId(int cellId, const int neighborIndex);

  double posWeight(const CalHit& calHit, const GlobalMethodsClass::WeightingMethod_t method) const;

  double posWeightTrueCluster(const CalHit& calHit, const double cellEngy,
                              GlobalMethodsClass::WeightingMethod_t method) const;

  LCCluster calculateEngyPosCM(const VInt& cellIdV, const MapIntCalHit& calHitsCellId,
                               const GlobalMethodsClass::WeightingMethod_t method);

  void calculateEngyPosCM_EngyV(const VInt& cellIdV, const VDouble& cellEngyV, MapIntCalHit const& calHitsCellId,
                                MapIntLCCluster& clusterCM, int clusterId,
                                GlobalMethodsClass::WeightingMethod_t method);

  void updateEngyPosCM(const CalHit& calHit, LCCluster& clusterCM);

  int checkClusterMergeCM(int clusterId1, int clusterId2, MapIntVInt const& clusterIdToCellId,
                          MapIntCalHit const& calHitsCellId, double distanceAroundCM, double percentOfEngyAroungCM,
                          GlobalMethodsClass::WeightingMethod_t method);

  double getEngyInMoliereFraction(MapIntCalHit const& calHitsCellId, VInt const& clusterIdToCellId,
                                  LCCluster const& clusterCM, double moliereFraction);

  double getEngyInMoliereFraction(MapIntCalHit const& calHitsCellId, VInt const& clusterIdToCellId,
                                  LCCluster const& clusterCM, double moliereFraction, MapIntInt& flag);

  // void dumpClusters( MapIntLCCluster const& clusterCM );

  std::string printClusters(const int armNow, MapIntMapIntLCCluster const& superClusterCM) const;
  std::string printClusters(MapIntLCCluster const& superClusterCM) const;
};

#endif // K4RECO_LUMICALCLUSTERER_H
