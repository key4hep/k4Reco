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
#include "LCCluster.h"

#include <DDSegmentation/BitFieldCoder.h>

#include <Gaudi/Algorithm.h>

#include <edm4hep/CalorimeterHitCollection.h>
#include <edm4hep/Cluster.h>
#include <edm4hep/ReconstructedParticle.h>
#include <edm4hep/SimCalorimeterHitCollection.h>
#include <edm4hep/Vector3f.h>

#include <array>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <variant>

enum RETVAL {
  NOK = 0,
  OK = 1,
};

class LumiCalClustererClass {

public:
  enum WeightingMethod_t { LogMethod = -1, EnergyMethod = 1 };

  enum Coordinate_t { COTheta, COPhi, COZ, COR, COP, COA };

  // Constructor
  LumiCalClustererClass(const Gaudi::Algorithm* alg);
  LumiCalClustererClass(LumiCalClustererClass const& rhs) = delete;
  LumiCalClustererClass& operator=(LumiCalClustererClass const& rhs) = delete;
  LumiCalClustererClass(LumiCalClustererClass&& rhs) = delete;
  LumiCalClustererClass& operator=(LumiCalClustererClass&& rhs) = default;
  ~LumiCalClustererClass() = default;

  // initialization routine - Called at the begining of the job.
  void init(const std::map<std::string, std::variant<int, float, std::string>>& _lcalRecoPars);

  /// set the cutOnFiducialVolume flag
  void setCutOnFiducialVolume(bool cutFlag) { m_cutOnFiducialVolume = cutFlag; }

  void createDecoder(const std::string& decoderString);

  // main actions in each event -Called for every event - the working horse.
  std::pair<RETVAL, edm4hep::CalorimeterHitCollection> processEvent(const edm4hep::SimCalorimeterHitCollection& col);

  MapIntMapIntVInt m_superClusterIdToCellId;
  MapIntMapIntVDouble m_superClusterIdToCellEngy;
  MapIntMapIntLCCluster m_superClusterIdClusterInfo;

  // Methods from GlobalMethodsClass
  void setConstants(const std::map<std::string, std::variant<int, float, std::string>>& _lcalRecoPars);
  WeightingMethod_t getMethod(const std::string& methodName) const;
  double toGev(const double valNow) const;
  static void cellIdZPR(const int cellId, int& cellZ, int& cellPhi, int& cellR, int& arm);
  static int cellIdZPR(const int cellZ, const int cellPhi, const int cellR, const int arm);
  static int cellIdZPR(const int cellId, const Coordinate_t ZPR);
  static double posWeight(const double cellEngy, const double totEngy, const WeightingMethod_t method,
                          const double logWeightConstNow);
  inline double getCalibrationFactor() const { return m_signalToGeV; }
  std::array<double, 3> rotateToLumiCal(const edm4hep::Vector3f& glob) const;
  std::tuple<std::optional<edm4hep::MutableCluster>, std::optional<edm4hep::MutableReconstructedParticle>>
  getLCIOObjects(const LCCluster& thisClusterInfo, const double minClusterEnergy, const bool cutOnFiducialVolume,
                 const edm4hep::CalorimeterHitCollection& calohits) const;

  // Direct parameter access
  double m_beamCrossingAngle{0.0};
  double m_zStart{0.0};
  double m_zEnd{0.0};
  double m_rMin{0.0};
  double m_rMax{0.0};
  int m_numCellsR{0};
  int m_numCellsPhi{0};
  int m_numCellsZ{0};
  double m_rCellLength{0.0};
  double m_rCellOffset{0.0};
  double m_phiCellLength{0.0};
  double m_zLayerPhiOffset{0.0};
  double m_thetaMin{0.0};
  double m_thetaMax{0.0};
  double m_logWeightConstant{0.0};
  double m_moliereRadius{0.0};
  double m_minSeparationDist{0.0};
  double m_elementsPercentInShowerPeakLayer{0.03};
  int m_numOfNearNeighbor{6};
  int m_clusterMinNumHits{15};
  double m_minHitEnergy{5 * 1e-6};
  double m_minClusterEngyGeV{0.0};
  double m_middleEnergyHitBoundFrac{0.01};
  std::string m_weightingMethod{"LogMethod"};
  double m_signalToGeV{1.0};

private:
  // Processor Parameters
  double m_hitMinEnergy{5 * 1e-6};

  // global variables
  size_t m_maxLayerToAnalyse{0};
  double m_logWeightConst{0.0};
  int m_nNearNeighbor{6};
  int m_cellRMax{0};
  int m_cellPhiMax{0};
  WeightingMethod_t m_methodCM{LogMethod};
  double m_minSeparationDistance{0.0};

  MapIntDouble m_totEngyArm;
  MapIntInt m_numHitsInArm;
  //  VInt _armsToCluster;

  std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder> m_mydecoder{};

  bool m_cutOnFiducialVolume{false};

  const Gaudi::Algorithm* m_alg;

  // From GlobalMethodsClass
  double m_backwardRotationPhi{0.0};
  std::map<int, double> m_armCosAngle{};
  std::map<int, double> m_armSinAngle{};

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
                            MapIntVector3d const& virtualClusterCM);

  int virtualCMPeakLayersFix(MapIntCalHit const& calHitsCellId, MapIntInt& cellIdToClusterId,
                             MapIntVInt& clusterIdToCellId, MapIntLCCluster& clusterCM,
                             MapIntVector3d virtualClusterCM);

  int buildSuperClusters(MapIntCalHit& calHitsCellIdGlobal, VMapIntCalHit const& calHitsCellId,
                         VMapIntVInt const& clusterIdToCellId, VMapIntLCCluster const& clusterCM,
                         VMapIntVector3d const& virtualClusterCM, MapIntInt& cellIdToSuperClusterId,
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

  double posWeight(const CalHit& calHit, const WeightingMethod_t method) const;

  double posWeightTrueCluster(const CalHit& calHit, const double cellEngy, WeightingMethod_t method) const;

  LCCluster calculateEngyPosCM(const VInt& cellIdV, const MapIntCalHit& calHitsCellId, const WeightingMethod_t method);

  void calculateEngyPosCM_EngyV(const VInt& cellIdV, const VDouble& cellEngyV, MapIntCalHit const& calHitsCellId,
                                MapIntLCCluster& clusterCM, int clusterId, WeightingMethod_t method);

  void updateEngyPosCM(const CalHit& calHit, LCCluster& clusterCM);

  int checkClusterMergeCM(int clusterId1, int clusterId2, MapIntVInt const& clusterIdToCellId,
                          MapIntCalHit const& calHitsCellId, double distanceAroundCM, double percentOfEngyAroungCM,
                          WeightingMethod_t method);

  double getEngyInMoliereFraction(MapIntCalHit const& calHitsCellId, VInt const& clusterIdToCellId,
                                  LCCluster const& clusterCM, double moliereFraction);

  double getEngyInMoliereFraction(MapIntCalHit const& calHitsCellId, VInt const& clusterIdToCellId,
                                  LCCluster const& clusterCM, double moliereFraction, MapIntInt& flag);

  // void dumpClusters( MapIntLCCluster const& clusterCM );

  std::string printClusters(const int armNow, MapIntMapIntLCCluster const& superClusterCM) const;
  std::string printClusters(MapIntLCCluster const& superClusterCM) const;

  // Private methods from GlobalMethodsClass
  edm4hep::Vector3f rotateToGlobal(const edm4hep::Vector3f& loc) const;
  void printAllParameters() const;
  bool setGeometryDD4hep();
};

#endif // K4RECO_LUMICALCLUSTERER_H
