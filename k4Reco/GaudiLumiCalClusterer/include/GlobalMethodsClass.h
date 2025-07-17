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
#ifndef K4RECO_GLOBALMETHODSCLASS_H
#define K4RECO_GLOBALMETHODSCLASS_H

#include <edm4hep/Cluster.h>
#include <edm4hep/ReconstructedParticle.h>
#include <edm4hep/Vector3f.h>

#include <array>
#include <map>
#include <optional>
#include <string>
#include <tuple>
#include <variant>

class LCCluster;

class GlobalMethodsClass {
public:
  GlobalMethodsClass();
  GlobalMethodsClass(const GlobalMethodsClass& rhs) = delete;
  GlobalMethodsClass& operator=(const GlobalMethodsClass& rhs) = default;
  GlobalMethodsClass(GlobalMethodsClass&& rhs) = default;
  GlobalMethodsClass& operator=(GlobalMethodsClass&& rhs) = default;
  ~GlobalMethodsClass() = default;

  enum WeightingMethod_t { LogMethod = -1, EnergyMethod = 1 };

  enum Parameter_t {
    ZStart,
    ZEnd,
    RMin,
    RMax,
    NumCellsR,
    NumCellsPhi,
    NumCellsZ,
    RCellLength,
    RCellOffset,
    PhiCellLength,
    PhiCellOffset,
    ZLayerThickness,
    ZLayerPhiOffset,
    ZLayerZOffset,
    ThetaMin,
    ThetaMax,
    LogWeightConstant,
    MoliereRadius,
    MinSeparationDist,
    ElementsPercentInShowerPeakLayer,
    NumOfNearNeighbor,
    ClusterMinNumHits,
    MinHitEnergy,
    MinClusterEngyGeV,
    MiddleEnergyHitBoundFrac,
    WeightingMethod,
    Signal_to_GeV,
    BeamCrossingAngle,
    LumiInColName,
    BetaGamma,
    Gamma
  };

  std::map<Parameter_t, int> m_globalParamI;
  std::map<Parameter_t, double> m_globalParamD;
  std::map<Parameter_t, std::string> m_globalParamS;

  enum Coordinate_t { COTheta, COPhi, COZ, COR, COP, COA };

  static std::string getParameterName(Parameter_t par);

  void setConstants(const std::map<std::string, std::variant<int, float, std::string>>& _lcalRecoPars);
  WeightingMethod_t getMethod(const std::string& methodName) const;

  double toGev(const double valNow) const;

  static void cellIdZPR(const int cellId, int& cellZ, int& cellPhi, int& cellR, int& arm);
  static int cellIdZPR(const int cellZ, const int cellPhi, const int cellR, const int arm);
  static int cellIdZPR(const int cellId, const Coordinate_t ZPR);

  void initializeAdditionalParameters();

  static double posWeight(const double cellEngy, const double totEngy,
                          const GlobalMethodsClass::WeightingMethod_t method, const double logWeightConstNow);

  inline double getCalibrationFactor() const { return m_globalParamD.at(Signal_to_GeV); }

  std::array<double, 3> rotateToLumiCal(const edm4hep::Vector3f& glob) const;

  std::tuple<std::optional<edm4hep::MutableCluster>, std::optional<edm4hep::MutableReconstructedParticle>>
  getLCIOObjects(const LCCluster& thisClusterInfo, const double minClusterEnergy, const bool cutOnFiducialVolume,
                 const edm4hep::CalorimeterHitCollection& calohits) const;

private:
  double m_backwardRotationPhi;

  edm4hep::Vector3f rotateToGlobal(const edm4hep::Vector3f& loc) const;

  void printAllParameters() const;

  bool setGeometryDD4hep();

  // LumiCal rotations angles ( local->Global )
  std::map<int, double> m_armCosAngle{};
  std::map<int, double> m_armSinAngle{};
};

#endif
