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
#ifndef k4RECO_GLOBALMETHODSCLASS_H
#define k4RECO_GLOBALMETHODSCLASS_H 1

#include <edm4hep/Cluster.h>
#include <edm4hep/ReconstructedParticle.h>
#include <edm4hep/Vector3f.h>

#include <map>
#include <string>
#include <tuple>
#include <variant>

class LCCluster;
class TGeoHMatrix;

class GlobalMethodsClass {
public:
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
  //

  enum Coordinate_t { COTheta, COPhi, COZ, COR, COP, COA };

  static std::string GetParameterName(Parameter_t par);

  typedef std::map<Parameter_t, int> ParametersInt;
  typedef std::map<Parameter_t, double> ParametersDouble;
  typedef std::map<Parameter_t, std::string> ParametersString;

  GlobalMethodsClass();
  GlobalMethodsClass(const GlobalMethodsClass& rhs) = delete;
  GlobalMethodsClass& operator=(const GlobalMethodsClass& rhs) = default;
  GlobalMethodsClass(GlobalMethodsClass&& rhs) = default;
  GlobalMethodsClass& operator=(GlobalMethodsClass&& rhs) = default;
  ~GlobalMethodsClass() = default;

  void SetConstants(const std::map<std::string, std::variant<int, float, std::string>>& _lcalRecoPars);
  WeightingMethod_t getMethod(const std::string& methodName) const;

  WeightingMethod_t method = LogMethod;

  double _backwardRotationPhi;
  ParametersInt GlobalParamI;
  ParametersDouble GlobalParamD;
  ParametersString GlobalParamS;

  double toSignal(const double valNow) const;
  double toGev(const double valNow) const;

  void ThetaPhiCell(const int cellId, std::map<GlobalMethodsClass::Coordinate_t, double>& thetaPhiCell) const;

  static void CellIdZPR(const int cellId, int& cellZ, int& cellPhi, int& cellR, int& arm);
  static int CellIdZPR(const int cellZ, const int cellPhi, const int cellR, const int arm);
  static int CellIdZPR(const int cellId, const Coordinate_t ZPR);

  void PrintAllParameters() const;
  void initializeAdditionalParameters();

  static double posWeight(const double cellEngy, const double totEngy,
                          const GlobalMethodsClass::WeightingMethod_t method, const double logWeightConstNow);

  inline double getCalibrationFactor() const { return GlobalParamD.at(Signal_to_GeV); }

  template <class T, class U>
  inline void rotateToGlobal(const T* loc, U* glob) const;
  void rotateToLumiCal(const edm4hep::Vector3f& glob, double* loc) const;

  std::tuple<std::optional<edm4hep::MutableCluster>, std::optional<edm4hep::MutableReconstructedParticle>>
  getLCIOObjects(const LCCluster& thisClusterInfo, const double minClusterEnergy, const bool cutOnFiducialVolume,
                 const edm4hep::CalorimeterHitCollection& calohits) const;

private:
  bool SetGeometryDD4HEP();

  // LumiCal rotations angles ( local->Global )
  std::map<int, double> m_armCosAngle{};
  std::map<int, double> m_armSinAngle{};
};

template <class T, class U>
void GlobalMethodsClass::rotateToGlobal(const T* loc, U* glob) const {
  const int armNow = (loc[2] < 0) ? -1 : 1;
  glob[0] = +m_armCosAngle.at(armNow) * loc[0] + m_armSinAngle.at(armNow) * loc[2];
  glob[1] = loc[1];
  glob[2] = -m_armSinAngle.at(armNow) * loc[0] + m_armCosAngle.at(armNow) * loc[2];
}

#endif
