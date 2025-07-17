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
#include "GlobalMethodsClass.h"
#include "LCCluster.h"
#include "LumiCalHit.h"

#include <DD4hep/Alignments.h>
#include <DD4hep/DD4hepUnits.h>
#include <DD4hep/DetElement.h>
#include <DD4hep/Detector.h>
#include <DD4hep/Objects.h>
#include <DD4hep/Readout.h>
#include <DD4hep/Segmentations.h>
#include <DDRec/DetectorData.h>
#include <DDSegmentation/Segmentation.h>
#include <DDSegmentation/SegmentationParameter.h>

#include <TGeoMatrix.h>

#include <edm4hep/CalorimeterHitCollection.h>
#include <edm4hep/MutableCluster.h>
#include <edm4hep/MutableReconstructedParticle.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

GlobalMethodsClass::GlobalMethodsClass()
    : m_globalParamI(), m_globalParamD(), m_globalParamS(), m_backwardRotationPhi(0.0) {}

/* --------------------------------------------------------------------------
   (1):	return a cellId for a given Z (layer), R (cylinder) and Phi (sector)
   (2):	return Z (layer), R (cylinder) and Phi (sector) for a given cellId
   -------------------------------------------------------------------------- */

#define SHIFT_I_32Fcal 0  // I = 10 bits  ( ring )
#define SHIFT_J_32Fcal 10 // J = 10 bits  ( sector)
#define SHIFT_K_32Fcal 20 // K = 10 bits  ( layer )
#define SHIFT_S_32Fcal 30 // S =  2 bits  ( side/arm )

#define MASK_I_32Fcal (unsigned int)0x000003FF
#define MASK_J_32Fcal (unsigned int)0x000FFC00
#define MASK_K_32Fcal (unsigned int)0x3FF00000
#define MASK_S_32Fcal (unsigned int)0xC0000000

std::array<double, 3> GlobalMethodsClass::rotateToLumiCal(const edm4hep::Vector3f& glob) const {
  const int armNow = (glob[2] < 0) ? -1 : 1;
  std::array<double, 3> loc;
  loc[0] = +m_armCosAngle.at(armNow) * glob[0] - m_armSinAngle.at(armNow) * glob[2];
  loc[1] = glob[1];
  loc[2] = +m_armSinAngle.at(armNow) * glob[0] + m_armCosAngle.at(armNow) * glob[2];
  return loc;
}

edm4hep::Vector3f GlobalMethodsClass::rotateToGlobal(const edm4hep::Vector3f& loc) const {
  const int armNow = (loc[2] < 0) ? -1 : 1;
  edm4hep::Vector3f glob;
  glob.x = +m_armCosAngle.at(armNow) * loc[0] + m_armSinAngle.at(armNow) * loc[2];
  glob.y = loc[1];
  glob.z = -m_armSinAngle.at(armNow) * loc[0] + m_armCosAngle.at(armNow) * loc[2];
  return glob;
}

int GlobalMethodsClass::cellIdZPR(const int cellZ, const int cellPhi, const int cellR, const int arm) {

  int cellId = 0;
  int side = (arm < 0) ? 0 : arm;
  cellId = (((side << SHIFT_S_32Fcal) & MASK_S_32Fcal) | ((cellR << SHIFT_I_32Fcal) & MASK_I_32Fcal) |
            ((cellPhi << SHIFT_J_32Fcal) & MASK_J_32Fcal) | ((cellZ << SHIFT_K_32Fcal) & MASK_K_32Fcal));
  return cellId;
}

void GlobalMethodsClass::cellIdZPR(const int cellID, int& cellZ, int& cellPhi, int& cellR, int& arm) {

  // compute Z,Phi,R indices according to the cellId

  cellR = ((((unsigned int)cellID) & MASK_I_32Fcal) >> SHIFT_I_32Fcal);
  cellPhi = ((((unsigned int)cellID) & MASK_J_32Fcal) >> SHIFT_J_32Fcal);
  cellZ = ((((unsigned int)cellID) & MASK_K_32Fcal) >> SHIFT_K_32Fcal);
  arm = ((((unsigned int)cellID) & MASK_S_32Fcal) >> SHIFT_S_32Fcal);
}

int GlobalMethodsClass::cellIdZPR(const int cellId, const GlobalMethodsClass::Coordinate_t ZPR) {

  int cellZ, cellPhi, cellR, arm;
  cellIdZPR(cellId, cellZ, cellPhi, cellR, arm);
  arm = (arm == 0) ? -1 : 1;
  if (ZPR == GlobalMethodsClass::COZ)
    return cellZ;
  else if (ZPR == GlobalMethodsClass::COR)
    return cellR;
  else if (ZPR == GlobalMethodsClass::COP)
    return cellPhi;
  else if (ZPR == GlobalMethodsClass::COA)
    return arm;

  return 0;
}

void GlobalMethodsClass::setConstants(
    const std::map<std::string, std::variant<int, float, std::string>>& _lcalRecoPars) {

  setGeometryDD4hep();

  //------------------------------------------------------------------------
  // Processor Parameters
  // Clustering/Reco parameters
  //(BP) layer relative phi offset - must go sometimes to GEAR params
  const std::string parname = "ZLayerPhiOffset";
  double val = std::get<float>(_lcalRecoPars.at(parname));

  // check units just in case ( convert to rad as needed )
  val = (val <= m_globalParamD[PhiCellLength]) ? val : val * M_PI / 180.;
  m_globalParamD[ZLayerPhiOffset] = val;

  m_globalParamD[Signal_to_GeV] = 1. / std::get<float>(_lcalRecoPars.at("EnergyCalibConst"));

  // logarithmic constant for position reconstruction
  m_globalParamD[LogWeightConstant] = std::get<float>(_lcalRecoPars.at("LogWeigthConstant"));

  m_globalParamD[MinHitEnergy] = toGev(std::get<float>(_lcalRecoPars.at("MinHitEnergy")));
  m_globalParamD[MiddleEnergyHitBoundFrac] = std::get<float>(_lcalRecoPars.at("MiddleEnergyHitBoundFrac"));
  m_globalParamD[ElementsPercentInShowerPeakLayer] =
      std::get<float>(_lcalRecoPars.at("ElementsPercentInShowerPeakLayer"));
  m_globalParamI[ClusterMinNumHits] = std::get<int>(_lcalRecoPars.at("ClusterMinNumHits"));

  // Moliere radius of LumiCal [mm]
  m_globalParamD[MoliereRadius] = std::get<float>(_lcalRecoPars.at("MoliereRadius"));

  // Geometrical fiducial volume of LumiCal - minimal and maximal polar angles [rad]
  // (BP) Note, this in local LumiCal Reference System ( crossing angle not accounted )
  // quite large - conservative, further reco-particles selection can be done later if desired
  m_globalParamD[ThetaMin] = (m_globalParamD[RMin] + m_globalParamD[MoliereRadius]) / m_globalParamD[ZEnd];
  m_globalParamD[ThetaMax] = (m_globalParamD[RMax] - m_globalParamD[MoliereRadius]) / m_globalParamD[ZStart];

  // minimal separation distance between any pair of clusters [mm]
  m_globalParamD[MinSeparationDist] = m_globalParamD[MoliereRadius];

  // minimal energy of a single cluster
  m_globalParamD[MinClusterEngyGeV] = std::get<float>(_lcalRecoPars.at("MinClusterEngy"));

  // hits positions weighting method
  m_globalParamS[WeightingMethod] = std::get<std::string>(_lcalRecoPars.at("WeightingMethod"));
  m_globalParamD[ElementsPercentInShowerPeakLayer] =
      std::get<float>(_lcalRecoPars.at("ElementsPercentInShowerPeakLayer"));
  m_globalParamI[NumOfNearNeighbor] = std::get<int>(_lcalRecoPars.at("NumOfNearNeighbor"));

  initializeAdditionalParameters();
}

double GlobalMethodsClass::toSignal(const double value) const { return value * getCalibrationFactor(); }

double GlobalMethodsClass::toGev(const double value) const { return value / getCalibrationFactor(); }

void GlobalMethodsClass::thetaPhiCell(const int cellId,
                                      std::map<GlobalMethodsClass::Coordinate_t, double>& thetaPhiCell) const {
  // compute Z,Phi,R coordinates according to the cellId
  // returned Phi is in the range (-M_PI, M_PI )

  int cellIdZ, cellIdPhi, cellIdR, arm;
  cellIdZPR(cellId, cellIdZ, cellIdPhi, cellIdR, arm);

  // theta
  double rCell =
      m_globalParamD.at(RMin) + (cellIdR + 0.5) * m_globalParamD.at(RCellLength) - m_globalParamD.at(RCellOffset);
  double zCell = fabs(m_globalParamD.at(ZStart)) + m_globalParamD.at(ZLayerThickness) * (cellIdZ) +
                 m_globalParamD.at(ZLayerZOffset);
  double thetaCell = atan(rCell / zCell);

  // phi
  //(BP) use phiCell size and account for possible layers relative offset/stagger
  // double phiCell   = 2 * M_PI * (double(cellIdPhi) + .5) / double(m_globalParamI[NumCellsPhi]) + double( cellIdZ % 2
  // )
  // * m_globalParamD[;
  double phiCell = (double(cellIdPhi)) * m_globalParamD.at(PhiCellLength) +
                   double((cellIdZ) % 2) * m_globalParamD.at(ZLayerPhiOffset) + m_globalParamD.at(PhiCellOffset);
  phiCell = (phiCell > M_PI) ? phiCell - 2. * M_PI : phiCell;
  // fill output container
  thetaPhiCell[GlobalMethodsClass::COTheta] = thetaCell;
  thetaPhiCell[GlobalMethodsClass::COPhi] = phiCell;
  thetaPhiCell[GlobalMethodsClass::COR] = rCell;
  thetaPhiCell[GlobalMethodsClass::COZ] = zCell;
  return;
}

std::string GlobalMethodsClass::getParameterName(Parameter_t par) {

  switch (par) {
  case ZStart:
    return "ZStart";
  case ZEnd:
    return "ZEnd";
  case RMin:
    return "RMin";
  case RMax:
    return "RMax";
  case NumCellsR:
    return "NumCellsR";
  case NumCellsPhi:
    return "NumCellsPhi";
  case NumCellsZ:
    return "NumCellsZ";
  case RCellLength:
    return "RCellLength";
  case RCellOffset:
    return "RCellOffset";
  case PhiCellLength:
    return "PhiCellLength";
  case PhiCellOffset:
    return "PhiCellOffset";
  case ZLayerThickness:
    return "ZLayerThickness";
  case ZLayerPhiOffset:
    return "ZLayerPhiOffset";
  case ZLayerZOffset:
    return "ZLayerZOffset";
  case ThetaMin:
    return "ThetaMin";
  case ThetaMax:
    return "ThetaMax";
  case LogWeightConstant:
    return "LogWeightConstant";
  case MoliereRadius:
    return "MoliereRadius";
  case MinSeparationDist:
    return "MinSeparationDist";
  case ElementsPercentInShowerPeakLayer:
    return "ElementsPercentInShowerPeakLayer";
  case NumOfNearNeighbor:
    return "NumOfNearNeighbor";
  case ClusterMinNumHits:
    return "ClusterMinNumHits";
  case MinHitEnergy:
    return "MinHitEnergy";
  case MinClusterEngyGeV:
    return "MinClusterEngyGeV";
  case MiddleEnergyHitBoundFrac:
    return "MiddleEnergyHitBoundFrac";
  case WeightingMethod:
    return "WeightingMethod";
  case Signal_to_GeV:
    return "Signal_to_GeV";
  case BeamCrossingAngle:
    return "BeamCrossingAngle";
  case LumiInColName:
    return "LumiInColName";
  case BetaGamma:
    return "BetaGamma";
  case Gamma:
    return "Gamma";
  default:
    return "Unknown Parameter";
  }
}

void GlobalMethodsClass::printAllParameters() const {
  // info() << "------------------------------------------------------------------" << std::endl;
  // info() << "********* LumiCalReco Parameters set in GlobalMethodClass ********" << std::endl;

  // for (ParametersInt::const_iterator it = m_globalParamI.begin();it != m_globalParamI.end() ;++it) {
  //   info() << " - (int)     " << getParameterName(it->first) << "  =  " << it->second<< std::endl;
  // }

  // for (ParametersDouble::const_iterator it = m_globalParamD.begin();it != m_globalParamD.end() ;++it) {
  //   info() << " - (double)  " << getParameterName(it->first) << "  =  " << it->second<< std::endl;
  // }

  // for (ParametersString::const_iterator it = m_globalParamS.begin();it != m_globalParamS.end() ;++it) {
  //   info() << " - (string)  " << getParameterName(it->first) << "  =  " << it->second<< std::endl;
  // }

  // info() << "---------------------------------------------------------------" << std::endl;
}

bool GlobalMethodsClass::setGeometryDD4hep() {

  dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();

  if (theDetector.detectors().count("LumiCal") == 0)
    return false;

  dd4hep::DetElement lumical(theDetector.detector("LumiCal"));
  dd4hep::Segmentation readout(theDetector.readout("LumiCalCollection").segmentation());

  // info() << "Segmentation Type" << readout.type() << std::endl;
  // info() << "FieldDef: " << readout.segmentation()->fieldDescription() << std::endl;

  const dd4hep::rec::LayeredCalorimeterData* theExtension = lumical.extension<dd4hep::rec::LayeredCalorimeterData>();
  const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& layers = theExtension->layers;

  m_globalParamD[RMin] = theExtension->extent[0] / dd4hep::mm;
  m_globalParamD[RMax] = theExtension->extent[1] / dd4hep::mm;

  // starting/end position [mm]
  m_globalParamD[ZStart] = theExtension->extent[2] / dd4hep::mm;
  m_globalParamD[ZEnd] = theExtension->extent[3] / dd4hep::mm;

  // cell division numbers
  typedef dd4hep::DDSegmentation::TypedSegmentationParameter<double> ParDou;
  ParDou* rPar = dynamic_cast<ParDou*>(readout.segmentation()->parameter("grid_size_r"));
  ParDou* rOff = dynamic_cast<ParDou*>(readout.segmentation()->parameter("offset_r"));
  ParDou* pPar = dynamic_cast<ParDou*>(readout.segmentation()->parameter("grid_size_phi"));
  ParDou* pOff = dynamic_cast<ParDou*>(readout.segmentation()->parameter("offset_phi"));

  if (!rPar || !pPar) {
    throw std::runtime_error("Could not obtain parameters from segmentation");
  }

  m_globalParamD[RCellLength] = rPar->typedValue() / dd4hep::mm;
  m_globalParamD[RCellOffset] =
      0.5 * m_globalParamD[RCellLength] + m_globalParamD[RMin] - rOff->typedValue() / dd4hep::mm;
  m_globalParamD[PhiCellLength] = pPar->typedValue() / dd4hep::radian;
  m_globalParamD[PhiCellOffset] = pOff->typedValue() / dd4hep::radian;
  m_globalParamI[NumCellsR] = (int)((m_globalParamD[RMax] - m_globalParamD[RMin]) / m_globalParamD[RCellLength]);
  m_globalParamI[NumCellsPhi] = (int)(2.0 * M_PI / m_globalParamD[PhiCellLength] + 0.5);
  m_globalParamI[NumCellsZ] = layers.size();

  // beam crossing angle ( convert to rad )
  dd4hep::DetElement::Children children = lumical.children();

  if (children.empty()) {
    throw std::runtime_error("Cannot obtain crossing angle from this LumiCal, update lcgeo?");
  }

  for (auto it = children.begin(); it != children.end(); ++it) {
    dd4hep::Position loc(0.0, 0.0, 0.0);
    dd4hep::Position glob(0.0, 0.0, 0.0);
    it->second.nominal().localToWorld(loc, glob);
    m_globalParamD[BeamCrossingAngle] = 2.0 * fabs(atan(glob.x() / glob.z()) / dd4hep::rad);
    if (glob.z() > 0.0) {
    } else {
      const auto& backwardCalo = &it->second.nominal().worldTransformation();

      // get phi rotation from global to local transformation
      TGeoHMatrix* tempMat = (TGeoHMatrix*)backwardCalo->Clone();
      double nulltr[] = {0.0, 0.0, 0.0};
      // undo backward and crossing angle rotation
      tempMat->SetTranslation(nulltr);
      // root matrices need degrees as argument
      tempMat->RotateY(m_globalParamD[BeamCrossingAngle] / 2.0 * 180 / M_PI);
      tempMat->RotateY(-180.0);
      double local[] = {0.0, 1.0, 0.0};
      double global[] = {0.0, 0.0, 0.0};

      tempMat->LocalToMaster(local, global);

      m_backwardRotationPhi = atan2(local[1], local[0]) - atan2(global[1], global[0]);
      if (m_backwardRotationPhi < M_PI)
        m_backwardRotationPhi += 2 * M_PI;

      delete tempMat;
    }
  }

  // layer thickness
  m_globalParamD[ZLayerThickness] = (layers[0].inner_thickness + layers[0].outer_thickness) / dd4hep::mm;
  m_globalParamD[ZLayerZOffset] = (layers[0].inner_thickness) / dd4hep::mm;

  // successfully created geometry from DD4hep
  return true;
}

GlobalMethodsClass::WeightingMethod_t GlobalMethodsClass::getMethod(const std::string& methodName) const {
  if (methodName == "LogMethod") {
    return LogMethod;
  }
  if (methodName == "EnergyMethod") {
    return EnergyMethod;
  }
  throw std::runtime_error("Unkown weighting method");
}

double GlobalMethodsClass::posWeight(const double cellEngy, const double totEngy,
                                     const GlobalMethodsClass::WeightingMethod_t method,
                                     const double logWeightConstNow) {
  if (method == GlobalMethodsClass::EnergyMethod)
    return cellEngy;
  // ???????? DECIDE/FIX - improve the log weight constants ????????
  if (method == GlobalMethodsClass::LogMethod) {
    return std::max(0.0, log(cellEngy / totEngy) + logWeightConstNow);
  }
  return -1;
}

void GlobalMethodsClass::initializeAdditionalParameters() {
  // Lorentz boost params
  const double beta = tan(m_globalParamD[BeamCrossingAngle] / 2.0);
  m_globalParamD[BetaGamma] = beta;
  m_globalParamD[Gamma] = sqrt(1. + beta * beta);

  m_armCosAngle[-1] = cos(-m_globalParamD[GlobalMethodsClass::BeamCrossingAngle] / 2.);
  m_armCosAngle[1] = cos(m_globalParamD[GlobalMethodsClass::BeamCrossingAngle] / 2.);
  m_armSinAngle[-1] = sin(-m_globalParamD[GlobalMethodsClass::BeamCrossingAngle] / 2.);
  m_armSinAngle[1] = sin(m_globalParamD[GlobalMethodsClass::BeamCrossingAngle] / 2.);
}

std::tuple<std::optional<edm4hep::MutableCluster>, std::optional<edm4hep::MutableReconstructedParticle>>
GlobalMethodsClass::getLCIOObjects(const LCCluster& thisClusterInfo, const double minClusterEnergy,
                                   const bool cutOnFiducialVolume,
                                   const edm4hep::CalorimeterHitCollection& calohits) const {
  double ThetaMid =
      (m_globalParamD.at(GlobalMethodsClass::ThetaMin) + m_globalParamD.at(GlobalMethodsClass::ThetaMax)) / 2.;
  double ThetaTol =
      (m_globalParamD.at(GlobalMethodsClass::ThetaMax) - m_globalParamD.at(GlobalMethodsClass::ThetaMin)) / 2.;

  const double clusterEnergy = thisClusterInfo.getE();
  if (clusterEnergy < minClusterEnergy)
    return std::make_tuple(std::nullopt, std::nullopt);

  if (cutOnFiducialVolume) {
    const double clusterTheta = thisClusterInfo.getTheta();
    if (fabs(clusterTheta - ThetaMid) > ThetaTol)
      return std::make_tuple(std::nullopt, std::nullopt);
  }

  auto cluster = edm4hep::MutableCluster();
  cluster.setEnergy(clusterEnergy);

  auto particle = edm4hep::MutableReconstructedParticle();
  particle.setMass(0.0);
  particle.setCharge(0.0);
  particle.setEnergy(clusterEnergy);
  particle.addToClusters(cluster);

  const edm4hep::Vector3f locPos{float(thisClusterInfo.getX()), float(thisClusterInfo.getY()),
                                 float(thisClusterInfo.getZ())};
  edm4hep::Vector3f gP = rotateToGlobal(locPos);
  cluster.setPosition(gP);

  const float norm = clusterEnergy / std::sqrt(gP[0] * gP[0] + gP[1] * gP[1] + gP[2] * gP[2]);
  const float clusterMomentum[3] = {float(gP[0] * norm), float(gP[1] * norm), float(gP[2] * norm)};
  particle.setMomentum(clusterMomentum);
  for (const auto& lumicalHit : thisClusterInfo.getCaloHits()) {
    for (const auto hitIndex : lumicalHit->getHits()) {
      cluster.addToHits(calohits.at(hitIndex));
    }
  }
  // cluster.subdetectorEnergies().resize(6);
  // cluster.subdetectorEnergies()[3] = clusterEnergy;
  // LCAL_INDEX=3 in DDPFOCreator.hh
  for (const auto val : {0., 0., 0., clusterEnergy, 0., 0.}) {
    cluster.addToSubdetectorEnergies(val);
  }

  return std::make_tuple(cluster, particle);
}
