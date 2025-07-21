/*
 * Copyright (c) 2020-2024 Key4hep-Project.
 *
 * This file is part of Key4hep.
 * See https://key4hep-doc/ for further info.
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
#include "LumiCalClusterer.h"
#include "LumiCalHit.h"

#include <Gaudi/Algorithm.h>

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
#include <edm4hep/SimCalorimeterHitCollection.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

const int SHIFT_I_32Fcal = 0;  // I = 10 bits  ( ring )
const int SHIFT_J_32Fcal = 10; // J = 10 bits  ( sector)
const int SHIFT_K_32Fcal = 20; // K = 10 bits  ( layer )
const int SHIFT_S_32Fcal = 30; // S =  2 bits  ( side/arm )

const unsigned int MASK_I_32Fcal = 0x000003FF; // 10 bits
const unsigned int MASK_J_32Fcal = 0x000FFC00; // 10 bits
const unsigned int MASK_K_32Fcal = 0x3FF00000; // 10 bits
const unsigned int MASK_S_32Fcal = 0xC0000000; // 2 bits

LumiCalClustererClass::LumiCalClustererClass(const Gaudi::Algorithm* alg) : m_alg(alg) {}

void LumiCalClustererClass::createDecoder(const std::string& decoderString) {
  m_mydecoder = std::make_unique<dd4hep::DDSegmentation::BitFieldCoder>(decoderString);
}

void LumiCalClustererClass::init(const std::map<std::string, std::variant<int, float, std::string>>& _lcalRecoPars) {
  setConstants(_lcalRecoPars);

  if (std::get<std::string>(_lcalRecoPars.at("WeightingMethod")) == "LogMethod") {
    m_methodCM = LogMethod;
  } else if (std::get<std::string>(_lcalRecoPars.at("WeightingMethod")) == "EnergyMethod") {
    m_methodCM = EnergyMethod;
  } else {
    throw std::runtime_error("Unknown weighting method: " + std::get<std::string>(_lcalRecoPars.at("WeightingMethod")));
  }

  m_hitMinEnergy = m_minHitEnergy;
  m_nNearNeighbor = m_numOfNearNeighbor;

  m_logWeightConst = m_logWeightConstant;

  m_minSeparationDistance = m_minSeparationDist;

  m_maxLayerToAnalyse = m_numCellsZ;
  m_cellRMax = m_numCellsR;
  m_cellPhiMax = m_numCellsPhi;
}

/* --------------------------------------------------------------------------
   (1):	return a cellId for a given Z (layer), R (cylinder) and Phi (sector)
   (2):	return Z (layer), R (cylinder) and Phi (sector) for a given cellId
   -------------------------------------------------------------------------- */

std::array<double, 3> LumiCalClustererClass::rotateToLumiCal(const edm4hep::Vector3f& glob) const {
  const int armNow = (glob[2] < 0) ? -1 : 1;
  std::array<double, 3> loc;
  loc[0] = +m_armCosAngle.at(armNow) * glob[0] - m_armSinAngle.at(armNow) * glob[2];
  loc[1] = glob[1];
  loc[2] = +m_armSinAngle.at(armNow) * glob[0] + m_armCosAngle.at(armNow) * glob[2];
  return loc;
}

edm4hep::Vector3f LumiCalClustererClass::rotateToGlobal(const edm4hep::Vector3f& loc) const {
  const int armNow = (loc[2] < 0) ? -1 : 1;
  edm4hep::Vector3f glob;
  glob.x = +m_armCosAngle.at(armNow) * loc[0] + m_armSinAngle.at(armNow) * loc[2];
  glob.y = loc[1];
  glob.z = -m_armSinAngle.at(armNow) * loc[0] + m_armCosAngle.at(armNow) * loc[2];
  return glob;
}

int LumiCalClustererClass::cellIdZPR(const int cellZ, const int cellPhi, const int cellR, const int arm) {

  int cellId = 0;
  int side = (arm < 0) ? 0 : arm;
  cellId = (((side << SHIFT_S_32Fcal) & MASK_S_32Fcal) | ((cellR << SHIFT_I_32Fcal) & MASK_I_32Fcal) |
            ((cellPhi << SHIFT_J_32Fcal) & MASK_J_32Fcal) | ((cellZ << SHIFT_K_32Fcal) & MASK_K_32Fcal));
  return cellId;
}

void LumiCalClustererClass::cellIdZPR(const int cellID, int& cellZ, int& cellPhi, int& cellR, int& arm) {

  // compute Z,Phi,R indices according to the cellId

  cellR = ((((unsigned int)cellID) & MASK_I_32Fcal) >> SHIFT_I_32Fcal);
  cellPhi = ((((unsigned int)cellID) & MASK_J_32Fcal) >> SHIFT_J_32Fcal);
  cellZ = ((((unsigned int)cellID) & MASK_K_32Fcal) >> SHIFT_K_32Fcal);
  arm = ((((unsigned int)cellID) & MASK_S_32Fcal) >> SHIFT_S_32Fcal);
}

int LumiCalClustererClass::cellIdZPR(const int cellId, const Coordinate_t ZPR) {

  int cellZ, cellPhi, cellR, arm;
  cellIdZPR(cellId, cellZ, cellPhi, cellR, arm);
  arm = (arm == 0) ? -1 : 1;
  if (ZPR == COZ)
    return cellZ;
  else if (ZPR == COR)
    return cellR;
  else if (ZPR == COP)
    return cellPhi;
  else if (ZPR == COA)
    return arm;

  return 0;
}

void LumiCalClustererClass::setConstants(
    const std::map<std::string, std::variant<int, float, std::string>>& _lcalRecoPars) {

  setGeometryDD4hep();

  //------------------------------------------------------------------------
  // Processor Parameters
  // Clustering/Reco parameters
  //(BP) layer relative phi offset - must go sometimes to GEAR params
  const std::string parname = "ZLayerPhiOffset";
  double val = std::get<float>(_lcalRecoPars.at(parname));

  // check units just in case ( convert to rad as needed )
  val = (val <= m_phiCellLength) ? val : val * M_PI / 180.;
  m_zLayerPhiOffset = val;

  m_signalToGeV = 1. / std::get<float>(_lcalRecoPars.at("EnergyCalibConst"));

  // logarithmic constant for position reconstruction
  m_logWeightConstant = std::get<float>(_lcalRecoPars.at("LogWeigthConstant"));

  m_minHitEnergy = toGev(std::get<float>(_lcalRecoPars.at("MinHitEnergy")));
  m_middleEnergyHitBoundFrac = std::get<float>(_lcalRecoPars.at("MiddleEnergyHitBoundFrac"));
  m_elementsPercentInShowerPeakLayer = std::get<float>(_lcalRecoPars.at("ElementsPercentInShowerPeakLayer"));
  m_clusterMinNumHits = std::get<int>(_lcalRecoPars.at("ClusterMinNumHits"));

  // Moliere radius of LumiCal [mm]
  m_moliereRadius = std::get<float>(_lcalRecoPars.at("MoliereRadius"));

  // Geometrical fiducial volume of LumiCal - minimal and maximal polar angles [rad]
  // (BP) Note, this in local LumiCal Reference System ( crossing angle not accounted )
  // quite large - conservative, further reco-particles selection can be done later if desired
  m_thetaMin = (m_rMin + m_moliereRadius) / m_zEnd;
  m_thetaMax = (m_rMax - m_moliereRadius) / m_zStart;

  // minimal separation distance between any pair of clusters [mm]
  m_minSeparationDist = m_moliereRadius;

  // minimal energy of a single cluster
  m_minClusterEngyGeV = std::get<float>(_lcalRecoPars.at("MinClusterEngy"));

  // hits positions weighting method
  m_elementsPercentInShowerPeakLayer = std::get<float>(_lcalRecoPars.at("ElementsPercentInShowerPeakLayer"));
  m_numOfNearNeighbor = std::get<int>(_lcalRecoPars.at("NumOfNearNeighbor"));

  m_armCosAngle[-1] = cos(-m_beamCrossingAngle / 2.);
  m_armCosAngle[1] = cos(m_beamCrossingAngle / 2.);
  m_armSinAngle[-1] = sin(-m_beamCrossingAngle / 2.);
  m_armSinAngle[1] = sin(m_beamCrossingAngle / 2.);
}

double LumiCalClustererClass::toGev(const double value) const { return value / getCalibrationFactor(); }

void LumiCalClustererClass::printAllParameters() const {
  // info() << "------------------------------------------------------------------" << std::endl;
  // info() << "********* LumiCalReco Parameters set in LumiCalClustererClass ********" << std::endl;
  // info() << "---------------------------------------------------------------" << std::endl;
}

bool LumiCalClustererClass::setGeometryDD4hep() {

  const dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();

  if (theDetector.detectors().count("LumiCal") == 0)
    return false;

  const dd4hep::DetElement lumical(theDetector.detector("LumiCal"));
  const dd4hep::Segmentation readout(theDetector.readout("LumiCalCollection").segmentation());

  // info() << "Segmentation Type" << readout.type() << std::endl;
  // info() << "FieldDef: " << readout.segmentation()->fieldDescription() << std::endl;

  const dd4hep::rec::LayeredCalorimeterData* theExtension = lumical.extension<dd4hep::rec::LayeredCalorimeterData>();
  const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& layers = theExtension->layers;

  m_rMin = theExtension->extent[0] / dd4hep::mm;
  m_rMax = theExtension->extent[1] / dd4hep::mm;

  // starting/end position [mm]
  m_zStart = theExtension->extent[2] / dd4hep::mm;
  m_zEnd = theExtension->extent[3] / dd4hep::mm;

  // cell division numbers
  typedef dd4hep::DDSegmentation::TypedSegmentationParameter<double> ParDou;
  const ParDou* rPar = dynamic_cast<ParDou*>(readout.segmentation()->parameter("grid_size_r"));
  const ParDou* rOff = dynamic_cast<ParDou*>(readout.segmentation()->parameter("offset_r"));
  const ParDou* pPar = dynamic_cast<ParDou*>(readout.segmentation()->parameter("grid_size_phi"));
  // const ParDou* pOff = dynamic_cast<ParDou*>(readout.segmentation()->parameter("offset_phi"));

  if (!rPar || !pPar) {
    throw std::runtime_error("Could not obtain parameters from segmentation");
  }

  m_rCellLength = rPar->typedValue() / dd4hep::mm;
  m_rCellOffset = 0.5 * m_rCellLength + m_rMin - rOff->typedValue() / dd4hep::mm;
  m_phiCellLength = pPar->typedValue() / dd4hep::radian;
  m_numCellsR = (int)((m_rMax - m_rMin) / m_rCellLength);
  m_numCellsPhi = (int)(2.0 * M_PI / m_phiCellLength + 0.5);
  m_numCellsZ = layers.size();

  // beam crossing angle ( convert to rad )
  const dd4hep::DetElement::Children children = lumical.children();

  if (children.empty()) {
    throw std::runtime_error("Cannot obtain crossing angle from this LumiCal, update lcgeo?");
  }

  for (auto it = children.begin(); it != children.end(); ++it) {
    dd4hep::Position loc(0.0, 0.0, 0.0);
    dd4hep::Position glob(0.0, 0.0, 0.0);
    it->second.nominal().localToWorld(loc, glob);
    m_beamCrossingAngle = 2.0 * fabs(atan(glob.x() / glob.z()) / dd4hep::rad);
    if (glob.z() > 0.0) {
    } else {
      const auto& backwardCalo = &it->second.nominal().worldTransformation();

      // get phi rotation from global to local transformation
      TGeoHMatrix* tempMat = (TGeoHMatrix*)backwardCalo->Clone();
      double nulltr[] = {0.0, 0.0, 0.0};
      // undo backward and crossing angle rotation
      tempMat->SetTranslation(nulltr);
      // root matrices need degrees as argument
      tempMat->RotateY(m_beamCrossingAngle / 2.0 * 180 / M_PI);
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

  // successfully created geometry from DD4hep
  return true;
}

double LumiCalClustererClass::posWeight(const double cellEngy, const double totEngy, const WeightingMethod_t method,
                                        const double logWeightConstNow) {
  if (method == EnergyMethod)
    return cellEngy;
  // ???????? DECIDE/FIX - improve the log weight constants ????????
  if (method == LogMethod) {
    return std::max(0.0, log(cellEngy / totEngy) + logWeightConstNow);
  }
  return -1;
}

std::tuple<std::optional<edm4hep::MutableCluster>, std::optional<edm4hep::MutableReconstructedParticle>>
LumiCalClustererClass::getLCIOObjects(const LCCluster& thisClusterInfo, const double minClusterEnergy,
                                      const bool cutOnFiducialVolume,
                                      const edm4hep::CalorimeterHitCollection& calohits) const {
  double ThetaMid = (m_thetaMin + m_thetaMax) / 2.;
  double ThetaTol = (m_thetaMax - m_thetaMin) / 2.;

  const double clusterEnergy = thisClusterInfo.getE();
  if (clusterEnergy < minClusterEnergy)
    return std::make_tuple(std::nullopt, std::nullopt);

  if (cutOnFiducialVolume && std::abs(thisClusterInfo.getTheta() - ThetaMid) > ThetaTol) {
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

/* ============================================================================
   main actions in each event:
   ========================================================================= */
std::optional<edm4hep::CalorimeterHitCollection>
LumiCalClustererClass::processEvent(const edm4hep::SimCalorimeterHitCollection& col) {
  // increment / initialize global variables
  m_totEngyArm[-1] = m_totEngyArm[1] = 0.;
  m_numHitsInArm[-1] = m_numHitsInArm[1] = 0;

  MapIntMapIntVCalHit calHits;
  MapIntMapIntLCCluster superClusterCM;

  MapIntMapIntCalHit calHitsCellIdGlobal;

  m_superClusterIdToCellId.clear();
  m_superClusterIdToCellEngy.clear();
  m_superClusterIdClusterInfo.clear();

  auto [retval, calohits] = getCalHits(col, calHits);
  if (!retval)
    return {};

  for (const int armNow : {-1, 1}) {
    m_alg->debug() << endmsg << "ARM = " << armNow << " : " << endmsg << endmsg;
    /* --------------------------------------------------------------------------
       Construct clusters for each arm
       -------------------------------------------------------------------------- */
    m_alg->debug() << "\tRun LumiCalClustererClass::buildClusters()" << endmsg;
    m_alg->debug() << "\tEnergy deposit: " << m_totEngyArm[armNow] << "\tNumber of hits: " << m_numHitsInArm[armNow]
                   << endmsg;
    buildClusters(calHits[armNow], calHitsCellIdGlobal[armNow], m_superClusterIdToCellId[armNow],
                  m_superClusterIdToCellEngy[armNow], superClusterCM[armNow], armNow);

    /* --------------------------------------------------------------------------
       Merge superClusters according the minDistance and minEngy rules
       -------------------------------------------------------------------------- */
    m_alg->debug() << "\tRun LumiCalClustererClass::clusterMerger()" << endmsg;
    m_alg->debug() << printClusters(armNow, superClusterCM);
    clusterMerger(m_superClusterIdToCellEngy[armNow], m_superClusterIdToCellId[armNow], superClusterCM[armNow],
                  calHitsCellIdGlobal[armNow]);

    /* --------------------------------------------------------------------------
       Perform fiducial volume cuts
       -------------------------------------------------------------------------- */
    m_alg->debug() << "\tRun LumiCalClustererClass::fiducialVolumeCuts()" << endmsg;
    m_alg->debug() << printClusters(armNow, superClusterCM);
    fiducialVolumeCuts(m_superClusterIdToCellId[armNow], m_superClusterIdToCellEngy[armNow], superClusterCM[armNow]);

    /* --------------------------------------------------------------------------
       Perform energy correction for inter-mixed superClusters
       -------------------------------------------------------------------------- */
#if _CLUSTER_MIXING_ENERGY_CORRECTIONS == 1
    if (superClusterCM[armNow].size() == 2) {
      m_alg->debug() << "Run LumiCalClustererClass::energyCorrections()" << endmsg;
      m_alg->debug() << printClusters(armNow, superClusterCM);
      energyCorrections(m_superClusterIdToCellId[armNow], m_superClusterIdToCellEngy[armNow], superClusterCM[armNow],
                        calHitsCellIdGlobal[armNow]);

      m_alg->debug() << "After LumiCalClustererClass::energyCorrections()" << endmsg;
      m_alg->debug() << printClusters(armNow, superClusterCM);
    }
#endif
  }

  /* --------------------------------------------------------------------------
     verbosity
     -------------------------------------------------------------------------- */
  m_alg->debug() << "Final clusters:" << endmsg;
  for (int armNow = -1; armNow < 2; armNow += 2) {
    for (MapIntLCCluster::iterator superClusterCMIterator = superClusterCM[armNow].begin();
         superClusterCMIterator != superClusterCM[armNow].end(); ++superClusterCMIterator) {
      const int superClusterId = superClusterCMIterator->first;

      m_alg->debug() << "  Arm:" << std::setw(4) << armNow << "  Id:" << std::setw(4) << superClusterId
                     << superClusterCMIterator->second << endmsg;
      m_superClusterIdClusterInfo[armNow][superClusterId] = superClusterCMIterator->second;
    }
  }

  return std::move(calohits);
}
