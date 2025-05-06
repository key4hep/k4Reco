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
#include "LumiCalClusterer.h"

#include <Gaudi/Algorithm.h>

#include <edm4hep/SimCalorimeterHitCollection.h>

#include <iomanip>

LumiCalClustererClass::LumiCalClustererClass(const Gaudi::Algorithm* alg)
    : m_superClusterIdToCellId(), m_superClusterIdToCellEngy(), m_superClusterIdClusterInfo(), m_clusterMinNumHits(15),
      m_hitMinEnergy(5 * 1e-6),
      // global variables
      m_numEventsPerTree(0), m_resetRootTrees(0), m_maxLayerToAnalyse(0), m_zFirstLayer(0), m_zLayerThickness(0.0),
      m_zLayerPhiOffset(0.0), m_rMin(0.0), m_rMax(0.0), m_rCellLength(0.0), m_phiCellLength(0.0),
      m_beamCrossingAngle(0.), m_elementsPercentInShowerPeakLayer(0.03), m_logWeightConst(0.0), m_nNearNeighbor(6),
      m_cellRMax(0), m_cellPhiMax(0), m_middleEnergyHitBoundFrac(0.01), m_methodCM(GlobalMethodsClass::LogMethod),
      m_moliereRadius(), m_thetaContainmentBounds(), m_minSeparationDistance(), m_minClusterEngyGeV(), m_totEngyArm(),
      m_numHitsInArm(), m_gmc(), m_alg(alg) {}

void LumiCalClustererClass::createDecoder(const std::string& decoderString) {
  m_mydecoder = std::make_unique<dd4hep::DDSegmentation::BitFieldCoder>(decoderString);
}

/* ============================================================================
   initial action before first event analysis starts:
   Called at the begining of the job before anything is read.
   ========================================================================= */
void LumiCalClustererClass::init(GlobalMethodsClass const& gmc) {

  m_gmc = gmc;
  /* --------------------------------------------------------------------------
     constants specific to this class
  _armsToCluster.clear();
  _armsToCluster.push_back(-1);
  _armsToCluster.push_back(1);
     -------------------------------------------------------------------------- */
  m_methodCM = gmc.getMethod(gmc.GlobalParamS.at(GlobalMethodsClass::WeightingMethod)); // GlobalMethodsClass::LogMethod
  m_clusterMinNumHits = gmc.GlobalParamI.at(GlobalMethodsClass::ClusterMinNumHits);     // = 15
  m_hitMinEnergy = gmc.GlobalParamD.at(GlobalMethodsClass::MinHitEnergy);               // = 5e-6
  m_zLayerThickness = gmc.GlobalParamD.at(GlobalMethodsClass::ZLayerThickness);         // = 4.5
  m_zLayerPhiOffset = gmc.GlobalParamD.at(GlobalMethodsClass::ZLayerPhiOffset);         // = 3.75 [deg]
  m_elementsPercentInShowerPeakLayer =
      gmc.GlobalParamD.at(GlobalMethodsClass::ElementsPercentInShowerPeakLayer); // = 0.03  //APS 0.04;
  m_nNearNeighbor =
      gmc.GlobalParamI.at(GlobalMethodsClass::NumOfNearNeighbor); // = 6; // number of near neighbors to consider
  m_beamCrossingAngle = gmc.GlobalParamD.at(GlobalMethodsClass::BeamCrossingAngle) / 2.;

  // the minimal energy to take into account in the initial clustering pass is
  // defined as m_middleEnergyHitBoundFrac of the minimal energy that is taken into
  // account when computing weighted averages in the log' weighting method
  m_middleEnergyHitBoundFrac = gmc.GlobalParamD.at(GlobalMethodsClass::MiddleEnergyHitBoundFrac); // =.01;

  /* --------------------------------------------------------------------------
     constants set by: GlobalMethodsClass
     -------------------------------------------------------------------------- */
  m_logWeightConst = gmc.GlobalParamD.at(GlobalMethodsClass::LogWeightConstant);
  m_moliereRadius = gmc.GlobalParamD.at(GlobalMethodsClass::MoliereRadius);

  // minimal separation distance and energy (of either cluster) to affect a merge
  m_minSeparationDistance = gmc.GlobalParamD.at(GlobalMethodsClass::MinSeparationDist);
  m_minClusterEngyGeV = gmc.GlobalParamD.at(GlobalMethodsClass::MinClusterEngyGeV);

  m_thetaContainmentBounds[0] = gmc.GlobalParamD.at(GlobalMethodsClass::ThetaMin);
  m_thetaContainmentBounds[1] = gmc.GlobalParamD.at(GlobalMethodsClass::ThetaMax);

  m_maxLayerToAnalyse = gmc.GlobalParamI.at(GlobalMethodsClass::NumCellsZ);
  m_cellRMax = gmc.GlobalParamI.at(GlobalMethodsClass::NumCellsR);
  m_cellPhiMax = gmc.GlobalParamI.at(GlobalMethodsClass::NumCellsPhi);

  m_zFirstLayer = gmc.GlobalParamD.at(GlobalMethodsClass::ZStart);
  m_rMin = gmc.GlobalParamD.at(GlobalMethodsClass::RMin);
  m_rMax = gmc.GlobalParamD.at(GlobalMethodsClass::RMax);

  m_rCellLength = gmc.GlobalParamD.at(GlobalMethodsClass::RCellLength);
  m_phiCellLength = gmc.GlobalParamD.at(GlobalMethodsClass::PhiCellLength);

  /* --------------------------------------------------------------------------
     Print out Parameters
     -------------------------------------------------------------------------- */
#if _GENERAL_CLUSTERER_DEBUG == 1
  info() << std::endl << "Global parameters for LumiCalClustererClass:" << std::endl;
  info() << " m_cellRMax: " << m_cellRMax << std::endl
         << " m_cellPhiMax: " << m_cellPhiMax << std::endl
         << " m_zFirstLayer: " << m_zFirstLayer << std::endl
         << " m_zLayerThickness: " << m_zLayerThickness << std::endl
         << " m_zLayerPhiOffset[deg]: " << m_zLayerPhiOffset * 180. / M_PI << std::endl
         << " m_rMin: " << m_rMin << std::endl
         << " m_rMax: " << m_rMax << std::endl
         << " m_rCellLength [mm]: " << m_rCellLength << std::endl
         << " m_phiCellLength [rad]:" << m_phiCellLength << std::endl
         << " m_methodCM: " << m_methodCM << std::endl
         << " m_logWeightConst: " << m_logWeightConst << std::endl
         << " m_elementsPercentInShowerPeakLayer: " << m_elementsPercentInShowerPeakLayer << std::endl
         << " m_moliereRadius: " << m_moliereRadius << std::endl
         << " m_minSeparationDistance: " << m_minSeparationDistance << std::endl
         << " m_minClusterEngy - GeV: " << m_minClusterEngyGeV << std::endl
         << " m_hitMinEnergy: " << m_hitMinEnergy << std::endl
         << " m_thetaContainmentBounds[0]: " << m_thetaContainmentBounds[0] << std::endl
         << " m_thetaContainmentBounds[1]: " << m_thetaContainmentBounds[1] << std::endl
         << " m_middleEnergyHitBoundFrac: " << m_middleEnergyHitBoundFrac << std::endl
         << " Clustering Options : " << std::endl
         << "           _CLUSTER_MIDDLE_RANGE_ENGY_HITS    " << _CLUSTER_MIDDLE_RANGE_ENGY_HITS << std::endl
         << "           _MOLIERE_RADIUS_CORRECTIONS        " << _MOLIERE_RADIUS_CORRECTIONS << std::endl
         << "           _CLUSTER_MIXING_ENERGY_CORRECTIONS " << _CLUSTER_MIXING_ENERGY_CORRECTIONS << std::endl
         << std::endl;
#endif
}

/* ============================================================================
   main actions in each event:
   ========================================================================= */
std::pair<RETVAL, edm4hep::CalorimeterHitCollection>
LumiCalClustererClass::processEvent(edm4hep::SimCalorimeterHitCollection& col) {
  // increment / initialize global variables
  m_totEngyArm[-1] = m_totEngyArm[1] = 0.;
  m_numHitsInArm[-1] = m_numHitsInArm[1] = 0;

  MapIntMapIntVCalHit calHits;
  MapIntMapIntLCCluster superClusterCM;

  MapIntMapIntCalHit calHitsCellIdGlobal;

  m_superClusterIdToCellId.clear();
  m_superClusterIdToCellEngy.clear();
  m_superClusterIdClusterInfo.clear();

  /* --------------------------------------------------------------------------
     Loop over all hits in the LCCollection and write the hits into std::vectors
     of IMPL::CalorimeterHitImpl. Hits are split in two std::vectors, one for each arm
     of LumiCal.
     -------------------------------------------------------------------------- */
  auto [retval, calohits] = getCalHits(col, calHits);
  if (!retval)
    return {RETVAL::NOK, std::move(calohits)};

  /* --------------------------------------------------------------------------
     cccccccccccccc
     --------------------------------------------------------------------------
  const int numArmsToCluster = _armsToCluster.size();
    int armNow = _armsToCluster[armToClusterNow];
 */
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

  return {RETVAL::OK, std::move(calohits)};
}
