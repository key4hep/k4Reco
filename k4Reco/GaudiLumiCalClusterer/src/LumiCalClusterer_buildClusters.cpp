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

#include "Global.h"
#include "LCCluster.h"
#include "LumiCalClusterer.h"
#include "LumiCalHit.h"

#include <TF1.h>
#include <TH1.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

int LumiCalClustererClass::buildClusters(const MapIntVCalHit& calHits, MapIntCalHit& calHitsCellIdGlobal,
                                         MapIntVInt& superClusterIdToCellId, MapIntVDouble& superClusterIdToCellEngy,
                                         MapIntLCCluster& superClusterCM, const int detectorArm) {
  int maxEngyLayerN(-1);
  double maxEngyLayer;
  int numSuperClusters;

  VMapIntCalHit calHitsCellId(m_maxLayerToAnalyse), calHitsSmallEngyCellId(m_maxLayerToAnalyse);

  std::vector<std::map<int, int>> cellIdToClusterId(m_maxLayerToAnalyse + 1);

  std::vector<std::map<int, std::vector<int>>> clusterIdToCellId(m_maxLayerToAnalyse + 1);

  std::vector<std::map<int, LCCluster>> clusterCM(m_maxLayerToAnalyse + 1);
  std::map<int, LCCluster>::iterator clusterCMIterator;

  std::map<int, int> cellIdToSuperClusterId;
  std::map<int, int>::iterator cellIdToSuperClusterIdIterator;

  std::map<int, LCCluster>::iterator superClusterCMIterator;

  std::vector<int> initialClusterControlVar(4), isShowerPeakLayer(m_maxLayerToAnalyse);

  std::map<int, std::vector<LCCluster>> engyPosCMLayer;
  std::map<int, std::vector<LCCluster>>::iterator engyPosCMLayerIterator;

  std::map<int, std::vector<int>> thisLayer;

  std::vector<std::map<int, edm4hep::Vector3d>> virtualClusterCM(m_maxLayerToAnalyse);

  std::map<int, TH1F> xLineFitCM, yLineFitCM;
  std::vector<std::vector<double>> fitParamX, fitParamY;

  std::map<int, double> layerToPosX, layerToPosY, layerToEngy;
  std::map<int, double>::iterator layerToPosXYIterator;

  /* --------------------------------------------------------------------------
     determine the total energy of the hits in the arm
     -------------------------------------------------------------------------- */
#if _GENERAL_CLUSTERER_DEBUG == 1 || _CLUSTER_BUILD_DEBUG == 1
  std::string detectorArmName;
  if (detectorArm > 0)
    detectorArmName = "positive detector arm";
  if (detectorArm < 0)
    detectorArmName = "negative detector arm";
  m_alg->debug() << "************************ buildClusters Arm " << detectorArm
                 << " *****************************************\n";
  m_alg->debug() << "\tTotal " << detectorArmName << " energy =  " << m_totEngyArm[detectorArm] << endmsg << endmsg;
#endif

  /* --------------------------------------------------------------------------
    1.  set the <middleEnergyHitBound> minimal energy to take into account in the initial  clustering pass
    2.  determine the <minNumElementsInShowerPeakLayer>  min number of hits that makeup a showerPeak layer and
        flag the layers that make the cut.
     -------------------------------------------------------------------------- */

  const double middleEnergyHitBound =
      exp(-1 * m_logWeightConst) * m_totEngyArm[detectorArm] * m_middleEnergyHitBoundFrac;
  const int minNumElementsInShowerPeakLayer = int(m_numHitsInArm[detectorArm] * m_elementsPercentInShowerPeakLayer);

#if _CLUSTER_BUILD_DEBUG == 1
  for (MapIntVCalHit::const_iterator calHitsIt = calHits.begin(); calHitsIt != calHits.end(); ++calHitsIt) {
    m_alg->debug() << "Hits in layer " << std::setw(3) << calHitsIt->first << std::setw(6) << calHitsIt->second.size()
                   << endmsg;
  }

  m_alg->debug() << "Shower peak layers:" << endmsg;
  m_alg->debug() << "\t min # of hits, min energy/hit([signal]) : " << minNumElementsInShowerPeakLayer << "\t("
                 << middleEnergyHitBound << " , "
                 << ")" << endmsg << "\t layers chosen : \n";
#endif

  // 1. set ShowerPeak tags in isShowerPeakLayer vector
  // 2. fill calHitsCellId with layer tagged calhits, skip hits with energy below <minHitEnergy>
  // 3. for ShowerPeak layers optionally split hits according to value <middleEnergyHitBound>

  for (const auto& [layerNow, hitsVec] : calHits) {
    std::size_t numHitsInLayer = hitsVec.size();
    isShowerPeakLayer[layerNow] = ((int)numHitsInLayer > minNumElementsInShowerPeakLayer) ? 1 : 0;
#if _CLUSTER_BUILD_DEBUG == 1
    if (isShowerPeakLayer[layerNow] == 1)
      m_alg->debug() << "\t" << layerNow << "\t nhits(" << numHitsInLayer << ")\n";
#endif
    for (std::size_t j = 0; j < numHitsInLayer; j++) {
      int cellIdHit = hitsVec[j]->getCellID0();
      double cellEngy = hitsVec[j]->getEnergy();
      if (cellEngy >= m_minHitEnergy) {
#if _CLUSTER_MIDDLE_RANGE_ENGY_HITS == 1
        /* split hits in ShowerPeakLayer into two sets one with hit energy below
         * and the other above middleEnergyHitBound
         */
        if (isShowerPeakLayer[layerNow]) {
          if (cellEngy <= middleEnergyHitBound) {
            calHitsSmallEngyCellId[layerNow][cellIdHit] = hitsVec[j];
          } else {
            calHitsCellId[layerNow][cellIdHit] = hitsVec[j];
          }
        } else {
          calHitsCellId[layerNow][cellIdHit] = hitsVec[j];
        }
#else
        // all hits assigned to one set
        calHitsCellId[layerNow][cellIdHit] = hitsVec[j];
#endif
      }
    }
  }

  /* --------------------------------------------------------------------------
     fill calHitsCellId with the layer tagged cal hits. for showerPeak layers
     separate cal hits with energy above/below the middleEnergyHitBound.
     in any case only choose hits with energy above the m_minHitEnergy cut
     -------------------------------------------------------------------------- */
  /*(BP)
  for (MapIntVCalHit::const_iterator calHitsIt = calHits.begin(); calHitsIt!=calHits.end(); ++calHitsIt) {
    //  for(size_t layerNow = 0; layerNow < m_maxLayerToAnalyse; layerNow++) {
    for(size_t j=0; j<calHitsIt->second.size(); j++){
      int       cellIdHit = (int)calHitsIt->second[j]->getCellID0();
      double    cellEngy = (double)calHitsIt->second[j]->getEnergy();
      const int layerNow = calHitsIt->first;
      if(cellEngy < m_minHitEnergy) continue;
      //(BP) Bug ? No matter what is value of the option: _CLUSTER_MIDDLE_RANGE_ENGY_HITS
        -  calHitsCellId contains all hits
      if(cellEngy > middleEnergyHitBound) calHitsCellId[layerNow][cellIdHit] = calHitsIt->second[j];

#if _CLUSTER_MIDDLE_RANGE_ENGY_HITS == 1
      if(isShowerPeakLayer[layerNow] == 1) {
        if(cellEngy <= middleEnergyHitBound) calHitsSmallEngyCellId[layerNow][cellIdHit] = calHitsIt->second[j];
      } else {
        calHitsCellId[layerNow][cellIdHit] = calHitsIt->second[j];
      }
#endif
    }
  }
  */

  /* --------------------------------------------------------------------------
     form initial clusters for the shower-peak layers
     -------------------------------------------------------------------------- */
  m_alg->debug() << "run initialClusterBuild and initialLowEngyClusterBuild:" << endmsg;

  // set the control vector for the initialClusterBuild clustering options
  initialClusterControlVar[0] = 1; // mergeOneHitClusters
  initialClusterControlVar[1] = 1; // mergeSmallToLargeClusters
  initialClusterControlVar[2] = 1; // mergeLargeToSmallClusters
  initialClusterControlVar[3] = 1; // forceMergeSmallToLargeClusters

  for (size_t layerNow = 0; layerNow < m_maxLayerToAnalyse; layerNow++) {
    if (isShowerPeakLayer[layerNow] == 1) {
      // run the initial clustering algorithm for the high energy hits
      m_alg->debug() << "\t layer " << layerNow << endmsg;
      // streamlog_message(
      //     DEBUG2, std::stringstream p; for (auto const& callHitCellID : calHitsCellId[layerNow]) {
      //       int cellId = callHitCellID.first;
      //       const auto& calHit = callHitCellID.second;
      //       const double* pos = calHit->getPosition();
      //       p << "\t\t CellId, pos(x,y,z), signal energy [MeV]: " << cellId << "\t (" << pos[0] << ", " << pos[1]
      //         << ", " << pos[2] << "), " << std::fixed << std::setprecision(3) << 1000. * (calHit->getEnergy())
      //         << endmsg;
      //     } p << endmsg;
      //     , p.str();)
      initialClusterBuild(calHitsCellId[layerNow],     // <--
                          cellIdToClusterId[layerNow], // -->
                          clusterIdToCellId[layerNow], // -->
                          clusterCM[layerNow],         // -->
                          initialClusterControlVar);   // <--

#if _CLUSTER_MIDDLE_RANGE_ENGY_HITS == 1
      // cluster the low energy hits
      initialLowEngyClusterBuild(calHitsSmallEngyCellId[layerNow], calHitsCellId[layerNow], cellIdToClusterId[layerNow],
                                 clusterIdToCellId[layerNow], clusterCM[layerNow]);
#endif
      // store max number of hits in ShowerPeakLayer

#if _CLUSTER_BUILD_DEBUG == 1
      dumpClusters(clusterCM[layerNow]);
#endif
    }
  }

  /* --------------------------------------------------------------------------
     check how many global clusters there are
     find the most frequent value of numClusters within ShowerPeakLayer set
     -------------------------------------------------------------------------- */
  MapIntInt numClustersCounter;
  // find the number of clusters in the majority of layers
  m_alg->debug() << "Searching number of clusters in the majority of layers..." << endmsg;
  for (size_t layerNow = 0; layerNow < m_maxLayerToAnalyse; layerNow++) {
    if (isShowerPeakLayer[layerNow] == 1) {
      const int numClusters = clusterCM[layerNow].size();
      m_alg->debug() << "\t -> layer " << layerNow << "\t global clusters " << numClusters << endmsg;
      numClustersCounter[numClusters]++;
    }
  }

  auto maxCluster = std::ranges::max_element(numClustersCounter, {}, &std::pair<const int, int>::second);
  int numClustersMajority = maxCluster->first;

  m_alg->debug() << "\t -> Found that there are " << maxCluster->second << " ShowerPeakLayers \n"
                 << "\t with " << numClustersMajority << " global clusters" << endmsg << endmsg;
  numClustersCounter.clear();

  /* --------------------------------------------------------------------------
     choose only layers which have a numClustersMajority number of clusters,
     and add their CM to engyPosCMLayer[clusterNow] the decision to choose a
     certain clusterNow is made according to the distance of the given
     cluster CM from an averaged CM
     -------------------------------------------------------------------------- */

  // find the layer with the highest energy which has numClustersMajority clusters
  maxEngyLayer = 0.;
  for (size_t layerNow = 0; layerNow < m_maxLayerToAnalyse; layerNow++) {
    if (isShowerPeakLayer[layerNow] == 1) {
      clusterCMIterator = clusterCM[layerNow].begin();
      const int numClusters = clusterCM[layerNow].size();

      double engyLayerNow = 0.;
      if (numClusters != numClustersMajority)
        continue;
      for (int clusterNow1 = 0; clusterNow1 < numClustersMajority; clusterNow1++, clusterCMIterator++) {
        int clusterId = (int)(*clusterCMIterator).first;

        engyLayerNow += clusterCM[layerNow][clusterId].getE();
      }

      if (maxEngyLayer < engyLayerNow) {
        maxEngyLayer = engyLayerNow;
        maxEngyLayerN = layerNow;
      }
    }
  }
  // for the layer with the most energy which has numClustersMajority clusters,
  // initialize the averageCM vector
  // (BP) FIX: Implicit assumptions is that ShowerPeakLayer is also maxEngyLayer, may not be true
  //      but anyway we do not need loop , frome above loop we know layer number

  std::vector<LCCluster> avrgCM;
  /*
  for(size_t layerNow = 0; layerNow < m_maxLayerToAnalyse; layerNow++) {
    if(isShowerPeakLayer[layerNow] == 1) {
      clusterCMIterator = clusterCM[layerNow].begin();
      const int numClusters       = clusterCM[layerNow].size();
      if( (numClusters != numClustersMajority) || (layerNow != maxEngyLayerN)) continue;

      for(int clusterNow = 0; clusterNow < numClustersMajority; clusterNow++, clusterCMIterator++){
      int clusterId = (int)(*clusterCMIterator).first;

      // store the CM energy/position vector for each cluster
      engyPosCMLayer[clusterNow].push_back( clusterCM[layerNow][clusterId] );
      thisLayer[clusterNow].push_back( layerNow );

      // initialize the averaged CM position vector
      avrgCM.push_back( clusterCM[layerNow][clusterId] );
      }
      break; // only do it in isShowerPeakLayer==1
    }// if isShowerPeakLayer
  } //check all Layers
  */
  if (maxEngyLayerN + 1 > 0) {
    clusterCMIterator = clusterCM[maxEngyLayerN].begin();
    for (int clusterNow = 0; clusterNow < numClustersMajority; clusterNow++, clusterCMIterator++) {
      int clusterId = (int)(*clusterCMIterator).first;
      // store the CM energy/position vector for each cluster
      engyPosCMLayer[clusterNow].push_back(clusterCM[maxEngyLayerN][clusterId]);
      thisLayer[clusterNow].push_back(maxEngyLayerN);
      // initialize the averaged CM position vector
      avrgCM.push_back(clusterCM[maxEngyLayerN][clusterId]);
    }
  }

  // for all layers in ShowerPeak except the layer with the most energy which has numClustersMajority clusters,
  // update the averageCM vector

  for (size_t layerNow = 0; layerNow < m_maxLayerToAnalyse; layerNow++) {
    if (isShowerPeakLayer[layerNow] == 1) {
      clusterCMIterator = clusterCM[layerNow].begin();
      int numClustersNow = int(clusterCM[layerNow].size());

      if ((numClustersNow != numClustersMajority) || (static_cast<int>(layerNow) == maxEngyLayerN))
        continue;

      for (int clusterNow1 = 0; clusterNow1 < numClustersMajority; clusterNow1++, clusterCMIterator++) {
        LCCluster const& thisCluster = clusterCM[layerNow][clusterCMIterator->first];
        const double CM1[2] = {thisCluster.getX(), thisCluster.getY()};

        std::map<int, double> weightedDistanceV;
        // compare the position of the CM to the averaged CM positions
        for (int clusterNow2 = 0; clusterNow2 < numClustersMajority; clusterNow2++) {
          const double distanceCM =
              std::hypot(CM1[0] - avrgCM[clusterNow2].getX(), CM1[1] - avrgCM[clusterNow2].getY());
          weightedDistanceV[clusterNow2] = (distanceCM > 0) ? 1. / distanceCM : 1e10;
        }

        auto closestCluster = std::ranges::max_element(weightedDistanceV, {}, &std::pair<const int, double>::second);
        // add the CM to the right vector
        engyPosCMLayer[closestCluster->first].push_back(thisCluster);
        thisLayer[closestCluster->first].push_back(layerNow);

        // update the multi-layer CM position
        // APS: BUGFIX This used to have the CM2 from the clusterNow2 loop above, instead of closestCluster
        avrgCM[closestCluster->first].addToEnergy(thisCluster.getE());
#pragma message("(BP) temporary fix, need modify cluster CM method ")
        //
        double wt_closest = avrgCM[closestCluster->first].getWeight();
        double wt_this = thisCluster.getWeight();
        double new_wt = wt_closest + wt_this;
        double new_x = CM1[0] * wt_this + (avrgCM[closestCluster->first].getX()) * wt_closest;
        double new_y = CM1[1] * wt_this + (avrgCM[closestCluster->first].getY()) * wt_closest;
        new_x /= new_wt;
        new_y /= new_wt;
        avrgCM[closestCluster->first].setX(new_x);
        avrgCM[closestCluster->first].setY(new_y);
        avrgCM[closestCluster->first].setWeight(new_wt);

        //
        //	avrgCM[closestCluster->first].setX( (CM1[0]+(avrgCM[closestCluster->first].getX()))/2.);
        //	avrgCM[closestCluster->first].setY( (CM1[1]+(avrgCM[closestCluster->first].getY()))/2.);

      } // for all clusters
    } // if isShowerPeakLayer
  } // for all layers

  /* --------------------------------------------------------------------------
     fit a stright line through each cluster from engyPosCMLayer. results
     (line parametrizations) are stored in fitParamX(Y)
     engyPosCMLayer contains now info for ShowerPeakLayer only (?!)
     -------------------------------------------------------------------------- */
#if _CLUSTER_BUILD_DEBUG == 1
  m_alg->debug() << "Fit lines through the averaged CM" << endmsg << endmsg;
  m_alg->debug() << "Fit Param should be this size: " << engyPosCMLayer.size() << endmsg;
#endif
  auto sFF = std::make_unique<TF1>("fitFunc", [](double* x, double* p) { return p[0] + p[1] * x[0]; }, -3000, -2000, 2);
  TF1& fitFunc(*sFF);

  //  for(size_t clusterNow=0; clusterNow < engyPosCMLayer.size(); clusterNow++, engyPosCMLayerIterator++) {
  for (engyPosCMLayerIterator = engyPosCMLayer.begin(); engyPosCMLayerIterator != engyPosCMLayer.end();
       engyPosCMLayerIterator++) {
    int clusterId = (int)(*engyPosCMLayerIterator).first;
    int clusterNow = engyPosCMLayerIterator->first;

    if (engyPosCMLayer[clusterId].size() < 3) {
#if _CLUSTER_BUILD_DEBUG == 1
      m_alg->debug() << "\t decrease the global cluster number by 1" << endmsg << endmsg;
#endif
      numClustersMajority--;
      // clusterNow--;//APS
      continue;
    }

#if _CLUSTER_BUILD_DEBUG == 1
    m_alg->debug() << "clusterId " << clusterId << endmsg;
#endif

    std::string hisName = "_xLineFitCM_Cluster";
    std::stringstream clusterNum;
    clusterNum << clusterNow;
    hisName += clusterNum.str();
    xLineFitCM[clusterNow] = TH1F(hisName.c_str(), hisName.c_str(), m_maxLayerToAnalyse * 10, 0, m_maxLayerToAnalyse);

    hisName = "_yLineFitCM_Cluster";
    hisName += clusterNum.str();
    yLineFitCM[clusterNow] = TH1F(hisName.c_str(), hisName.c_str(), m_maxLayerToAnalyse * 10, 0, m_maxLayerToAnalyse);

    /* --------------------------------------------------------------------------
       since more than on cluster may have choosen the same averagedCM in a
       given layer, some layers may have more than one entry in the engyPosCMLayer
       map. therefore an averaging is performed for each layer
       -------------------------------------------------------------------------- */

    // initialize the layer-averaging position/energy maps
    for (std::size_t layerN = 0; layerN < engyPosCMLayer[clusterId].size(); layerN++) {
      const int layerNow = thisLayer[clusterId][layerN];
      layerToPosX[layerNow] = layerToPosY[layerNow] = layerToEngy[layerNow] = 0.;
    }

    // fill the maps with energy-weighted positions
    for (std::size_t layerN = 0; layerN < engyPosCMLayer[clusterId].size(); layerN++) {
      const int layerNow = thisLayer[clusterId][layerN];
      /*(BP) FIX:
       *  it should be weighted with WeightingMethod
       */
      double weightNow = engyPosCMLayer[clusterNow][layerN].getWeight();
      layerToPosX[layerNow] += engyPosCMLayer[clusterNow][layerN].getX() * weightNow;
      layerToPosY[layerNow] += engyPosCMLayer[clusterNow][layerN].getY() * weightNow;
      layerToEngy[layerNow] += weightNow;
      /*
      layerToPosX[layerNow] += engyPosCMLayer[clusterNow][layerN].getX() * engyPosCMLayer[clusterNow][layerN].getE();
      layerToPosY[layerNow] += engyPosCMLayer[clusterNow][layerN].getY() * engyPosCMLayer[clusterNow][layerN].getE();
      layerToEngy[layerNow] += engyPosCMLayer[clusterNow][layerN].getE();
      */
    }

    // fill histograms of x(z) and y(z) of the CM positions
    layerToPosXYIterator = layerToPosX.begin();
    for (std::size_t layerN = 0; layerN < layerToPosX.size(); layerN++, layerToPosXYIterator++) {
      const int layerNow = (int)(*layerToPosXYIterator).first;

      // get back to units of position
      layerToPosX[layerNow] /= layerToEngy[layerNow];
      layerToPosY[layerNow] /= layerToEngy[layerNow];

      xLineFitCM[clusterNow].Fill(layerNow, layerToPosX[layerNow]);
      yLineFitCM[clusterNow].Fill(layerNow, layerToPosY[layerNow]);

#if _CLUSTER_BUILD_DEBUG == 1
      m_alg->debug() << "\tlayer , avPos(x,y) : " << std::setw(3) << layerNow << " (" << std::setw(6)
                     << layerToPosX[layerNow] << " , " << std::setw(6) << layerToPosY[layerNow] << ")" << endmsg;
#endif
    }

    // fit a straight line for each histogram, and store the fit results
    xLineFitCM[clusterNow].Fit("fitFunc", "+CQ0");
    fitParamX.push_back(std::vector<double>(2, 0.0));
    fitParamX.back()[0] = fitFunc.GetParameter(0);
    fitParamX.back()[1] = fitFunc.GetParameter(1);

#if _CLUSTER_BUILD_DEBUG == 1
    m_alg->debug() << "\t -> xFitPar 0,1:  " << fitFunc.GetParameter(0) << " (+-) " << fitFunc.GetParError(0)
                   << " \t,\t " << fitFunc.GetParameter(1) << " (+-) " << fitFunc.GetParError(1) << endmsg;
#endif

    yLineFitCM[clusterNow].Fit("fitFunc", "+CQ0");
    fitParamY.push_back(std::vector<double>(2, 0.0));
    fitParamY.back()[0] = fitFunc.GetParameter(0);
    fitParamY.back()[1] = fitFunc.GetParameter(1);

#if _CLUSTER_BUILD_DEBUG == 1
    m_alg->debug() << "\t -> yFitPar 0,1:  " << fitFunc.GetParameter(0) << " (+-) " << fitFunc.GetParError(0)
                   << " \t,\t " << fitFunc.GetParameter(1) << " (+-) " << fitFunc.GetParError(1) << endmsg << endmsg;
#endif

    // cleanUp
    layerToPosX.clear();
    layerToPosY.clear();
    layerToEngy.clear();
  }
  // cleanUp
  xLineFitCM.clear();
  yLineFitCM.clear();

  /* --------------------------------------------------------------------------
     extrapolate the CM positions in all layers (and form clusters in the
     non shower-peak layers)
     -------------------------------------------------------------------------- */
#if _CLUSTER_BUILD_DEBUG == 1
  m_alg->debug() << "Extrapolate virtual cluster CMs" << endmsg << endmsg;
#endif

  // fill virtual cluster CM vectors for all the layers
  for (size_t layerNow = 0; layerNow < m_maxLayerToAnalyse; layerNow++) {

    if (calHitsCellId[layerNow].empty())
      continue;

    for (int clusterNow = 0; clusterNow < numClustersMajority; clusterNow++) {
      int maxLayerToRaiseVirtualClusterSize = int(0.75 * m_maxLayerToAnalyse);
      double fitPar0, fitPar1, hitLayerRatio;

      if (fitParamX[clusterNow].empty())
        continue; // APS

      edm4hep::Vector3d virtualClusterCMV;

      // extrapolated x/y positions
      virtualClusterCMV.x = fitParamX[clusterNow][0] + fitParamX[clusterNow][1] * layerNow; // x position
      virtualClusterCMV.y = fitParamY[clusterNow][0] + fitParamY[clusterNow][1] * layerNow; // y position

      // ???????? DECIDE/FIX - incorparate the parameters given here better in the code ????????
      // ???????? DECIDE/FIX - consider a different middle layer for the else condition ????????
      // extrapolated cluster radius around CM position
#pragma message("WARNING: Fix these parameters")
      if (avrgCM[clusterNow].getE() > 1) {
        fitPar0 = 236.7;
        fitPar1 = 9.11;
        hitLayerRatio = 22 / 2618.;
      } else {
        fitPar0 = 226.5;
        fitPar1 = 10.3;
        hitLayerRatio = 22 / 2570.;
      }

      if (static_cast<int>(layerNow) < maxLayerToRaiseVirtualClusterSize)
        virtualClusterCMV.z = exp((fitPar0 + fitPar1 * layerNow) * hitLayerRatio);
      else
        virtualClusterCMV.z = exp((fitPar0 + fitPar1 * maxLayerToRaiseVirtualClusterSize) * hitLayerRatio);

      // The numbers above for fitPar0/1 were derived for a detector with moliereRadius=18.2
      // they must, therefore, they must be corrected for according to the m_moliereRadius used now
      virtualClusterCM[layerNow][clusterNow] = virtualClusterCMV;
    }

    // form clusters for the non shower-peak layers in the non shower-peak layers only.
    if (isShowerPeakLayer[layerNow] == 0) {
      try {
        virtualCMClusterBuild(calHitsCellId[layerNow], cellIdToClusterId[layerNow], clusterIdToCellId[layerNow],
                              clusterCM[layerNow], virtualClusterCM[layerNow]);
      } catch (std::out_of_range& e) {
        m_alg->debug() << "exception" << endmsg;
        m_alg->debug() << e.what() << endmsg;
        throw;
      }
#if _CLUSTER_BUILD_DEBUG == 1
      m_alg->debug() << "\tbuild cluster around a virtual CM in layer " << layerNow << endmsg;
#endif
    }

    // only do something if there are less real clusters than virtual clusters in the layer
    else {
      int numRealClusters = clusterIdToCellId[layerNow].size();
      int numVirtualClusters = virtualClusterCM[layerNow].size();
      if (numRealClusters < numVirtualClusters) {
#if _VIRTUALCLUSTER_BUILD_DEBUG == 1
        m_alg->debug() << "\tin layer " << layerNow << " there are real/virtual clusters: " << numRealClusters
                       << "  ,  " << numVirtualClusters << endmsg;
#endif

        virtualCMPeakLayersFix(calHitsCellId[layerNow], cellIdToClusterId[layerNow], clusterIdToCellId[layerNow],
                               clusterCM[layerNow], virtualClusterCM[layerNow]);

#if _CLUSTER_BUILD_DEBUG == 1
        m_alg->debug() << "\tre-cluster in layer " << layerNow << endmsg;
#endif
      }
    }
  }
  // cleanUp
  avrgCM.clear();
  fitParamX.clear();
  fitParamY.clear();

  /* --------------------------------------------------------------------------
     merge all the clusters from the different layers into superClusters.
     the number of final clusters will be numClustersMajority (the majority
     found in the shower-peak layers), and each cluster will be merged into
     a superCluster according to its distance from it (stored as virtual clusters
     in virtualClusterCM).
     -------------------------------------------------------------------------- */

  m_alg->debug() << endmsg << "Building superClusters" << endmsg << endmsg;
  m_alg->debug() << printClusters(superClusterCM);
  int buildSuperClustersFlag =
      buildSuperClusters(calHitsCellIdGlobal, calHitsCellId, clusterIdToCellId, clusterCM, virtualClusterCM,
                         cellIdToSuperClusterId, superClusterIdToCellId, superClusterCM);

  if (buildSuperClustersFlag == 0)
    return 0;

  // cleanUp
  virtualClusterCM.clear();

#if _MOLIERE_RADIUS_CORRECTIONS == 1

  m_alg->debug() << endmsg << "RUN engyInMoliereCorrections() ..." << endmsg;
  m_alg->debug() << printClusters(superClusterCM);

  int engyInMoliereFlag = engyInMoliereCorrections(
      calHitsCellIdGlobal, calHits, calHitsCellId, (clusterIdToCellId), (clusterCM), (cellIdToClusterId),
      (cellIdToSuperClusterId), (superClusterIdToCellId), superClusterCM, middleEnergyHitBound, detectorArm);

  if (engyInMoliereFlag == 0) {
    m_alg->debug() << "Ran engyInMoliereCorrections ... not successful " << endmsg;
    m_alg->debug() << printClusters(superClusterCM);

    return 0;
  }
  m_alg->debug() << "Ran engyInMoliereCorrections ... successful " << endmsg;
  m_alg->debug() << printClusters(superClusterCM);

#endif // #if _MOLIERE_RADIUS_CORRECTIONS == 1

  /* --------------------------------------------------------------------------
     --------------------------------------------------------------------------
     NOTE:
     --------------------------------------------------------------------------
     FROM THIS POINT ON ENERGY OF HITS BELONGING TO A SUPERCLUSTER
     MUST BE ACCESSED BY: superClusterIdToCellEngy
     AND NOT BY:          calHitsCellIdGlobal[cellIdHit]->getEnergy()
     AS THE ENERGY OF SINGLE HITS MAY BE shared BY SEVERAL CLUSTERS !
     --------------------------------------------------------------------------
     -------------------------------------------------------------------------- */

  /* --------------------------------------------------------------------------
     fill the superClusterIdToCellEngy container
     -------------------------------------------------------------------------- */
  superClusterCMIterator = superClusterCM.begin();
  numSuperClusters = superClusterCM.size();
  for (int superClusterNow = 0; superClusterNow < numSuperClusters; superClusterNow++, superClusterCMIterator++) {
    const int superClusterId = (int)(*superClusterCMIterator).first;

    const int numElementsInSuperCluster = superClusterIdToCellId[superClusterId].size();
    for (int cellNow = 0; cellNow < numElementsInSuperCluster; cellNow++) {
      int cellIdHit = superClusterIdToCellId[superClusterId][cellNow];

      double engyHit = calHitsCellIdGlobal.at(cellIdHit)->getEnergy();

      superClusterIdToCellEngy[superClusterId].push_back(engyHit);
    }
  }

  /* --------------------------------------------------------------------------
     verbosity
     -------------------------------------------------------------------------- */

  // streamlog_message(
  //     DEBUG5, std::stringstream p; p << "\tCreated SuperClusters:" << endmsg;
  //     for (auto const& superCluster : superClusterCM) {
  //       const int superClusterId = superCluster.first;
  //       auto const& cluster = superCluster.second;
  //       p << "\t Id " << superClusterId << cluster << "  \t energy " << cluster.getEnergy() << "     \t pos(x,y) =  (
  //       "
  //         << cluster.getX() << " , " << cluster.getY() << " )"
  //         << "     \t pos(theta,phi) =  ( " << cluster.getTheta() << " , " << cluster.getPhi() << " )" << endmsg;
  //     };
  //     , p.str(););

  return 1;
}
