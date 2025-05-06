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
#include "SuperTrueClusterWeights.h"

#include <algorithm>
#include <map>
#include <utility>
#include <vector>

void LumiCalClustererClass::clusterMerger(MapIntVDouble& clusterIdToCellEngy, MapIntVInt& clusterIdToCellId,
                                          MapIntLCCluster& clusterCM, MapIntCalHit& calHitsCellIdGlobal) {
  int clusterId, clusterId1, clusterId2;
  int cellId;
  double engyNow;

  std::map<int, std::vector<int>>::iterator clusterIdToCellIdIterator;

  std::vector<int> allClusterIds, cellIdV;
  std::vector<double> cellEngyV;

  std::map<int, double> idToCellEngy;
  std::map<int, double>::iterator idToCellEngyIterator;

  std::vector<SuperTrueClusterWeights*> clusterPairWeightsV;

  while (1) {
    /* --------------------------------------------------------------------------
       make a std::vector of all of the cluster Ids, and than write the pair-wise
       properties for all pair combinations
       -------------------------------------------------------------------------- */
    allClusterIds.clear();
    clusterIdToCellIdIterator = clusterIdToCellId.begin();

    for (size_t clusterNow = 0; clusterNow < clusterIdToCellId.size(); clusterNow++, clusterIdToCellIdIterator++) {

      clusterId = (int)(*clusterIdToCellIdIterator).first; // Id of cluster
      allClusterIds.push_back(clusterId);
    }

    for (size_t clusterNow1 = 0; clusterNow1 < allClusterIds.size(); clusterNow1++) {
      clusterId1 = allClusterIds[clusterNow1];
      //(BP) avoid creation of "reverse order" pairs
      //     for(int clusterNow2=0; clusterNow2<numClusters2; clusterNow2++) {
      for (size_t clusterNow2 = clusterNow1 + 1; clusterNow2 < allClusterIds.size(); clusterNow2++) {
        clusterId2 = allClusterIds[clusterNow2];

        if (clusterId1 == clusterId2)
          continue;

        SuperTrueClusterWeights* clusterPairWeightsNow =
            new SuperTrueClusterWeights(clusterId1, clusterId2, clusterCM[clusterId1], clusterCM[clusterId2]);
        /* BP. in the following only distance seems to be important
         * so make it simpler
        clusterPairWeightsNow->setWeight( "minEngyDistance",
                                          m_minSeparationDistance,
                                          m_minClusterEngySignal );

        if(clusterPairWeightsNow->weight < 0) {         // at least one condition is not true
          clusterPairWeightsNow->setWeight("distance");
          if( clusterPairWeightsNow->weight <= minSeparationDistance ){
            clusterPairWeightsV.push_back(clusterPairWeightsNow);
          }else{
            delete clusterPairWeightsNow;
          }
        } else {                                        // both conditions are false
          delete clusterPairWeightsNow;
        }
        */
        clusterPairWeightsNow->setWeight("distance");
        if (clusterPairWeightsNow->weight < m_minSeparationDistance)
          clusterPairWeightsV.push_back(clusterPairWeightsNow);
        else
          delete clusterPairWeightsNow;
      }
    }

    // if all the pairs made the cut than finish the loop
    if (clusterPairWeightsV.empty())
      break;

    // choose the pair with the shortest weight (distance)
    sort(clusterPairWeightsV.begin(), clusterPairWeightsV.end(), SuperTrueClusterWeights::Compare);

    clusterId1 = clusterPairWeightsV[0]->superClusterId;
    clusterId2 = clusterPairWeightsV[0]->trueClusterId;

    /* --------------------------------------------------------------------------
       go over all hits in the discarded cluster write a new hit-energy std::map
       -------------------------------------------------------------------------- */
    for (int clusterNow = 0; clusterNow < 2; clusterNow++) {
      if (clusterNow == 0)
        clusterId = clusterId1;
      if (clusterNow == 1)
        clusterId = clusterId2;

      for (size_t hitNow = 0; hitNow < clusterIdToCellId[clusterId].size(); hitNow++) {
        cellId = clusterIdToCellId[clusterId][hitNow];

        engyNow = clusterIdToCellEngy[clusterId][hitNow];
        idToCellEngy[cellId] += engyNow;
      }
    }

    // cleanUp both clusters
    clusterIdToCellId.erase(clusterId1);
    clusterCM.erase(clusterId1);
    clusterIdToCellEngy.erase(clusterId1);
    clusterIdToCellId.erase(clusterId2);
    clusterCM.erase(clusterId2);
    clusterIdToCellEngy.erase(clusterId2);

    /* --------------------------------------------------------------------------
       write the new merged cluster with the Id of one of the clusters
       -------------------------------------------------------------------------- */
    clusterId = clusterId1;
    idToCellEngyIterator = idToCellEngy.begin();
    for (size_t clusterNow = 0; clusterNow < idToCellEngy.size(); clusterNow++, idToCellEngyIterator++) {
      cellId = (int)(*idToCellEngyIterator).first;

      engyNow = idToCellEngy[cellId];

      clusterIdToCellId[clusterId].push_back(cellId);
      clusterIdToCellEngy[clusterId].push_back(engyNow);
    }

    /* --------------------------------------------------------------------------
       compute the total energy and center of mass of the new merged cluster
       -------------------------------------------------------------------------- */
    cellIdV = clusterIdToCellId[clusterId];
    cellEngyV = clusterIdToCellEngy[clusterId];

    // initialize the energy/position std::vector for new clusters only
    clusterCM[clusterId] = LCCluster();

    // calculate/update the energy/position of the CM
    calculateEngyPosCM_EngyV(cellIdV, cellEngyV, calHitsCellIdGlobal, clusterCM, clusterId, m_methodCM);

    // cleanUp
    cellIdV.clear();
    cellEngyV.clear();

    /* --------------------------------------------------------------------------
       cleanUp virtual memory allocations
       -------------------------------------------------------------------------- */
    for (size_t clusterNow1 = 0; clusterNow1 < clusterPairWeightsV.size(); clusterNow1++)
      delete clusterPairWeightsV[clusterNow1];

    clusterPairWeightsV.clear();
    idToCellEngy.clear();
  }

  /* --------------------------------------------------------------------------
     verbosity
     -------------------------------------------------------------------------- */
#if _GENERAL_CLUSTERER_DEBUG == 1
  if (not clusterIdToCellId.empty()) {
    m_alg->debug() << "Clusters:" << endmsg;
  } else {
    m_alg->debug() << "No Clusters on this side" << endmsg;
  }

  clusterIdToCellIdIterator = clusterIdToCellId.begin();
  numClusters = clusterIdToCellId.size();
  for (int clusterNow = 0; clusterNow < numClusters; clusterNow++, clusterIdToCellIdIterator++) {
    clusterId = (int)(*clusterIdToCellIdIterator).first;

    m_alg->debug() << "Id " << std::setw(4) << clusterId << clusterCM[clusterId] << endmsg;
  }

#endif

  // cleanUp
  allClusterIds.clear();
}
