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
#include "GlobalMethodsClass.h"
#include "LCCluster.h"
#include "LumiCalClusterer.h"
#include "LumiCalHit.h"
#include "ProjectionInfo.h"
#include "VirtualCluster.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <map>
#include <memory>
#include <utility>
#include <vector>
#include <ranges>

int LumiCalClustererClass::initialClusterBuild(const MapIntCalHit& calHitsCellId, MapIntInt& cellIdToClusterId,
                                               MapIntVInt& clusterIdToCellId, MapIntLCCluster& clusterCM,
                                               const VInt& controlVar) {
  /* --------------------------------------------------------------------------
     layer parameters
     -------------------------------------------------------------------------- */
  // (*). Small clusters are merged with large clusters if they are close. A weight is computed for each
  //	close cluster candidate, and the max-weighted cluster is chosen. This is done twice. First small
  //	clusters try to be merged into large ones, and then large clusters try to merge small ones into themselvs.
  // ----------------------------------------------------------------------------------------------------------------
  //	(1). definitions of 'small' and 'large' sizes
  const size_t numElementsInLayer(calHitsCellId.size());
  const int numElementsSmallClusterToMerge(0.10 * numElementsInLayer); // 5;  // max number of elements in small cluster
  const int numElementsLargeClusterToMerge(0.15 *
                                           numElementsInLayer); // 10; // min number of elements in a large cluster

  //	(2). decide to merge if the two clusters are close to each other, and if the merged cluster's CM
  //	isnt shifted to an area where there are few cal hits (small percentage of CM energy)
  const double mergeScanDistanceFromLargeClusterCM(
      m_moliereRadius); // max distance of small cluster from CM of large one
  const double distanceToCollectEngyAroundCM(.5 * m_moliereRadius); // distance from the CM to collect energy
  const double percentOfEngyAroungCM(0.4);                          // percentage of cluster CM energy

  //	(3). the weight for merging two clusters is a Power() function of the energy of the large cluster and the
  //	distance to its CM. two different weight functions for the two pases.
  int smallClusterEngyCMPower = 1;
  int smallClusterDistanceCMPower = -3;

  int largeClusterEngyCMPower = 1;
  int largeClusterDistanceCMPower = -1;

  // (*). Remaining small clusters are force-merged into large clusters. The choise of who to merge with is done through
  //	weights.
  // -------------------------------------------------------------------------------------------------------------------
  //	(1). definitions of 'small' and 'large' sizes
  const unsigned int numElementsSmallClusterToMergeForced =
      .08 * numElementsInLayer; // 5;  // max number of elements in small cluster

  //	(2). the weight for merging two clusters is a Power() function of the energy of the
  //		large cluster and the distance to its CM.
  int smallClusterEngyCMPowerForced = 0;
  int smallClusterDistanceCMPowerForced = -1;

  // control options for enabeling the different rutines (1 is for active)
  int mergeOneHitClusters = controlVar.at(0);
  int mergeSmallToLargeClusters = controlVar.at(1);
  int mergeLargeToSmallClusters = controlVar.at(2);
  int forceMergeSmallToLargeClusters = controlVar.at(3);

  /* --------------------------------------------------------------------------
     maps for................
     -------------------------------------------------------------------------- */
  // map param: (1). cellId of cal hit , (2). Id of cluster the cal hit belongs to
  // std::map < int , int >		cellIdToClusterId;

  // map param: (1). Id of cluster the cal hit belongs to , (2). std::vector of cellIds of cal hit
  // std::map < int , std::vector<int> >		clusterIdToCellId;
  // std::map < int , std::vector<int> > :: iterator	clusterIdToCellIdIterator  ;

  // std::vector for holding the Ids of cells/clusters
  std::vector<int> clusterIdV;

  // counter for cluster Id
  int clusterId = 0;

  /* --------------------------------------------------------------------------
     maps for storing for each cal hit the highest-energy near neighbor, and for
     each cal hit to keep track of the other cal hits for which it is the highest
     energy nearest neighbor.
     -------------------------------------------------------------------------- */
  // map param: (1). cellId of cal hit , (2). cellId of highest energy near neighbor
  std::map<int, int> isConnectedToNeighbor;
  // map param: (1). cellId of cal hit , (2). cellIds of neighbors which are connected to the cal hit
  std::map<int, std::vector<int>> neighborsConectedToMe;

  // copy hits in this layer to a cal hit std::vector
  VecCalHit calHitsLayer;

  for (const auto& [cellIdHit, calHit] : calHitsCellId) {
    calHitsLayer.push_back(calHit);
    // initialization
    isConnectedToNeighbor[cellIdHit] = 0;
    cellIdToClusterId[cellIdHit] = 0;
  }

  // sort acording to energy in ascending order (lowest energy is first)
  std::ranges::sort(calHitsLayer, [](const auto& a, const auto& b) { return a->getEnergy() < b->getEnergy(); });

  /* --------------------------------------------------------------------------
     connect each cal hit to it's highest-energy near neighbor.
     -------------------------------------------------------------------------- */
  for (const auto& calHit : calHitsLayer) {

    int cellIdHit = calHit->getCellID0();

    // go on to next cal hit if this hit has already been registered
    if (isConnectedToNeighbor[cellIdHit] != 0)
      continue;

    bool neighborFound = false;
    do {
      double maxEngyNeighbor = 0.;
      double engyCalHit = calHitsCellId.at(cellIdHit)->getEnergy();

      // get the Ids of the neighbor
      for (int neighborIndex = 0; neighborIndex < m_nNearNeighbor; neighborIndex++) {
        // find cellId of neighbor
        const int cellIdNeighbor = getNeighborId(cellIdHit, neighborIndex);
        if (cellIdNeighbor == 0)
          continue;
        // if the neighbor has a cal hit...
        const auto neighborIt = calHitsCellId.find(cellIdNeighbor);
        if (neighborIt == calHitsCellId.end()) {
          continue;
        }
        const auto& neighbor = neighborIt->second;
        double engyNeighbor = neighbor->getEnergy();
        if ((maxEngyNeighbor < engyNeighbor) && (engyNeighbor >= engyCalHit)) {
          // register the neighbor at the cal hit
          isConnectedToNeighbor[cellIdHit] = cellIdNeighbor;
          // update highest-energy counter
          maxEngyNeighbor = engyNeighbor;
        }

      } // For all neighbor

      if (maxEngyNeighbor > 0) {
        // check if the neighbor has already been registered
        const int cellIdNeighbor = isConnectedToNeighbor[cellIdHit];

        // modify the conditional control variable
        if (isConnectedToNeighbor[cellIdNeighbor] != 0)
          neighborFound = true;

        // register the cal hit at the neighbor
        neighborsConectedToMe[cellIdNeighbor].push_back(cellIdHit);
        // in the next iteration work with the cal hit's highest-energy near neighbor
        cellIdHit = isConnectedToNeighbor[cellIdHit];
      } else {
        neighborFound = true;
      }
    } while (!neighborFound); // while not neighborFound
  } // for all hits in layer

  /* --------------------------------------------------------------------------
     create clusters from each connected bunch of cal hits.
     -------------------------------------------------------------------------- */
  // sort according to energy in descending order (highest energy is first)
  std::ranges::sort(calHitsLayer, [](const auto& a, const auto& b) { return a->getEnergy() > b->getEnergy(); });

  for (const auto& calHit : calHitsLayer) {

    std::vector<int> neighborFoundId;

    int cellIdHit = calHit->getCellID0();

    // if the cal hit has already been registered in a cluster continue to the next one
    if (cellIdToClusterId[cellIdHit] > 0)
      continue;

    // add the hit to a new cluster
    clusterId++;
    cellIdToClusterId[cellIdHit] = clusterId;
    clusterIdToCellId[clusterId].push_back(cellIdHit);

    bool neighborFound = true;
    while (neighborFound) {
      // get the Ids of the connected neighbor
      for (const auto cellIdNeighbor : neighborsConectedToMe[cellIdHit]) {

        // register the neighbor in the cluster
        cellIdToClusterId[cellIdNeighbor] = clusterId;
        clusterIdToCellId[clusterId].push_back(cellIdNeighbor);

        // store the connected cal hit Id - it will be registered in a cluster next...
        neighborFoundId.push_back(cellIdNeighbor);
      }

      // in case of a loop in connections, make sure that if a cell has already been
      // considered, then it wont be reanalyzed again...
      neighborsConectedToMe[cellIdHit].clear();

      // every neighbor found is taken off from the neighborFoundId std::vector
      // once this std::vector is empty, end the while(neighborFound == 1) loop

      if (neighborFoundId.size() == 0) {
        neighborFound = false;
      } else {
        cellIdHit = neighborFoundId.back();
        neighborFoundId.pop_back();
      }
    }
  }
  // cleanUp
  isConnectedToNeighbor.clear();
  neighborsConectedToMe.clear();

  /* --------------------------------------------------------------------------
     merge clusters that have one hit only to the nearest cluster which is
     a near neighbor (choose the neighbor with the highest energy).
     -------------------------------------------------------------------------- */
  if (mergeOneHitClusters == 1) {
    // go over all cellIds which have been registered in a cluster (first parameter of cellIdToClusterId).
    std::map<int, int>::iterator cellIdToClusterIdIterator = cellIdToClusterId.begin(), cEnd = cellIdToClusterId.end();
    for (; cellIdToClusterIdIterator != cEnd; ++cellIdToClusterIdIterator) {
      const int cellIdHit = cellIdToClusterIdIterator->first;     // cellId of cal hit
      const int clusterIdHit = cellIdToClusterIdIterator->second; // cluster Id of cal hit

      if (clusterIdToCellId[clusterIdHit].size() == 1) {
        int maxEngyNeighborId = 0;
        double maxEngyNeighbor = 0.;

        // get the Ids of the neighbor
        for (int neighborIndex = 0; neighborIndex < m_nNearNeighbor; neighborIndex++) {
          // find cellId of neighbor
          const int cellIdNeighbor = getNeighborId(cellIdHit, neighborIndex);
          if (cellIdNeighbor == 0)
            continue;

          // if the neighbor has a cal hit...
          if (calHitsCellId.find(cellIdNeighbor) != calHitsCellId.end()) {
            double engyNeighbor = calHitsCellId.at(cellIdNeighbor)->getEnergy();

            // find the neighbor with the highest energy
            if (maxEngyNeighbor < engyNeighbor) {
              // store the cellId of the neighbor
              maxEngyNeighborId = cellIdNeighbor;
              // update highest-energy counter
              maxEngyNeighbor = engyNeighbor;
            }
          }
        }

        if (maxEngyNeighbor > 0) {
          // connect to the neighbor with the highest energy
          int clusterIdNeighbor = cellIdToClusterId[maxEngyNeighborId];
          cellIdToClusterId[cellIdHit] = clusterIdNeighbor;
          clusterIdToCellId[clusterIdNeighbor].push_back(cellIdHit);

          // cleanUp the discarded cluster
          clusterIdToCellId.erase(clusterIdHit);
        }
      }
    }
  } // mergeOneHitClusters

  /* --------------------------------------------------------------------------
   *     compute total energy and center of mass of each cluster
   *     create clusterCM map from scratch
     -------------------------------------------------------------------------- */
  clusterCM.clear();

  for (const auto& [clsID, cellIds] : clusterIdToCellId) {
    clusterCM[clsID] = calculateEngyPosCM(cellIds, calHitsCellId, m_methodCM);
  }
  // printMapVector(__func__,__LINE__,clusterCM);

  /* --------------------------------------------------------------------------
     merge clusters with close CM positions -
     (1). first go over small clusters and merge with large clusters.
     -------------------------------------------------------------------------- */
  if (mergeSmallToLargeClusters == 1) {
    // sort the clusterIds according to clusterCM energy in ascending order (lowest energy is first)
    VVDouble clusterIdEngyV2;

    for (const auto& [clsID, cellIds] : clusterIdToCellId) {
      clusterIdEngyV2.emplace_back(std::vector<double>(2, 0.0));
      clusterIdEngyV2.back()[0] = clusterCM[clsID].getE();
      clusterIdEngyV2.back()[1] = clsID;
    }

    std::ranges::sort(clusterIdEngyV2, [](const auto& a, const auto& b) {
      return a[0] < b[0]; // sort by energy in ascending order
    });

    // copy the Ids that are now in order to a std::vector
    for (const auto& id : clusterIdEngyV2)
      clusterIdV.push_back((int)id[1]);

    // cleanUp
    clusterIdEngyV2.clear();

    // merge close neighboring clusters
    while (clusterIdV.size() > 1) {

      std::map<int, double> weightedDistanceV;

      // start with the highest-energy cluster
      if ((int)clusterIdToCellId[clusterIdV[0]].size() > numElementsSmallClusterToMerge) {
        clusterIdV.erase(clusterIdV.begin());
        continue;
      }

      // compute distance to neighboring clusters
      // for (std::size_t j = 1; j < clusterIdV.size(); j++) {
      for (std::size_t j = 1; j < clusterIdV.size(); j++) {
        clusterId = clusterIdV[j];

        // const double distanceCM =
        // distance2D(clusterCM[clusterIdV[0]].getPosition(),clusterCM[clusterId].getPosition());
        const double distanceCM = std::hypot(clusterCM[clusterIdV[0]].getX() - clusterCM[clusterId].getX(),
                                             clusterCM[clusterIdV[0]].getY() - clusterCM[clusterId].getY());

        int considerCloseCluster = 0;
        // choose close clusters
        if ((int)clusterIdToCellId[clusterId].size() >= numElementsLargeClusterToMerge)
          if (distanceCM < mergeScanDistanceFromLargeClusterCM)
            considerCloseCluster =
                checkClusterMergeCM(clusterIdV[0], clusterIdV[j], clusterIdToCellId, calHitsCellId,
                                    distanceToCollectEngyAroundCM, percentOfEngyAroungCM, m_methodCM);

        // if all the conditions regarding the close cluster were satisfied, add
        // it to the list of possible merging partner, one of which will be chosen next
        if (considerCloseCluster == 1) {
          // store the weight of the close cluster
          const double engyCM = clusterCM[clusterId].getE();
          double closeClusterWeight =
              pow(engyCM, smallClusterEngyCMPower) * pow(distanceCM, smallClusterDistanceCMPower);
          weightedDistanceV[clusterId] = closeClusterWeight;
        }
      }

      if (!weightedDistanceV.empty()) {
        // find the close cluster with the highest weightedDistance
        auto maxElement = std::ranges::max_element(weightedDistanceV, {}, &std::pair<const int, double>::second);
        // merge the two clusters
        const int clusterId1 = clusterIdV[0];
        const int clusterId2 = maxElement->first;

        for (const auto cellIdHit : clusterIdToCellId[clusterId1]) {
          // add hit from clusterIdToCellId with clusterId to one with maxWDClusterId
          clusterIdToCellId[clusterId2].push_back(cellIdHit);
          cellIdToClusterId[cellIdHit] = clusterId2;

          //  update the totalEnergy counter and CM position of the cluster
          updateEngyPosCM(calHitsCellId.at(cellIdHit), clusterCM[clusterId2]);
        }

        // cleanUp
        clusterIdToCellId.erase(clusterId1);
        clusterCM.erase(clusterId1);
      }

      // erase the merged cluster Id from the clusterIdV std::vector.
      // if there are no close neighbors than positionInVec still has the initial zero
      // value, and so the initial cluster is erased from clusterIdV.
      clusterIdV.erase(clusterIdV.begin());
    }
    // cleanUp
    clusterIdV.clear();
  }

  /* --------------------------------------------------------------------------
     merge clusters with close CM positions
     (2). go over all clusters and try to merge.
     -------------------------------------------------------------------------- */
  if (mergeLargeToSmallClusters == 1) {
    // sort the clusterIds according to clusterCM energy in
    // descending order (highest energy is first)

    VVDouble clusterIdEngyV2;

    for (const auto& [clsID, cellIds] : clusterIdToCellId) {
      clusterIdEngyV2.push_back(std::vector<double>(2, 0.0));
      clusterIdEngyV2.back()[0] = clusterCM[clsID].getE();
      clusterIdEngyV2.back()[1] = clsID;
    }

    std::ranges::sort(clusterIdEngyV2, [](const auto& a, const auto& b) {
      return a[0] > b[0]; // sort by energy in descending order
    });

    // copy the Ids that are now in order to a std::vector
    for (const auto& id : clusterIdEngyV2)
      clusterIdV.push_back((int)id[1]);

    // cleanUp
    clusterIdEngyV2.clear();

    // merge close neighboring clusters
    while (clusterIdV.size() > 1) {

      std::map<int, double> weightedDistanceV;

      // start with the highest-energy cluster
      clusterId = clusterIdV[0];
      LCCluster const& highestCluster = clusterCM[clusterId];
      double CM1[2] = {highestCluster.getX(), highestCluster.getY()};

      // compute distance to neighboring clusters
      for (std::size_t j = 1; j < clusterIdV.size(); j++) {

        const int clusterIdJ = clusterIdV[j];
        LCCluster const& thisCluster = clusterCM[clusterIdJ];
        double CM2[2] = {thisCluster.getX(), thisCluster.getY()};

        const double distanceCM = std::hypot(CM1[0] - CM2[0], CM1[1] - CM2[1]);

        int considerCloseCluster = 0;
        if (distanceCM < mergeScanDistanceFromLargeClusterCM)
          considerCloseCluster = checkClusterMergeCM(clusterIdV[0], clusterIdV[j], clusterIdToCellId, calHitsCellId,
                                                     distanceToCollectEngyAroundCM, percentOfEngyAroungCM, m_methodCM);

        // if all the conditions regarding the close cluster were satisfied, add
        // it to the list of possible merging partner, one of which will be chosen next
        if (considerCloseCluster == 1) {
          // store the weight of the close cluster
          const double engyCM = thisCluster.getE();
          assert(distanceCM >= 0 && engyCM > 0); // APS:: assert >= now so 0.0 passes

          double closeClusterWeight =
              pow(engyCM, largeClusterEngyCMPower) * pow(distanceCM, largeClusterDistanceCMPower);
          weightedDistanceV[clusterIdJ] = closeClusterWeight;
        }
      } // for all other clusters

      if (not weightedDistanceV.empty()) {
        // find the close cluster with the highest weightedDistance
        auto maxElement = std::ranges::max_element(weightedDistanceV, {}, &std::pair<const int, double>::second);
        // merge the two clusters
        const int clusterId1 = maxElement->first;
        const int clusterId2 = clusterIdV[0];

        for (std::size_t j = 0; j < clusterIdToCellId[clusterId1].size(); j++) {
          // add hit from clusterIdToCellId with clusterId to one with maxWDClusterId
          const int cellIdHit = clusterIdToCellId[clusterId1][j];
          clusterIdToCellId[clusterId2].push_back(cellIdHit);
          cellIdToClusterId[cellIdHit] = clusterId2;

          //  update the totalEnergy counter and CM position of the cluster
          updateEngyPosCM(calHitsCellId.at(cellIdHit), clusterCM[clusterId2]);
        }

        // cleanUp
        clusterIdToCellId.erase(clusterId1);
        clusterCM.erase(clusterId1);

        // erase the merged cluster Id from the clusterIdV std::vector.
        // if there are no close neighbors than positionInVec still has the initial zero
        // value, and so the initial cluster is erased from clusterIdV.
        clusterIdV.erase(std::find(clusterIdV.begin(), clusterIdV.end(), clusterId1));
      } else {
        clusterIdV.erase(clusterIdV.begin());
      }

    } // while
  }

  /* --------------------------------------------------------------------------
     merge small clusters with large ones
     -------------------------------------------------------------------------- */
  if (forceMergeSmallToLargeClusters == 1) {
    std::vector<int> smallClusterIdV, largeClusterIdV;

    // sort clusters into two std::vectors according to the number of elements in each cluster
    for (const auto& [loopClusterId, cellIds] : clusterIdToCellId) {
      // ???????? DECIDE/FIX - improve this number ????????
      if (cellIds.size() <= numElementsSmallClusterToMergeForced)
        smallClusterIdV.push_back(loopClusterId);
      else
        largeClusterIdV.push_back(loopClusterId);
    }

    for (size_t j = 0; j < smallClusterIdV.size(); j++) {

      std::map<int, double> weightedDistanceV;

      // start with the highest-energy cluster
      double CM1[2] = {clusterCM[smallClusterIdV[j]].getX(), clusterCM[smallClusterIdV[j]].getY()};

      // compute distance to neighboring clusters
      for (size_t k = 0; k < largeClusterIdV.size(); k++) {
        const int largeClusterId = largeClusterIdV[k];
        double CM2[2] = {clusterCM[largeClusterId].getX(), clusterCM[largeClusterId].getY()};

        const double distanceCM = std::hypot(CM1[0] - CM2[0], CM1[1] - CM2[1]);

        // if all the conditions regarding the close cluster were satisfied, add
        // it to the list of possible merging partner, one of which will be chosen next

        // store the weight of the close cluster
        const double engyCM = clusterCM[largeClusterId].getE();

        assert(distanceCM > 0 && engyCM > 0);

        double closeClusterWeight =
            pow(engyCM, smallClusterEngyCMPowerForced) * pow(distanceCM, smallClusterDistanceCMPowerForced);
        weightedDistanceV[largeClusterId] = closeClusterWeight;
      }

      if (!weightedDistanceV.empty()) {
        // find the close cluster with the highest weightedDistance
        auto maxElement = std::ranges::max_element(weightedDistanceV, {}, &std::pair<const int, double>::second);
        // merge the two clusters
        const int clusterId1 = smallClusterIdV[j];
        const int clusterId2 = maxElement->first;

        for (std::size_t k = 0; k < clusterIdToCellId[clusterId1].size(); k++) {
          // add hit from clusterIdToCellId with clusterId to one with maxWDClusterId
          const int cellIdHit = clusterIdToCellId[clusterId1][k];
          clusterIdToCellId[clusterId2].push_back(cellIdHit);
          cellIdToClusterId[cellIdHit] = clusterId2;

          //  update the totalEnergy counter and CM position of the cluster
          updateEngyPosCM(calHitsCellId.at(cellIdHit), clusterCM[clusterId2]);
        }

        // cleanUp
        clusterIdToCellId.erase(clusterId1);
        clusterCM.erase(clusterId1);
      }
    }
  }

  return 1;
}

/* =========================================================================
   LumiCalClustererClass :: initialLowEngyClusterBuild
   ============================================================================
   (1). Description:
   --------------------------------
   -     merge the unclustered cal hits with the existing clusters
   ============================================================================ */

int LumiCalClustererClass::initialLowEngyClusterBuild(MapIntCalHit const& calHitsSmallEngyCellId,
                                                      MapIntCalHit& calHitsCellId, MapIntInt& cellIdToClusterId,
                                                      MapIntVInt& clusterIdToCellId, MapIntLCCluster& clusterCM) {
  /* --------------------------------------------------------------------------
     merge the unclustered cal hits with the existing clusters
     -------------------------------------------------------------------------- */
  for (const auto& [cellIdHit, thisHit] : calHitsSmallEngyCellId) {
    // add the small energy hits that have now been clustred to the cal hit list at calHitsCellId
    calHitsCellId[cellIdHit] = thisHit;
    // position of the cal hit
    double CM1[2] = {thisHit->getPosition()[0], thisHit->getPosition()[1]};

    std::map<int, double> weightedDistanceV;
    // compute weight for the cal hit and each cluster
    for (const auto& [clusterId, cellIds] : clusterIdToCellId) {
      LCCluster const& thisCluster = clusterCM[clusterId];
      double CM2[2] = {thisCluster.getX(), thisCluster.getY()};
      const double distanceCM = std::hypot(CM1[0] - CM2[0], CM1[1] - CM2[1]);

      if (distanceCM > 0)
        weightedDistanceV[clusterId] = 1. / distanceCM;
      else
        weightedDistanceV[clusterId] = 1e10;
    }

    // decide to which cluster to merge the hit according to a proper weight
    auto closestCluster = std::ranges::max_element(weightedDistanceV, {}, &std::pair<const int, double>::second);

    // add the hit to the cluster with the max weight
    if (closestCluster != weightedDistanceV.end()) {
      // add the hit to the chosen cluster
      clusterIdToCellId[closestCluster->first].push_back(cellIdHit);
      cellIdToClusterId[cellIdHit] = closestCluster->first;
      //  update the totalEnergy counter and CM position of the cluster
      updateEngyPosCM(calHitsCellId[cellIdHit], clusterCM[closestCluster->first]);
    }
  }

  return 1;
}

int LumiCalClustererClass::virtualCMClusterBuild(MapIntCalHit const& calHitsCellId, MapIntInt& cellIdToClusterId,
                                                 MapIntVInt& clusterIdToCellId, MapIntLCCluster& clusterCM,
                                                 MapIntVirtualCluster const& virtualClusterCM) {
  std::vector<int> unClusteredCellId;

#if _VIRTUALCLUSTER_BUILD_DEBUG == 1
  cout << endl << coutBlue << "Hits inside virtualRadius of virtualClusters" << coutDefault << endl << endl;
#endif

  /* --------------------------------------------------------------------------
     form clusters by gathering hits inside the virtual cluster radius
     -------------------------------------------------------------------------- */
  for (const auto& [cellIdHit, thisHit] : calHitsCellId) {
    double CM1[2] = {thisHit->getPosition()[0], thisHit->getPosition()[1]};

    // compute the distance of the cal hit from the virtual cluster CMs, and keep score of
    // the clusters that the hit is in range of
    std::map<int, double> weightedDistanceV;

    for (const auto& [virtualClusterId, virtualCluster] : virtualClusterCM) {
      double CM2[2] = {virtualCluster.getX(), virtualCluster.getY()};
      const double distanceCM = std::hypot(CM1[0] - CM2[0], CM1[1] - CM2[1]);

      if (distanceCM <= virtualClusterCM.at(virtualClusterId).getZ()) {
        weightedDistanceV[virtualClusterId] = (distanceCM > 0) ? 1. / distanceCM : 1e10;
      }

#if _VIRTUALCLUSTER_BUILD_DEBUG == 1
      cout << "\tclusterId: " << virtualClusterId << " \t hit(x,y): " << CM1[0] << " \t " << CM1[1]
           << " \t virtualCM(x,y):" << CM2[0] << " \t " << CM2[1] << endl
           << "\t\tvirtualRadius: " << virtualClusterCM[virtualClusterId][0] << " \t distance: " << distanceCM << endl
           << endl;
#endif
    }

    // decide to which cluster to merge the hit according to a proper weight
    auto closestCluster = std::ranges::max_element(weightedDistanceV, {}, &std::pair<const int, double>::second);

    // create clusters from the hits which made the cut, and keep score of those which didnt
    if (closestCluster != weightedDistanceV.end()) {
      // add the hit to the chosen cluster
      clusterIdToCellId[closestCluster->first].push_back(cellIdHit);
      cellIdToClusterId[cellIdHit] = closestCluster->first;

#if _VIRTUALCLUSTER_BUILD_DEBUG == 1
      std::cout << "\tnum possible -  " << weightedDistanceV.size() << " \t chosen - " << closestCluster->first
                << coutDefault << endl
                << endl;
#endif

    } else {
      unClusteredCellId.push_back(cellIdHit);

#if _VIRTUALCLUSTER_BUILD_DEBUG == 1
      cout << "\tno cluster is within range ..." << coutDefault << endl << endl;
#endif
    }
  } // for all Calo Hits in this map

  /* --------------------------------------------------------------------------
     compute total energy and center of mass of each cluster
     -------------------------------------------------------------------------- */
  clusterCM.clear();

  for (const auto& [clusterId, cellIds] : clusterIdToCellId) {
    // calculate the energy/position of the CM
    clusterCM[clusterId] = calculateEngyPosCM(cellIds, calHitsCellId, m_methodCM);
  }

  /* --------------------------------------------------------------------------
     if no clusters were created, make them into temporary clusters with
     zero energy (and no cal hits), so that hits will have a chance to merge with them
     -------------------------------------------------------------------------- */
  for (const auto& [virtualClusterId, virtualCluster] : virtualClusterCM) {
    if (clusterIdToCellId[virtualClusterId].empty()) {
      clusterCM[virtualClusterId] = LCCluster(virtualClusterCM.at(virtualClusterId));
    }
  }

  /* --------------------------------------------------------------------------
     merge unclustered cal hits with the existing clusters
     -------------------------------------------------------------------------- */
  for (size_t hitNow = 0; hitNow < unClusteredCellId.size(); hitNow++) {

    const int cellIdHit = unClusteredCellId[hitNow];

    std::map<int, double> weightedDistanceV;

    // position of the cal hit
    const auto& thisHit = calHitsCellId.at(cellIdHit);
    double CM1[2] = {thisHit->getPosition()[0], thisHit->getPosition()[1]};

    // compute weight for the cal hit and each cluster
    for (MapIntVInt::const_iterator clusterIdToCellIdIterator = clusterIdToCellId.begin();
         clusterIdToCellIdIterator != clusterIdToCellId.end(); ++clusterIdToCellIdIterator) {
      LCCluster const& thisCluster = clusterCM[clusterIdToCellIdIterator->first];
      double CM2[2] = {thisCluster.getX(), thisCluster.getY()};
      const double distanceCM = std::hypot(CM1[0] - CM2[0], CM1[1] - CM2[1]);
      weightedDistanceV[clusterIdToCellIdIterator->first] = (distanceCM > 0) ? 1. / distanceCM : 1e10;
    }

    // decide to which cluster to merge the hit according to a proper weight
    auto closestCluster = std::ranges::max_element(weightedDistanceV, {}, &std::pair<const int, double>::second);

    // add the hit to the cluster with the max weight
    if (closestCluster != weightedDistanceV.end()) {
      // add the hit to the chosen cluster
      clusterIdToCellId[closestCluster->first].push_back(cellIdHit);
      cellIdToClusterId[cellIdHit] = closestCluster->first;
      //  update the totalEnergy counter and CM position of the cluster
      updateEngyPosCM(calHitsCellId.at(cellIdHit), clusterCM[closestCluster->first]);
    }
  } // for all unclustered hits

  return 1;
}

int LumiCalClustererClass::virtualCMPeakLayersFix(MapIntCalHit const& calHitsCellId, MapIntInt& cellIdToClusterId,
                                                  MapIntVInt& clusterIdToCellId, MapIntLCCluster& clusterCM,
                                                  MapIntVirtualCluster virtualClusterCM) {
  // general variables
  std::vector<std::vector<double>> unClusteredCellId;

  std::map<int, int> virtualToRealClusterId, oldVirtualClusterCM, oldClusterCM, oldToNewVirtualClusterIds;

  std::map<int, double>::iterator weightedDistanceVIterator;

  /* --------------------------------------------------------------------------
     make sure that all the virtual cluster Ids are different than the
     existing real cluster Ids.
     -------------------------------------------------------------------------- */
  for (const auto& [virtualClusterId, _] : virtualClusterCM) {
    oldVirtualClusterCM[virtualClusterId] = 1;
  }

  for (const auto& [clusterId, _] : clusterIdToCellId) {
    oldClusterCM[clusterId] = 1;
  }

  for (const auto& [virtualClusterId, _] : virtualClusterCM) {
    int newVirtualClusterId = virtualClusterId;
    int loopFlag = 1;
    while (loopFlag == 1) {
      if ((oldClusterCM[newVirtualClusterId] == 1) || (oldVirtualClusterCM[newVirtualClusterId] == 1))
        newVirtualClusterId++;
      else {
        loopFlag = 0;
        oldVirtualClusterCM[newVirtualClusterId] = 1;
      }
    }
    oldToNewVirtualClusterIds[virtualClusterId] = newVirtualClusterId;
  }

  for (const auto& [virtualClusterId, newVirtualClusterId] : oldToNewVirtualClusterIds) {
    virtualClusterCM[newVirtualClusterId] = virtualClusterCM[virtualClusterId];
    virtualClusterCM.erase(virtualClusterId);

    // initialization for a later stage
    virtualToRealClusterId[newVirtualClusterId] = 0;
  }

  oldToNewVirtualClusterIds.clear();
  oldVirtualClusterCM.clear();
  oldClusterCM.clear();

  /* --------------------------------------------------------------------------
     associate the real clusters with virtual clusters and thus find the
     virtual clusters which dont have a real cluster of their own.
     -------------------------------------------------------------------------- */

  for (const auto& [clusterId, _] : clusterIdToCellId) {
    std::map<int, double> weightedDistanceV;

    double CM1[2] = {clusterCM[clusterId].getX(), clusterCM[clusterId].getY()};
    for (const auto& [virtualClusterId, virtualCluster] : virtualClusterCM) {
      const double distanceCM = std::hypot(CM1[0] - virtualCluster.getX(), CM1[1] - virtualCluster.getY());
      weightedDistanceV[virtualClusterId] = (distanceCM > 0) ? 1. / distanceCM : 1e10;
    }

    // decide which virtualCluster to associate with the real cluster
    auto closestCluster = std::ranges::max_element(weightedDistanceV, {}, &std::pair<const int, double>::second);

#if _VIRTUALCLUSTER_BUILD_DEBUG == 1
    cout << "cluster " << clusterId << " is associated with virtualCluster " << closestCluster->first << endl;
#endif

    virtualToRealClusterId[closestCluster->first] = 1;
  }

  /* --------------------------------------------------------------------------
     go over all of the hits, and change the cluster id of the hits inside
     the empty virtual clusters to that of the virtual cluster
     -------------------------------------------------------------------------- */

  for (const auto& [cellIdHit, thisHit] : calHitsCellId) {
    std::map<int, double> weightedDistanceV;

    // compute the distance of the cal hit from the virtual cluster CMs
    for (const auto& [virtualClusterId, virtualCluster] : virtualClusterCM) {
      if (virtualToRealClusterId[virtualClusterId] == 1)
        continue;

      const double distanceCM = std::hypot(thisHit->getPosition()[0] - virtualCluster.getX(),
                                           thisHit->getPosition()[1] - virtualCluster.getY());
      if (distanceCM <= virtualCluster.getZ()) {
        weightedDistanceV[virtualClusterId] = (distanceCM > 0) ? 1. / distanceCM : 1e10;
      }
    }

    if (weightedDistanceV.empty())
      continue;
    // decide to which virtualCluster to merge the hit according to a proper weight
    auto closestCluster = std::ranges::max_element(weightedDistanceV, {}, &std::pair<const int, double>::second);

#if _VIRTUALCLUSTER_BUILD_DEBUG == 1
    cout << " - num possible -  " << weightedDistanceV.size() << " \t id of chosen - (" << closestCluster->first
         << ") \t virtualClusterId/hit(x,y): "
         << " \t " << virtualClusterCM[closestCluster->first][1] << " \t " << virtualClusterCM[closestCluster->first][2]
         << endl;
#endif

    // add the hit to the chosen cluster
    clusterIdToCellId[closestCluster->first].push_back(cellIdHit);
    cellIdToClusterId[cellIdHit] = closestCluster->first;
  }

  /* --------------------------------------------------------------------------
     initialize and rewrite the clusterIdToCellId[] std::vectors for all clusters
     -------------------------------------------------------------------------- */

  for (auto& [clusterId, cellIds] : clusterIdToCellId) {
    cellIds.clear();
  }

  for (const auto& [cellIdHit, _] : calHitsCellId) {
    const int clusterId = cellIdToClusterId[cellIdHit];
    clusterIdToCellId[clusterId].push_back(cellIdHit);
  }

  /* --------------------------------------------------------------------------
     compute total energy and center of mass of each new cluster
     -------------------------------------------------------------------------- */
  clusterCM.clear();
  for (const auto& [clusterId, cellIdV] : clusterIdToCellId) {
    if (cellIdV.empty())
      continue;
    clusterCM[clusterId] = calculateEngyPosCM(cellIdV, calHitsCellId, m_methodCM);
  }

  return 1;
}

int LumiCalClustererClass::buildSuperClusters(MapIntCalHit& calHitsCellIdGlobal,
                                              VMapIntCalHit const& calHitsCellIdLayer,
                                              VMapIntVInt const& clusterIdToCellId, VMapIntLCCluster const& clusterCM,
                                              VMapIntVirtualCluster const& virtualClusterCM,
                                              MapIntInt& cellIdToSuperClusterId, MapIntVInt& superClusterIdToCellId,
                                              MapIntLCCluster& superClusterCM) {
/* --------------------------------------------------------------------------
     merge all layer-clusters into global superClusters according to distance
     from the virtual clusters' CM
     -------------------------------------------------------------------------- */
#if _VIRTUALCLUSTER_BUILD_DEBUG == 1
  cout << coutBlue << "Assign each cluster to a virtual cluster:" << coutDefault << endl;
#endif

  for (size_t layerNow = 0; layerNow < static_cast<size_t>(m_maxLayerToAnalyse); ++layerNow) {
    std::vector<int> reClusterHits;
    std::map<int, LCCluster> const& clusterInLayer = clusterCM.at(layerNow);

    for (std::map<int, LCCluster>::const_iterator clIter = clusterInLayer.begin(); clIter != clusterInLayer.end();
         ++clIter) {

      LCCluster const& thisCluster = clIter->second;
      std::map<int, double> weightedDistanceV;
      double CM1[2] = {thisCluster.getX(), thisCluster.getY()};

#if _VIRTUALCLUSTER_BUILD_DEBUG == 1
      const int clusterId = (int)(clIter->first);
      cout << " - layer " << layerNow << " \t real cluster id,x,y,engy : " << clusterId << " \t " << CM1[0] << " \t "
           << CM1[1] << " \t " << clusterCM[layerNow][clusterId][0] << endl;
#endif

      for (std::map<int, VirtualCluster>::const_iterator virtualClusterCMIterator =
               virtualClusterCM.at(layerNow).begin();
           virtualClusterCMIterator != virtualClusterCM.at(layerNow).end(); ++virtualClusterCMIterator) {
        const double CM2[2] = {virtualClusterCMIterator->second.getX(), virtualClusterCMIterator->second.getY()};
        const double distanceCM = std::hypot(CM1[0] - CM2[0], CM1[1] - CM2[1]);

#if _VIRTUALCLUSTER_BUILD_DEBUG == 1
        cout << "\t\t virtual cluster id,x,y : " << virtualClusterCMIterator->first << " \t " << CM2[0] << " \t "
             << CM2[1] << " \t distanceCM : " << distanceCM << endl;
#endif

        weightedDistanceV[virtualClusterCMIterator->first] = (distanceCM > 0) ? 1. / distanceCM : 1e10;
      }

      auto closestCluster = std::ranges::max_element(weightedDistanceV, {}, &std::pair<const int, double>::second);

      // make sure that all clusters have been added to superclusters
      if (closestCluster == weightedDistanceV.end())
        return 0;

      /* --------------------------------------------------------------------------
         in the case where there are more than 1 virtualClusters and the real
         cluster is not close to any virtualCM, than instead of assigning the whole
         real cluster to a virtualCluster, the real cluster will be dismanteled
         and each hit will be assigned to a virtualCluster separatly
         -------------------------------------------------------------------------- */
      const double distanceCM = 1. / closestCluster->second;
      if (distanceCM > m_moliereRadius && virtualClusterCM[layerNow].size() > 1) {
        reClusterHits.push_back(clIter->first); // ID of the Cluster

#if _VIRTUALCLUSTER_BUILD_DEBUG == 1
        cout << "\t\t no virtual cluster is chosen... " << coutDefault << endl << endl;
        ;
#endif
        continue;
      }

#if _VIRTUALCLUSTER_BUILD_DEBUG == 1
      cout << "\t\t chosen virtual cluster " << closestCluster->first << coutDefault << endl << endl;
      ;
#endif

      // go over all hits in the cluster and add to the chosen superCluster
      std::vector<int> const& cellIdsInThisCluster = clusterIdToCellId.at(layerNow).at(clIter->first);
      for (size_t j = 0; j < cellIdsInThisCluster.size(); ++j) {
        // add hit from clusterIdToCellId with clusterId to one with maxWDClusterId
        const int cellIdHit = cellIdsInThisCluster[j];
        superClusterIdToCellId[closestCluster->first].push_back(cellIdHit);
        cellIdToSuperClusterId[cellIdHit] = closestCluster->first;
      }
    } // over all clusters in this layer

    /* --------------------------------------------------------------------------
       real clusters far away from virtualCMs are dismanteled and each hit is
       assigned to a virtualCluster separatly according to distance
       -------------------------------------------------------------------------- */
    for (const auto& clusterId : reClusterHits) {

#if _VIRTUALCLUSTER_BUILD_DEBUG == 1
      cout << coutBlue << "\t - dismantel cluster " << clusterId << " and assign each hit to a virtualCluster "
           << coutDefault << endl;
#endif
      std::vector<int> const& cellIdsInThisCluster = clusterIdToCellId.at(layerNow).at(clusterId);

      for (const auto& cellIdHit : cellIdsInThisCluster) {
        // add hit from clusterIdToCellId with clusterId to one with maxWDClusterId

        std::map<int, double> weightedDistanceV;

        // position of the cal hit
        const auto& thisHit = calHitsCellIdLayer.at(layerNow).at(cellIdHit);
        double CM1[2] = {thisHit->getPosition()[0], thisHit->getPosition()[1]};

#if _VIRTUALCLUSTER_BUILD_DEBUG == 1
        cout << "\t cellId " << cellIdHit << " \t x,y : " << CM1[0] << " \t " << CM1[1] << endl;
#endif

        for (const auto& [virtualClusterId, virtualCluster] : virtualClusterCM[layerNow]) {
          const double CM2[2] = {virtualCluster.getX(), virtualCluster.getY()};
          const double distanceCM = std::hypot(CM1[0] - CM2[0], CM1[1] - CM2[1]);

#if _VIRTUALCLUSTER_BUILD_DEBUG == 1
          cout << "\t\t virtual cluster id,x,y : " << virtualClusterId << virtualCluster
               << " \t distanceCM : " << distanceCM << endl;
#endif

          weightedDistanceV[virtualClusterId] = (distanceCM > 0) ? 1. / distanceCM : 1e10;
        }

        // decide to which superCluster to merge the cluster according to a proper weight
        auto closestCluster = std::ranges::max_element(weightedDistanceV, {}, &std::pair<const int, double>::second);

        // make sure that all clusters have been added to superclusters
        assert(closestCluster != weightedDistanceV.end());

#if _VIRTUALCLUSTER_BUILD_DEBUG == 1
        cout << "\t\t chosen virtual cluster " << closestCluster->first << coutDefault << endl << endl;
        ;
#endif

        superClusterIdToCellId[closestCluster->first].push_back(cellIdHit);
        cellIdToSuperClusterId[cellIdHit] = closestCluster->first;
      }
    }

    // write all the layer cal hits from calHitsCellIdLayer to the global calHitsCellId map
    MapIntCalHit const& layerHits = calHitsCellIdLayer.at(layerNow);
    for (const auto& [cellIdHit, hitPtr] : layerHits) {
      calHitsCellIdGlobal[cellIdHit] = hitPtr;
    }
  }

  /* --------------------------------------------------------------------------
     compute total energy and center of mass of each new superCluster
     -------------------------------------------------------------------------- */
  superClusterCM.clear();

  for (const auto& [superClusterId, cellIdV] : superClusterIdToCellId) {
    // calculate/update the energy/position of the CM
    superClusterCM[superClusterId] = calculateEngyPosCM(cellIdV, calHitsCellIdGlobal, m_methodCM);
  }

#if _CLUSTER_BUILD_DEBUG == 1
  m_alg->debug() << " - superClusters:" << endmsg;

  for (MapIntVInt::const_iterator superClusterIdToCellIdIterator = superClusterIdToCellId.begin();
       superClusterIdToCellIdIterator != superClusterIdToCellId.end(); ++superClusterIdToCellIdIterator) {
    const int superClusterId = (int)(*superClusterIdToCellIdIterator).first; // Id of cluster

    m_alg->debug() << "\t Id " << superClusterId << "  \t  energy " << superClusterCM[superClusterId].getE()
                   << "     \t pos(x,y) =  ( " << superClusterCM[superClusterId].getX() << " , "
                   << superClusterCM[superClusterId].getY() << " )"
                   << "     \t pos(theta,phi) =  ( " << superClusterCM[superClusterId].getTheta() << " , "
                   << superClusterCM[superClusterId].getPhi() << " )" << endmsg;
  }
#endif

  return 1;
}

int LumiCalClustererClass::engyInMoliereCorrections(MapIntCalHit const& calHitsCellIdGlobal,
                                                    MapIntVCalHit const& calHits,
                                                    VMapIntCalHit const&, //   calHitsCellIdLayer,
                                                    VMapIntVInt& clusterIdToCellId, VMapIntLCCluster& clusterCM,
                                                    VMapIntInt& cellIdToClusterId, MapIntInt& cellIdToSuperClusterId,
                                                    MapIntVInt& superClusterIdToCellId, MapIntLCCluster& superClusterCM,
                                                    double middleEnergyHitBound, int /* detectorArm */) {
  // ???????? DECIDE/FIX - incorparate the parameter given here better in the code ????????
  int engyHitBoundMultiply = 1;
  double molRadPercentage = 0.4;
  double engyPercentInMolFrac = 0.8;
  double minEngyPercentInMol = 0.5;
  double midEngyPercentInMol = 0.8;
  double baseEngyPercentInMol = 0.9;

  // general variables
  MapIntCalHit calHitsCellIdProjection, calHitsCellIdProjectionFull;

  int rejectFlag;
  double superClusterMolRatio = 0., superClusterMolRatio_Tmp = 0., projectionClusterMolRatio = 0.;

  std::map<int, double> superClusterEngyInMoliere;
  std::map<int, ProjectionInfo> calHitsProjection, calHitsProjectionFull;
  std::map<int, int> projectionFlag;
  std::vector<int> initialClusterControlVar(4), superClusterRejected, superClusterAccepted;

  m_alg->debug() << "Before Moliere corrections:" << endmsg;

  double totEngyArmAboveMin = 0.;
  double totEngyInAllMol = 0.;

  /* --------------------------------------------------------------------------
     find the percentage of energy for each cluster within m_moliereRadius
     -------------------------------------------------------------------------- */
  for (const auto& [superClusterId, superCluster] : superClusterCM) {
    superClusterEngyInMoliere[superClusterId] =
        getEngyInMoliereFraction(calHitsCellIdGlobal, superClusterIdToCellId[superClusterId], superCluster, 1.);

    totEngyInAllMol += superClusterEngyInMoliere[superClusterId];
    totEngyArmAboveMin += superCluster.getE();

    m_alg->debug() << "\t Id " << superClusterId << "  \t energy " << superCluster.getE() << "     \t pos(x,y) =  ( "
                   << superCluster.getX() << " , " << superCluster.getY() << " )"
                   << "     \t pos(theta,phi) =  ( " << superCluster.getTheta() << " , " << superCluster.getPhi()
                   << " )" << endmsg
                   << "\t\t engy in m_moliereRadius  \t=   " << superClusterEngyInMoliere[superClusterId]
                   << " \t -> totEngy in Moliere = \t "
                   << superClusterEngyInMoliere[superClusterId] / superCluster.getE() << " %" << endmsg;
  }

  superClusterMolRatio = totEngyInAllMol / totEngyArmAboveMin;

  m_alg->debug() << "\t (-) Tot engy in arm = " << totEngyArmAboveMin
                 << " , tot engy in all super clusters in MolRad = " << totEngyInAllMol << "  -> their ratio = \t "
                 << superClusterMolRatio << endmsg;

  /* --------------------------------------------------------------------------
     set a global variable that warns about an unstable clustering, in
     the case that not enough energy is found within m_moliereRadius around
     the superClusters:
     (A). reject if below baseEngyPercentInMol of the energy is in moliere
     for all clusters combined
     (B). reject if below midEngyPercentInMol of the energy is in moliere
     for each cluster seperatly
     (*). also store the id's of clusters with very below minEngyPercentInMol
     energy in moliere, for force merging later, in case the profile
     clusters arent accepted.
     -------------------------------------------------------------------------- */
  rejectFlag = (superClusterMolRatio < baseEngyPercentInMol) ? 1 : 0;

  for (const auto& [superClusterId, superCluster] : superClusterCM) {
    const double engyPercentInMol = superClusterEngyInMoliere[superClusterId] / superCluster.getE();
    if (engyPercentInMol < midEngyPercentInMol)
      rejectFlag = 1;
    if (engyPercentInMol < minEngyPercentInMol)
      superClusterRejected.push_back(superClusterId);
    else
      superClusterAccepted.push_back(superClusterId);
  }

  /* --------------------------------------------------------------------------
     in the case that not enough energy is found within m_moliereRadius around
     each one of the superClusters, try to build projectionClusters.
     a projection layer is created out of all the hits of the real layers.
     clusters are built from high energy hits and the energy distribution is
     measured around the projectionCluster CMs.
     -------------------------------------------------------------------------- */
  if (rejectFlag == 1) {

    /* --------------------------------------------------------------------------
       sum up the energy for each Phi/R cell for all Z layers
       -------------------------------------------------------------------------- */
    for (const auto& [layerId, hitsVec] : calHits) {
      for (const auto& thisCalHit : hitsVec) {
        const int cellIdHit = thisCalHit->getCellID0();

        double cellEngy = thisCalHit->getEnergy();
        /// APS: This encoding needs to be fixed, now using
        // get the new cellId for the projection
        int cellIdHitZ = m_maxLayerToAnalyse + 1;
        int cellIdHitPhi = GlobalMethodsClass::cellIdZPR(cellIdHit, GlobalMethodsClass::COP);
        int cellIdHitR = GlobalMethodsClass::cellIdZPR(cellIdHit, GlobalMethodsClass::COR);
        int cellIdHitArm = GlobalMethodsClass::cellIdZPR(cellIdHit, GlobalMethodsClass::COA);
        int cellIdProjection = GlobalMethodsClass::cellIdZPR(cellIdHitZ, cellIdHitPhi, cellIdHitR, cellIdHitArm);

        // the original hit's layer number is stored in (the previously unused) CellID1
        cellIdHitZ = (cellIdHit >> 0) & (int)((1 << 10) - 1);
        // only high energy projection hits will be considered in the first stage
        if (cellEngy > middleEnergyHitBound * engyHitBoundMultiply) {
          ProjectionInfo& thisProjection = calHitsProjection[cellIdProjection];
          if (thisProjection.newObject) {
            thisProjection = ProjectionInfo(thisCalHit, cellIdHitZ);
          } else {
            thisProjection.addHit(thisCalHit);
          }
        }
        // store all hits above m_hitMinEnergy energy
        if (cellEngy > m_hitMinEnergy) {
          ProjectionInfo& thisProjection = calHitsProjectionFull[cellIdProjection];
          if (thisProjection.newObject) {
            thisProjection = ProjectionInfo(thisCalHit, cellIdHitZ);
          } else {
            thisProjection.addHit(thisCalHit);
          }
        }
      }
    }

    /* --------------------------------------------------------------------------
       input the results into new cal hit objects
       -------------------------------------------------------------------------- */
    for (const auto& [cellIdProjection, projection] : calHitsProjection) {
      auto calHitNew = std::make_shared<LumiCalHit>(cellIdProjection, projection);

      calHitsCellIdProjection[cellIdProjection] = std::move(calHitNew);
    }

    /* --------------------------------------------------------------------------
       build clusters out of the projection hits
       -------------------------------------------------------------------------- */
    // set the control std::vector for the initialClusterBuild clustering options
    initialClusterControlVar[0] = 1; // mergeOneHitClusters
    initialClusterControlVar[1] = 1; // mergeSmallToLargeClusters
    initialClusterControlVar[2] = 1; // mergeLargeToSmallClusters
    initialClusterControlVar[3] = 1; // forceMergeSmallToLargeClusters

    initialClusterBuild(calHitsCellIdProjection, cellIdToClusterId[m_maxLayerToAnalyse],
                        clusterIdToCellId[m_maxLayerToAnalyse], clusterCM[m_maxLayerToAnalyse],
                        initialClusterControlVar);

    /* --------------------------------------------------------------------------
       find the percentage of energy for each cluster within m_moliereRadius
       -------------------------------------------------------------------------- */
    m_alg->debug() << "Projection clusters (initial clustering with all hits):" << endmsg;

    int engyInMoliereFlag = 0;
    for (const auto& [clusterId, cluster] : clusterCM[m_maxLayerToAnalyse]) {
      const double thisProjectionClusterEngyInMoliere =
          getEngyInMoliereFraction(calHitsCellIdProjection, clusterIdToCellId[m_maxLayerToAnalyse][clusterId],
                                   cluster, molRadPercentage);

      const double engyPercentInMol = thisProjectionClusterEngyInMoliere / cluster.getE();

      if (engyPercentInMol < engyPercentInMolFrac)
        engyInMoliereFlag = 1;

      m_alg->debug() << "\tprojection " << clusterId << "  at (x,y) = ("
                     << clusterCM[m_maxLayerToAnalyse][clusterId].getX() << " , "
                     << clusterCM[m_maxLayerToAnalyse][clusterId].getY() << ")   \t engy in (" << molRadPercentage
                     << " * m_moliereRadius)  =   " << thisProjectionClusterEngyInMoliere << " \t-> % totEngy = \t "
                     << engyPercentInMol << endmsg;
    }

    /* --------------------------------------------------------------------------
       optimise the projection clusters (if needed).
       stop when there is a percentage of engyPercentInMolFrac total energy within
       a percentage molRadPercentage od moliere radius (each iteration the total
       energy decreases due to a higher cut on minHitEnergy).
       -------------------------------------------------------------------------- */
    if (engyInMoliereFlag == 1) {
      m_alg->debug() << "Optimizing ...  -  ignore hits with energy below    ";
    }

    while (engyInMoliereFlag == 1) {

      engyInMoliereFlag = 0;
      engyHitBoundMultiply += 30;

      m_alg->debug() << "  ->  " << middleEnergyHitBound * engyHitBoundMultiply;

      // remove hits that are of low energy
      std::vector<int> idsToErase;
      for (const auto& [cellIdProjection, hitPtr] : calHitsCellIdProjection) {
        double engyHit = (double)hitPtr->getEnergy();
        if (engyHit < middleEnergyHitBound * engyHitBoundMultiply)
          idsToErase.push_back(cellIdProjection);
      }

      int numIdsToErase = idsToErase.size();
      int numHitsRemaining = (int)calHitsCellIdProjection.size() - numIdsToErase;
      if (numHitsRemaining < 5) {
        m_alg->debug() << "  -- optimization of the projection clusters has failed ... --" << endmsg;
        break;
      }

      for (int hitNow = 0; hitNow < numIdsToErase; hitNow++) {
        int idsToEraseNow = idsToErase[hitNow];
        // erase entry from map
        calHitsCellIdProjection.erase(idsToEraseNow);
      }

      /* --------------------------------------------------------------------------
         build clusters out of the projection hits
         -------------------------------------------------------------------------- */
      // clean up the clustering results from the previous run
      cellIdToClusterId[m_maxLayerToAnalyse].clear();
      clusterIdToCellId[m_maxLayerToAnalyse].clear();
      clusterCM[m_maxLayerToAnalyse].clear();

      // set the control std::vector for the initialClusterBuild clustering options
      initialClusterControlVar[0] = 1; // mergeOneHitClusters
      initialClusterControlVar[1] = 1; // mergeSmallToLargeClusters
      initialClusterControlVar[2] = 1; // mergeLargeToSmallClusters
      initialClusterControlVar[3] = 1; // forceMergeSmallToLargeClusters

      initialClusterBuild(calHitsCellIdProjection, cellIdToClusterId[m_maxLayerToAnalyse],
                          clusterIdToCellId[m_maxLayerToAnalyse], clusterCM[m_maxLayerToAnalyse],
                          initialClusterControlVar);

      /* --------------------------------------------------------------------------
         find the percentage of energy for each cluster within m_moliereRadius
         -------------------------------------------------------------------------- */
      for (const auto& [clusterIdHit, cluster] : clusterCM[m_maxLayerToAnalyse]) {
        const double thisProjectionClusterEngyInMoliere =
            getEngyInMoliereFraction(calHitsCellIdProjection, clusterIdToCellId[m_maxLayerToAnalyse][clusterIdHit],
                                     cluster, molRadPercentage);

        double engyPercentInMol =
            thisProjectionClusterEngyInMoliere / clusterCM[m_maxLayerToAnalyse][clusterIdHit].getE();

        if (engyPercentInMol < engyPercentInMolFrac)
          engyInMoliereFlag = 1;
      }

    } // engyInMoliereFlag == 1

    /* --------------------------------------------------------------------------
       since the optimization of the projection clusters involved deleting of low
       energy cal hits, these need to be re-registered so as to calculate the
       total energy around the CM within m_moliereRadius.
       -------------------------------------------------------------------------- */
    for (const auto& [cellIdProjection, projection] : calHitsProjectionFull) {
      auto calHitNew = std::make_shared<LumiCalHit>(cellIdProjection, projection);

      calHitsCellIdProjectionFull[cellIdProjection] = std::move(calHitNew);
      // a flag map for avoiding double counting of clustered cells
      projectionFlag[cellIdProjection] = 0;
    }

    m_alg->debug() << "Projection clusters (all hits):" << endmsg;

    /* --------------------------------------------------------------------------
       find the percentage of energy for each cluster within m_moliereRadius
       -------------------------------------------------------------------------- */
    totEngyInAllMol = 0.;

    for (const auto& [clusterIdHit, cluster] : clusterCM[m_maxLayerToAnalyse]) {
      const double thisProjectionClusterEngyInMoliere =
          getEngyInMoliereFraction(calHitsCellIdProjectionFull, clusterIdToCellId[m_maxLayerToAnalyse][clusterIdHit],
                                   cluster, 1., projectionFlag);

      totEngyInAllMol += thisProjectionClusterEngyInMoliere;

      m_alg->debug() << "\tfull-Projection " << clusterIdHit << cluster
                     << ")   \t engy in m_moliereRadius  \t=   " << thisProjectionClusterEngyInMoliere
                     << " (possibly under-counted...)" << endmsg;
    }

    projectionClusterMolRatio = totEngyInAllMol / totEngyArmAboveMin;

    m_alg->debug() << "\t (-) Tot engy in arm = " << totEngyArmAboveMin
                   << " , tot engy in all projection clusters  in MolRad = " << totEngyInAllMol
                   << "  -> their ratio = \t " << projectionClusterMolRatio << endmsg;
  } // if rejectFlag == 1

  /* --------------------------------------------------------------------------
     if the projectionClusters (after possible optimization) have enough energy
     within m_moliereRadius around their CMs, and the superClusters don't,
     then dump the superClusters and re-cluster around the projectionClusters.
     this time a very rudementary clustering algorithm is used.
     -------------------------------------------------------------------------- */
  if (superClusterMolRatio < projectionClusterMolRatio) {
    MapIntVInt superClusterIdToCellId_Tmp;
    std::map<int, int> cellIdToSuperClusterId_Tmp;

    // create new clusters around the projection CMs
    for (const auto& [cellIdHit, thisHit] : calHitsCellIdGlobal) {
      std::map<int, double> weightedDistanceV;

      for (const auto& [clusterId, cluster] : clusterCM[m_maxLayerToAnalyse]) {
        const double distanceCM =
            std::hypot(thisHit->getPosition()[0] - cluster.getX(), thisHit->getPosition()[1] - cluster.getY());
        weightedDistanceV[clusterId] = (distanceCM > 0) ? 1. / distanceCM : 1e10;
      }

      // decide to which superCluster to merge the cluster according to a proper weight
      auto closestCluster = std::ranges::max_element(weightedDistanceV, {}, &std::pair<const int, double>::second);
      // make sure that all hits have have been added to clusters
      assert(closestCluster != weightedDistanceV.end());
      // add the hit to the chosen superCluster
      superClusterIdToCellId_Tmp[closestCluster->first].push_back(cellIdHit);
      cellIdToSuperClusterId_Tmp[cellIdHit] = closestCluster->first;
    }

    /* --------------------------------------------------------------------------
       compute total energy and center of mass of each new superCluster
       -------------------------------------------------------------------------- */
    MapIntLCCluster superClusterCM_Tmp;

    for (MapIntVInt::const_iterator superClusterIdToCellIdIterator = superClusterIdToCellId_Tmp.begin();
         superClusterIdToCellIdIterator != superClusterIdToCellId_Tmp.end(); ++superClusterIdToCellIdIterator) {
      // calculate/update the energy/position of the CM
      superClusterCM_Tmp[superClusterIdToCellIdIterator->first] =
          calculateEngyPosCM(superClusterIdToCellIdIterator->second, calHitsCellIdGlobal, m_methodCM);
    }

    /* --------------------------------------------------------------------------
       find the percentage of energy for each new superCluster within m_moliereRadius
       -------------------------------------------------------------------------- */
    totEngyArmAboveMin = 0.;
    totEngyInAllMol = 0.;

    m_alg->debug() << "New superClusters:" << endmsg;

    std::map<int, double> superClusterEngyInMoliere_Tmp;

    for (const auto& [superClusterId, superCluster] : superClusterCM_Tmp) {
      superClusterEngyInMoliere_Tmp[superClusterId] =
          getEngyInMoliereFraction(calHitsCellIdGlobal, superClusterIdToCellId_Tmp[superClusterId], superCluster, 1.);

      totEngyInAllMol += superClusterEngyInMoliere_Tmp[superClusterId];
      totEngyArmAboveMin += superCluster.getE();

      m_alg->debug() << "superCluster " << superClusterId << " \tat (x,y) = (" << superCluster.getX() << " , "
                     << superCluster.getY()
                     << ")   \t engy in m_moliereRadius  \t=   " << superClusterEngyInMoliere_Tmp[superClusterId]
                     << " \t-> % totEngy = \t " << superClusterEngyInMoliere_Tmp[superClusterId] / superCluster.getE()
                     << endmsg;
    }

    superClusterMolRatio_Tmp = totEngyInAllMol / totEngyArmAboveMin;

    m_alg->debug() << "\t (-) Tot engy in arm = " << totEngyArmAboveMin
                   << " , tot engy in all super clusters in MolRad = " << totEngyInAllMol << "  -> their ratio = \t "
                   << superClusterMolRatio_Tmp << endmsg;

    if (superClusterMolRatio < superClusterMolRatio_Tmp) {
      superClusterIdToCellId = superClusterIdToCellId_Tmp;
      cellIdToSuperClusterId = cellIdToSuperClusterId_Tmp;
      superClusterCM = superClusterCM_Tmp;
      superClusterEngyInMoliere = superClusterEngyInMoliere_Tmp;

      superClusterMolRatio = superClusterMolRatio_Tmp;

      rejectFlag = 0;
      m_alg->debug() << "\t -- ACCEPTED new superCluster(s) -- " << endmsg;

    } else {
      rejectFlag = 1;
      m_alg->debug() << "\t -- REJECTED new superCluster(s) -- " << endmsg;
    }

  } else {
    if (projectionClusterMolRatio > 0) {
      m_alg->debug() << " \t -- the projection is not an improvement on the SuperCluster(s) -- " << endmsg;
    }
    rejectFlag = 1;
  }

  /* --------------------------------------------------------------------------
     in the case that the projection cluster was not accepted, and some of the
     superClusters have a low energy content in their moliere radius, the
     merge these superClusters with the rest
     -------------------------------------------------------------------------- */
  if (rejectFlag == 1 && superClusterRejected.size() > 0 && superClusterAccepted.size() > 0) {

    for (const auto superClusterIdBad : superClusterRejected) {

      const std::vector<int>& cellIds = superClusterIdToCellId[superClusterIdBad];
      // Find the clusters closest to the previously rejected Hits
      for (const auto cellIdHit : cellIds) {
        std::map<int, double> weightedDistanceV;
        const auto& thisHit = calHitsCellIdGlobal.at(cellIdHit);

        for (const auto superClusterIdGood : superClusterAccepted) {
          const double distanceCM = std::hypot(thisHit->getPosition()[0] - superClusterCM[superClusterIdGood].getX(),
                                               thisHit->getPosition()[1] - superClusterCM[superClusterIdGood].getY());
          weightedDistanceV[superClusterIdGood] = (distanceCM > 0) ? 1. / distanceCM : 1e10;
        }

        // decide to which superCluster to merge the hit according to a proper weight
        auto closestCluster = std::ranges::max_element(weightedDistanceV, {}, &std::pair<const int, double>::second);
        // make sure that all hits have have been added to clusters
        assert(closestCluster != weightedDistanceV.end());

        // add the hit to the chosen superCluster
        superClusterIdToCellId[closestCluster->first].push_back(cellIdHit);
        cellIdToSuperClusterId[cellIdHit] = closestCluster->first;

        //  update the totalEnergy counter and CM position of the superCluster
        updateEngyPosCM(calHitsCellIdGlobal.at(cellIdHit), superClusterCM[closestCluster->first]);

      } // for elements in SuperCluster

      // cleanUp
      superClusterIdToCellId.erase(superClusterIdBad);
      superClusterCM.erase(superClusterIdBad);
    }

    // find the percentage of energy for each cluster within m_moliereRadius
    totEngyArmAboveMin = 0.;
    totEngyInAllMol = 0.;

#if _MOL_RAD_CORRECT_DEBUG == 1
    cout << endl << coutBlue << "Merged superCluster(s):" << coutDefault << endl;
#endif

    for (auto& [superClusterId, superCluster] : superClusterCM) {
      superClusterEngyInMoliere[superClusterId] =
          getEngyInMoliereFraction(calHitsCellIdGlobal, superClusterIdToCellId[superClusterId], superCluster, 1.);

      totEngyInAllMol += superClusterEngyInMoliere[superClusterId];
      totEngyArmAboveMin += superCluster.getE();

#if _MOL_RAD_CORRECT_DEBUG == 1
      double engyPercentInMol = superClusterEngyInMoliere[superClusterId] / superCluster.getE();
      cout << "superCluster " << superClusterId << " \tat (x,y) = (" << superCluster.getX() << " , "
           << superCluster.getY()
           << ")   \t engy in m_moliereRadius  \t=   " << superClusterEngyInMoliere[superClusterId]
           << " \t-> % totEngy = \t " << engyPercentInMol << coutDefault << endl;
#endif
    }

    superClusterMolRatio = totEngyInAllMol / totEngyArmAboveMin;

    m_alg->debug() << "\t (-) Tot engy in arm = " << totEngyArmAboveMin
                   << " , tot engy in all super clusters in MolRad = " << totEngyInAllMol << "  -> their ratio = \t "
                   << superClusterMolRatio << endmsg;
  }

  // cleanUp dynamically allocated memory
  calHitsCellIdProjection.clear();
  calHitsCellIdProjectionFull.clear();

  /* --------------------------------------------------------------------------
     re-compute total energy and center of mass of each superCluster (just in case...)
     -------------------------------------------------------------------------- */
  for (const auto& [superClusterId, cellIds] : superClusterIdToCellId) {
    superClusterCM[superClusterId] = calculateEngyPosCM(cellIds, calHitsCellIdGlobal, m_methodCM);
  }

  return 1;
}
