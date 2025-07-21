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

#include <cassert>
#include <cmath>
#include <iomanip>
#include <map>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

/* --------------------------------------------------------------------------
   compute the nearest neigbour cellId
   - in case of the next neighbor being outside the detector, return 0
   - get neighbors in Phi +/- 1 (neighborIndex < 2), and neighbors in R +/- n ( n = neighborIndex/2 )
   -------------------------------------------------------------------------- */
int LumiCalClustererClass::getNeighborId(int cellId, const int neighborIndex) {

  int cellZ, cellPhi, cellR, arm;
  // compute Z,Phi,R coordinates according to the cellId
  cellIdZPR(cellId, cellZ, cellPhi, cellR, arm);

  // change iRho cell index  according to the neighborIndex
  if (neighborIndex >= 2) {
    if (neighborIndex % 2)
      cellR -= neighborIndex / 2;
    else
      cellR += neighborIndex / 2;
  }

  // change iPhi cell index
  if (neighborIndex == 0) {
    cellPhi += 1;
  } else if (neighborIndex == 1) {
    cellPhi -= 1;
  }

  if (cellR == m_cellRMax || cellR < 0)
    return 0;

  if (cellPhi == m_cellPhiMax)
    cellPhi = 0;
  else if (cellPhi < 0)
    cellPhi = m_cellPhiMax - 1;
  // compute neighbor cellId according to the new Z,Phi,R coordinates
  cellId = cellIdZPR(cellZ, cellPhi, cellR, arm);

  return cellId;
}

/* --------------------------------------------------------------------------
   calculate weight for cluster CM according to different methods
   - provide cellEnergy to use
   -------------------------------------------------------------------------- */
double LumiCalClustererClass::posWeightTrueCluster(const CalHit& calHit, const double cellEngy,
                                                   const WeightingMethod_t method) const {
  int detectorArm = ((calHit->getPosition()[2] < 0) ? -1 : 1);
  return posWeight(cellEngy, m_totEngyArm.at(detectorArm), method, m_logWeightConst);
}

// /* --------------------------------------------------------------------------
//    calculate weight for cluster CM according to different methods
//    -------------------------------------------------------------------------- */
double LumiCalClustererClass::posWeight(const CalHit& calHit, const WeightingMethod_t method) const {
  return posWeightTrueCluster(calHit, calHit->getEnergy(), method);
}

// /* --------------------------------------------------------------------------
//    compute center of mass of each cluster
//    (3). calculate the map clusterCM from scratch
//    -------------------------------------------------------------------------- */
LCCluster LumiCalClustererClass::calculateEngyPosCM(const VInt& cellIdV, const MapIntCalHit& calHitsCellId,
                                                    WeightingMethod_t method) {

  double totEngy(0.0), xHit(0.0), yHit(0.0), zHit(0.0), thetaHit(0.0), weightSum(0.0);
  int loopFlag = 1;
  VecCalHit caloHits;
  while (loopFlag == 1) {
    caloHits.clear();
    for (VInt::const_iterator it = cellIdV.begin(); it != cellIdV.end(); ++it) {
      auto& calHit = calHitsCellId.at(*it);
      caloHits.push_back(calHit);

      const double weightHit = posWeight(calHit, method);
      weightSum += weightHit;

      const auto position = calHit->getPosition();

      xHit += position[0] * weightHit;
      yHit += position[1] * weightHit;
      zHit += position[2] * weightHit;
      totEngy += calHit->getEnergy();
    }
    if (weightSum > 0.) {
      xHit /= weightSum;
      yHit /= weightSum;
      zHit /= weightSum;
      thetaHit = atan(sqrt(xHit * xHit + yHit * yHit) / fabs(zHit));
      loopFlag = 0;

    } else {
      // initialize counters and recalculate with the Energy-weights method
      method = LumiCalClustererClass::EnergyMethod;
      totEngy = xHit = yHit = zHit = thetaHit = weightSum = 0.;
    }
  }

  return LCCluster(totEngy, xHit, yHit, zHit, weightSum, method, thetaHit, 0.0, caloHits);
}

// /* --------------------------------------------------------------------------
//    compute center of mass of each cluster
//    (3). calculate the map clusterCM from scratch
//    -------------------------------------------------------------------------- */
void LumiCalClustererClass::calculateEngyPosCM_EngyV(const VInt& cellIdV, const VDouble& cellEngyV,
                                                     const MapIntCalHit& calHitsCellId, MapIntLCCluster& clusterCM,
                                                     int clusterId, WeightingMethod_t method) {

  double totEngy(0.0), xHit(0.0), yHit(0.0), zHit(0.0), thetaHit(0.0), weightSum(0.0);
  int loopFlag = 1;
  VecCalHit caloHits;
  while (loopFlag == 1) {
    caloHits.clear();
    for (auto it = cellIdV.begin(); it != cellIdV.end(); ++it) {
      const size_t k = it - cellIdV.begin();
      const auto& calHit = calHitsCellId.at(*it);
      caloHits.push_back(calHit);
      const double weightHit = posWeightTrueCluster(calHit, cellEngyV[k], method);
      weightSum += weightHit;
      const auto position = calHit->getPosition();

      xHit += position[0] * weightHit;
      yHit += position[1] * weightHit;
      zHit += position[2] * weightHit;
      totEngy += cellEngyV[k];
    }
    if (weightSum > 0) {
      xHit /= weightSum;
      yHit /= weightSum;
      zHit /= weightSum;
      thetaHit = atan(sqrt(xHit * xHit + yHit * yHit) / fabs(zHit));
      //      thetaHit /= weightSum;
      loopFlag = 0;

    } else {
      // initialize counters and recalculate with the Energy-weights method
      method = LumiCalClustererClass::EnergyMethod;
      totEngy = xHit = yHit = zHit = thetaHit = weightSum = 0.;
    }
  }

  clusterCM[clusterId] = LCCluster(totEngy, xHit, yHit, zHit, weightSum, method, thetaHit, 0.0, caloHits);
}

// /* --------------------------------------------------------------------------
//    compute center of mass of each cluster
//    (2). update the map clusterCM with the new cal hit
//    -------------------------------------------------------------------------- */
void LumiCalClustererClass::updateEngyPosCM(const CalHit& calHit, LCCluster& clusterCM) {
  double engyHit = calHit->getEnergy();
  WeightingMethod_t method = clusterCM.getMethod() == WeightingMethod_t::LogMethod ? WeightingMethod_t::LogMethod
                                                                                   : WeightingMethod_t::EnergyMethod;

  clusterCM.addToEnergy(engyHit);

  double weightHit = posWeight(calHit, method);

  if (weightHit > 0) {
    double weightCM = clusterCM.getWeight();

    double xCM = clusterCM.getX() * weightCM;
    double yCM = clusterCM.getY() * weightCM;
    double zCM = clusterCM.getZ() * weightCM;

    double xHit = calHit->getPosition()[0] * weightHit;
    double yHit = calHit->getPosition()[1] * weightHit;
    double zHit = calHit->getPosition()[2] * weightHit;

    weightCM += weightHit;
    xCM += xHit;
    xCM /= weightCM;
    yCM += yHit;
    yCM /= weightCM;
    zCM += zHit;
    zCM /= weightCM;
    clusterCM.setPosition(xCM, yCM, zCM);
    clusterCM.setWeight(weightCM);
  }
}

// /* --------------------------------------------------------------------------
//    make sure that the CM of two merged clusters is where most of the merged
//    cluster's energy is deposited
//    -------------------------------------------------------------------------- */
int LumiCalClustererClass::checkClusterMergeCM(int clusterId1, int clusterId2, MapIntVInt const& clusterIdToCellId,
                                               MapIntCalHit const& calHitsCellId, double distanceAroundCM,
                                               double percentOfEngyAroungCM, WeightingMethod_t method) {

  // std::vector for holding the Ids of clusters
  VInt cellIdV;
  double totEngyAroundCM = 0.;

  // add to cellIdV hits from both clusters, and sum up each cluster's energy
  for (VInt::const_iterator cellIt = clusterIdToCellId.at(clusterId1).begin();
       cellIt != clusterIdToCellId.at(clusterId1).end(); ++cellIt) {
    cellIdV.push_back(*cellIt);
  }

  for (VInt::const_iterator cellIt = clusterIdToCellId.at(clusterId2).begin();
       cellIt != clusterIdToCellId.at(clusterId2).end(); ++cellIt) {
    cellIdV.push_back(*cellIt);
  }

  LCCluster engyPosCM(calculateEngyPosCM(cellIdV, calHitsCellId, method));
  double CM1[3] = {engyPosCM.getX(), engyPosCM.getY(), 0.0};
  const double engyCM = engyPosCM.getE();

  // check that the CM position has a significant amount of the energy around it
  for (auto cellIt = cellIdV.begin(); cellIt != cellIdV.end(); ++cellIt) {
    const auto& hit = calHitsCellId.at(*cellIt);
    const double distanceCM = std::hypot(CM1[0] - hit->getPosition()[0], CM1[1] - hit->getPosition()[1]);
    if (distanceCM < distanceAroundCM) {
      double engyHit = hit->getEnergy();
      totEngyAroundCM += engyHit;
    }
  }

  if (totEngyAroundCM > (engyCM * percentOfEngyAroungCM))
    return 1;
  else
    return 0;
}

// /* --------------------------------------------------------------------------
//    get the energy around a cluster CM within a distanceToScan raduis
//    -------------------------------------------------------------------------- */
double LumiCalClustererClass::getEngyInMoliereFraction(MapIntCalHit const& calHitsCellId,
                                                       VInt const&, // clusterIdToCellId,
                                                       LCCluster const& clusterCM, double moliereFraction) {

  const double distanceToScan = m_moliereRadius * moliereFraction;
  double engyAroundCM = 0.0;

  for (MapIntCalHit::const_iterator calHitsCellIdIterator = calHitsCellId.begin();
       calHitsCellIdIterator != calHitsCellId.end(); ++calHitsCellIdIterator) {
    const auto& calHit = calHitsCellIdIterator->second;
    const double distanceCM =
        std::hypot(clusterCM.getX() - calHit->getPosition()[0], clusterCM.getY() - calHit->getPosition()[1]);
    if (distanceCM < distanceToScan)
      engyAroundCM += calHit->getEnergy();
  }

  return engyAroundCM;
}

// // overloaded with different variables and functionality...
double LumiCalClustererClass::getEngyInMoliereFraction(MapIntCalHit const& calHitsCellId,
                                                       VInt const&, // clusterIdToCellId,
                                                       LCCluster const& clusterCM, double moliereFraction,
                                                       MapIntInt& flag) {

  const double distanceToScan = m_moliereRadius * moliereFraction;
  double engyAroundCM = 0.;

  for (MapIntCalHit::const_iterator calHitsCellIdIterator = calHitsCellId.begin();
       calHitsCellIdIterator != calHitsCellId.end(); ++calHitsCellIdIterator) {
    const auto& calHit = calHitsCellIdIterator->second;
    const double distanceCM =
        std::hypot(clusterCM.getX() - calHit->getPosition()[0], clusterCM.getY() - calHit->getPosition()[1]);

    const int cellIdHit = calHitsCellIdIterator->first;
    if (distanceCM < distanceToScan && flag[cellIdHit] == 0) {
      engyAroundCM += calHit->getEnergy();
      flag[cellIdHit] = 1;
    }
  }

  return engyAroundCM;
}

std::string LumiCalClustererClass::printClusters(const int armNow, const MapIntMapIntLCCluster& superClusterCM) const {
  std::stringstream p;
  for (const auto& superClusterCMIterator : superClusterCM.at(armNow)) {
    p << "  Arm:" << std::setw(4) << armNow << "  Id:" << std::setw(4) << superClusterCMIterator.first
      << superClusterCMIterator.second << std::endl;
  }
  return p.str();
}

std::string LumiCalClustererClass::printClusters(const MapIntLCCluster& superClusterCM) const {
  std::stringstream p;
  for (const auto& superClusterCMIterator : superClusterCM) {
    p << "  Id:" << std::setw(4) << superClusterCMIterator.first << superClusterCMIterator.second << std::endl;
  }
  return p.str();
}
