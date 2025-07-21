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
#include "LCCluster.h"
#include "LumiCalClusterer.h"

#include <cmath>
#include <map>
#include <utility>
#include <vector>

void LumiCalClustererClass::fiducialVolumeCuts(std::map<int, std::vector<int>>& superClusterIdToCellId,
                                               std::map<int, std::vector<double>>& superClusterIdToCellEngy,
                                               std::map<int, LCCluster>& superClusterCM) {
  if (!m_cutOnFiducialVolume) {
    m_alg->debug() << "Skipping fiducial volume cuts" << endmsg;
    return;
  }

  int superClusterId;
  double thetaSuperCluster;

  std::map<int, std::vector<int>>::iterator superClusterIdToCellIdIterator;

  std::vector<int> clusterIdToErase;

  /* --------------------------------------------------------------------------
     discard true/reconstructed clusters that are outside the fiducial volume
     -------------------------------------------------------------------------- */
  superClusterIdToCellIdIterator = superClusterIdToCellId.begin();
  for (size_t superClusterNow = 0; superClusterNow < superClusterIdToCellId.size();
       superClusterNow++, superClusterIdToCellIdIterator++) {
    superClusterId = (int)(*superClusterIdToCellIdIterator).first; // Id of cluster

    thetaSuperCluster = fabs(superClusterCM[superClusterId].getTheta());
    if (thetaSuperCluster < m_thetaMin || thetaSuperCluster > m_thetaMax)
      clusterIdToErase.push_back(superClusterId);
  }

  for (const auto superClusterId : clusterIdToErase) {
    superClusterIdToCellId.erase(superClusterId);
    superClusterIdToCellEngy.erase(superClusterId);
    superClusterCM.erase(superClusterId);
  }
  clusterIdToErase.clear();
}
