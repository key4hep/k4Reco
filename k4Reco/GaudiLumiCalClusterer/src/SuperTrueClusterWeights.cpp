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
#include "SuperTrueClusterWeights.h"
#include "LCCluster.h"

#include <algorithm>
#include <cmath>

SuperTrueClusterWeights::SuperTrueClusterWeights(int superClusterIdNow, int trueClusterIdNow,
                                                 LCCluster const& superClusterCM, LCCluster const& trueClusterCM)
    :

      m_superClusterId(superClusterIdNow), m_trueClusterId(trueClusterIdNow),
      m_distance(std::hypot(superClusterCM.getX() - trueClusterCM.getX(), superClusterCM.getY() - trueClusterCM.getY())),
      m_deltaEngy(fabs(superClusterCM.getE() - trueClusterCM.getE())),
      m_minEngy(std::min(superClusterCM.getE(), trueClusterCM.getE())), m_weight(-1) {}

void SuperTrueClusterWeights::setWeight(std::string weightMethod) {

  m_weight = (weightMethod == "distance") ? m_distance : m_deltaEngy;
}

void SuperTrueClusterWeights::setWeight(std::string weightMethod, double minSeparationDistance,
                                        double minClusterEngyGeV) {

  if (weightMethod == "minEngyDistance") {

    m_weight = (m_distance > minSeparationDistance && m_minEngy > minClusterEngyGeV) ? 1. : -1.;
  }
}
