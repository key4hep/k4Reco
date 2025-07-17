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
#ifndef SuperTrueClusterWeights_hh
#define SuperTrueClusterWeights_hh 1

#include <string>

class LCCluster;

/* --------------------------------------------------------------------------
   class and sort rule for computing weights required for assignement of
   reconstructed to true clusters
   -------------------------------------------------------------------------- */
class SuperTrueClusterWeights {

public:
  SuperTrueClusterWeights(int superClusterIdNow, int trueClusterIdNow, LCCluster const& superClusterCM,
                          LCCluster const& trueClusterCM);

  double distance2D(double* pos1, double* pos2);
  void setWeight(std::string weightMethod);
  void setWeight(std::string weightMethod, double minSeparationDistance, double minClusterEngyGeV);

  int m_superClusterId, m_trueClusterId;
  double m_distance, m_deltaEngy, m_minEngy, m_weight;

};

#endif // SuperTrueClusterWeights_hh
