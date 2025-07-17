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
#include "LumiCalHit.h"
#include "VirtualCluster.h"

#include <iomanip>

LCCluster::LCCluster(const VirtualCluster& vc) : m_position{vc.getX(), vc.getY(), vc.getZ()} {}

LCCluster::LCCluster(double energy, double x, double y, double z, double weight,
                     GlobalMethodsClass::WeightingMethod_t method, double theta, double phi,
                     VecCalHit const& caloHitVector)
    : m_position{x, y, z}, m_energy(energy), m_weight(weight), m_method(method), m_theta(theta), m_phi(phi),
      m_caloHits(caloHitVector) {
  CalculatePhi();
}

void LCCluster::clear() {
  m_energy = 0.0;
  m_position[0] = 0.0;
  m_position[1] = 0.0;
  m_position[2] = 0.0;
  m_weight = 0.0;
  m_method = GlobalMethodsClass::LogMethod;
  m_theta = 0.0;
  m_phi = 0.0;
  m_caloHits.clear();
}

std::ostream& operator<<(std::ostream& o, const LCCluster& rhs) {
  o << "  Energy " << std::setw(10) << rhs.m_energy << "  Method " << std::setw(4) << rhs.m_method << "  Weight "
    << std::setw(10) << rhs.m_weight << "  N Calo Hits " << std::setw(7) << rhs.m_caloHits.size()
    << "  pos(x,y,z) =  ( " << std::setw(10) << rhs.m_position[0] << " , " << std::setw(10) << rhs.m_position[1]
    << " , " << std::setw(10) << rhs.m_position[2] << " )"
    << "  pos(theta,phi) =  ( " << std::setw(10) << rhs.m_theta << " , " << std::setw(10) << rhs.m_phi << " )";
  return o;
}

/** recalculate the position of the cluster based on the theta and phi averages
 *
 * Resolution in Theta (R) is better than in RPhi so averaging theta gives better results
 */
void LCCluster::recalculatePositionFromHits(const GlobalMethodsClass& gmc) {
  const double logConstant(gmc.m_globalParamD.at(GlobalMethodsClass::LogWeightConstant));

  // re-set new cluster energy
  m_energy = 0.0;
  for (auto const& calHit : m_caloHits) {
    m_energy += calHit->getEnergy();
  }

  const double rMin = gmc.m_globalParamD.at(GlobalMethodsClass::RMin);

  double thetaTemp(0.0), weightsTemp(0.0), xTemp(0.0), yTemp(0.0), zTemp(0.0);
  for (auto const& calHit : m_caloHits) {
    const auto pos = calHit->getPosition();
    const double rCell = sqrt(pos[0] * pos[0] + pos[1] * pos[1]);
    // cell area scales with radius, reduce weight for cells at larger radii
    const double weight =
        GlobalMethodsClass::posWeight(calHit->getEnergy(), m_energy, m_method, logConstant) * rMin / rCell;

    if (not(weight > 0))
      continue;

    const double theta = atan(rCell / pos[2]);

    thetaTemp += theta * weight;

    xTemp += pos[0] * weight;
    yTemp += pos[1] * weight;
    zTemp += pos[2] * weight;

    weightsTemp += weight;
  }

  xTemp /= weightsTemp;
  yTemp /= weightsTemp;
  zTemp /= weightsTemp;

  thetaTemp /= weightsTemp;

  m_phi = atan2(yTemp, xTemp);
  m_theta = thetaTemp;

  if (m_theta < 0)
    m_theta += M_PI;
  const double sign = zTemp < 0 ? -1.0 : 1.0;

  const double r = sqrt(xTemp * xTemp + yTemp * yTemp + zTemp * zTemp);
  const double zStart = sign * gmc.m_globalParamD.at(GlobalMethodsClass::ZStart);

  m_position[0] = r * sin(m_theta) * cos(m_phi);
  m_position[1] = r * sin(m_theta) * sin(m_phi);
  m_position[2] = r * cos(m_theta);

  // cluster position at front-face of LumiCal
  m_rzStart = tan(m_theta) * zStart;
  m_positionAtFront[0] = m_rzStart * cos(m_phi);
  m_positionAtFront[1] = m_rzStart * sin(m_phi);
  m_positionAtFront[2] = zStart;
}
