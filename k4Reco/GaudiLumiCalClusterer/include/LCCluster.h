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
#ifndef K4RECO_LCCLUSTER_H
#define K4RECO_LCCLUSTER_H 1

#include "Global.h"
#include "GlobalMethodsClass.h"

#include <array>
#include <cmath>
#include <memory>
#include <ostream>

class VirtualCluster;

class LCCluster {

public:
  LCCluster() = default;
  explicit LCCluster(const VirtualCluster& vc);
  LCCluster(double energy, double x, double y, double z, double weight, GlobalMethodsClass::WeightingMethod_t method,
            double theta, double phi, VecCalHit const& caloHitVector);
  // LCCluster( const LCCluster& rhs ) = delete;
  // LCCluster& operator=( const LCCluster& rhs ) = delete;
  // LCCluster( LCCluster&& rhs ) = default;
  // LCCluster& operator=( LCCluster&& rhs ) = default;

  void clear();

  inline double getX() const { return m_position[0]; }
  inline double getY() const { return m_position[1]; }
  inline double getZ() const { return m_position[2]; }
  inline double getE() const { return m_energy; }
  inline double getWeight() const { return m_weight; }
  inline GlobalMethodsClass::WeightingMethod_t getMethod() const { return m_method; }
  inline double getEnergy() const { return getE(); }
  inline double getTheta() const { return m_theta; }
  inline double getPhi() const { return m_phi; }
  inline double getRZStart() const { return m_rzStart; }

  inline const std::array<double, 3>& getPosition() const { return m_position; }
  inline const std::array<double, 3>& getPositionAtFront() const { return m_positionAtFront; }

  inline void setX(double x) {
    m_position[0] = x;
    CalculatePhi();
  }
  inline void setY(double y) {
    m_position[1] = y;
    CalculatePhi();
  }
  inline void setZ(double z) {
    m_position[2] = z;
    CalculateTheta();
  }
  inline void setPosition(double x, double y, double z) {
    m_position[0] = x;
    m_position[1] = y;
    m_position[2] = z;
    CalculatePhi();
    CalculateTheta();
  }
  inline void setWeight(double w) { m_weight = w; }
  inline void setTheta(double t) { m_theta = t; }
  inline void setPhi(double p) { m_phi = p; }

  void addToEnergy(double E) { m_energy += E; }

  friend std::ostream& operator<<(std::ostream& o, const LCCluster& rhs);

  inline VecCalHit const& getCaloHits() const { return m_caloHits; }

  /// calculate the cluster position based on the caloHits associated to the cluster
  void recalculatePositionFromHits(const GlobalMethodsClass& gmc);

private:
  void CalculatePhi();
  void CalculateTheta();
  std::array<double, 3> m_position{0.0, 0.0, 0.0};
  std::array<double, 3> m_positionAtFront{0.0, 0.0, 0.0};
  double m_energy = 0.0, m_weight = 0.0;
  GlobalMethodsClass::WeightingMethod_t m_method = GlobalMethodsClass::LogMethod;
  double m_theta = 0.0, m_phi = 0.0, m_rzStart = 0.0;
  VecCalHit m_caloHits{};
};

inline void LCCluster::CalculatePhi() { m_phi = atan2(m_position[1], m_position[0]); }
inline void LCCluster::CalculateTheta() {
  m_theta = atan(sqrt(m_position[1] * m_position[1] + m_position[0] * m_position[0]) / fabs(m_position[2]));
}

using SLCCluster = std::shared_ptr<LCCluster>;

#endif // K4RECO_LCCLUSTER_H
