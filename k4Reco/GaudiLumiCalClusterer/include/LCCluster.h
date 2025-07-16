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
#ifndef LCCluster_hh
#define LCCluster_hh 1

#include "Global.h"
#include "GlobalMethodsClass.h"

#include <cmath>
#include <memory>
#include <ostream>

class VirtualCluster;

class LCCluster {

public:
  LCCluster();
  explicit LCCluster(const VirtualCluster& vc);
  LCCluster(double energy, double x, double y, double z, double weight, GlobalMethodsClass::WeightingMethod_t method,
            double theta, double phi, VecCalHit const& caloHitVector);
  // LCCluster( const LCCluster& rhs ) = delete;
  // LCCluster& operator=( const LCCluster& rhs ) = delete;
  // LCCluster( LCCluster&& rhs ) = default;
  // LCCluster& operator=( LCCluster&& rhs ) = default;

  void clear();

  inline double getX() const { return _position[0]; }
  inline double getY() const { return _position[1]; }
  inline double getZ() const { return _position[2]; }
  inline double getE() const { return _energy; }
  inline double getWeight() const { return _weight; }
  inline GlobalMethodsClass::WeightingMethod_t getMethod() const { return _method; }
  inline double getEnergy() const { return getE(); }
  inline double getTheta() const { return _theta; }
  inline double getPhi() const { return _phi; }
  inline double getRZStart() const { return _rzStart; }

  inline const double* getPosition() const { return _position; }
  inline const double* getPositionAtFront() const { return _positionAtFront; }

  inline void setX(double x) {
    _position[0] = x;
    CalculatePhi();
  }
  inline void setY(double y) {
    _position[1] = y;
    CalculatePhi();
  }
  inline void setZ(double z) {
    _position[2] = z;
    CalculateTheta();
  }
  inline void setPosition(double x, double y, double z) {
    _position[0] = x;
    _position[1] = y;
    _position[2] = z;
    CalculatePhi();
    CalculateTheta();
  }
  inline void setWeight(double w) { _weight = w; }
  inline void setTheta(double t) { _theta = t; }
  inline void setPhi(double p) { _phi = p; }

  void addToEnergy(double E) { _energy += E; }

  friend std::ostream& operator<<(std::ostream& o, const LCCluster& rhs);

  inline VecCalHit const& getCaloHits() const { return _caloHits; }

  /// calculate the cluster position based on the caloHits associated to the cluster
  void recalculatePositionFromHits(GlobalMethodsClass const& gmc);

private:
  void CalculatePhi();
  void CalculateTheta();
  double _position[3] = {0.0, 0.0, 0.0};
  double _positionAtFront[3] = {0.0, 0.0, 0.0};
  double _energy = 0.0, _weight = 0.0;
  GlobalMethodsClass::WeightingMethod_t _method = GlobalMethodsClass::LogMethod;
  double _theta = 0.0, _phi = 0.0, _rzStart = 0.0;
  VecCalHit _caloHits{};
};

inline void LCCluster::CalculatePhi() {
  _phi = atan2(_position[1], _position[0]);
  return;
}
inline void LCCluster::CalculateTheta() {
  _theta = atan(sqrt(_position[1] * _position[1] + _position[0] * _position[0]) / fabs(_position[2]));
  return;
}

using SLCCluster = std::shared_ptr<LCCluster>;

#endif // LCCluster_hh
