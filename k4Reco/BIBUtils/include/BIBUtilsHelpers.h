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
#ifndef K4RECO_BIBUTILSHELPERS_H
#define K4RECO_BIBUTILSHELPERS_H 1

#include <cmath>

namespace k4reco::bibutils {

/// Angle (in radians, in [0, pi]) between two 3-vectors. Reproduces the value
/// returned by TVector3::Angle used in the original Marlin processors without
/// pulling in the ROOT Physics library.
inline double angleBetween(double ax, double ay, double az, double bx, double by, double bz) {
  const double na = std::sqrt(ax * ax + ay * ay + az * az);
  const double nb = std::sqrt(bx * bx + by * by + bz * bz);
  if (na == 0. || nb == 0.) {
    return 0.;
  }
  double cosAngle = (ax * bx + ay * by + az * bz) / (na * nb);
  if (cosAngle > 1.) {
    cosAngle = 1.;
  } else if (cosAngle < -1.) {
    cosAngle = -1.;
  }
  return std::acos(cosAngle);
}

} // namespace k4reco::bibutils

#endif
