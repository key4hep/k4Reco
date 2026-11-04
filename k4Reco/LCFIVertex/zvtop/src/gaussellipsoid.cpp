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
#include "../include/gaussellipsoid.h"
#include "../../util/inc/matrix.h"
#include "../include/interactionpoint.h"
#include <cmath>
namespace vertex_lcfi {
namespace ZVTOP {
  GaussEllipsoid::GaussEllipsoid(InteractionPoint* IP) : _IP(IP) {}

  double GaussEllipsoid::valueAt(const Vector3& Point) const {
    // Normalize to origin
    Vector3 RelativePoint = Point - (_IP->position());
    // Calculate value of UNNORMALISED gaussian at point from covarience matrix
    // Lyons pp 60
    SymMatrix3x3 InvErr = _IP->inverseErrorMatrix();
    return exp(-0.5 * prec_inner_prod(RelativePoint, prec_prod(InvErr, RelativePoint)));
  }

  InteractionPoint* GaussEllipsoid::ip() { return _IP; }

} // namespace ZVTOP
} // namespace vertex_lcfi
