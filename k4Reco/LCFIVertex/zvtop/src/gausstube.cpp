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
#include "../include/gausstube.h"
#include "../../inc/track.h"
#include "../../inc/trackstate.h"
#include <math.h>

namespace vertex_lcfi {
namespace ZVTOP {
  GaussTube::GaussTube(Track* Track)
      : // We make a trackstate to be used by the tube, with the appropriate swimmer
        _TrackState(new TrackState(Track)) {}

  double GaussTube::valueAt(const Vector3& Point) const {
    // Calculate value of UNNORMALISED gaussian at point from covarience matrix
    boost::numeric::ublas::bounded_vector<double, 2> Residual;
    // XY Dist in 2D
    _TrackState->swimToStateNearestXY(Point);
    Residual(0) = _TrackState->xyDistanceTo(Point);
    // Z in 3D
    _TrackState->swimToStateNearest(Point);
    // The 3Ddist , 2Ddist and distance on z plane form a right triangle, convert to zaxis by dividing by sin theta
    // Check hypotenuse longest
    if (_TrackState->distanceTo(Point) < Residual(0))
      Residual(1) = 0;
    else
      Residual(1) = sqrt(_TrackState->distanceTo2(Point) - pow(Residual(0), 2)) *
                    sqrt(pow(_TrackState->parentTrack()->helixRep().tanLambda(), 2) + 1.0);

    // Value of tube = -0.5exp(res.inv(V).res) - Lyons pp 60
    boost::numeric::ublas::bounded_vector<double, 2> temp =
        prec_prod(_TrackState->inversePositionCovarMatrix(), Residual);
    // std::cout << exp(-0.5 * prec_inner_prod(Residual,temp)) << std::endl;
    return exp(-0.5 * prec_inner_prod(Residual, temp));
  }

  GaussTube::~GaussTube() {
    delete _TrackState;
    _TrackState = nullptr;
  }

} // namespace ZVTOP
} // namespace vertex_lcfi
