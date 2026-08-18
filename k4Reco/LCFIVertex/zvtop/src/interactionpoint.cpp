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
#include "../include/interactionpoint.h"
#include "../../util/inc/matrix.h"

namespace vertex_lcfi {
namespace ZVTOP {
  InteractionPoint::InteractionPoint() {}

  InteractionPoint::InteractionPoint(const Vector3& Position, const SymMatrix3x3& ErrorMatrix)
      : _Position(Position), _ErrorMatrix(ErrorMatrix), _InvErrorMatrix(InvertMatrix(ErrorMatrix)) {}

  double InteractionPoint::distanceTo(const Vector3& Point) const { return _Position.distanceTo(Point); }

  const Vector3& InteractionPoint::position() const { return _Position; }

  const SymMatrix3x3& InteractionPoint::errorMatrix() const { return _ErrorMatrix; }

  const Matrix3x3& InteractionPoint::inverseErrorMatrix() const { return _InvErrorMatrix; }

  double InteractionPoint::chi2(const Vector3& Point) const {
    Vector3 Residual = Point - (this->position());
    // std::cout << "Res: " << Residual << std::endl ;
    // std::cout << "Err: " << _InvErrorMatrix<< std::endl;
    // std::cout << "Err: " << this->inverseErrorMatrix()<< std::endl;
    // std::cout << "Chi2: " << prec_inner_prod(trans(Residual),prec_prod(this->inverseErrorMatrix(), Residual))<<
    // std::endl<< std::endl;
    return prec_inner_prod(trans(Residual), prec_prod(this->inverseErrorMatrix(), Residual));
  }
} // namespace ZVTOP
} // namespace vertex_lcfi
