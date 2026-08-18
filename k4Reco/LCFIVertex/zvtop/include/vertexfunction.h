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
#ifndef VERTEXFUNCTION_H
#define VERTEXFUNCTION_H

#include "../../util/inc/matrix.h"
#include "../../util/inc/vector3.h"
#include <vector>

namespace vertex_lcfi {
class TrackState;

namespace ZVTOP {
  using namespace vertex_lcfi::util;
  // Forward Declaration
  class VertexFunctionElement;
  class InteractionPoint;

  //! Vertex Function Interface
  /*!
  Pure virtual class interface class, cannot be instantiated.
   \author Ben Jeffery (b.jeffery1@physics.ox.ac.uk)
   \version 0.1
   \date    20/09/05
  */
  class VertexFunction {

  public:
    // Query Methods
    virtual double valueAt(const Vector3& Point) const = 0;
    virtual Matrix3x3 firstDervAt(const Vector3& Point) const = 0;
    virtual Matrix3x3 secondDervAt(const Vector3& Point) const = 0;
    virtual ~VertexFunction() {}
  };
} // namespace ZVTOP
} // namespace vertex_lcfi
#endif // VERTEXFUNCTION_H
