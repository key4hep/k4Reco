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
#ifndef VERTEXFUNCTIONELEMENT_H
#define VERTEXFUNCTIONELEMENT_H

#include "../../util/inc/vector3.h"

namespace vertex_lcfi {
namespace ZVTOP {
  using namespace vertex_lcfi::util;

  //! Vertex Fuction Element (Tubes, ellipse) Interface
  /*!
  Pure virtual class interface class, cannot be instantiated.
   \author Ben Jeffery (b.jeffery1@physics.ox.ac.uk)
   \version 0.1
   \date    20/09/05
  */
  class VertexFunctionElement {
  public:
    virtual double valueAt(const Vector3& Point) const = 0;
    virtual ~VertexFunctionElement() {}
  };
} // namespace ZVTOP
} // namespace vertex_lcfi

#endif // VERTEXFUNCTIONELEMENT_H
