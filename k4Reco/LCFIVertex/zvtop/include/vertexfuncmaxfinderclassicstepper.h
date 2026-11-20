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
#ifndef VERTEXFUNCMAXFINDERDUMBSTEPPER_H
#define VERTEXFUNCMAXFINDERDUMBSTEPPER_H

#include "../../util/inc/vector3.h"
#include "vertexfuncmaxfinder.h"

namespace vertex_lcfi {
namespace ZVTOP {
  class VertexFunction;

  //! Robust VertexFuncMaxFinder
  /*!
  Similar to original SLD implementaion, minimises along each axis in turn
  by stepping till minimum reached. Currently hard wired step size.
   \author Ben Jeffery (b.jeffery1@physics.ox.ac.uk)
   \version 0.1
   \date    20/09/05
  */
  class VertexFuncMaxFinderClassicStepper : public VertexFuncMaxFinder {
  public:
    VertexFuncMaxFinderClassicStepper();
    Vector3 findNearestMaximum(const Vector3& StartPoint, VertexFunction* VertexFunction);
    VertexFuncMaxFinderClassicStepper(const vertex_lcfi::ZVTOP::VertexFuncMaxFinderClassicStepper&) = delete;
    VertexFuncMaxFinderClassicStepper& operator=(const vertex_lcfi::ZVTOP::VertexFuncMaxFinderClassicStepper&) = delete;

  private:
    Vector3 _CurrentPos{};
    double _CurrentValue = 0;
    VertexFunction* _Function = nullptr;

    void _minimiseAlongAxis(const Vector3& Step);
  };
} // namespace ZVTOP
} // namespace vertex_lcfi
#endif // VERTEXFUNCMAXFINDERDUMBSTEPPER_H
