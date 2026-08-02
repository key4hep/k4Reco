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
#ifndef GAUSSELLIPSOID_H
#define GAUSSELLIPSOID_H

#include "../include/vertexfunctionelement.h"

#include "../../util/inc/vector3.h"

using namespace vertex_lcfi::util;

namespace vertex_lcfi {
namespace ZVTOP {
  // Forward Declarations
  class InteractionPoint;

  //! Gaussian ellipsoid component of the vertex function
  /*!
  Gaussian Ellipsoid in detector space (usually representing IP)
  The ellipsoid position and size are determined by the position and
  covariance matrix of a given InteractionPoint by:
  \f[
  V(\mathbf{r})=e^{-0.5rV^{-1}r^{T}}
  \f]
  where r is the vector from ip position to query point and V is the
  covariance matrix. (Both in x,y,z space)
  <br>Note this is a deliberatly unnormalised gaussian.
  <br>The interaction point is not modified at any point by this class
   \author Ben Jeffery (b.jeffery1@physics.ox.ac.uk)
   \version 0.1
   \date    20/09/05
  */
  class GaussEllipsoid : public VertexFunctionElement {
  public:
    //! Construct from an interaction point
    /*!
    \param IP Pointer to InteractionPoint to use
    */
    GaussEllipsoid(InteractionPoint* IP);
    GaussEllipsoid(const GaussEllipsoid&) = delete;
    GaussEllipsoid& operator=(const GaussEllipsoid&) = delete;
    //! Calculate the value of the ellipsoid at a point
    /*!
    \param Point Vector3 of the spatial point
    \return Value of ellipsoid at point.
    */
    double valueAt(const Vector3& Point) const;

    //! InteractionPoint object used
    /*!
    \return Pointer to InteractionPoint used by this instance
    */
    InteractionPoint* ip();

  private:
    InteractionPoint* _IP;
  };
} // namespace ZVTOP
} // namespace vertex_lcfi

#endif // GAUSSELLIPSOID_H
