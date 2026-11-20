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
#ifndef GAUSSTUBE_H
#define GAUSSTUBE_H

#include "../../util/inc/vector3.h"
#include "../include/vertexfunctionelement.h"

namespace vertex_lcfi {
class Track;
class TrackState;
namespace ZVTOP {
  //! Gaussian tube component of the vertex function
  /*!
  Gaussian Tube in detector space (usually representing a track)
  The track position and size are determined by the position and
  covariance matrix of a given Track by:
  \f[
  V(\mathbf{r})=e^{-0.5rV^{-1}r^{T}}
  \f]
  where r is the vector from to query point to the closest point on the track
  and V is the covariance matrix of the track. (Both in rPhi,z space)
  <br>Note this is a deliberatly unnormalised gaussian.
  <br>The Track is not modified at any point by this class
  <br>This guassian tube makes its own trackstate object from the track given
  at construction which it uses to perform the calculation.
  \author Ben Jeffery (b.jeffery1@physics.ox.ac.uk)
   \version 0.1
   \date    20/09/05
  */
  class GaussTube : public VertexFunctionElement {
  public:
    //! Construct from a Track, makes a trackstate for its own use.
    /*!
    As a trackstate is made from the track, note that any changes to the track will not
    propogate to the tube.
    \param Track Track that forms the guassian tube
    */
    GaussTube(Track* Track);
    GaussTube(const GaussTube&) = delete;
    GaussTube& operator=(const GaussTube&) = delete;

    //! Delete tube and trackstate used
    ~GaussTube();

    //! Calculate the value of the tube at point
    /*!
    Swims the track to the point of closest approach and calculate the tube value.
    \param Point Vector3 of the spacial point
    \return Value of tube at point
    */
    double valueAt(const Vector3& Point) const;

  private:
    TrackState* _TrackState = nullptr;
  };
} // namespace ZVTOP
} // namespace vertex_lcfi
#endif // GAUSSTUBE_H
