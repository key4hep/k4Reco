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
#ifndef VERTEXFINDERCLASSIC_H
#define VERTEXFINDERCLASSIC_H

#include "../../util/inc/vector3.h"
#include <list>
#include <vector>

using namespace vertex_lcfi::util;

namespace vertex_lcfi {
class Track;
//! Namespace containing ZVTOP Implementation
namespace ZVTOP {
  // Forward Declarations
  class CandidateVertex;
  class InteractionPoint;
  class VertexFunction;

  //! Vertex Finding object - classic ZVTOP
  /*!
  Very unfinished - A set of tracks that vertex finding is performed on,
  I originally was going to call this a jet, but then realised this may not always be the case.
  Subject to change as we think about how to interface ZVTOP
  Note new mwmory management
   \author Ben Jeffery (b.jeffery1@physics.ox.ac.uk)
   \version 0.1
   \date    20/09/05
  */
  class VertexFinderClassic {
  public:
    // Constructors NB remember algoritm parameters are set per vertexfinder
    VertexFinderClassic(const std::vector<Track*>& Tracks, InteractionPoint* IP, const Vector3& JetAxis,
                        double Kip = 1.0, double Kalpha = 5.0, double TwoProngCut = 10.0, double TrackTrimCut = 10.0,
                        double ResolverCutOff = 0.6);

    VertexFinderClassic(const vertex_lcfi::ZVTOP::VertexFinderClassic&) = delete;
    VertexFinderClassic& operator=(const vertex_lcfi::ZVTOP::VertexFinderClassic&) = delete;
    // Need to invaliate vertex result if these changed
    void addTrack(Track* const Track);
    void setIP(InteractionPoint* const IP);

    // returns true if track was in set and removed
    bool removeTrack(Track* const Track);
    bool clearIP();

    // run ZVRES!
    std::list<CandidateVertex*> findVertices();

  private:
    std::vector<CandidateVertex*> _removeOneTrackNoIPVertices(std::list<CandidateVertex*>* CVList);
    void _ifNoIPAddIP(std::list<CandidateVertex*>* CVList);

    std::vector<Track*> _TrackList{};
    InteractionPoint* _IP = nullptr;
    VertexFunction* _VF = nullptr;

    double _Kip = 0.0;
    double _Kalpha = 0.0;
    Vector3 _JetAxis{};
    double _TwoProngCut = 0.0;
    double _TrackTrimCut = 0.0;
    double _ResolverCutOff = 0.0;
  };
} // namespace ZVTOP
} // namespace vertex_lcfi
#endif // VERTEXFINDERCLASSIC_H
