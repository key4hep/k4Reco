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
#ifndef VERTEXFINDERGHOST_H
#define VERTEXFINDERGHOST_H

#include "../../util/inc/vector3.h"
#include <list>
#include <vector>

using vertex_lcfi::util::Vector3;

namespace vertex_lcfi {
class Track;

namespace ZVTOP {
  // Forward Declarations
  class CandidateVertex;
  class InteractionPoint;

  //! Vertex Finding object - classic ZVTOP
  /*!
  Perform the ghost track algorithm on a set of tracks
  \author Ben Jeffery (b.jeffery1@physics.ox.ac.uk)
   \version 0.1
   \date    10/03/06
  */
  class VertexFinderGhost {
  public:
    // Constructors NB remember algoritm parameters are set per vertexfinder
    VertexFinderGhost(const std::vector<Track*>& Tracks, InteractionPoint* IP);
    VertexFinderGhost(const vertex_lcfi::ZVTOP::VertexFinderGhost&) = delete;
    VertexFinderGhost& operator=(const vertex_lcfi::ZVTOP::VertexFinderGhost&) = delete;
    // Need to invaliate vertex result if these changed
    void addTrack(Track* const Track);
    void setIP(InteractionPoint* const IP);
    Vector3 seedDirection() const { return _SeedDirection; }
    Vector3& seedDirection() { return _SeedDirection; }
    double minimumProbability() const { return _MinimumProbability; }
    double& minimumProbability() { return _MinimumProbability; }
    double initialGhostWidth() const { return _InitialGhostWidth; }
    double& initialGhostWidth() { return _InitialGhostWidth; }
    double maxChi2Allowed() const { return _MaxChi2Allowed; }
    double& maxChi2Allowed() { return _MaxChi2Allowed; }

    // returns true if track was in set and removed
    bool removeTrack(Track* const Track);
    bool clearIP();

    Track* lastGhost() const;
    // run ZVGST!
    std::list<CandidateVertex*> findVertices();

  private:
    Vector3 _SeedDirection{};
    double _MinimumProbability = 0.0;
    double _InitialGhostWidth = 0.0;
    double _MaxChi2Allowed = 0.0;

    std::vector<Track*> _TrackList{};
    InteractionPoint* _IP = nullptr;
    Track* _LastGhost = nullptr;
  };
} // namespace ZVTOP
} // namespace vertex_lcfi
#endif // VERTEXFINDERGHOST_H
