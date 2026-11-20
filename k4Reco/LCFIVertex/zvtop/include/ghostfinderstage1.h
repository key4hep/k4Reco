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
#ifndef GHOSTFINDERSTAGE1_H
#define GHOSTFINDERSTAGE1_H

#include "../../inc/track.h"
#include "../../util/inc/vector3.h"
#include "../include/vertexfitterlsm.h"
#include <map>
#include <vector>

namespace vertex_lcfi {
class TrackState;
namespace ZVTOP {
  // Forward Declarations
  class CandidateVertex;
  class InteractionPoint;

  //! First stage of ghost track alorithm - ghost finder
  /*!
  From a given set of Tracks and an InteractionPoint find a track
  with a direction and width consistant with a decaying particle forming the
  jet.
  <br> Note that currently the ghost track always originates at the origin
  and ignores the position of the Interaction Point. This should be fixed in a future release.
  \todo Upgrade to movable IP
  \author Ben Jeffery (b.jeffery1@physics.ox.ac.uk)
   \version 0.1
   \date    12/12/06
  */
  class GhostFinderStage1 {
  public:
    //! Default Constructor
    /*!
    Creates a finder
    */
    GhostFinderStage1();

    //! Find the ghost track
    /*!
    Using a given initial width and direction the ghost track is swivelled in Phi and Theta to
    minimise the chi sqaured to the other tracks. The ghost track width is then inflated to be consistant
    with the set of tracks to the level of MaxChiAllowed. The minimisation is then repeated and the
    width adjusted again.
    <br> There are some details of the minimisation which are deatiled in the GhostTrack paper
    */
    Track* findGhost(double InitialWidth, double MaxChi2Allowed, const Vector3& JetDir,
                     const std::vector<Track*>& JetTracks, InteractionPoint* IP);

    // method that gives a value for chi2 at a point, this specific name
    // is used so that the function minimiser template can be used.
    double valueAt(std::vector<double> CurrentAngles);

  private:
    Track _makeGhost(std::vector<double> Angles, double Width);
    double _tanLambda(double theta);
    void _fillLZeroChis();
    double _findAdjustedWidth(TrackState* GhostTrack, double CurrentWidth, std::vector<TrackState*>& JetTracks,
                              double MaxChi2Allowed);

    double _CurrentWidth = 0.0;
    Vector3 _JetDir{};
    int _UseChiEquation = 0;
    std::vector<TrackState*> _JetTracks{};
    std::map<TrackState*, double> _ChiToLZero{};
    VertexFitterLSM _Fitter{};
  };
} // namespace ZVTOP
} // namespace vertex_lcfi
#endif // GHOSTFINDERSTAGE1_H
