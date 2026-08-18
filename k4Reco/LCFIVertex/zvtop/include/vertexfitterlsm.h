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
#ifndef VERTEXFITTERLSM_H
#define VERTEXFITTERLSM_H

#include "../../util/inc/matrix.h"
#include "../../util/inc/vector3.h"
#include "vertexfitter.h"
#include <map>
#include <vector>

namespace vertex_lcfi {
class TrackState;

namespace ZVTOP {
  using namespace vertex_lcfi::util;
  class CandidateVertex;
  class InteractionPoint;

  class VertexFitterLSM : public VertexFitter {
  public:
    VertexFitterLSM();
    ~VertexFitterLSM() {}
    VertexFitterLSM(const VertexFitterLSM&) = delete;
    VertexFitterLSM& operator=(const VertexFitterLSM&) = delete;
    // CandidateVertex fitVertex(const std::vector<TrackState*> & Tracks, InteractionPoint* IP, bool CalculateError);
    void fitVertex(const std::vector<TrackState*>& Tracks, InteractionPoint* IP, Vector3& Result);
    void fitVertex(const std::vector<TrackState*>& Tracks, InteractionPoint* IP, Vector3& Result,
                   double& ChiSquaredOfFit);
    void fitVertex(const std::vector<TrackState*>& Tracks, InteractionPoint* IP, Vector3& Result,
                   double& ChiSquaredOfFit, std::map<TrackState*, double>& ChiSquaredOfTrack, double& ChiSquaredOfIP);
    void fitVertex(const std::vector<TrackState*>& Tracks, InteractionPoint* IP, Vector3& Result,
                   Matrix3x3& ResultError, double& ChiSquaredOfFit, std::map<TrackState*, double>& ChiSquaredOfTrack,
                   double& ChiSquaredOfIP);
    // method that gives a value for chi2 at a point, this specific name
    // is used so that the function minimiser template can be used.
    double valueAt(const Vector3& point);
    double valueAt(const std::vector<double>& point);

    void setSeed(Vector3 Seed);
    void setInitialStep(double Step);

  private:
    std::vector<TrackState*> _trackStateList{}; // a copy of the trackStates being fitted
    InteractionPoint* _ip = nullptr;
    Vector3 _ManualSeed{};
    bool _UseManualSeed = false;
    double _InitialStep = 0.0;
    double _chi2Contribution(const Vector3& point, TrackState* pTrackState); // contribution from each individual track
    double _chi2Contribution(const Vector3& point,
                             InteractionPoint* pIP); // the contribution from the ip only (N.B. pIP could be NULL)
  };
} // namespace ZVTOP
} // namespace vertex_lcfi
#endif // VERTEXFITTERLSM_H
