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
#ifndef VERTEXFITTERKALMAN_H
#define VERTEXFITTERKALMAN_H

#include "../../inc/TState.h"
#include "vertexfitter.h"

namespace vertex_lcfi {
class TrackState;

namespace ZVTOP {
  class InteractionPoint;

  class VertexFitterKalman : public VertexFitter {
  public:
    VertexFitterKalman();
    ~VertexFitterKalman() {}

    void fitVertex(const std::vector<TrackState*>& Tracks, InteractionPoint* IP, Vector3& Result);

    void fitVertex(const std::vector<TrackState*>& Tracks, InteractionPoint* IP, Vector3& Result,
                   double& ChiSquaredOfFit);

    void fitVertex(const std::vector<TrackState*>& Tracks, InteractionPoint* IP, Vector3& Result,
                   double& ChiSquaredOfFit, std::map<TrackState*, double>& ChiSquaredOfTrack, double& ChiSquaredOfIP);

    void fitVertex(const std::vector<TrackState*>& Tracks, InteractionPoint* IP, Vector3& Result,
                   Matrix3x3& ResultError, double& ChiSquaredOfFit, std::map<TrackState*, double>& ChiSquaredOfTrack,
                   double& ChiSquaredOfIP);

    void fitVertex(const std::vector<TrackState*>& Tracks, InteractionPoint* IP, Vector3& Result,
                   Matrix3x3& ResultError, double& ChiSquaredOfFit);

    double getDeviationFromVertex(const TState* state, const double v[], const double Cv[]) const;

    void setSeed(Vector3 Seed);

    bool estimateVertex(double vtx[]);

    double robustMean(std::vector<double> vec);

  private:
    std::vector<TState> fStates{};
    std::vector<double> fChi2chain{};

    Vector3 m_manualSeed{};
    bool m_useManualSeed = false;

    double fP[6];
    double fC[21];
    double fVtxGuess[3];

    int fNDF = 0;
    double fChi2 = 0.0;
  };
} // namespace ZVTOP
} // namespace vertex_lcfi

#endif // VERTEXFITTERKALMAN_H
