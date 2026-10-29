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
#ifndef INCLUDE_TSTATE_H
#define INCLUDE_TSTATE_H 1

// Include files

#include "track.h"
#include "trackstate.h"

/** class TState TState.h include/TState.h
 *
 *
 *  author Tomas Lastovicka (LCFI)
 *  date   2007-07-26
 */

namespace vertex_lcfi {

class TState {

public:
  virtual ~TState();

  TState(TrackState* TrackState);
  TState(const TState&) = default;
  TState& operator=(const vertex_lcfi::TState&) = default;

  void GetMeasurement(const double xyz[], double m[], double V[]) const;
  void TransportBz(double dS, double P[], double C[]) const;

  double GetDStoPointBz(const double xyz[]) const;
  bool GetDStoTStateBz(const TState* p, double& DS, double& DS1) const;

  inline int charge() const { return fQ; }
  inline Track* track() const { return fParentTrack; }
  inline TrackState* trackState() const { return fParentState; }

protected:
  // Pointer to Track that created this state
  Track* fParentTrack = nullptr;

  // Pointer to TrackState that created this state
  TrackState* fParentState = nullptr;

  double fP[6];    //* Main particle parameters {X,Y,Z,Px,Py,Pz}
  double fC[21];   //* Low-triangle covariance matrix of fP
  int fQ = 0;      //* Particle charge
  double fB = 0.0; //* B-field (Bz only)

  double fCLight = 0.0;

private:
};
} // namespace vertex_lcfi

#endif // INCLUDE_TSTATE_H
