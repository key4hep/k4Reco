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

// Checks TrackHelix::getDistanceToPoint against trajectories obtained by
// integrating the Lorentz force numerically, so that the reference does not
// share any convention (helix centre, sense of rotation) with TrackHelix.

#include "TrackHelix.h"

#include <array>
#include <cmath>
#include <cstdio>
#include <initializer_list>

namespace {

using Vec3 = std::array<double, 3>;

constexpr double kBField = 3.57; // [T], along +z
// Speed of light in GeV / (T mm): dp/ds = charge * kC * (p_hat x B) with p in GeV, s in mm, B in T.
constexpr double kC = 0.299792458e-3;
constexpr double kMaxArcStep = 0.5; // [mm] integration step along the trajectory
constexpr double kTimeTolerance = 1e-6;
constexpr double kDistanceTolerance = 1e-5; // [mm]

int nFailures = 0;

Vec3 cross(const Vec3& a, const Vec3& b) {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

/// Equation of motion in the path-length parameter t of TrackHelix (dr/dt = p, so that
/// z = z_ref + t * pz): dp/dt = charge * kC * (p x B).
void derivatives(const Vec3& mom, double charge, Vec3& dPos, Vec3& dMom) {
  dPos = mom;
  dMom = cross(mom, {0., 0., kBField});
  for (auto& c : dMom) {
    c *= charge * kC;
  }
}

/// Propagates (pos, mom) by the path-length parameter t with a fourth-order Runge-Kutta.
Vec3 propagate(const Vec3& startPos, const Vec3& startMom, double charge, double t) {
  const double p = std::sqrt(startMom[0] * startMom[0] + startMom[1] * startMom[1] + startMom[2] * startMom[2]);
  const int nSteps = static_cast<int>(std::ceil(t * p / kMaxArcStep));
  const double h = t / nSteps;

  Vec3 pos = startPos;
  Vec3 mom = startMom;
  for (int iStep = 0; iStep < nSteps; ++iStep) {
    std::array<Vec3, 4> kPos, kMom;
    derivatives(mom, charge, kPos[0], kMom[0]);
    for (int k = 1; k < 4; ++k) {
      const double f = k < 3 ? 0.5 * h : h;
      Vec3 midMom;
      for (int i = 0; i < 3; ++i) {
        midMom[i] = mom[i] + f * kMom[k - 1][i];
      }
      derivatives(midMom, charge, kPos[k], kMom[k]);
    }
    for (int i = 0; i < 3; ++i) {
      pos[i] += h / 6. * (kPos[0][i] + 2. * kPos[1][i] + 2. * kPos[2][i] + kPos[3][i]);
      mom[i] += h / 6. * (kMom[0][i] + 2. * kMom[1][i] + 2. * kMom[2][i] + kMom[3][i]);
    }
  }
  return pos;
}

void check(const char* label, const Vec3& ref, const Vec3& mom, double charge, double t) {
  k4reco::bibutils::TrackHelix helix;
  helix.initialize(ref.data(), mom.data(), charge, kBField);

  Vec3 point = propagate(ref, mom, charge, t);
  double distance[3] = {0., 0., 0.};
  const double time = helix.getDistanceToPoint(point.data(), distance);

  if (std::fabs(time - t) > kTimeTolerance * std::fmax(1., t) || distance[2] > kDistanceTolerance) {
    std::printf("FAIL %s: charge=%+.0f mom=(%g, %g, %g) t=%g -> time=%g, d3D=%g\n", label, charge, mom[0], mom[1],
                mom[2], t, time, distance[2]);
    ++nFailures;
  }
}

} // namespace

int main() {
  const Vec3 origin{0., 0., 0.};
  const Vec3 displaced{12., -7., 30.};

  for (const double charge : {-1., 1.}) {
    // Transverse track starting at the origin, hit on the forward trajectory.
    check("transverse", origin, {10., 0., 0.}, charge, 10.);

    // Same with a small longitudinal component: must agree with the case above.
    check("almost transverse", origin, {10., 0., 1e-6}, charge, 10.);

    check("forward", origin, {10., 0., 5.}, charge, 10.);

    // Low-pT tracks scanned over the starting direction and over most of the
    // first turn, so that the arc crosses the -pi/pi boundary of atan2.
    const double pt = 0.5;
    const double turnTime = 2. * M_PI / (kC * kBField);
    for (int iDir = 0; iDir < 12; ++iDir) {
      const double phiMom = -M_PI + (iDir + 0.5) * M_PI / 6.;
      const Vec3 looperTransverse{pt * std::cos(phiMom), pt * std::sin(phiMom), 0.};
      const Vec3 looperForward{pt * std::cos(phiMom), pt * std::sin(phiMom), 0.3};
      for (int iFrac = 1; iFrac < 20; ++iFrac) {
        const double t = 0.05 * iFrac * turnTime;
        check("looper transverse", displaced, looperTransverse, charge, t);
        check("looper forward", displaced, looperForward, charge, t);
        check("looper forward, later turn", displaced, looperForward, charge, t + 2. * turnTime);
      }
    }
  }

  if (nFailures > 0) {
    std::printf("%d check(s) failed\n", nFailures);
    return 1;
  }
  std::printf("All TrackHelix checks passed\n");
  return 0;
}
