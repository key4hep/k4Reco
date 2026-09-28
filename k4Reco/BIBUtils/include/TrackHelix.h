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
#ifndef K4RECO_TRACKHELIX_H
#define K4RECO_TRACKHELIX_H 1

#include <cmath>

namespace k4reco::bibutils {

/**
 * Minimal, self-contained helix used by FilterConeHits, so the package does not
 * depend on MarlinUtil. It reproduces exactly the two operations the cone filter
 * needs from MarlinUtil's HelixClass (iLCSoft, BSD-3-Clause):
 *   - the distance from a space point to the helix (getDistanceToPoint), and
 *   - the parameter at which the helix crosses a barrel cylinder of a given
 *     radius (getPointOnCircle).
 * The same conventions are kept (positions in mm, momenta in GeV, B in Tesla
 * along +z, the "time" return value being the path-length parameter t such that
 * z = z_ref + t * pz).
 */
class TrackHelix {
public:
  /// Initialize from a reference point (production vertex), momentum, charge and
  /// the magnetic field Bz [Tesla]. Mirrors HelixClass::Initialize_VP.
  void initialize(const double pos[3], const double mom[3], double charge, double bField) {
    for (int i = 0; i < 3; ++i) {
      m_ref[i] = pos[i];
      m_mom[i] = mom[i];
    }
    m_charge = charge;
    m_pxy = std::sqrt(mom[0] * mom[0] + mom[1] * mom[1]);
    m_radius = m_pxy / (kFCT * bField);
    m_tanLambda = m_pxy > 0. ? mom[2] / m_pxy : 0.;
    const double phiMom = std::atan2(mom[1], mom[0]);
    m_xCentre = pos[0] + m_radius * std::cos(phiMom - kHalfPi * charge);
    m_yCentre = pos[1] + m_radius * std::sin(phiMom - kHalfPi * charge);
  }

  /// Distance from xPoint [mm] to the helix. Fills distance = {dXY, dZ, d3D} and
  /// returns the path-length parameter at the point of closest approach.
  double getDistanceToPoint(const double xPoint[3], double distance[3]) const {
    const double phi = std::atan2(xPoint[1] - m_yCentre, xPoint[0] - m_xCentre);
    const double phi0 = std::atan2(m_ref[1] - m_yCentre, m_ref[0] - m_xCentre);

    const double dx = m_xCentre - xPoint[0];
    const double dy = m_yCentre - xPoint[1];
    const double distXY = std::fabs(std::sqrt(dx * dx + dy * dy) - m_radius);

    int nCircles = 0;
    if (std::fabs(m_tanLambda * m_radius) > 1.0e-20) {
      const double xCircles = (phi0 - phi - m_charge * (xPoint[2] - m_ref[2]) / (m_tanLambda * m_radius)) / kTwoPi;
      nCircles = nearestInteger(xCircles);
    }

    const double dPhi = kTwoPi * static_cast<double>(nCircles) + phi - phi0;
    const double zOnHelix = m_ref[2] - m_charge * m_radius * m_tanLambda * dPhi;
    const double distZ = std::fabs(zOnHelix - xPoint[2]);

    double time;
    if (std::fabs(m_mom[2]) > 1.0e-20) {
      time = (zOnHelix - m_ref[2]) / m_mom[2];
    } else {
      time = m_charge * m_radius * dPhi / m_pxy;
    }

    distance[0] = distXY;
    distance[1] = distZ;
    distance[2] = std::sqrt(distXY * distXY + distZ * distZ);
    return time;
  }

  /// Path-length parameter at which the helix first crosses the cylinder of the
  /// given radius, measured from ref. Returns -1e20 if the helix never reaches
  /// (or always lies outside) the cylinder, matching HelixClass::getPointOnCircle.
  double getPointOnCircle(double cylinderRadius, const double ref[3]) const {
    const double distCenterToIP = std::sqrt(m_xCentre * m_xCentre + m_yCentre * m_yCentre);

    if ((distCenterToIP + m_radius) < cylinderRadius) {
      return -1.0e20;
    }
    if ((m_radius + cylinderRadius) < distCenterToIP) {
      return -1.0e20;
    }

    const double phiCentre = std::atan2(m_yCentre, m_xCentre);
    double phiStar = cylinderRadius * cylinderRadius + distCenterToIP * distCenterToIP - m_radius * m_radius;
    phiStar = 0.5 * phiStar / std::fmax(1.0e-20, cylinderRadius * distCenterToIP);
    if (phiStar > 1.0) {
      phiStar = 0.9999999;
    } else if (phiStar < -1.0) {
      phiStar = -0.9999999;
    }
    phiStar = std::acos(phiStar);

    const double xx1 = cylinderRadius * std::cos(phiCentre + phiStar);
    const double yy1 = cylinderRadius * std::sin(phiCentre + phiStar);
    const double xx2 = cylinderRadius * std::cos(phiCentre - phiStar);
    const double yy2 = cylinderRadius * std::sin(phiCentre - phiStar);

    const double phi1 = std::atan2(yy1 - m_yCentre, xx1 - m_xCentre);
    const double phi2 = std::atan2(yy2 - m_yCentre, xx2 - m_xCentre);
    const double phi0 = std::atan2(ref[1] - m_yCentre, ref[0] - m_xCentre);

    double dphi1 = phi1 - phi0;
    double dphi2 = phi2 - phi0;
    if (dphi1 < 0 && m_charge < 0) {
      dphi1 += kTwoPi;
    } else if (dphi1 > 0 && m_charge > 0) {
      dphi1 -= kTwoPi;
    }
    if (dphi2 < 0 && m_charge < 0) {
      dphi2 += kTwoPi;
    } else if (dphi2 > 0 && m_charge > 0) {
      dphi2 -= kTwoPi;
    }

    const double tt1 = -m_charge * dphi1 * m_radius / m_pxy;
    const double tt2 = -m_charge * dphi2 * m_radius / m_pxy;
    return tt1 < tt2 ? tt1 : tt2;
  }

private:
  /// Integer closest to x (same tie-breaking as the original HelixClass).
  static int nearestInteger(double x) {
    int n1;
    if (x >= 0.) {
      n1 = static_cast<int>(x);
    } else {
      n1 = static_cast<int>(x) - 1;
    }
    const int n2 = n1 + 1;
    return std::fabs(n1 - x) < std::fabs(n2 - x) ? n1 : n2;
  }

  // 0.1 * 2.99792458 GeV / (T m) expressed for positions in mm: radius[mm] = pT[GeV] / (kFCT * B[T]).
  static constexpr double kFCT = 2.99792458e-4;
  static constexpr double kHalfPi = 0.5 * M_PI;
  static constexpr double kTwoPi = 2.0 * M_PI;

  double m_ref[3]{};
  double m_mom[3]{};
  double m_charge{0.};
  double m_pxy{0.};
  double m_radius{0.};
  double m_tanLambda{0.};
  double m_xCentre{0.};
  double m_yCentre{0.};
};

} // namespace k4reco::bibutils

#endif
