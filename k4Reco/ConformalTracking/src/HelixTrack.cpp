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
#include "HelixTrack.h"

// Remove when the warnings are fixed in KalTest
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-copy"
#include <kaltest/THelicalTrack.h>
#pragma GCC diagnostic pop

#include <TVector3.h>
#include <cmath>

// defines if s of the helix increases in the direction of x2 to x3
bool HelixTrack::forwards = true;

HelixTrack::HelixTrack(const edm4hep::Vector3d& x1, const edm4hep::Vector3d& x2, const edm4hep::Vector3d& x3, double Bz,
                       bool direction) {
  // Make a KalTest THelicalTrack
  TVector3 p1(x1[0], x1[1], x1[2]);
  TVector3 p2(x2[0], x2[1], x2[2]);
  TVector3 p3(x3[0], x3[1], x3[2]);

  THelicalTrack helicalTrack(p1, p2, p3, Bz, direction);

  // Set the track parameters and convert from the KalTest system to the lcio system

  m_phi0 = toBaseRange(helicalTrack.GetPhi0() + M_PI / 2.);
  m_omega = 1. / helicalTrack.GetRho();
  m_z0 = helicalTrack.GetDz();
  m_d0 = -helicalTrack.GetDrho();
  m_tanLambda = helicalTrack.GetTanLambda();

  m_ref_point_x = helicalTrack.GetPivot().X();
  m_ref_point_y = helicalTrack.GetPivot().Y();
  m_ref_point_z = helicalTrack.GetPivot().Z();
}

double HelixTrack::moveRefPoint(double x, double y, double z) {
  const double radius = 1.0 / m_omega;

  const double sinPhi0 = sin(m_phi0);
  const double cosPhi0 = cos(m_phi0);

  const double deltaX = x - m_ref_point_x;
  const double deltaY = y - m_ref_point_y;

  double phi0Prime = atan2(sinPhi0 - (deltaX / (radius - m_d0)), cosPhi0 + (deltaY / (radius - m_d0)));

  while (phi0Prime < 0)
    phi0Prime += 2.0 * M_PI;
  while (phi0Prime >= 2.0 * M_PI)
    phi0Prime -= 2.0 * M_PI;

  const double d0Prime = m_d0 + deltaX * sinPhi0 - deltaY * cosPhi0 +
                         ((deltaX * cosPhi0 + deltaY * sinPhi0) * tan((phi0Prime - m_phi0) / 2.0));

  // In order to have terms which behave well as Omega->0 we make use of deltaX and deltaY to replace sin( phi0Prime -
  // phi0 ) and cos( phi0Prime - phi0 )

  const double sinDeltaPhi = (-m_omega / (1.0 - (m_omega * d0Prime))) * (deltaX * cosPhi0 + deltaY * sinPhi0);

  const double cosDeltaPhi = 1.0 + (m_omega * m_omega / (2.0 * (1.0 - m_omega * d0Prime))) *
                                       (d0Prime * d0Prime - (deltaX + m_d0 * sinPhi0) * (deltaX + m_d0 * sinPhi0) -
                                        (deltaY - m_d0 * cosPhi0) * (deltaY - m_d0 * cosPhi0));

  const double s = atan2(-sinDeltaPhi, cosDeltaPhi) / m_omega;

  const double z0Prime = m_ref_point_z - z + m_z0 + m_tanLambda * s;

  phi0Prime = toBaseRange(phi0Prime);

  m_d0 = d0Prime;
  m_phi0 = phi0Prime;
  m_z0 = z0Prime;

  m_ref_point_x = x;
  m_ref_point_y = y;
  m_ref_point_z = z;

  return (s / radius);
}
