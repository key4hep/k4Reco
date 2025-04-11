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
#ifndef K4RECO_HELIXTRACK_H
#define K4RECO_HELIXTRACK_H

#include <edm4hep/Vector3d.h>

#include <cmath>

class HelixTrack {
public:
  HelixTrack(double ref_point_x, double ref_point_y, double ref_point_z, double d0, double z0, double phi0,
             double omega, double tanLambda)
      : m_ref_point_x(ref_point_x), m_ref_point_y(ref_point_y), m_ref_point_z(ref_point_z), m_d0(d0), m_z0(z0),
        m_phi0(phi0), m_omega(omega), m_tanLambda(tanLambda) {
    m_phi0 = toBaseRange(m_phi0);
  }

  HelixTrack(const edm4hep::Vector3d& x1, const edm4hep::Vector3d& x2, const edm4hep::Vector3d& x3, double Bz,
             bool direction);

  HelixTrack(const HelixTrack&) = delete;
  HelixTrack& operator=(const HelixTrack&) = delete;
  HelixTrack(HelixTrack&&) = delete;
  HelixTrack& operator=(HelixTrack&&) = delete;
  ~HelixTrack() = default;

  double moveRefPoint(double x, double y, double z);

  double getRefPointX() const { return m_ref_point_x; }
  double getRefPointY() const { return m_ref_point_y; }
  double getRefPointZ() const { return m_ref_point_z; }
  double getD0() const { return m_d0; }
  double getZ0() const { return m_z0; }
  double getPhi0() const { return m_phi0; }
  double getOmega() const { return m_omega; }
  double getTanLambda() const { return m_tanLambda; }

  // defines if s of the helix increases in the direction of x2 to x3
  static bool forwards;

private:
  double m_ref_point_x = 0.0;
  double m_ref_point_y = 0.0;
  double m_ref_point_z = 0.0;
  double m_d0 = 0.0;
  double m_z0 = 0.0;
  double m_phi0 = 0.0;
  double m_omega = 0.0;
  double m_tanLambda = 0.0;

  /** helper function to restrict the range of the azimuthal angle to ]-pi,pi]*/
  inline double toBaseRange(double phi) const {
    phi = std::fmod(phi + M_PI, 2. * M_PI);
    if (phi < 0)
      phi += M_PI;
    else
      phi -= M_PI;
    return phi;
  }
};

#endif
