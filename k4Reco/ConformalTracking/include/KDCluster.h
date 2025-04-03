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
#ifndef K4RECO_KDCLUSTERS_H
#define K4RECO_KDCLUSTERS_H

#include <edm4hep/TrackerHitPlane.h>

#include <cmath>
#include <memory>
#include <vector>

// ------------------------------------------------------------------------------------
// The KDCluster class is a simple hit class used in the Cellular Automaton tracking.
// They are designed to be a lightweight object containing all of the information
// needed for tracking in conformal space, to be ordered in a binary tree (KDtree,
// hence the name). The hits contain their co-ordinates in conformal space, both in
// cartesian (u-v) and polar (r-theta) notation. Cartesian positions allow fast nearest
// neighbour searches, while the polar positions allow both nearest neighbour searches
// in theta (avoiding sector definitions) and directed tracking to flow inside-out or
// outside-in. They additionally contain minimal detector information (id, layer and side)
// ------------------------------------------------------------------------------------

class KDCluster {
public:
  // Constructors, main initialisation is with tracker hit
  KDCluster() = default;
  KDCluster(const edm4hep::TrackerHitPlane& hit, bool endcap, bool forward = false)
      : m_x(hit.getPosition()[0]), m_y(hit.getPosition()[1]),
        m_z(hit.getPosition()[2]), // Store the (unaltered) z position
        m_endcap(endcap), m_forward(forward) {
    // Calculate conformal position in cartesian co-ordinates
    const double radius2 = (m_x * m_x + m_y * m_y);
    const double radius2Inv = 1. / radius2;
    const double radius = sqrt(radius2);
    m_u = m_x * radius2Inv;
    m_v = m_y * radius2Inv;
    // Note the position in polar co-ordinates
    m_r = 1. / radius;
    m_theta = atan2(m_v, m_u) + M_PI;
    m_radius = radius;
    // Get the error in the conformal (uv) plane
    // This is the xy error projected. Unfortunately, the
    // dU is not always aligned with the xy plane, it might
    // be dV. Check and take the smallest
    if (hit.getDv() < m_error) {
      m_error = hit.getDv();
      m_errorZ = hit.getDu();
    } else {
      m_error = hit.getDu();
      m_errorZ = hit.getDv();
    }

    const double sinTheta = sin(m_theta);
    const double cosTheta = cos(m_theta);
    if (endcap) {
      m_errorU = (m_error * std::abs(sinTheta) + m_errorZ * std::abs(cosTheta)) * (m_r * m_r);
      m_errorV = (m_error * std::abs(cosTheta) + m_errorZ * std::abs(sinTheta)) * (m_r * m_r);
      m_errorX = m_error * sinTheta;
      m_errorY = m_error * cosTheta;

      // Need to set endcap error in z!
      m_errorZ = 0.25;

    } else {
      m_errorU = m_error * m_r * m_r * sinTheta;
      m_errorV = m_error * m_r * m_r * cosTheta;
      m_errorX = m_error * sinTheta;
      m_errorY = m_error * cosTheta;
    }
  }
  KDCluster(const KDCluster&) = delete;
  KDCluster& operator=(const KDCluster&) = delete;
  KDCluster(KDCluster&&) = default;
  KDCluster& operator=(KDCluster&&) = default;
  ~KDCluster() = default;

  double getX() const { return m_x; }
  double getY() const { return m_y; }
  double getU() const { return m_u; }
  double getV() const { return m_v; }
  double getR() const { return m_r; }
  double getRadius() const { return m_radius; } // "real" radius in xy
  double getTheta() const { return m_theta; }
  double getZ() const { return m_z; }
  float getError() const { return m_error; }
  float getErrorX() const { return m_errorX; }
  float getErrorY() const { return m_errorY; }
  float getErrorU() const { return m_errorU; }
  float getErrorV() const { return m_errorV; }
  float getErrorZ() const { return m_errorZ; }
  bool used() const { return m_used; }

  void setU(double u) { m_u = u; }
  void setV(double v) { m_v = v; }
  void setR(double r) { m_r = r; }
  void setTheta(double theta) { m_theta = theta; }
  void setZ(double z) { m_z = z; }
  void setError(float error) { m_error = error; }
  void used(bool used) { m_used = used; }

  // Subdetector information
  void setDetectorInfo(long subdet, long side, long layer, long module, long sensor) {
    m_subdet = subdet;
    m_side = side;
    m_layer = layer;
    m_module = module;
    m_sensor = sensor;
  }
  // Used long because this is what's used in DD4hep
  long getSubdetector() const { return m_subdet; }
  long getSide() const { return m_side; }
  long getLayer() const { return m_layer; }
  long getModule() const { return m_module; }
  long getSensor() const { return m_sensor; }
  bool forward() const { return m_forward; }
  bool endcap() const { return m_endcap; }

  // Check if another hit is on the same detecting layer
  bool sameLayer(const std::shared_ptr<KDCluster> kdhit) const {
    if (kdhit->getSubdetector() == m_subdet && kdhit->getSide() == m_side && kdhit->getLayer() == m_layer)
      return true;
    return false;
  }

  // Check if another hit is on the same sensor of the same detecting layer
  bool sameSensor(const std::shared_ptr<KDCluster> kdhit) const {
    return kdhit->getLayer() == m_layer && kdhit->getSubdetector() == m_subdet && kdhit->getSide() == m_side &&
           kdhit->getModule() == m_module && kdhit->getSensor() == m_sensor;
  }

private:
  // Each hit contains the conformal co-ordinates in cartesian
  // and polar notation, plus the subdetector information
  double m_x = 0.0;
  double m_y = 0.0;
  double m_u = 0.0;
  double m_v = 0.0;
  double m_r = 0.0;
  double m_radius = 0.0; // "real" radius in xy
  double m_z = 0.0;
  float m_error = 0.0;
  float m_errorX = 0.0;
  float m_errorY = 0.0;
  float m_errorU = 0.0;
  float m_errorV = 0.0;
  float m_errorZ = 0.0;
  double m_theta = 0.0;
  long m_subdet = 0;
  long m_side = 0;
  long m_layer = 0;
  long m_module = 0;
  long m_sensor = 0;
  bool m_used = false;
  bool m_endcap = false;
  bool m_forward = false;
};

typedef std::shared_ptr<KDCluster> SKDCluster;
typedef std::vector<SKDCluster> SharedKDClusters;

#endif
