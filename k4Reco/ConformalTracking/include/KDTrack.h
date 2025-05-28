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
#ifndef K4RECO_KDTRACK_H
#define K4RECO_KDTRACK_H

#include "KDCluster.h"
#include "Parameters.h"

#include <memory>
#include <stdexcept>
#include <vector>

class TH2F;

// ------------------------------------------------------------------------------------
// The KDTrack class is a simple track class designed to allow fast linear fitting in
// conformal space. The class holds a vector of KDCluster objects (the hits attached to
// the track) and a psuedo-kalman filter object.
// ------------------------------------------------------------------------------------

class KalmanTrack;

class KDTrack {
public:
  //--- Constructor and destructor
  KDTrack(Parameters const& par);

  KDTrack(const KDTrack&) = default;
  KDTrack& operator=(const KDTrack&) = delete;
  KDTrack(KDTrack&&) = default;
  KDTrack& operator=(KDTrack&&) = default;
  ~KDTrack() = default;

  //--- Functions to add and remove clusters
  void add(SKDCluster cluster) { m_clusters.push_back(cluster); }
  void insert(SKDCluster cluster) { m_clusters.insert(m_clusters.begin(), cluster); }
  void remove(int clusterN) {
    if (clusterN < 0 || clusterN >= static_cast<int>(m_clusters.size())) {
      throw std::out_of_range("KDTrack::remove: clusterN out of range");
    }
    m_clusters.erase(m_clusters.begin() + clusterN);
  }

  //--- Fit functions
  double calculateChi2();
  double calculateChi2SZ(TH2F* histo = NULL, bool debug = false);
  void linearRegression(bool highPTfit = false);
  void linearRegressionConformal(bool debug = false);
  double sinc(double) const;
  void FillDistribution(TH2F*);

  //--- Functions to set and return member variables

  // UV fit parameters
  double intercept() const { return m_intercept; }
  double gradient() const { return m_gradient; }
  double quadratic() const { return m_quadratic; }
  double chi2() const { return m_chi2; }
  double chi2ndof() const { return m_chi2ndof; }
  bool rotated() const { return m_rotated; }

  // SZ fit parameters
  double interceptZS() const { return m_interceptZS; }
  double gradientZS() const { return m_gradientZS; }
  double chi2ZS() const { return m_chi2ZS; }
  double chi2ndofZS() const { return m_chi2ndofZS; }

  // Clusters and kalman track pointer
  SharedKDClusters clusters() const { return m_clusters; }

  double pt() const { return m_pT; }

  //--- Each KDTrack contains parameters for the two separate fits,
  //--- along with the errors and list of clusters used

  // UV fit parameters
  double m_gradient = 0.0;
  double m_gradientError = 0.0;
  double m_intercept = 0.0;
  double m_interceptError = 0.0;
  double m_quadratic = 0.0;
  double m_chi2 = 0.0;
  double m_chi2ndof = 0.0;
  bool m_rotated = false;

  // SZ fit parameters
  double m_gradientZS = 0.0;
  double m_interceptZS = 0.0;
  double m_chi2ZS = 0.0;
  double m_chi2ndofZS = 0.0;
  bool m_rotatedSZ = 0.0;
  bool fillFit = false;

  // Clusters and kalman track pointer
  double m_pT = 0.0;
  SharedKDClusters m_clusters{};
  bool m_kalmanFitForward = true;
};

typedef std::vector<std::unique_ptr<KDTrack>> UniqueKDTracks;
typedef std::unique_ptr<KDTrack> UKDTrack;

#endif
