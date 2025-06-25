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
#ifndef K4RECO_GAUDIDDKALTESTTRACK_H
#define K4RECO_GAUDIDDKALTESTTRACK_H

#include "GaudiDDKalTest.h"

#include <edm4hep/TrackState.h>
#include <edm4hep/TrackerHitPlane.h>

#include <TMatrixD.h>

#include <cmath>
#include <memory>

class DDVMeasLayer;
class DDVTrackHit;
class THelicalTrack;
class TKalTrack;
class TKalTrackSite;

namespace Gaudi {
class Algorithm;
}

class TObjArray;

const int modeClosest = 0;

/** Implementation of the IMarlinTrack interface, using KalTest and KalDet to provide
 *  the needed functionality for a Kalman Filter.
 */

class GaudiDDKalTestTrack {
public:
  GaudiDDKalTestTrack(
      const Gaudi::Algorithm* thisAlg, GaudiDDKalTest* ktest,
      std::shared_ptr<std::map<const edm4hep::TrackerHit*, DDVTrackHit*>> edm4hep_hits_to_kaltest_hits = nullptr);

public:
  GaudiDDKalTestTrack(const GaudiDDKalTestTrack&) = delete;
  GaudiDDKalTestTrack& operator=(const GaudiDDKalTestTrack&) = delete;
  GaudiDDKalTestTrack(GaudiDDKalTestTrack&&) = delete;
  GaudiDDKalTestTrack& operator=(GaudiDDKalTestTrack&&) = delete;
  ~GaudiDDKalTestTrack();

  /** add hit to track - the hits have to be added ordered in time ( i.e. typically outgoing )
   *  this order will define the direction of the energy loss used in the fit
   */
  int addHit(const edm4hep::TrackerHit* hit);

  /** add hit to track - the hits have to be added ordered in time ( i.e. typically outgoing )
   *  this order will define the direction of the energy loss used in the fit
   */
  int addHit(const edm4hep::TrackerHit* trkhit, const DDVMeasLayer* ml);

  /** add hit to track - the hits have to be added ordered in time ( i.e. typically outgoing )
   *  this order will define the direction of the energy loss used in the fit
   */
  int addHit(const edm4hep::TrackerHit* trkhit, DDVTrackHit* kalhit, const DDVMeasLayer* ml);

  /** initialise the fit with a track state
   *  the fit direction has to be specified using IMarlinTrack::backward or IMarlinTrack::forward.
   *  this is the order that will be used in the fit().
   *  it is the users responsibility that the track state is consistent with the order
   *  of the hits used in addHit() ( i.e. the direction of energy loss )
   *  Note: the bfield_z is not taken from the argument but from the first hit
   *  should consider changing the interface ...
   */
  int initialise(const edm4hep::TrackState& ts, bool fitDirection);

  /** perform the fit of all current hits, returns error code ( IMarlinTrack::success if no error ) .
   *  the fit will be performed  in the order specified at initialise() wrt the order used in addHit(), i.e.
   *  IMarlinTrack::backward implies fitting from the outside to the inside for tracks comming from the IP.
   */
  int fit(double maxChi2Increment = DBL_MAX);

  /** smooth track states from the last filtered hit back to the measurement site associated with the given hit
   */
  int smooth(const edm4hep::TrackerHit* hit);

  /** update the current fit using the supplied hit, return code via int. Provides the Chi2 increment to the fit from
   * adding the hit via reference. the given hit will not be added if chi2increment > maxChi2Increment.
   */
  int addAndFit(const edm4hep::TrackerHit* hit, double& chi2increment, double maxChi2Increment = DBL_MAX);

  /** update the current fit using the supplied hit, return code via int. Provides the Chi2 increment to the fit from
   * adding the hit via reference. the given hit will not be added if chi2increment > maxChi2Increment.
   */
  int addAndFit(DDVTrackHit* kalhit, double& chi2increment, TKalTrackSite*& site, double maxChi2Increment = DBL_MAX);

  /** get track state at measurement associated with the given hit, returning TrackState, chi2 and ndf via reference
   */
  int getTrackState(const edm4hep::TrackerHit* hit, edm4hep::TrackState& ts, double& chi2, int& ndf) const;

  /** get the list of hits included in the fit, together with the chi2 contributions of the hits.
   *  Pointers to the hits together with their chi2 contribution will be filled into a vector of
   *  pairs consitining of the pointer as the first part of the pair and the chi2 contribution as
   *  the second.
   */
  const std::vector<std::pair<const edm4hep::TrackerHit*, double>>& getHitsInFit() const;

  /** get the list of hits which have been rejected by from the fit due to the a chi2 increment greater than threshold,
   *  Pointers to the hits together with their chi2 contribution will be filled into a vector of
   *  pairs consitining of the pointer as the first part of the pair and the chi2 contribution as
   *  the second.
   */
  const std::vector<std::pair<const edm4hep::TrackerHit*, double>>& getOutliers() const;

  /** get the current number of degrees of freedom for the fit.
   */
  int getNDF() const;

  /** get TrackeHit at which fit became constrained, i.e. ndf >= 0
   */
  const edm4hep::TrackerHit* getTrackerHitAtPositiveNDF() const;

  /** propagate the fit at the measurement site associated with the given hit, to the point of closest approach to the
   * given point, returning TrackState, chi2 and ndf via reference
   */
  int propagate(const edm4hep::Vector3d& point, const edm4hep::TrackerHit* hit, edm4hep::TrackState& ts, double& chi2,
                int& ndf);

  /** propagate the fit at the provided measurement site, to the point of closest approach to the given point,
   *  returning TrackState, chi2 and ndf via reference
   */
  int propagate(const edm4hep::Vector3d& point, const TKalTrackSite& site, edm4hep::TrackState& ts, double& chi2,
                int& ndf, const DDVMeasLayer* ml = nullptr);

  /** propagate the fit at the measurement site associated with the given hit, to numbered sensitive layer,
   *  returning TrackState, chi2, ndf and integer ID of sensitive detector element via reference
   */
  int propagateToLayer(int layerID, const edm4hep::TrackerHit* hit, edm4hep::TrackState& ts, double& chi2, int& ndf,
                       int& detElementID, int mode = modeClosest);

  /** propagate the fit at the measurement site, to numbered sensitive layer,
   *  returning TrackState, chi2, ndf and integer ID of sensitive detector element via reference
   */
  int propagateToLayer(int layerID, const TKalTrackSite& site, edm4hep::TrackState& ts, double& chi2, int& ndf,
                       int& detElementID, int mode = modeClosest);

  /** extrapolate the fit at the measurement site, to numbered sensitive layer,
   *  returning intersection point in global coordinates and integer ID of the intersected sensitive detector element
   * via reference
   */
  int intersectionWithLayer(int layerID, const TKalTrackSite& site, edm4hep::Vector3d& point, int& detElementID,
                            const DDVMeasLayer*& ml, int mode = modeClosest);

  /** extrapolate the fit at the measurement site, to sensitive detector elements contained in the std::vector,
   *  and return intersection point in global coordinates via reference
   */
  int findIntersection(std::vector<DDVMeasLayer const*>& meas_modules, const TKalTrackSite& site,
                       edm4hep::Vector3d& point, int& detElementID, const DDVMeasLayer*& ml, int mode = modeClosest);

  /** extrapolate the fit at the measurement site, to the DDVMeasLayer,
   *  and return intersection point in global coordinates via reference
   */
  int findIntersection(const DDVMeasLayer& meas_module, const TKalTrackSite& site, edm4hep::Vector3d& point,
                       double& dphi, int& detElementIDconst, int mode = modeClosest);

  // //** end of memeber functions from IMarlinTrack interface

  /** fill LCIO Track State with parameters from helix and cov matrix
   */
  void ToLCIOTrackState(const TKalTrackSite& site, edm4hep::TrackState& ts, double& chi2, int& ndf) const;

  /** fill LCIO Track State with parameters from helix and cov matrix
   */
  void ToLCIOTrackState(const THelicalTrack& helix, const TMatrixD& cov, edm4hep::TrackState& ts, double& chi2,
                        int& ndf) const;

  /** get the measurement site associated with the given lcio TrackerHit trkhit
   */
  int getSiteFromLCIOHit(const edm4hep::TrackerHit* trkhit, TKalTrackSite*& site) const;

  /** helper function to restrict the range of the azimuthal angle to ]-pi,pi]*/
  inline double toBaseRange(double phi) const {
    phi = std::fmod(phi + M_PI, 2. * M_PI);
    if (phi < 0)
      phi += M_PI;
    else
      phi -= M_PI;
    return phi;
  }

  std::unique_ptr<TKalTrack> m_kaltrack{}; // unique ptr to be able to forward declare
  TObjArray* m_kalhits = nullptr;
  GaudiDDKalTest* m_ktest = nullptr;
  // used to store whether initial track state has been supplied or created
  bool m_initialised = false;
  // used to store the fit direction supplied to intialise
  bool m_fitDirection = false;
  // used to store whether smoothing has been performed
  bool m_smoothed = false;
  const edm4hep::TrackerHit* m_trackHitAtPositiveNDF = nullptr;
  int m_hitIndexAtPositiveNDF = -1;

  // map to store relation between lcio hits and measurement sites
  std::map<const edm4hep::TrackerHit*, TKalTrackSite*> m_hit_used_for_sites{};

  // map to store relation between lcio hits kaltest hits
  std::map<DDVTrackHit*, const edm4hep::TrackerHit*> m_kaltest_hits_to_edm4hep_hits{};
  std::shared_ptr<std::map<const edm4hep::TrackerHit*, DDVTrackHit*>> m_edm4hep_hits_to_kaltest_hits{};

  // vector to store lcio hits rejected for measurement sites
  std::vector<const edm4hep::TrackerHit*> m_hit_not_used_for_sites{};

  // vector to store the chi-sqaure increment for measurement sites
  std::vector<std::pair<const edm4hep::TrackerHit*, double>> m_hit_chi2_values{};

  // vector to store the chi-sqaure increment for measurement sites
  std::vector<std::pair<const edm4hep::TrackerHit*, double>> m_outlier_chi2_values{};

  // Originally not present in MarlinDDKalTest, this is needed to be able to use
  // logging from Gaudi
  const Gaudi::Algorithm* m_thisAlg;
};

#endif
