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
#include "GaudiDDKalTestTrack.h"

#include <IMPL/TrackerHitPlaneImpl.h>

#include <kaltest/TKalFilterCond.h>
#include <kaltest/TKalTrack.h>
#include <kaltest/TKalTrackSite.h>

// Remove when the warnings are fixed in DDKalTest
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
#include <DDKalTest/DDCylinderHit.h>
#include <DDKalTest/DDCylinderMeasLayer.h>
#include <DDKalTest/DDPlanarHit.h>
#include <DDKalTest/DDVMeasLayer.h>
#include <DDKalTest/DDVTrackHit.h>
#pragma GCC diagnostic pop

#include <edm4hep/TrackerHitPlane.h>

#include <Gaudi/Algorithm.h>

#include <stdexcept>

const int site_fails_chi2_cut = 3;
const int no_intersection = 4;

/** Helper class for defining a filter condition based on the delta chi2 in the AddAndFilter step.
 */
class KalTrackFilter : public TKalFilterCond {
public:
  /** C'tor - takes as optional argument the maximum allowed delta chi2 for adding the hit (in IsAccepted() )
   */
  KalTrackFilter(double maxDeltaChi2 = DBL_MAX) : m_maxDeltaChi2(maxDeltaChi2), m_passed_last_filter_step(true) {}

  virtual Bool_t IsAccepted(const TKalTrackSite& site) override {
    double deltaChi2 = site.GetDeltaChi2();

    // m_thisAlg->debug() << " KalTrackFilter::IsAccepted called  !  deltaChi2 = " << std::scientific << deltaChi2
    //                    << " m_maxDeltaChi2 = " << m_maxDeltaChi2 << endmsg;

    m_passed_last_filter_step = deltaChi2 < m_maxDeltaChi2;

    return (m_passed_last_filter_step);
  }

  void resetFilterStatus() { m_passed_last_filter_step = true; }
  bool passedLastFilterStep() const { return m_passed_last_filter_step; }

protected:
  double m_maxDeltaChi2;
  bool m_passed_last_filter_step;
};

// //---------------------------------------------------------------------------------------------------------------
GaudiDDKalTestTrack::GaudiDDKalTestTrack(
    const Gaudi::Algorithm* algorithm, GaudiDDKalTest* ktest,
    std::shared_ptr<std::map<const edm4hep::TrackerHitPlane*, DDVTrackHit*>> edm4hep_hits_to_kaltest_hits)
    : m_ktest(ktest), m_edm4hep_hits_to_kaltest_hits(edm4hep_hits_to_kaltest_hits), m_thisAlg(algorithm) {
  m_kaltrack.SetOwner();

  m_kalhits = new TObjArray();
  m_kalhits->SetOwner();

  m_initialised = false;
  m_fitDirection = false;
  m_smoothed = false;

  m_trackHitAtPositiveNDF = nullptr;
  m_hitIndexAtPositiveNDF = 0;

  if (!m_edm4hep_hits_to_kaltest_hits) {
    m_edm4hep_hits_to_kaltest_hits = std::make_shared<std::map<const edm4hep::TrackerHitPlane*, DDVTrackHit*>>();
  }
}

// GaudiDDKalTestTrack::~GaudiDDKalTestTrack() {
//   delete m_kalhits;
// }

int GaudiDDKalTestTrack::addHit(const edm4hep::TrackerHitPlane* trkhit) {
  return this->addHit(trkhit, m_ktest->findMeasLayer(*trkhit));
}

int GaudiDDKalTestTrack::addHit(const edm4hep::TrackerHitPlane* trkhit, const DDVMeasLayer* ml) {
  m_thisAlg->debug() << "GaudiDDKalTestTrack::addHit: trkhit = " << trkhit->id() << " addr: " << trkhit
                     << " ml = " << ml << endmsg;

  if (trkhit && ml) {
    // TODO: a LCIO hit has to be created because it's needed downstream
    auto hit = IMPL::TrackerHitPlaneImpl();
    double pos[3] = {trkhit->getPosition()[0], trkhit->getPosition()[1], trkhit->getPosition()[2]};
    hit.setPosition(pos);
    hit.setCellID0(trkhit->getCellID());
    hit.setdU(trkhit->getDu());
    hit.setdV(trkhit->getDv());

    // Not needed
    // hit.setCovMatrix(lcioCov);
    // static_cast<IMPL::TrackerHitImpl*>(static_cast<EVENT::TrackerHit*>(hit))->setCovMatrix(lcioCov);
    // hit.setQuality(trkhit->getQuality());
    // hit.setType(trkhit->getType());
    // hit.setEDep(trkhit->getEDep());
    // hit.setEDepError(trkhit->getEDepError());
    // float u[2] = {trkhit->getU()[0], trkhit->getU()[1]};
    // hit.setU(u);
    // float v[2] = {trkhit->getV()[0], trkhit->getV()[1]};
    // hit.setV(v);
    // hit.setTime(trkhit->getTime());

    auto* kalhit = ml->ConvertLCIOTrkHit(&hit);
    return this->addHit(trkhit, kalhit, ml);
  } else {
    m_thisAlg->warning() << " GaudiDDKalTestTrack::addHit - bad inputs " << trkhit << " ml : " << ml << endmsg;
    return 1;
  }
}

int GaudiDDKalTestTrack::addHit(const edm4hep::TrackerHitPlane* trkhit, DDVTrackHit* kalhit, const DDVMeasLayer* ml) {
  if (kalhit && ml) {
    m_kalhits->Add(kalhit);                             // Add hit and set surface found
    m_kaltest_hits_to_edm4hep_hits[kalhit] = trkhit;    // add hit to map relating lcio and kaltest hits
    (*m_edm4hep_hits_to_kaltest_hits)[trkhit] = kalhit; // add hit to map relating lcio and kaltest hits
  } else {
    delete kalhit;
    return 1;
  }

  m_thisAlg->debug() << "GaudiDDKalTestTrack::addHit: hit added "
                     << "number of hits for track = " << m_kalhits->GetEntries() << endmsg;

  return 0;
}

int GaudiDDKalTestTrack::initialise(const edm4hep::TrackState& ts, bool fitDirection) {
  if (m_kalhits->GetEntries() == 0) {
    m_thisAlg->error() << "<<<<<< GaudiDDKalTestTrack::Initialise: Number of Hits is Zero. Cannot Initialise >>>>>>>"
                       << endmsg;
    return 1;
  }

  // SJA:FIXME: check here if the track is already initialised, and for now don't allow it to be re-initialised
  //            if the track is going to be re-initialised then we would need to do it directly on the first site
  if (m_initialised) {
    throw std::runtime_error("Track fit already initialised");
  }

  m_thisAlg->debug() << "GaudiDDKalTestTrack::initialise using TrackState: track parameters used for init : "
                     << "\t D0 " << ts.D0 << "\t Phi :" << ts.phi << "\t Omega " << ts.omega << "\t Z0 " << ts.Z0
                     << "\t tan(Lambda) " << ts.tanLambda

                     << "\t pivot : [" << ts.referencePoint[0] << ", " << ts.referencePoint[1] << ", "
                     << ts.referencePoint[2] << " - r: "
                     << std::sqrt(ts.referencePoint[0] * ts.referencePoint[0] +
                                  ts.referencePoint[1] * ts.referencePoint[1])
                     << "]" << endmsg;

  m_fitDirection = fitDirection;

  // get Bz from first hit
  TVTrackHit& h1 = *dynamic_cast<TVTrackHit*>(m_kalhits->At(0));
  double Bz = h1.GetBfield();

  // for GeV, Tesla, R in mm
  double alpha = Bz * 2.99792458E-4;

  double kappa = (Bz == 0.0 ? DBL_MAX : ts.omega / alpha);

  THelicalTrack helix(-ts.D0, toBaseRange(ts.phi - M_PI / 2.), kappa, ts.Z0, ts.tanLambda, ts.referencePoint[0],
                      ts.referencePoint[1], ts.referencePoint[2], Bz);

  TMatrixD cov(kSdim, kSdim);

  cov(0, 0) = ts.covMatrix[0];          //   d0, d0
  cov(0, 1) = -ts.covMatrix[1];         //   d0, phi
  cov(0, 2) = -ts.covMatrix[3] / alpha; //   d0, kappa
  cov(0, 3) = -ts.covMatrix[6];         //   d0, z0
  cov(0, 4) = -ts.covMatrix[10];        //   d0, tanl

  cov(1, 0) = -ts.covMatrix[1];        //   phi, d0
  cov(1, 1) = ts.covMatrix[2];         //   phi, phi
  cov(1, 2) = ts.covMatrix[4] / alpha; //   phi, kappa
  cov(1, 3) = ts.covMatrix[7];         //   phi, z0
  cov(1, 4) = ts.covMatrix[11];        //   tanl, phi

  cov(2, 0) = -ts.covMatrix[3] / alpha;          //   kappa, d0
  cov(2, 1) = ts.covMatrix[4] / alpha;           //   kappa, phi
  cov(2, 2) = ts.covMatrix[5] / (alpha * alpha); //   kappa, kappa
  cov(2, 3) = ts.covMatrix[8] / alpha;           //   kappa, z0
  cov(2, 4) = ts.covMatrix[12] / alpha;          //   kappa, tanl

  cov(3, 0) = -ts.covMatrix[6];        //   z0, d0
  cov(3, 1) = ts.covMatrix[7];         //   z0, phi
  cov(3, 2) = ts.covMatrix[8] / alpha; //   z0, kappa
  cov(3, 3) = ts.covMatrix[9];         //   z0, z0
  cov(3, 4) = ts.covMatrix[13];        //   z0, tanl

  cov(4, 0) = -ts.covMatrix[10];        //   tanl, d0
  cov(4, 1) = ts.covMatrix[11];         //   tanl, phi
  cov(4, 2) = ts.covMatrix[12] / alpha; //   tanl, kappa
  cov(4, 3) = ts.covMatrix[13];         //   tanl, z0
  cov(4, 4) = ts.covMatrix[14];         //   tanl, tanl

  //    cov.Print();

  // move the helix to either the position of the last hit or the first depending on initalise_at_end

  // default case initalise_at_end
  int index = m_kalhits->GetEntries() - 1;
  // or initialise at start
  if (m_fitDirection == true) {
    index = 0;
  }

  TVTrackHit* kalhit = dynamic_cast<TVTrackHit*>(m_kalhits->At(index));

  double dphi;

  TVector3 initial_pivot;

  // Leave the pivot at the origin for a 1-dim hit
  if (kalhit->GetDimension() > 1) {
    initial_pivot = kalhit->GetMeasLayer().HitToXv(*kalhit);
  } else {
    initial_pivot = TVector3(0.0, 0.0, 0.0);
  }

  // ---------------------------
  //  Create an initial start site for the track using the  hit
  // ---------------------------
  // set up a dummy hit needed to create initial site

  TVTrackHit* pDummyHit = nullptr;

  if ((pDummyHit = dynamic_cast<DDCylinderHit*>(kalhit))) {
    pDummyHit = (new DDCylinderHit(*static_cast<DDCylinderHit*>(kalhit)));

  } else if ((pDummyHit = dynamic_cast<DDPlanarHit*>(kalhit))) {
    pDummyHit = (new DDPlanarHit(*static_cast<DDPlanarHit*>(kalhit)));
    //}
    // else if ( (pDummyHit = dynamic_cast<DDPlanarStripHit *>( kalhit )) ) {
    // pDummyHit = (new DDPlanarStripHit(*static_cast<DDPlanarStripHit*>( kalhit )));
    if (pDummyHit->GetDimension() == 1) {
      const TVMeasLayer* ml = &pDummyHit->GetMeasLayer();

      const TVSurface* surf = dynamic_cast<const TVSurface*>(ml);

      if (surf) {
        double phi;

        surf->CalcXingPointWith(helix, initial_pivot, phi);
        m_thisAlg->debug() << "  GaudiDDKalTestTrack::initialise - CalcXingPointWith called for 1d hit ... " << endmsg;

      } else {
        m_thisAlg->debug() << "<<<<<<<<< GaudiDDKalTestTrack::initialise: dynamic_cast failed for TVSurface  >>>>>>>"
                           << endmsg;
        return 1;
      }
    }
  } else {
    m_thisAlg->warning() << "<<<<<<<<< GaudiDDKalTestTrack::initialise: dynamic_cast failed for hit type >>>>>>>"
                         << endmsg;
    return 1;
  }

  TVTrackHit& dummyHit = *pDummyHit;

  // SJA:FIXME: this constants should go in a header file
  //  give the dummy hit huge errors so that it does not contribute to the fit
  dummyHit(0, 1) = 1.e16; // give a huge error to d

  if (dummyHit.GetDimension() > 1)
    dummyHit(1, 1) = 1.e16; // give a huge error to z

  // use dummy hit to create initial site
  TKalTrackSite& initialSite = *new TKalTrackSite(dummyHit);

  initialSite.SetHitOwner(); // site owns hit
  initialSite.SetOwner();    // site owns states

  // ---------------------------
  //  Set up initial track state
  // ---------------------------

  helix.MoveTo(initial_pivot, dphi, nullptr, &cov);

  static TKalMatrix initialState(kSdim, 1);
  initialState(0, 0) = helix.GetDrho();      // d0
  initialState(1, 0) = helix.GetPhi0();      // phi0
  initialState(2, 0) = helix.GetKappa();     // kappa
  initialState(3, 0) = helix.GetDz();        // dz
  initialState(4, 0) = helix.GetTanLambda(); // tan(lambda)
  if (kSdim == 6)
    initialState(5, 0) = 0.; // t0

  // make sure that the pivot is in the right place
  initialSite.SetPivot(initial_pivot);

  // ---------------------------
  //  Set up initial Covariance Matrix
  // ---------------------------

  TKalMatrix covK(5, 5);
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      covK[i][j] = cov[i][j];
    }
  }
  if (kSdim == 6)
    covK(5, 5) = 1.e6; // t0

  //    covK.Print();

  // Add initial states to the site
  initialSite.Add(new TKalTrackState(initialState, covK, initialSite, TVKalSite::kPredicted));
  initialSite.Add(new TKalTrackState(initialState, covK, initialSite, TVKalSite::kFiltered));

  // add the initial site to the track: that is, give the track initial parameters and covariance
  // matrix at the starting measurement layer
  m_kaltrack.Add(&initialSite);

  m_initialised = true;

  return 0;
}

int GaudiDDKalTestTrack::addAndFit(DDVTrackHit* kalhit, double& chi2increment, TKalTrackSite*& site,
                                   double maxChi2Increment) {
  m_thisAlg->debug() << "GaudiDDKalTestTrack::addAndFit called : maxChi2Increment = " << std::scientific
                     << maxChi2Increment << endmsg;

  if (!m_initialised) {
    throw std::runtime_error("Track fit not initialised");
  }

  const auto* ml = dynamic_cast<const DDVMeasLayer*>(&(kalhit->GetMeasLayer()));

  if (m_thisAlg->msgLevel(MSG::DEBUG)) {
    m_thisAlg->debug() << "Kaltrack::addAndFit :  add site to track at index : " << ml->GetIndex() << " for type "
                       << ml->GetName();
    m_thisAlg->debug() << " with CellIDs:";
    for (unsigned int i = 0; i < ml->getNCellIDs(); ++i) {
      m_thisAlg->debug() << " : " << ml->getCellIDs()[i];
    }
    m_thisAlg->debug() << endmsg;
  }

  TKalTrackSite* temp_site = new TKalTrackSite(*kalhit); // create new site for this hit

  KalTrackFilter filter(maxChi2Increment);
  filter.resetFilterStatus();

  temp_site->SetFilterCond(&filter);

  // this is the only point at which a hit is actually filtered
  // and it is here that we can get the GetDeltaChi2 vs the maxChi2Increment
  // it will always be possible to get the delta chi2 so long as we have a link to the sites ...
  // although calling smooth will natrually update delta chi2.

  if (!m_kaltrack.AddAndFilter(*temp_site)) {
    chi2increment = temp_site->GetDeltaChi2();
    // get the measurement layer of the current hit
    TVector3 pos = ml->HitToXv(*kalhit);
    m_thisAlg->debug() << "Kaltrack::addAndFit : site discarded! at index : " << ml->GetIndex() << " for type "
                       << ml->GetName() << " chi2increment = " << chi2increment
                       << " maxChi2Increment = " << maxChi2Increment << " x = " << pos.x() << " y = " << pos.y()
                       << " z = " << pos.z() << " with CellIDs: " << endmsg;

    for (unsigned int i = 0; i < (dynamic_cast<const DDVMeasLayer*>(&(kalhit->GetMeasLayer()))->getNCellIDs()); ++i) {
      m_thisAlg->debug() << " CellID = "
                         << dynamic_cast<const DDVMeasLayer*>(&(kalhit->GetMeasLayer()))->getCellIDs()[i] << endmsg;
    }

    delete temp_site; // delete site if filter step failed

    // and this also works..
    m_thisAlg->debug() << " addAndFit : Site passed Chi2 filter step ? " << filter.passedLastFilterStep() << endmsg;

    if (filter.passedLastFilterStep() == false) {
      return site_fails_chi2_cut;
    } else {
      return 1;
    }
  }

  site = temp_site;
  chi2increment = site->GetDeltaChi2();

  return 0;
}

int GaudiDDKalTestTrack::addAndFit(const edm4hep::TrackerHitPlane* trkhit, double& chi2increment,
                                   double maxChi2Increment) {
  if (!trkhit) {
    m_thisAlg->debug() << "GaudiDDKalTestTrack::addAndFit( EVENT::TrackerHit* trkhit, double& chi2increment, "
                          "double maxChi2Increment): trkhit == 0"
                       << endmsg;
    return 1;
  }

  const DDVMeasLayer* ml = m_ktest->findMeasLayer(*trkhit);

  if (!ml) {
    // fg: not sure if ml should ever be 0 - but it seems to happen,
    //     if point is not on surface and more than one surface exists ...

    m_thisAlg->warning() << ">>>>>>>>>>>  no measurment layer found for trkhit cellid0 : " << trkhit->getCellID()
                         << " at " << edm4hep::Vector3d(trkhit->getPosition()) << endmsg;

    return 1;
  }

  auto* kalhit = (*m_edm4hep_hits_to_kaltest_hits)[trkhit];

  if (!kalhit) { // fg: ml->ConvertLCIOTrkHit returns 0 if hit not on surface !!!
    return 1;
  }

  TKalTrackSite* site = nullptr;
  int error_code = this->addAndFit(kalhit, chi2increment, site, maxChi2Increment);

  if (error_code != 0) {
    delete kalhit;

    // if the hit fails for any reason other than the Chi2 cut record the Chi2 contibution as DBL_MAX
    if (error_code != site_fails_chi2_cut) {
      chi2increment = DBL_MAX;
    }

    m_outlier_chi2_values.push_back(std::make_pair(trkhit, chi2increment));

    m_thisAlg->debug() << ">>>>>>>>>>>  addAndFit Number of Outliers : " << m_outlier_chi2_values.size() << endmsg;

    return error_code;
  } else {
    this->addHit(trkhit, kalhit, ml);
    m_hit_used_for_sites[trkhit] = site;
    m_hit_chi2_values.push_back(std::make_pair(trkhit, chi2increment));
  }

  // set the values for the point at which the fit becomes constained
  if (!m_trackHitAtPositiveNDF && m_kaltrack.GetNDF() >= 0) {
    m_trackHitAtPositiveNDF = trkhit;
    m_hitIndexAtPositiveNDF = m_kaltrack.IndexOf(site);

    m_thisAlg->debug() << ">>>>>>>>>>>  Fit is now constrained at : " << trkhit->getCellID() << " pos "
                       << edm4hep::Vector3d(trkhit->getPosition()) << " trkhit = " << m_trackHitAtPositiveNDF
                       << " index of kalhit = " << m_hitIndexAtPositiveNDF << " NDF = " << m_kaltrack.GetNDF()
                       << endmsg;
  }

  return 0;
}

int GaudiDDKalTestTrack::fit(double maxChi2Increment) {
  // SJA:FIXME: what do we do about calling fit after we have already added hits and filtered
  // I guess this would created new sites when addAndFit is called
  // one option would be to remove the sites
  // need to check where the sites are stored ...  probably in the KalTrackSystem
  //
  m_thisAlg->debug() << "GaudiDDKalTestTrack::fit() called " << endmsg;

  if (!m_initialised) {
    throw std::runtime_error("Track fit not initialised");
  }

  // ---------------------------
  //  Prepare hit iterrator for adding hits to kaltrack
  // ---------------------------

  TIter next(m_kalhits, m_fitDirection);

  // ---------------------------
  //  Start Kalman Filter
  // ---------------------------

  DDVTrackHit* kalhit = nullptr;

  while ((kalhit = dynamic_cast<DDVTrackHit*>(next()))) {
    double chi2increment;
    TKalTrackSite* site = nullptr;
    int error_code = this->addAndFit(kalhit, chi2increment, site, maxChi2Increment);

    const edm4hep::TrackerHitPlane* trkhit = m_kaltest_hits_to_edm4hep_hits[kalhit];

    if (error_code == 0) { // add trkhit to map associating trkhits and sites
      m_hit_used_for_sites[trkhit] = site;
      m_hit_chi2_values.push_back(std::make_pair(trkhit, chi2increment));

      // set the values for the point at which the fit becomes constained
      if (!m_trackHitAtPositiveNDF && m_kaltrack.GetNDF() >= 0) {
        m_trackHitAtPositiveNDF = trkhit;
        m_hitIndexAtPositiveNDF = m_kaltrack.IndexOf(site);

        m_thisAlg->debug() << ">>>>>>>>>>>  Fit is now constrained at : " << trkhit->getCellID() << " pos "
                           << trkhit->getPosition() << " trkhit = " << m_trackHitAtPositiveNDF
                           << " index of kalhit = " << m_hitIndexAtPositiveNDF << " NDF = " << m_kaltrack.GetNDF()
                           << endmsg;
      }

    } else { // hit rejected by the filter, so store in the list of rejected hits

      // if the hit fails for any reason other than the Chi2 cut record the Chi2 contibution as DBL_MAX
      if (error_code != site_fails_chi2_cut) {
        chi2increment = DBL_MAX;
      }

      m_outlier_chi2_values.push_back(std::make_pair(trkhit, chi2increment));
      m_thisAlg->debug() << ">>>>>>>>>>>  fit(): Number of Outliers : " << m_outlier_chi2_values.size() << endmsg;

      m_hit_not_used_for_sites.push_back(trkhit);
    }

  } // end of Kalman filter

  // if (m_ktest->getOption(MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing)) {
  //   m_thisAlg->debug() << "Perform Smoothing for All Previous Measurement Sites " << endmsg;
  //   int error = this->smooth();

  //   if (error != 0)
  //     return error;
  // }

  // return m_hit_used_for_sites.empty() == false ? success : all_sites_fail_fit ;
  if (m_hit_used_for_sites.empty() == false) {
    return 0;
  } else {
    return 1;
  }
}

/** smooth track states from the last filtered hit back to the measurement site associated with the given hit
 */
int GaudiDDKalTestTrack::smooth(const edm4hep::TrackerHitPlane* trkhit) {
  m_thisAlg->debug() << "GaudiDDKalTestTrack::smooth( EVENT::TrackerHit* " << trkhit << "  ) " << endmsg;

  if (!trkhit) {
    return 1;
  }

  std::map<EVENT::TrackerHit*, TKalTrackSite*>::const_iterator it;

  TKalTrackSite* site = nullptr;
  int error_code = getSiteFromLCIOHit(trkhit, site);

  if (error_code != 0)
    return error_code;

  int index = m_kaltrack.IndexOf(site);

  m_kaltrack.SmoothBackTo(index);

  m_smoothed = true;

  return 0;
}

int GaudiDDKalTestTrack::getTrackState(const edm4hep::TrackerHitPlane* trkhit, edm4hep::TrackState& ts, double& chi2,
                                       int& ndf) const {
  m_thisAlg->debug()
      << "GaudiDDKalTestTrack::getTrackState( EVENT::TrackerHit* trkhit, IMPL::TrackStateImpl& ts ) using hit: "
      << trkhit << " with cellID0 = " << trkhit->getCellID() << endmsg;

  TKalTrackSite* site = nullptr;
  int error_code = getSiteFromLCIOHit(trkhit, site);

  if (error_code != 0)
    return error_code;

  m_thisAlg->debug() << "GaudiDDKalTestTrack::getTrackState: site " << site << endmsg;

  this->ToLCIOTrackState(*site, ts, chi2, ndf);

  return 0;
}

const std::vector<std::pair<const edm4hep::TrackerHitPlane*, double>>& GaudiDDKalTestTrack::getHitsInFit() const {
  return m_hit_chi2_values;

  // this needs more thought. What about when the hits are added using addAndFit?

  // need to check the order so that we can return the list ordered in time
  // as they will be added to m_hit_chi2_values in the order of fitting
  // not in the order of time
  //
  //    if( m_fitDirection == IMarlinTrack::backward ){
  //      std::reverse_copy( m_hit_chi2_values.begin() , m_hit_chi2_values.end() , std::back_inserter(  hits  )  ) ;
  //    } else {
  //      std::copy( m_hit_chi2_values.begin() , m_hit_chi2_values.end() , std::back_inserter(  hits  )  ) ;
  //    }
}

const std::vector<std::pair<const edm4hep::TrackerHitPlane*, double>>& GaudiDDKalTestTrack::getOutliers() const {
  return m_outlier_chi2_values;
  // this needs more thought. What about when the hits are added using addAndFit?
  //    // need to check the order so that we can return the list ordered in time
  //    // as they will be added to m_hit_chi2_values in the order of fitting
  //    // not in the order of time
  //
  //    if( m_fitDirection == IMarlinTrack::backward ){
  //      std::reverse_copy( m_outlier_chi2_values.begin() , m_outlier_chi2_values.end() , std::back_inserter(  hits  )
  //      ) ;
  //    } else {
  //      std::copy( m_outlier_chi2_values.begin() , m_outlier_chi2_values.end() , std::back_inserter(  hits  )  ) ;
  //    }
}

int GaudiDDKalTestTrack::getNDF() const {
  if (!m_initialised) {
    throw std::runtime_error("GaudiDDKalTestTrack is not initialised");
  }
  return const_cast<TKalTrack&>(m_kaltrack).GetNDF();
}

const edm4hep::TrackerHitPlane* GaudiDDKalTestTrack::getTrackerHitAtPositiveNDF() const {
  return m_trackHitAtPositiveNDF;
}

int GaudiDDKalTestTrack::propagate(const edm4hep::Vector3d& point, const edm4hep::TrackerHitPlane* trkhit,
                                   edm4hep::TrackState& ts, double& chi2, int& ndf) {
  TKalTrackSite* site = nullptr;
  int error_code = getSiteFromLCIOHit(trkhit, site);

  if (error_code != 0)
    return error_code;

  // check if the point is inside the beampipe
  // SJA:FIXME: really we should also check if the PCA to the point is also less than R

  const DDVMeasLayer* ml = m_ktest->getIPLayer();

  if (ml && std::sqrt(point.x * point.x + point.y * point.y + point.z * point.z) > m_ktest->getIPLayer()->GetR())
    ml = nullptr;

  return this->propagate(point, *site, ts, chi2, ndf, ml);
}

int GaudiDDKalTestTrack::propagate(const edm4hep::Vector3d& point, const TKalTrackSite& site, edm4hep::TrackState& ts,
                                   double& chi2, int& ndf, const DDVMeasLayer* ml) {
  m_thisAlg->debug() << "GaudiDDKalTestTrack::propagate( const Vector3D& point, const TKalTrackSite& site, "
                        "IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) called "
                     << endmsg;

  const TVector3 tpoint(point.x, point.y, point.z);

  TKalTrackState& trkState = static_cast<TKalTrackState&>(site.GetCurState()); // this segfaults if no hits are present

  THelicalTrack helix = trkState.GetHelix();
  double dPhi = 0.0;

  Int_t sdim = trkState.GetDimension(); // dimensions of the track state, it will be 5 or 6
  TKalMatrix sv(sdim, 1);

  TKalMatrix F(sdim, sdim); // propagator matrix to be returned by transport function
  F.UnitMatrix();           // set the propagator matrix to the unit matrix

  TKalMatrix Q(sdim, sdim); // noise matrix to be returned by transport function
  Q.Zero();
  TVector3 x0; // intersection point to be returned by transport

  TMatrixD c0(trkState.GetCovMat());

  // the last layer crossed by the track before point
  if (!ml) {
    ml = m_ktest->getLastMeasLayer(helix, tpoint);
  }

  if (ml) {
    m_ktest->m_det.Transport(site, *ml, x0, sv, F, Q); // transport to last layer cross before point

    // given that we are sure to have intersected the layer ml as this was provided via getLastMeasLayer, x0 will lie on
    // the layer this could be checked with the method isOnSurface so F will be the propagation matrix from the current
    // location to the last surface and Q will be the noise matrix up to this point

    TKalMatrix Ft = TKalMatrix(TMatrixD::kTransposed, F);
    c0 = F * c0 * Ft + Q; // update covaraince matrix and add the MS assosiated with moving to tvml

    helix.MoveTo(x0, dPhi, nullptr, nullptr); // move the helix to tvml

  } else { // the current site is at the last surface before the point to propagate to
    ml = dynamic_cast<const DDVMeasLayer*>(&(site.GetHit().GetMeasLayer()));
  }

  // get whether the track is incomming or outgoing at the last surface
  const TVSurface* sfp = dynamic_cast<const TVSurface*>(ml); // last surface

  TMatrixD dxdphi = helix.CalcDxDphi(0);                      // tangent vector at last surface
  TVector3 dxdphiv(dxdphi(0, 0), dxdphi(1, 0), dxdphi(2, 0)); // convert matirix diagonal to vector
  //    Double_t cpa = helix.GetKappa();                              // get pt

  Bool_t isout = -dPhi * dxdphiv.Dot(sfp->GetOutwardNormal(x0)) < 0
                     ? kTRUE
                     : kFALSE; // out-going or in-coming at the destination surface

  // now move to the point
  TKalMatrix DF(sdim, sdim);
  DF.UnitMatrix();
  helix.MoveTo(tpoint, dPhi, &DF, nullptr); // move helix to desired point, and get propagator matrix

  TKalMatrix Qms(sdim, sdim);
  ml->CalcQms(isout, helix, dPhi, Qms); // calculate MS for the final step through the present material

  TKalMatrix DFt = TKalMatrix(TMatrixD::kTransposed, DF);
  c0 = DF * c0 * DFt + Qms; // update the covariance matrix

  this->ToLCIOTrackState(helix, c0, ts, chi2, ndf);

  return 0;
}

int GaudiDDKalTestTrack::propagateToLayer(int layerID, const edm4hep::TrackerHitPlane* trkhit, edm4hep::TrackState& ts,
                                          double& chi2, int& ndf, int& detElementID, int mode) {
  TKalTrackSite* site = nullptr;
  int error_code = getSiteFromLCIOHit(trkhit, site);

  if (error_code != 0)
    return error_code;

  return this->propagateToLayer(layerID, *site, ts, chi2, ndf, detElementID, mode);
}

int GaudiDDKalTestTrack::propagateToLayer(int layerID, const TKalTrackSite& site, edm4hep::TrackState& ts, double& chi2,
                                          int& ndf, int& detElementID, int mode) {
  m_thisAlg->debug() << "GaudiDDKalTestTrack::propagateToLayer( int layerID, const TKalTrackSite& site, "
                        "IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID ) called "
                     << endmsg;

  edm4hep::Vector3d crossing_point;
  const DDVMeasLayer* ml = nullptr;

  int error_code = this->intersectionWithLayer(layerID, site, crossing_point, detElementID, ml, mode);

  if (error_code != 0)
    return error_code;

  return this->propagate(crossing_point, site, ts, chi2, ndf, ml);
}

int GaudiDDKalTestTrack::intersectionWithLayer(int layerID, const TKalTrackSite& site, edm4hep::Vector3d& point,
                                               int& detElementID, const DDVMeasLayer*& ml, int mode) {
  m_thisAlg->debug() << "GaudiDDKalTestTrack::intersectionWithLayer( int layerID, const TKalTrackSite& site, "
                        "Vector3D& point, int& detElementID, int mode) called layerID = "
                     << layerID << endmsg;

  std::vector<const DDVMeasLayer*> meas_modules = m_ktest->getSensitiveMeasurementModulesForLayer(layerID);

  if (meas_modules.size() == 0) {
    m_thisAlg->debug() << "GaudiDDKalTestTrack::intersectionWithLayer layer id unknown: layerID = " << layerID
                       << endmsg;
    return no_intersection;
  }

  //  int index_of_intersected;
  int error_code = this->findIntersection(meas_modules, site, point, detElementID, ml, mode);

  if (error_code == 0) {
    m_thisAlg->debug() << "GaudiDDKalTestTrack::intersectionWithLayer intersection with layerID = " << layerID
                       << ": at x = " << point.x << " y = " << point.y << " z = " << point.z
                       << " detElementID = " << detElementID << " " << detElementID << endmsg;

  } else if (error_code == no_intersection) {
    ml = nullptr;
    m_thisAlg->debug() << "GaudiDDKalTestTrack::intersectionWithLayer No intersection with layerID = " << layerID << " "
                       << layerID << endmsg;
  }

  return error_code;
}

int GaudiDDKalTestTrack::findIntersection(const DDVMeasLayer& meas_module, const TKalTrackSite& site,
                                          edm4hep::Vector3d& point, double& dphi, int& detElementID, int mode) {
  TKalTrackState& trkState = static_cast<TKalTrackState&>(site.GetCurState());
  THelicalTrack helix = trkState.GetHelix();
  TVector3 xto; // reference point at destination to be returned by CalcXinPointWith
  int crossing_exist = meas_module.getIntersectionAndCellID(helix, xto, dphi, detElementID, mode);

  m_thisAlg->debug() << "GaudiDDKalTestTrack::intersectionWithLayer crossing_exist = " << crossing_exist << " dphi "
                     << dphi << " with detElementIDs: " << detElementID;
  m_thisAlg->debug() << endmsg;

  if (crossing_exist == 0) {
    return no_intersection;
  } else {
    point = {xto.X(), xto.Y(), xto.Z()};
  }

  return 0;
}

int GaudiDDKalTestTrack::findIntersection(std::vector<DDVMeasLayer const*>& meas_modules, const TKalTrackSite& site,
                                          edm4hep::Vector3d& point, int& detElementID, const DDVMeasLayer*& ml,
                                          int mode) {
  double dphi_min = DBL_MAX; // use to store the min deflection angle found so that can avoid the crossing on the far
                             // side of the layer
  bool surf_found = false;

  for (size_t i = 0; i < meas_modules.size(); ++i) {
    double dphi = 0;
    // need to send a temporary point as we may get the crossing point with the layer on the oposite side of the layer
    edm4hep::Vector3d point_temp;

    int temp_detElementID;

    int error_code = findIntersection(*meas_modules[i], site, point_temp, dphi, temp_detElementID, mode);

    if (error_code == 0) {
      // make sure we get the next crossing
      if (std::abs(dphi) < dphi_min) {
        dphi_min = std::abs(dphi);
        surf_found = true;
        ml = meas_modules[i];
        detElementID = temp_detElementID;
        point = point_temp;
      }

    } else if (error_code != no_intersection) { // in which case error_code is an error rather than simply a lack of
                                                // intersection, so return

      return error_code;
    }
  }

  // check if the surface was found and return accordingly
  if (surf_found) {
    return 0;
  } else {
    return no_intersection;
  }
}

void GaudiDDKalTestTrack::ToLCIOTrackState(const THelicalTrack& helix, const TMatrixD& cov, edm4hep::TrackState& ts,
                                           double& chi2, int& ndf) const {
  chi2 = const_cast<TKalTrack&>(m_kaltrack).GetChi2();
  ndf = const_cast<TKalTrack&>(m_kaltrack).GetNDF();

  //============== convert parameters to LCIO convention ====

  //  this is for incomming tracks ...
  // To preserve the results from iLCSoft, omega has to be kept as a double
  // storing it directly in ts.omega loses precision
  double omega = 1. / helix.GetRho();

  ts.D0 = -helix.GetDrho();
  ts.phi = toBaseRange(helix.GetPhi0() + M_PI / 2.); // fi0  - M_PI/2.  ) ;
  ts.omega = omega;
  ts.Z0 = helix.GetDz();
  ts.tanLambda = helix.GetTanLambda();

  double alpha = omega / helix.GetKappa(); // conversion factor for omega (1/R) to kappa (1/Pt)

  ts.covMatrix[0] = cov(0, 0); //   d0,   d0

  ts.covMatrix[1] = -cov(1, 0); //   phi0, d0
  ts.covMatrix[2] = cov(1, 1);  //   phi0, phi

  ts.covMatrix[3] = -cov(2, 0) * alpha;        //   omega, d0
  ts.covMatrix[4] = cov(2, 1) * alpha;         //   omega, phi
  ts.covMatrix[5] = cov(2, 2) * alpha * alpha; //   omega, omega

  ts.covMatrix[6] = -cov(3, 0);        //   z0  , d0
  ts.covMatrix[7] = cov(3, 1);         //   z0  , phi
  ts.covMatrix[8] = cov(3, 2) * alpha; //   z0  , omega
  ts.covMatrix[9] = cov(3, 3);         //   z0  , z0

  ts.covMatrix[10] = -cov(4, 0);        //   tanl, d0
  ts.covMatrix[11] = cov(4, 1);         //   tanl, phi
  ts.covMatrix[12] = cov(4, 2) * alpha; //   tanl, omega
  ts.covMatrix[13] = cov(4, 3);         //   tanl, z0
  ts.covMatrix[14] = cov(4, 4);         //   tanl, tanl

  const auto& pivot = helix.GetPivot();
  ts.referencePoint = {static_cast<float>(pivot.X()), static_cast<float>(pivot.Y()), static_cast<float>(pivot.Z())};

  m_thisAlg->debug() << " kaltest track parameters: "
                     << " chi2/ndf " << chi2 / ndf << " chi2 " << chi2 << " ndf " << ndf << " prob "
                     << TMath::Prob(chi2, ndf) << endmsg

                     << "\t D0 " << ts.D0 << "[+/-" << sqrt(ts.covMatrix[0]) << "] "
                     << "\t Phi :" << ts.phi << "[+/-" << sqrt(ts.covMatrix[2]) << "] "
                     << "\t Omega " << omega << "[+/-" << sqrt(ts.covMatrix[5]) << "] "
                     << "\t Z0 " << ts.Z0 << "[+/-" << sqrt(ts.covMatrix[9]) << "] "
                     << "\t tan(Lambda) " << ts.tanLambda << "[+/-" << sqrt(ts.covMatrix[14]) << "] "

                     << "\t pivot : [" << pivot[0] << ", " << pivot[1] << ", " << pivot[2]
                     << " - r: " << std::sqrt(pivot[0] * pivot[0] + pivot[1] * pivot[1]) << "]" << endmsg;
}

void GaudiDDKalTestTrack::ToLCIOTrackState(const TKalTrackSite& site, edm4hep::TrackState& ts, double& chi2,
                                           int& ndf) const {
  TKalTrackState& trkState =
      static_cast<TKalTrackState&>(site.GetCurState()); // GetCutState will return the last added state to this site
                                                        // Assuming everything has proceeded as expected
                                                        // this will be Predicted -> Filtered -> Smoothed

  // streamlog_out( DEBUG3 ) << " GaudiDDKalTestTrack::ToLCIOTrackState : " << endmsg ;
  // trkState.DebugPrint() ;

  const TMatrixD& c0(trkState.GetCovMat());

  this->ToLCIOTrackState(trkState.GetHelix(), c0, ts, chi2, ndf);
}

int GaudiDDKalTestTrack::getSiteFromLCIOHit(const edm4hep::TrackerHitPlane* trkhit, TKalTrackSite*& site) const {
  const auto it = m_hit_used_for_sites.find(trkhit);

  if (it == m_hit_used_for_sites.end()) { // hit not associated with any site

    bool found = false;

    for (size_t i = 0; i < m_hit_not_used_for_sites.size(); ++i) {
      trkhit = m_hit_not_used_for_sites[i];
      if (trkhit) {
        found = true;
        break;
      }
    }

    if (found) {
      m_thisAlg->debug() << "GaudiDDKalTestTrack::getSiteFromLCIOHit: hit was rejected during filtering" << endmsg;
      return 1;
    } else {
      m_thisAlg->debug() << "GaudiDDKalTestTrack::getSiteFromLCIOHit: hit " << trkhit << " not in list of supplied hits"
                         << endmsg;
      return 1;
    }
  }

  site = it->second;

  m_thisAlg->debug() << "GaudiDDKalTestTrack::getSiteFromLCIOHit: site " << site << " found for hit " << trkhit
                     << endmsg;
  return 0;
}
