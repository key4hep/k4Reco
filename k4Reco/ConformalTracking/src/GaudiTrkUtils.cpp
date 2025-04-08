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
#include "GaudiTrkUtils.h"

#include "GaudiDDKalTest.h"
#include "GaudiDDKalTestTrack.h"

#include "HelixTrack.h"

#include <edm4hep/MutableTrack.h>
#include <edm4hep/Track.h>
#include <edm4hep/TrackState.h>
#include <edm4hep/TrackerHitPlane.h>
#include <edm4hep/Vector3d.h>
#include <edm4hep/utils/vector_utils.h>

#include <DD4hep/BitFieldCoder.h>

#include <bitset>
#include <memory>
#include <vector>

// From MarlinTrkUtils
const int MIN_NDF = 6;
const int no_intersection = 4;

void setStreamlogOutputLevel(const Gaudi::Algorithm* thisAlg, streamlog::logscope* scope) {
  MSG::Level outputLevel = thisAlg->msgLevel();
  switch (outputLevel) {
  case MSG::ERROR:
    scope->setLevel<streamlog::ERROR>();
    break;
  case MSG::WARNING:
    scope->setLevel<streamlog::WARNING>();
    break;
  case MSG::INFO:
    scope->setLevel<streamlog::MESSAGE>();
    break;
  case MSG::DEBUG:
    scope->setLevel<streamlog::DEBUG>();
    break;
  default:
    scope->setLevel<streamlog::MESSAGE>();
    break;
  }
}

int GaudiTrkUtils::createFinalisedLCIOTrack(GaudiDDKalTestTrack& marlinTrk,
                                            const std::vector<const edm4hep::TrackerHitPlane*>& hit_list,
                                            edm4hep::MutableTrack& track, bool fit_direction,
                                            const edm4hep::CovMatrix6f& initial_cov_for_prefit, float bfield_z,
                                            double maxChi2Increment) {
  if (hit_list.empty())
    return 1;

  int return_error = 0;
  edm4hep::TrackState pre_fit(0, 0, 0, 0, 0, 0, 0, {}, initial_cov_for_prefit);
  return_error = createPrefit(hit_list, pre_fit, bfield_z);
  // m_thisAlg->info() << " **** createFinalisedLCIOTrack - created pre-fit: " <<  toString( &pre_fit )  << endmsg ;

  // use prefit parameters to produce Finalised track
  if (return_error == 0) {
    return_error = createFinalisedLCIOTrack(marlinTrk, hit_list, track, fit_direction, pre_fit, maxChi2Increment);
  } else {
    m_thisAlg->debug() << "MarlinTrk::createFinalisedLCIOTrack : Prefit failed error = " << return_error << endmsg;
  }

  return return_error;
}

int GaudiTrkUtils::createFinalisedLCIOTrack(GaudiDDKalTestTrack& marlinTrk,
                                            const std::vector<const edm4hep::TrackerHitPlane*>& hit_list,
                                            edm4hep::MutableTrack& track, bool fit_direction,
                                            edm4hep::TrackState& pre_fit, double maxChi2Increment) {
  if (hit_list.empty())
    return 1;

  int fit_status = createFit(hit_list, marlinTrk, pre_fit, fit_direction, maxChi2Increment);

  if (fit_status != 0) {
    m_thisAlg->info() << "MarlinTrk::createFinalisedLCIOTrack fit failed: fit_status = " << fit_status << endmsg;
    return 1;
  }

  int error = finaliseLCIOTrack(marlinTrk, track, hit_list, fit_direction);
  return error;
}

int GaudiTrkUtils::createPrefit(const std::vector<const edm4hep::TrackerHitPlane*>& hit_list,
                                edm4hep::TrackState& pre_fit, float bfield_z) {
  if (hit_list.empty())
    return 1;

  // loop over all the hits and create a list consisting only 2D hits

  std::vector<const edm4hep::TrackerHitPlane*> twoD_hits;

  for (unsigned ihit = 0; ihit < hit_list.size(); ++ihit) {
    // check if this a space point or 2D hit
    // TODO: check if bitset works
    if (std::bitset<32>(hit_list[ihit]->getType())[29] == false) {
      // then add to the list
      twoD_hits.push_back(hit_list[ihit]);
    }
  }

  ///////////////////////////////////////////////////////
  // check that there are enough 2-D hits to create a helix
  ///////////////////////////////////////////////////////

  if (twoD_hits.size() < 3) { // no chance to initialise print warning and return
    m_thisAlg->warning() << "MarlinTrk::createFinalisedLCIOTrack Cannot create helix from less than 3 2-D hits"
                         << endmsg;
    // TODO
    // return IMarlinTrack::bad_intputs;
    return 1;
  }

  ///////////////////////////////////////////////////////
  // make a helix from 3 hits to get a trackstate
  ///////////////////////////////////////////////////////

  // SJA:FIXME: this may not be the optimal 3 hits to take in certain cases where the 3 hits are not well spread over
  // the track length
  const edm4hep::Vector3d x1 = twoD_hits[0]->getPosition();
  const edm4hep::Vector3d x2 = twoD_hits[twoD_hits.size() / 2]->getPosition();
  const edm4hep::Vector3d x3 = twoD_hits.back()->getPosition();

  HelixTrack helixTrack(x1, x2, x3, bfield_z, HelixTrack::forwards);

  helixTrack.moveRefPoint(0.0, 0.0, 0.0);

  const float referencePoint[3] = {float(helixTrack.getRefPointX()), float(helixTrack.getRefPointY()),
                                   float(helixTrack.getRefPointZ())};

  pre_fit.D0 = helixTrack.getD0();
  pre_fit.phi = helixTrack.getPhi0();
  pre_fit.omega = helixTrack.getOmega();
  pre_fit.Z0 = helixTrack.getZ0();
  pre_fit.tanLambda = helixTrack.getTanLambda();

  pre_fit.referencePoint = referencePoint;

  return 0;
}

int GaudiTrkUtils::createFit(const std::vector<const edm4hep::TrackerHitPlane*>& hit_list,
                             GaudiDDKalTestTrack& marlinTrk, edm4hep::TrackState& pre_fit, bool fit_direction,
                             double maxChi2Increment) {
  if (hit_list.empty())
    return 1;

  int return_error = 0;

  ///////////////////////////////////////////////////////
  // add hits to IMarlinTrk
  ///////////////////////////////////////////////////////

  //  start by trying to add the hits to the track we want to finally use.
  m_thisAlg->debug() << "MarlinTrk::createFit Start Fit: AddHits: number of hits to fit " << hit_list.size() << endmsg;

  std::vector<const edm4hep::TrackerHitPlane*> added_hits;
  unsigned int ndof_added = 0;

  for (auto it = hit_list.begin(); it != hit_list.end(); ++it) {
    const edm4hep::TrackerHitPlane* trkHit = *it;
    bool isSuccessful = false;

    // TODO: check if bitset works
    if (std::bitset<32>(trkHit->getType())[30]) { // it is a composite spacepoint

      // TODO:
      // //Split it up and add both hits to the MarlinTrk
      // const EVENT::LCObjectVec rawObjects = trkHit->getRawHits();

      // for (unsigned k = 0; k < rawObjects.size(); k++) {
      //   EVENT::TrackerHit* rawHit = dynamic_cast<EVENT::TrackerHit*>(rawObjects[k]);

      //   if (marlinTrk->addHit(rawHit) == IMarlinTrack::success) {
      //     isSuccessful = true;  //if at least one hit from the spacepoint gets added
      //     ++ndof_added;
      //     streamlog_out(DEBUG4) << "MarlinTrk::createFit ndof_added = " << ndof_added << std::endl;
      //   }
      // }
    } else { // normal non composite hit

      if (marlinTrk.addHit(trkHit) == 0) {
        isSuccessful = true;
        ndof_added += 2;
        m_thisAlg->debug() << "MarlinTrk::createFit ndof_added = " << ndof_added << endmsg;
      }
    }

    if (isSuccessful) {
      added_hits.push_back(trkHit);
    } else {
      m_thisAlg->debug() << "Hit " << it - hit_list.begin() << " Dropped " << endmsg;
    }
  }

  if (ndof_added < MIN_NDF) {
    m_thisAlg->debug() << "MarlinTrk::createFit : Cannot fit less with less than " << MIN_NDF
                       << " degrees of freedom. Number of hits =  " << added_hits.size() << " ndof = " << ndof_added
                       << endmsg;
    return 1;
  }

  ///////////////////////////////////////////////////////
  // set the initial track parameters
  ///////////////////////////////////////////////////////

  return_error = marlinTrk.initialise(pre_fit, fit_direction); // IMarlinTrack::backward ) ;

  if (return_error != 0) {
    m_thisAlg->debug() << "MarlinTrk::createFit Initialisation of track fit failed with error : " << return_error
                       << endmsg;
    return return_error;
  }

  return marlinTrk.fit(maxChi2Increment);
}

int GaudiTrkUtils::finaliseLCIOTrack(GaudiDDKalTestTrack& marlintrk, edm4hep::MutableTrack& track,
                                     const std::vector<const edm4hep::TrackerHitPlane*>& hit_list, bool fit_direction) {
  int ndf = 0;
  double chi2 = -DBL_MAX;

  // First check NDF to see if it make any sense to continue.
  // The track will be dropped if the NDF is less than 0
  ndf = marlintrk.getNDF();

  if (ndf < 0) {
    m_thisAlg->warning()
        << "MarlinTrk::finaliseLCIOTrack: number of degrees of freedom less than 0 track dropped : NDF = " << ndf
        << endmsg;
    return 1;
  } else {
    m_thisAlg->debug() << "MarlinTrk::finaliseLCIOTrack: NDF = " << ndf << endmsg;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // get the list of hits used in the fit
  // add these to the track, add spacepoints as long as at least on strip hit is used.
  ////////////////////////////////////////////////////////////////////////////////////////////////////////

  std::vector<const edm4hep::TrackerHitPlane*> used_hits;

  const auto& hits_in_fit = marlintrk.getHitsInFit();
  const auto& outliers = marlintrk.getOutliers();

  ///////////////////////////////////////////////
  // now loop over the hits provided for fitting
  // we do this so that the hits are added in the
  // order in which they have been fitted
  ///////////////////////////////////////////////

  for (unsigned ihit = 0; ihit < hit_list.size(); ++ihit) {
    const edm4hep::TrackerHitPlane* trkHit = hit_list[ihit];
    // TODO
    // if( UTIL::BitSet32( trkHit.getType() )[ UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT ]   ){ //it is a composite
    // spacepoint
    if (std::bitset<32>(trkHit->getType())[30]) { // it is a composite spacepoint
      // // TODO:
      throw std::runtime_error("MarlinTrk::finaliseLCIOTrack: composite spacepoints not yet supported");
      // // const auto rawObjects = trkHit.getRawHits();
      // const auto rawObjects = std::vector<edm4hep::TrackerHit>{};
      // for (unsigned k = 0; k < rawObjects.size(); k++) {
      //   // TODO
      //   // edm4hep::TrackerHit* rawHit = dynamic_cast< edm4hep::TrackerHit* >( rawObjects[k] );
      //   edm4hep::TrackerHit* rawHit;
      //   bool                 is_outlier = false;
      //   // here we loop over outliers as this will be faster than looping over the used hits
      //   for (unsigned ohit = 0; ohit < outliers.size(); ++ohit) {
      //     if (rawHit == outliers[ohit].first) {
      //       is_outlier = true;
      //       break;  // break out of loop over outliers
      //     }
      //   }

      //   if (is_outlier == false) {
      //     // TODO
      //     // used_hits.push_back(hit_list[ihit]);
      //     // track.addHit(used_hits.back());
      //     break;  // break out of loop over rawObjects
      //   }
      // }
    } else {
      bool is_outlier = false;
      // here we loop over outliers as this will be faster than looping over the used hits
      for (unsigned ohit = 0; ohit < outliers.size(); ++ohit) {
        // TODO
        if (trkHit == outliers[ohit].first) {
          is_outlier = true;
          break; // break out of loop over outliers
        }
      }
      if (is_outlier == false) {
        used_hits.push_back(trkHit);
        track.addToTrackerHits(*used_hits.back());
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // We now need to find out at which point the fit is constrained
  // and therefore be able to provide well formed (pos. def.) cov. matrices
  ///////////////////////////////////////////////////////////////////////////

  edm4hep::TrackState trkStateAtFirstHit;
  const edm4hep::TrackerHitPlane* firstHit =
      fit_direction == false ? hits_in_fit.back().first : hits_in_fit.front().first;
  const edm4hep::TrackerHitPlane* lastHit =
      fit_direction == false ? hits_in_fit.front().first : hits_in_fit.back().first;
  const edm4hep::TrackerHitPlane* last_constrained_hit = marlintrk.getTrackerHitAtPositiveNDF();

  // m_thisAlg->debug() << "MarlinTrk::finaliseLCIOTrack: firstHit : " << toString( firstHit )
  //     		  << " lastHit:                                " << toString( lastHit )
  //     		  << " last constrained hit:                   " << toString( last_constrained_hit )
  //     		  << " fit direction is forward : " << fit_direction << endmsg ;

  int return_error = marlintrk.smooth(last_constrained_hit);

  m_thisAlg->debug() << "MarlinTrk::finaliseLCIOTrack: return_code for smoothing to last constrained hit "
                     << last_constrained_hit << " = " << return_error << " NDF = " << ndf << endmsg;

  if (return_error != 0) {
    return return_error;
  }

  ///////////////////////////////////////////////////////
  // first create trackstate at IP
  ///////////////////////////////////////////////////////
  const edm4hep::Vector3d point(0., 0., 0.); // nominal IP

  edm4hep::TrackState trkStateIP;

  // m_thisAlg->debug() << "MarlinTrk::finaliseLCIOTrack: finalised kaltest track  : "
  //     		  << marlintrk.toString() << endmsg ;

  ///////////////////////////////////////////////////////
  // make sure that the track state can be propagated to the IP
  ///////////////////////////////////////////////////////

  // TODO
  // bool usingAidaTT = ( trksystem->name() == "AidaTT" ) ;
  bool usingAidaTT = false;

  // if we fitted backwards, the firstHit is the last one used in the fit and we simply propagate to the IP:

  if (fit_direction == false || usingAidaTT) {
    return_error = marlintrk.propagate(point, firstHit, trkStateIP, chi2, ndf);
  } else {
    // if we fitted forward, we start from the last_constrained hit
    // and then add the last inner hits with a Kalman step ...

    // create a temporary IMarlinTrack

    // TODO
    auto mTrk = std::make_shared<GaudiDDKalTestTrack>(m_thisAlg, const_cast<GaudiDDKalTest*>(&m_ddkaltest),
                                                      marlintrk.m_edm4hep_hits_to_kaltest_hits);
    // mTrk->m_edm4hep_hits_to_kaltest_hits = ;

    edm4hep::TrackState ts;

    double chi2Tmp = 0;
    int ndfTmp = 0;
    return_error = marlintrk.getTrackState(last_constrained_hit, ts, chi2, ndf);

    // m_thisAlg->debug()  << "  MarlinTrk::finaliseLCIOTrack:--  TrackState at last constrained hit : " << endmsg
    //   			<< toString( &ts )    << endmsg ;

    // need to add a dummy hit to the track
    mTrk->addHit(last_constrained_hit);
    mTrk->initialise(ts, fit_direction);

    auto hI = hits_in_fit.rbegin();

    while ((*hI).first != last_constrained_hit) {
      // m_thisAlg->debug()  << "  MarlinTrk::finaliseLCIOTrack:--  hit in reverse_iterator : "  << endmsg
      // 			  << toString( (*hI).first ) << endmsg ;
      ++hI;
    }

    ++hI;

    while (hI != hits_in_fit.rend()) {
      const edm4hep::TrackerHitPlane* h = (*hI).first;

      double deltaChi;
      double maxChi2Increment = 1e10; // ???

      int addHit = mTrk->addAndFit(h, deltaChi, maxChi2Increment);

      // m_thisAlg->debug() << " MarlinTrk::finaliseLCIOTrack: hit " << toString( h )
      // 			<< "  added : " << MarlinTrk::errorCode( addHit )
      // 			<< "  deltaChi2: " << deltaChi
      // 			<< endmsg ;

      // TODO
      if (addHit != 0) {
        m_thisAlg->error() << " ****  MarlinTrk::finaliseLCIOTrack:  could not add inner hit to track !!! " << endmsg;
      }

      ++hI;

    } //------------------------------------

    // m_thisAlg->debug() << "MarlinTrk::finaliseLCIOTrack: temporary kaltest track for track state at the IP: "
    //       	    <<  mTrk->toString() << endmsg ;

    // now propagate the temporary track to the IP
    return_error = mTrk->propagate(point, firstHit, trkStateIP, chi2Tmp, ndfTmp);

    // m_thisAlg->debug() << " ***  MarlinTrk::finaliseLCIOTrack: - propagated temporary track fromfirst hit to IP : "
    // <<  toString( trkStateIP ) << endmsg ;
  }

  // TODO
  if (return_error != 0) {
    m_thisAlg->debug() << "MarlinTrk::finaliseLCIOTrack: return_code for propagation = " << return_error
                       << " NDF = " << ndf << endmsg;

    return return_error;
  }

  trkStateIP.location = edm4hep::TrackState::AtIP;
  track.addToTrackStates(trkStateIP);
  track.setChi2(chi2);
  track.setNdf(ndf);

  // set the track states at the first and last hits
  // @ first hit

  m_thisAlg->debug() << "  >>>>>>>>>>>MarlinTrk::finaliseLCIOTrack:  create TrackState AtFirstHit" << endmsg;

  // TODO
  return_error = marlintrk.getTrackState(firstHit, trkStateAtFirstHit, chi2, ndf);

  // TODO
  if (return_error == 0) {
    trkStateAtFirstHit.location = edm4hep::TrackState::AtFirstHit;
    track.addToTrackStates(trkStateAtFirstHit);
  } else {
    m_thisAlg->warning()
        << "  >>>>>>>>>>>MarlinTrk::finaliseLCIOTrack:  MarlinTrk::finaliseLCIOTrack:  could not get TrackState "
           "at First Hit "
        << firstHit << endmsg;
  }

  // @ last hit

  m_thisAlg->debug() << "  >>>>>>>>>>> MarlinTrk::finaliseLCIOTrack: create TrackState AtLastHit : using trkhit "
                     << last_constrained_hit << endmsg;

  edm4hep::Vector3d last_hit_pos(lastHit->getPosition());

  edm4hep::TrackState trkStateAtLastHit;

  return_error = marlintrk.propagate(last_hit_pos, last_constrained_hit, trkStateAtLastHit, chi2, ndf);

  if (return_error == 0) {
    trkStateAtLastHit.location = edm4hep::TrackState::AtLastHit;
    track.addToTrackStates(trkStateAtLastHit);
  } else {
    m_thisAlg->debug() << "  >>>>>>>>>>> MarlinTrk::finaliseLCIOTrack:  could not get TrackState at Last Hit "
                       << last_constrained_hit << endmsg;
  }

  // set the track state at Calo Face
  edm4hep::TrackState trkStateCalo;
  // TODO
  bool tanL_is_positive = trkStateIP.tanLambda > 0;

  // TODO
  return_error = createTrackStateAtCaloFace(marlintrk, trkStateCalo, last_constrained_hit, tanL_is_positive);

  if (return_error == 0) {
    trkStateCalo.location = edm4hep::TrackState::AtCalorimeter;
    track.addToTrackStates(trkStateCalo);
  } else {
    m_thisAlg->debug() << "  >>>>>>>>>>> MarlinTrk::finaliseLCIOTrack:  could not get TrackState at Calo Face "
                       << endmsg;
    // FIXME: ignore track state at Calo face for debugging new tracking ...
    //  THIS IS ALSO PART OF THE PREVIOUS TODO
  }

  // This branch is never taken in the original MarlinTrkUtils code
  // as this function is always called without the last two arguments
  // meaning they default to nullptr

  // else {
  //   track.addToTrackStates(atLastHit);
  //   track.addToTrackStates(atCaloFace);
  // }

  return return_error;
}

void GaudiTrkUtils::addHitNumbersToTrack(std::vector<int32_t>& subdetectorHitNumbers,
                                         const std::vector<const edm4hep::TrackerHitPlane*>& hit_list, bool hits_in_fit,
                                         const dd4hep::DDSegmentation::BitFieldCoder& cellID_encoder) const {
  // Because in EDM4hep for vector members only hits can be added, we need the whole
  // vector before starting to assign
  std::map<int64_t, int32_t> hitNumbers;

  for (const auto& hit : hit_list) {
    // cellID_encoder.setValue(hit->getCellID());
    // int detID = cellID_encoder[UTIL::LCTrackerCellID::subdet()];
    int64_t detID = cellID_encoder.get(hit->getCellID(), 0);
    ++hitNumbers[detID];
  }

  int offset = 2;
  if (!hits_in_fit) { // all hit atributed by patrec
    offset = 1;
  }

  // this assumes that there is no tracker with an index larger than the ecal ...

  // subdetectorHitNumbers.resize(2 * lcio::ILDDetID::ECAL);
  subdetectorHitNumbers.resize(2 * 20);

  // for (std::map<int, int>::iterator it = hitNumbers.begin(); it != hitNumbers.end(); ++it) {
  for (const auto& [detIndex, hitNumber] : hitNumbers) {
    // track.subdetectorHitNumbers().at(2 * detIndex - offset) = it->second;
    subdetectorHitNumbers.at(2 * detIndex - offset) = hitNumber;
  }
}

int GaudiTrkUtils::createTrackStateAtCaloFace(GaudiDDKalTestTrack& marlintrk, edm4hep::TrackState& trkStateCalo,
                                              const edm4hep::TrackerHitPlane* trkhit, bool tanL_is_positive) {
  int return_error = 0;
  int return_error_barrel = 0;
  int return_error_endcap = 0;

  double chi2 = -DBL_MAX;
  int ndf = 0;

  std::string cellIDEncodingString = m_geoSvc->constantAsString(m_encodingStringVariable);
  dd4hep::DDSegmentation::BitFieldCoder encoder(cellIDEncodingString);
  std::uint64_t cellID = 0;

  // ================== need to get the correct ID(s) for the calorimeter face  ============================

  // TODO: ILDDet specific,
  // unsigned ecal_barrel_face_ID = lcio::ILDDetID::ECAL;
  // unsigned ecal_barrel_face_ID = 20;
  // unsigned ecal_endcap_face_ID = lcio::ILDDetID::ECAL_ENDCAP;
  unsigned ecal_endcap_face_ID = 29;

  //=========================================================================================================

  // subdet was in the original, corresponds to index 0
  // encoder[lcio::LCTrackerCellID::side()]   = lcio::ILDDetID::barrel;
  encoder.set(cellID, 0, 20);
  encoder.set(cellID, "side", 0);
  encoder.set(cellID, "layer", 0);

  int detElementID = 0;

  edm4hep::TrackState tsBarrel;
  edm4hep::TrackState tsEndcap;

  // propagate to the barrel layer
  return_error_barrel =
      marlintrk.propagateToLayer(encoder.lowWord(cellID), trkhit, tsBarrel, chi2, ndf, detElementID, 1);

  // propagate to the endcap layer
  // subdet was in the original, corresponds to index 0
  encoder.set(cellID, 0, ecal_endcap_face_ID);
  if (tanL_is_positive) {
    // encoder[lcio::LCTrackerCellID::side()] = lcio::ILDDetID::fwd;
    encoder.set(cellID, "side", 1);
  } else {
    // encoder[lcio::LCTrackerCellID::side()] = lcio::ILDDetID::bwd;
    encoder.set(cellID, "side", -1);
  }
  return_error_endcap =
      marlintrk.propagateToLayer(encoder.lowWord(cellID), trkhit, tsEndcap, chi2, ndf, detElementID, 1);

  // // check which is the right intersection / closer to the trkhit
  if (return_error_barrel == no_intersection) {
    // if barrel fails just return ts at the Endcap if exists
    return_error = return_error_endcap;
    trkStateCalo = tsEndcap;
  } else if (return_error_endcap == no_intersection) {
    // if barrel succeeded and endcap fails return ts at the barrel
    return_error = return_error_barrel;
    trkStateCalo = tsBarrel;
  } else {
    // this means both barrel and endcap have intersections. Return closest to the tracker hit
    edm4hep::Vector3d hitPos(trkhit->getPosition());
    edm4hep::Vector3d barrelPos = {tsBarrel.referencePoint[0], tsBarrel.referencePoint[1], tsBarrel.referencePoint[2]};
    edm4hep::Vector3d endcapPos = {tsEndcap.referencePoint[0], tsEndcap.referencePoint[1], tsEndcap.referencePoint[2]};
    double dToBarrel = edm4hep::utils::magnitude(hitPos - barrelPos);
    double dToendcap = edm4hep::utils::magnitude(hitPos - endcapPos);

    if (dToBarrel < dToendcap) {
      return_error = return_error_barrel;
      trkStateCalo = tsBarrel;
    } else {
      return_error = return_error_endcap;
      trkStateCalo = tsEndcap;
    }
  }

  // bd: d0 and z0 of the track state at the calorimeter must be 0 by definition for all tracks.
  //  info() << "  >>>>>>>>>>> createTrackStateAtCaloFace : setting d0 and z0 to 0. for track state at calorimeter : "
  //         << toString(trkStateCalo) << endmsg;

  trkStateCalo.D0 = 0.;
  trkStateCalo.Z0 = 0.;

  if (return_error != 0) {
    m_thisAlg->info()
        << "  >>>>>>>>>>> createTrackStateAtCaloFace :  could not get TrackState at Calo Face: return_error = "
        << return_error << endmsg;
  }

  return return_error;
}
