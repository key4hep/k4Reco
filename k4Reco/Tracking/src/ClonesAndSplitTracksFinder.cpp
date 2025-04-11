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
#include "ClonesAndSplitTracksFinder.h"
// #include "HitsSorterAndDebugger.h"

#include <streamlog/logstream.h>
#include <streamlog/streamlog.h>

#include <edm4hep/MutableTrack.h>
#include <edm4hep/Track.h>
#include <edm4hep/TrackCollection.h>

#include <GaudiKernel/ISvcLocator.h>

#include <algorithm>
#include <map>
#include <string>
#include <utility>
#include <vector>

ClonesAndSplitTracksFinder::ClonesAndSplitTracksFinder(const std::string& name, ISvcLocator* svcLoc)
    : Transformer(name, svcLoc,
                  {
                      KeyValues("InputTrackCollectionName", {"SiTracks"}),
                  },
                  {
                      KeyValues("OutputTrackCollectionName", {"SiTracksMerged"}),
                  }) {}

StatusCode ClonesAndSplitTracksFinder::initialize() {
  // Setting the streamlog output is necessary to avoid lots of overhead.
  // Otherwise it would be equivalent to running with every debug message
  // being computed
  streamlog::out.init(std::cout, "");
  streamlog::logscope* scope = new streamlog::logscope(streamlog::out);
  scope->setLevel<streamlog::MESSAGE0>();

  // // usually a good idea to
  // printParameters();

  // _trksystem = MarlinTrk::Factory::createMarlinTrkSystem("DDKalTest", nullptr, "");

  // _magneticField = MarlinUtil::getBzAtOrigin();
  // ///////////////////////////////

  // _encoder = std::make_shared<UTIL::BitField64>(lcio::LCTrackerCellID::encoding_string());

  // if (not _trksystem) {
  //   throw EVENT::Exception("Cannot initialize MarlinTrkSystem of Type: DDKalTest");
  // }

  // _trksystem->setOption(MarlinTrk::IMarlinTrkSystem::CFG::useQMS, _MSOn);
  // _trksystem->setOption(MarlinTrk::IMarlinTrkSystem::CFG::usedEdx, _ElossOn);
  // _trksystem->setOption(MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing, _SmoothOn);
  // _trksystem->init();

  // // Put default values for track fitting
  // m_initialTrackError_d0 = 1.e6;
  // m_initialTrackError_phi0 = 1.e2;
  // m_initialTrackError_omega = 1.e-4;
  // m_initialTrackError_z0 = 1.e6;
  // m_initialTrackError_tanL = 1.e2;
  // m_maxChi2perHit = 1.e2;

  // _n_run = 0;
  return StatusCode::SUCCESS;
}

edm4hep::TrackCollection ClonesAndSplitTracksFinder::operator()(const edm4hep::TrackCollection& input_track_col) const {
  // set the correct configuration for the tracking system for this event
  // TODO:
  // MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::useQMS> mson(_trksystem, _MSOn);
  // MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::usedEdx> elosson(_trksystem, _ElossOn);
  // MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing> smoothon(_trksystem, _SmoothOn);

  const size_t nTracks = input_track_col.size();
  debug() << " >> ClonesAndSplitTracksFinder starts with " << nTracks << " tracks." << endmsg;

  // establish the track collection that will be created
  auto trackVec = edm4hep::TrackCollection();
  // _encoder->reset();

  //------------
  // FIRST STEP: REMOVE CLONES
  //------------

  edm4hep::TrackCollection tracksWithoutClones = removeClones(input_track_col);
  const size_t ntracksWithoutClones = tracksWithoutClones.size();
  debug() << " >> ClonesAndSplitTracksFinder found " << ntracksWithoutClones << " tracks without clones." << endmsg;

  if (m_mergeSplitTracks && ntracksWithoutClones > 1) {
    throw std::runtime_error("Setting mergeSplitTracks to true is not yet implemented.");
    // debug() << " Try to merge tracks ..." << endmsg;

    //------------
    // SECOND STEP: MERGE TRACKS
    //------------

    // TODO:
    // mergeSplitTracks(trackVec, input_track_col, tracksWithoutClones);
  } else {
    debug() << " Not even try to merge tracks ..." << endmsg;
    trackVec = std::move(tracksWithoutClones);
  }

  return trackVec;
}

// Function to check if two KDtracks contain several hits in common
size_t ClonesAndSplitTracksFinder::overlappingHits(const edm4hep::Track& track1, const edm4hep::Track& track2) const {
  size_t nHitsInCommon = 0;
  const auto& trackVec1 = track1.getTrackerHits();
  const auto& trackVec2 = track2.getTrackerHits();
  for (size_t hit = 0; hit < trackVec1.size(); hit++) {
    if (std::find(trackVec2.begin(), trackVec2.end(), trackVec1.at(hit)) != trackVec2.end())
      nHitsInCommon++;
  }
  return nHitsInCommon;
}

// void ClonesAndSplitTracksFinder::fromTrackToTrackImpl(const edm4hep::Track* track, TrackImpl*& trackFinal) {
//   const TrackState* ts_atOther = 0;
//   ts_atOther = track->getTrackState(TrackState::AtOther);
//   if (ts_atOther)
//     trackFinal->addTrackState(new TrackStateImpl(*ts_atOther));
//   const TrackState* ts_atIP = 0;
//   ts_atIP = track->getTrackState(TrackState::AtIP);
//   if (ts_atIP)
//     trackFinal->addTrackState(new TrackStateImpl(*ts_atIP));
//   const TrackState* ts_atFirstHit = 0;
//   ts_atFirstHit = track->getTrackState(TrackState::AtFirstHit);
//   if (ts_atFirstHit)
//     trackFinal->addTrackState(new TrackStateImpl(*ts_atFirstHit));
//   const TrackState* ts_atLastHit = 0;
//   ts_atLastHit = track->getTrackState(TrackState::AtLastHit);
//   if (ts_atLastHit)
//     trackFinal->addTrackState(new TrackStateImpl(*ts_atLastHit));
//   const TrackState* ts_atCalorimeter = 0;
//   ts_atCalorimeter = track->getTrackState(TrackState::AtCalorimeter);
//   if (ts_atCalorimeter)
//     trackFinal->addTrackState(new TrackStateImpl(*ts_atCalorimeter));
//   const TrackState* ts_atVertex = 0;
//   ts_atVertex = track->getTrackState(TrackState::AtVertex);
//   if (ts_atVertex)
//     trackFinal->addTrackState(new TrackStateImpl(*ts_atVertex));
//   const TrackState* ts_atLastLocation = 0;
//   ts_atLastLocation = track->getTrackState(TrackState::LastLocation);
//   if (ts_atLastLocation)
//     trackFinal->addTrackState(new TrackStateImpl(*ts_atLastLocation));

//   for (UInt_t i = 0; i < (track->getTrackerHits()).size(); i++) {
//     trackFinal->addHit(track->getTrackerHits().at(i));
//   }
//   trackFinal->setRadiusOfInnermostHit(track->getRadiusOfInnermostHit());
//   trackFinal->setChi2(track->getChi2());
//   trackFinal->setNdf(track->getNdf());
//   trackFinal->setdEdx(track->getdEdx());
//   trackFinal->setdEdxError(track->getdEdxError());
// }

edm4hep::TrackCollection
ClonesAndSplitTracksFinder::removeClones(const edm4hep::TrackCollection& input_track_col) const {
  edm4hep::TrackCollection tracksWithoutClones;
  debug() << "ClonesAndSplitTracksFinder::removeClones " << endmsg;

  // loop over the input tracks

  std::multimap<size_t, std::pair<size_t, edm4hep::Track>> candidateClones;

  for (size_t iTrack = 0; iTrack < input_track_col.size(); ++iTrack) {
    int countClones = 0;
    const edm4hep::Track& track1 = input_track_col.at(iTrack);

    for (size_t jTrack = iTrack + 1; jTrack < input_track_col.size(); ++jTrack) {
      const edm4hep::Track& track2 = input_track_col.at(jTrack);

      if (track1 != track2) {
        const size_t nOverlappingHits = overlappingHits(track1, track2);
        if (nOverlappingHits >= 2) { // clones
          countClones++;
          edm4hep::Track bestTrack = bestInClones(track1, track2, nOverlappingHits);
          candidateClones.emplace(iTrack, std::make_pair(jTrack, bestTrack));
        } else {
          continue;
        }
      }

    } // end second track loop

    if (countClones == 0) {
      tracksWithoutClones.push_back(track1.clone());
    }

  } // end first track loop

  filterClonesAndMergedTracks(candidateClones, input_track_col, tracksWithoutClones, true);
  return tracksWithoutClones;
}

// void ClonesAndSplitTracksFinder::mergeSplitTracks(std::unique_ptr<LCCollectionVec>& trackVec,
//                                                   LCCollection*& input_track_col,
//                                                   EVENT::TrackVec& tracksWithoutClones) {
//   debug() << "ClonesAndSplitTracksFinder::mergeSplitTracks " << endmsg;

//   std::multimap<int, std::pair<int, edm4hep::Track*>> mergingCandidates;
//   std::set<int> iter_duplicates;

//   for (UInt_t iTrack = 0; iTrack < tracksWithoutClones.size(); ++iTrack) {
//     int countMergingPartners = 0;
//     bool toBeMerged = false;
//     edm4hep::Track* track_i = static_cast<edm4hep::Track*>(tracksWithoutClones.at(iTrack));

//     double pt_i = 0.3 * _magneticField / (fabs(track_i->getOmega()) * 1000.);
//     double theta_i = (M_PI / 2 - atan(track_i->getTanLambda())) * 180. / M_PI;
//     double phi_i = track_i->getPhi() * 180. / M_PI;

//     // Merge only tracks with min pt
//     // Try to avoid merging loopers for now
//     if (pt_i < _minPt) {
//       debug() << " Track #" << iTrack << ": pt = " << pt_i << ", theta = " << theta_i
//                             << ", phi = " << phi_i << endmsg;
//       debug() << " Track #" << iTrack << " does not fulfil min pt requirement." << endmsg;
//       debug() << " TRACK STORED" << endmsg;

//       TrackImpl* trackFinal = new TrackImpl;
//       fromTrackToTrackImpl(track_i, trackFinal);
//       trackVec->addElement(trackFinal);
//       continue;
//     }

//     for (UInt_t jTrack = iTrack + 1; jTrack < tracksWithoutClones.size(); ++jTrack) {
//       edm4hep::Track* track_j = static_cast<edm4hep::Track*>(tracksWithoutClones.at(jTrack));
//       bool isCloseInTheta = false, isCloseInPhi = false, isCloseInPt = false;

//       if (track_j != track_i) {
//         double pt_j = 0.3 * _magneticField / (fabs(track_j->getOmega() * 1000.));
//         double theta_j = (M_PI / 2 - atan(track_j->getTanLambda())) * 180. / M_PI;
//         double phi_j = track_j->getPhi() * 180. / M_PI;
//         debug() << " Track #" << iTrack << ": pt = " << pt_i << ", theta = " << theta_i
//                               << ", phi = " << phi_i << endmsg;
//         debug() << " Track #" << jTrack << ": pt = " << pt_j << ", theta = " << theta_j
//                               << ", phi = " << phi_j << endmsg;

//         if (pt_j < _minPt) {
//           debug() << " Track #" << jTrack << " does not fulfil min pt requirement. Skip. " << endmsg;
//           continue;
//         }

//         double significanceTanLambda = calculateSignificanceTanLambda(track_i, track_j);
//         double significancePhi = calculateSignificancePhi(track_i, track_j);
//         double significancePt = calculateSignificancePt(track_i, track_j);

//         debug() << " -> tanLambda significance = " << significanceTanLambda << " with cut at "
//                               << _maxSignificanceTheta << endmsg;
//         if (significanceTanLambda < _maxSignificanceTheta) {
//           isCloseInTheta = true;
//           debug() << " Tracks are close in theta " << endmsg;
//         }

//         debug() << " -> phi significance = " << significancePhi << " with cut at " << _maxSignificancePhi
//                               << endmsg;
//         if (significancePhi < _maxSignificancePhi) {
//           isCloseInPhi = true;
//           debug() << " Tracks are close in phi " << endmsg;
//         }

//         debug() << " -> pt significance = " << significancePt << " with cut at " << _maxSignificancePt
//                               << endmsg;
//         if (significancePt < _maxSignificancePt) {
//           isCloseInPt = true;
//           debug() << " Tracks are close in pt  " << endmsg;
//         }

//         if (streamlog::out.write<streamlog::DEBUG5>()) {
//           debug() << " Track #" << iTrack << ": " << endmsg;
//           printHits(track_i);
//           debug() << " Track #" << jTrack << ": " << endmsg;
//           printHits(track_j);
//         }

//         toBeMerged = isCloseInTheta && isCloseInPhi && isCloseInPt;

//         if (toBeMerged) { // merging, refitting, storing in a container of mergingCandidates (multimap <*track1,
//                           // pair<*track2,*trackMerged>>)
//           EVENT::edm4hep::Track* lcioTrkPtr = nullptr;
//           mergeAndFit(track_i, track_j, lcioTrkPtr);
//           if (not lcioTrkPtr) {
//             continue;
//           }
//           mergingCandidates.insert(make_pair(iTrack, make_pair(jTrack, lcioTrkPtr)));
//           countMergingPartners++;
//           iter_duplicates.insert(iTrack);
//           iter_duplicates.insert(jTrack);
//         } else { // no merging conditions met
//           continue;
//         }
//       }

//     } // end loop on jTracks

//     // Track was already found as duplicate
//     const bool is_in = iter_duplicates.find(iTrack) != iter_duplicates.end();
//     if (countMergingPartners == 0 && !is_in) { // if track_i has no merging partner, store it in the output vec
//       debug() << " Track #" << iTrack << " has no merging partners, so TRACK STORED." << endmsg;

//       TrackImpl* trackFinal = new TrackImpl;
//       fromTrackToTrackImpl(track_i, trackFinal);
//       trackVec->addElement(trackFinal);

//     } else {
//       debug() << " TRACK NOT STORED" << endmsg;
//     }
//     if (countMergingPartners != 0)
//       debug() << " possible merging partners for track #" << iTrack << " are = " << countMergingPartners
//                             << endmsg;

//   } // end loop on iTracks

//   EVENT::TrackVec finalTracks;
//   filterClonesAndMergedTracks(mergingCandidates, input_track_col, finalTracks, false);

//   for (UInt_t iTrk = 0; iTrk < finalTracks.size(); iTrk++) {
//     debug() << " TRACK STORED" << endmsg;
//     TrackImpl* trackFinal = new TrackImpl;
//     fromTrackToTrackImpl(finalTracks.at(iTrk), trackFinal);
//     trackVec->addElement(trackFinal);
//   }
// }

// double ClonesAndSplitTracksFinder::calculateSignificancePt(const edm4hep::Track* first, const edm4hep::Track* second)
// {
//   float omegaFirst = first->getOmega();
//   float omegaSecond = second->getOmega();

//   double ptFirst = 0.3 * _magneticField / (fabs(first->getOmega() * 1000.));
//   double ptSecond = 0.3 * _magneticField / (fabs(second->getOmega() * 1000.));

//   const float sigmaPOverPFirst = sqrt(first->getCovMatrix()[5]) / fabs(omegaFirst);
//   const float sigmaPOverPSecond = sqrt(second->getCovMatrix()[5]) / fabs(omegaSecond);
//   const float sigmaPtFirst = ptFirst * sigmaPOverPFirst;
//   const float sigmaPtSecond = ptSecond * sigmaPOverPSecond;

//   const double significance = calculateSignificance(ptFirst, ptSecond, sigmaPtFirst, sigmaPtSecond);

//   return significance;
// }

// double ClonesAndSplitTracksFinder::calculateSignificancePhi(const edm4hep::Track* first, const edm4hep::Track*
// second) {
//   float phiFirst = first->getPhi();
//   float phiSecond = second->getPhi();
//   float deltaPhi = (M_PI - std::abs(std::abs(phiFirst - phiSecond) - M_PI));

//   const float sigmaPhiFirst = sqrt(first->getCovMatrix()[2]);
//   const float sigmaPhiSecond = sqrt(second->getCovMatrix()[2]);

//   const double significance = calculateSignificance(deltaPhi, 0.0, sigmaPhiFirst, sigmaPhiSecond);

//   return significance;
// }

// double ClonesAndSplitTracksFinder::calculateSignificanceTanLambda(const edm4hep::Track* first, const edm4hep::Track*
// second) {
//   float tanLambdaFirst = first->getTanLambda();
//   float tanLambdaSecond = second->getTanLambda();

//   const float sigmaTanLambdaFirst = sqrt(first->getCovMatrix()[14]);
//   const float sigmaTanLambdaSecond = sqrt(second->getCovMatrix()[14]);

//   const double significance =
//       calculateSignificance(tanLambdaFirst, tanLambdaSecond, sigmaTanLambdaFirst, sigmaTanLambdaSecond);

//   return significance;
// }

// double ClonesAndSplitTracksFinder::calculateSignificance(const double firstPar, const double secondPar,
//                                                          const double firstPar_sigma, const double secondPar_sigma) {
//   const float delta = fabs(firstPar - secondPar);
//   const float sigmaDelta = sqrt(firstPar_sigma * firstPar_sigma + secondPar_sigma * secondPar_sigma);

//   return delta / sigmaDelta;
// }

void ClonesAndSplitTracksFinder::filterClonesAndMergedTracks(
    std::multimap<size_t, std::pair<size_t, edm4hep::Track>>& candidates, const edm4hep::TrackCollection& inputTracks,
    edm4hep::TrackCollection& trackVecFinal, bool clones) const {
  // TODO: Fix this vector
  std::vector<podio::RelationRange<edm4hep::TrackerHit>> savedHitVec;

  for (const auto& iter : candidates) {
    size_t track_a_id = iter.first;
    size_t track_b_id = iter.second.first;
    edm4hep::Track track_final = iter.second.second;
    size_t countConnections = candidates.count(track_a_id);
    bool multiConnection = (countConnections > 1);

    if (!multiConnection) { // if only 1 connection

      if (clones) { // clones: compare the track pointers
        auto it_trk = std::find(trackVecFinal.begin(), trackVecFinal.end(), track_final);
        if (it_trk != trackVecFinal.end()) { // if the track is already there, do nothing
          continue;
        }
        trackVecFinal.push_back(track_final.clone());
      } else { // mergeable tracks: compare the sets of tracker hits

        const auto& track_final_hits = track_final.getTrackerHits();
        bool toBeSaved = true;

        for (const auto& hitsVec : savedHitVec) {
          if (equal(hitsVec.begin(), hitsVec.end(), track_final_hits.begin())) {
            toBeSaved = false;
            break;
          }
        }

        if (toBeSaved) {
          savedHitVec.push_back(track_final_hits);
          trackVecFinal.push_back(track_final.clone());
        }
      }

    } else { // if more than 1 connection, clones and mergeable tracks have to be treated a little different

      if (clones) { // clones

        // look at the elements with equal range. If their bestTrack is the same, store it (if not already in). If their
        // bestTrack is different, don't store it
        auto ret = candidates.equal_range(
            track_a_id); // a std::pair of iterators on the multimap [
                         // std::pair<std::multimap<edm4hep::Track*,std::pair<edm4hep::Track*,edm4hep::Track*>>::iterator,
                         // std::multimap<edm4hep::Track*,std::pair<edm4hep::Track*,edm4hep::Track*>>::iterator> ]
        edm4hep::TrackCollection bestTracksMultiConnections;
        for (auto it = ret.first; it != ret.second; ++it) {
          edm4hep::Track track_best = it->second.second;
          bestTracksMultiConnections.push_back(track_best);
        }
        if (std::adjacent_find(bestTracksMultiConnections.begin(), bestTracksMultiConnections.end(),
                               std::not_equal_to<edm4hep::Track>()) ==
            bestTracksMultiConnections.end()) { // one best track with the same track key
          auto it_trk = std::find(trackVecFinal.begin(), trackVecFinal.end(), bestTracksMultiConnections.at(0));
          if (it_trk != trackVecFinal.end()) { // if the track is already there, do nothing
            continue;
          }
          trackVecFinal.push_back(bestTracksMultiConnections.at(0));

        } else { // multiple best tracks with the same track key
          continue;
        }

      } // end of clones

      else { // mergeable tracks -- at the moment they are all stored (very rare anyways)
        const edm4hep::Track& track_a = inputTracks.at(track_a_id);
        const edm4hep::Track& track_b = inputTracks.at(track_b_id);

        // auto trk1 = std::find(trackVecFinal.begin(), trackVecFinal.end(), track_a);
        auto trk1 = std::find(trackVecFinal.begin(), trackVecFinal.end(), inputTracks.at(track_a_id));

        if (trk1 != trackVecFinal.end()) { // if the track1 is already there
          continue;
        }
        // otherwise store the two tracks
        trackVecFinal.push_back(track_a);
        trackVecFinal.push_back(track_b);

      } // end of mergeable tracks
    }
  }
}

// void ClonesAndSplitTracksFinder::mergeAndFit(edm4hep::Track* track_i, edm4hep::Track* track_j, edm4hep::Track*&
// lcioTrkPtr) {
//   debug() << "ClonesAndSplitTracksFinder::mergeAndFit " << endmsg;
//   EVENT::TrackerHitVec trkHits_i = track_i->getTrackerHits();
//   EVENT::TrackerHitVec trkHits_j = track_j->getTrackerHits();

//   EVENT::TrackerHitVec trkHits;
//   for (UInt_t iHits = 0; iHits < trkHits_i.size(); iHits++) {
//     trkHits.push_back(trkHits_i.at(iHits));
//   }
//   // Remove common hits while filling for the second track
//   for (UInt_t jHits = 0; jHits < trkHits_j.size(); jHits++) {
//     if (std::find(trkHits.begin(), trkHits.end(), trkHits_j.at(jHits)) != trkHits.end()) {
//       debug() << " This hit is already in the track" << endmsg;
//       continue;
//     } else {
//       trkHits.push_back(trkHits_j.at(jHits));
//     }
//   }
//   std::sort(trkHits.begin(), trkHits.end(), sort_by_radius);
//   if (streamlog::out.write<streamlog::DEBUG5>()) {
//     debug() << " Hits in track to be merged: " << endmsg;
//     printHits(trkHits);
//   }

//   auto mergedTrack = std::unique_ptr<TrackImpl>(new TrackImpl);

//   auto marlin_trk = std::unique_ptr<MarlinTrk::IMarlinTrack>(_trksystem->createTrack());

//   // Make an initial covariance matrix with very broad default values
//   EVENT::FloatVec covMatrix(15, 0);          // Size 15, filled with 0s
//   covMatrix[0] = (_initialTrackError_d0);    // sigma_d0^2
//   covMatrix[2] = (_initialTrackError_phi0);  // sigma_phi0^2
//   covMatrix[5] = (_initialTrackError_omega); // sigma_omega^2
//   covMatrix[9] = (_initialTrackError_z0);    // sigma_z0^2
//   covMatrix[14] = (_initialTrackError_tanL); // sigma_tanl^2

//   const bool direction = _extrapolateForward ? MarlinTrk::IMarlinTrack::forward : MarlinTrk::IMarlinTrack::backward;

//   int fit_status = MarlinTrk::createFinalisedLCIOTrack(marlin_trk.get(), trkHits, mergedTrack.get(), direction,
//                                                        covMatrix, _magneticField, _maxChi2perHit);

//   if (fit_status != 0) {
//     streamlog_out(DEBUG4) << "Fit failed with error status " << fit_status << endmsg;
//     return;
//   }
//   debug() << " >> Fit not failed ! " << endmsg;

//   // fit finished - get hits in the fit
//   std::vector<std::pair<EVENT::TrackerHit*, double>> hits_in_fit;
//   std::vector<std::pair<EVENT::TrackerHit*, double>> outliers;

//   // remember the hits are ordered in the order in which they were fitted

//   marlin_trk->getHitsInFit(hits_in_fit);
//   if (hits_in_fit.size() < 3) {
//     streamlog_out(DEBUG4) << "Less than 3 hits in fit: Track discarded. Number of hits =  " << trkHits.size()
//                           << endmsg;
//     return;
//   }

//   marlin_trk->getOutliers(outliers);

//   std::vector<TrackerHit*> all_hits;
//   all_hits.reserve(hits_in_fit.size() + outliers.size());

//   for (unsigned ihit = 0; ihit < hits_in_fit.size(); ++ihit) {
//     all_hits.push_back(hits_in_fit[ihit].first);
//   }

//   for (unsigned ihit = 0; ihit < outliers.size(); ++ihit) {
//     all_hits.push_back(outliers[ihit].first);
//   }

//   UTIL::BitField64 encoder2(lcio::LCTrackerCellID::encoding_string());
//   encoder2.reset(); // reset to 0
//   MarlinTrk::addHitNumbersToTrack(mergedTrack.get(), all_hits, false, encoder2);
//   MarlinTrk::addHitNumbersToTrack(mergedTrack.get(), hits_in_fit, true, encoder2);

//   if (streamlog::out.write<streamlog::DEBUG5>()) {
//     debug() << " Merged track : " << endmsg;
//     printHits(&*mergedTrack);
//   }

//   lcioTrkPtr = mergedTrack.release();
// }

edm4hep::Track ClonesAndSplitTracksFinder::bestInClones(const edm4hep::Track& track_a, const edm4hep::Track& track_b,
                                                        size_t nOverlappingHits) const {
  // This function compares two tracks which have a certain number of overlapping hits and returns the best track
  // The best track is chosen based on length (in terms of number of hits) and chi2/ndf requirements
  // In general, the longest track is preferred. When clones have same length, the one with best chi2/ndf is chosen

  edm4hep::Track bestTrack;

  const auto& trackerHit_a = track_a.getTrackerHits();
  const auto& trackerHit_b = track_b.getTrackerHits();

  size_t trackerHit_a_size = trackerHit_a.size();
  size_t trackerHit_b_size = trackerHit_b.size();

  float b_chi2 = track_b.getChi2() / static_cast<float>(track_b.getNdf());
  float a_chi2 = track_a.getChi2() / static_cast<float>(track_a.getNdf());

  if (nOverlappingHits == trackerHit_a_size) { // if the second track is the first track + segment
    bestTrack = track_b;
  } else if (nOverlappingHits == trackerHit_b_size) { // if the second track is a subtrack of the first track
    bestTrack = track_a;
  } else if (trackerHit_b_size == trackerHit_a_size) { // if the two tracks have the same length
    if (b_chi2 <= a_chi2) {
      bestTrack = track_b;
    } else {
      bestTrack = track_a;
    }
  } else if (trackerHit_b_size > trackerHit_a_size) { // if the second track is longer
    bestTrack = track_b;
  } else if (trackerHit_b_size < trackerHit_a_size) { // if the second track is shorter
    bestTrack = track_a;
  }
  return bestTrack;
}
