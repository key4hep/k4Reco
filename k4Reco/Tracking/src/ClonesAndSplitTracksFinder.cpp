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

  auto trackVec = edm4hep::TrackCollection();

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
    if (std::ranges::find(trackVec2, trackVec1.at(hit)) != trackVec2.end())
      nHitsInCommon++;
  }
  return nHitsInCommon;
}

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
        auto it_trk = std::ranges::find(trackVecFinal, track_final);
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
          auto it_trk = std::ranges::find(trackVecFinal, bestTracksMultiConnections.at(0));
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

        auto trk1 = std::ranges::find(trackVecFinal, inputTracks.at(track_a_id));

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
