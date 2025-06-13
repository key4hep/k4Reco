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
#include "RefitFinal.h"

#include "GaudiDDKalTestTrack.h"
#include "GaudiTrkUtils.h"

#include <streamlog/logstream.h>
#include <streamlog/streamlog.h>

#include <DD4hep/BitFieldCoder.h>

#include <edm4hep/TrackCollection.h>
#include <edm4hep/TrackMCParticleLinkCollection.h>

#include "k4FWCore/Transformer.h"

#include <string>

RefitFinal::RefitFinal(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(name, svcLoc,
                       {
                           KeyValues("InputTrackCollectionName", {"TruthTracks"}),
                           KeyValues("InputRelationCollectionName", {"SiTrackRelations"}),
                       },
                       {
                           KeyValues("OutputTrackCollectionName", {"RefittedTracks"}),
                           KeyValues("OutputRelationCollectionName", {"RefittedRelation"}),
                       }) {}

StatusCode RefitFinal::initialize() {
  // Setting the streamlog output is necessary to avoid lots of overhead.
  // Otherwise it would be equivalent to running with every debug message
  // being computed
  streamlog::out.init(std::cout, "");
  streamlog::logscope* scope = new streamlog::logscope(streamlog::out);
  setStreamlogOutputLevel(this, scope);

  // usually a good idea to
  // printParameters();

  // _trksystem = MarlinTrk::Factory::createMarlinTrkSystem("DDKalTest", nullptr, "");

  m_geoSvc = serviceLocator()->service(m_geoSvcName);
  if (!m_geoSvc) {
    error() << "Unable to retrieve GeoSvc" << endmsg;
    return StatusCode::FAILURE;
  }
  std::string cellIDEncodingString = m_geoSvc->constantAsString(m_encodingStringVariable.value());
  m_encoder = dd4hep::DDSegmentation::BitFieldCoder(cellIDEncodingString);

  m_ddkaltest.init();
  m_ddkaltest.setEncoder(m_encoder);

  // ///////////////////////////////

  // _encoder = std::make_shared<UTIL::BitField64>(lcio::LCTrackerCellID::encoding_string());

  // if (not _trksystem) {
  //   throw EVENT::Exception("Cannot initialize MarlinTrkSystem of Type: DDKalTest");
  // }

  // _trksystem->setOption(MarlinTrk::IMarlinTrkSystem::CFG::useQMS, _MSOn);
  // _trksystem->setOption(MarlinTrk::IMarlinTrkSystem::CFG::usedEdx, _ElossOn);
  // _trksystem->setOption(MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing, _SmoothOn);
  // _trksystem->init();

  // _n_run = 0;

  return StatusCode::SUCCESS;
}

std::tuple<edm4hep::TrackCollection, edm4hep::TrackMCParticleLinkCollection>
RefitFinal::operator()(const edm4hep::TrackCollection& input_track_col,
                       const std::vector<const edm4hep::TrackMCParticleLinkCollection*>& input_rel_col) const {
  // // set the correct configuration for the tracking system for this event
  // MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::useQMS> mson(_trksystem, _MSOn);
  // MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::usedEdx> elosson(_trksystem, _ElossOn);
  // MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing> smoothon(_trksystem, _SmoothOn);

  // establish the track collection that will be created
  edm4hep::TrackCollection trackVec;
  edm4hep::TrackMCParticleLinkCollection trackRelationCollection;

  if (input_rel_col.size() == 0) {
    debug() << "No input relation collection, not creating one either" << endmsg;
  }

  const size_t nTracks = input_track_col.size();

  debug() << " Number of Tracks " << nTracks << endmsg;

  // loop over the input tracks and refit
  for (size_t iTrack = 0; iTrack < nTracks; ++iTrack) {
    const auto& track = input_track_col.at(iTrack);
    const auto& trkHits = track.getTrackerHits();
    // TODO: Don't create a copy because finaliseLCIOTrack and marlin_trk need pointers
    std::vector<const edm4hep::TrackerHit*> trkHitsPtr;
    trkHitsPtr.reserve(trkHits.size());
    for (const auto& hit : trkHits) {
      trkHitsPtr.push_back(&hit);
    }

    auto marlin_trk = GaudiDDKalTestTrack(this, const_cast<GaudiDDKalTest*>(&m_ddkaltest));

    debug() << "---- track n = " << iTrack << "  n hits = " << trkHits.size() << endmsg;

    for (const auto& ptr : trkHitsPtr) {
      marlin_trk.addHit(ptr);
    }

    int init_status = FitInit2(track, marlin_trk);

    if (init_status != 0) {
      continue;
    }

    // debug() << "Refit: Trackstate after initialisation\n" << marlin_trk.toString() << endmsg;

    debug() << "track initialised " << endmsg;

    int fit_status = marlin_trk.fit();

    // debug() << "RefitHit: Trackstate after fit()\n" << marlin_trk.toString() << endmsg;

    if (fit_status != 0) {
      continue;
    }

    edm4hep::MutableTrack lcio_trk = trackVec.create();

    GaudiTrkUtils trkUtils(static_cast<const Gaudi::Algorithm*>(this), m_ddkaltest, m_geoSvc,
                           m_encodingStringVariable.value());

    int return_code = trkUtils.finaliseLCIOTrack(marlin_trk, lcio_trk, trkHitsPtr, true);

    if (return_code != 0) {
      debug() << "finaliseLCIOTrack failed" << endmsg;
      continue;
    }

    // debug() << " *** created finalized LCIO track - return code " << return_code << std::endl << *lcio_trk << endmsg;

    // fit finished - get hits in the fit

    // remember the hits are ordered in the order in which they were fitted

    const auto hits_in_fit = marlin_trk.getHitsInFit();

    if (int(hits_in_fit.size()) < m_minClustersOnTrackAfterFit) {
      debug() << "Less than " << m_minClustersOnTrackAfterFit
              << " hits in fit: Track "
                 "Discarded. Number of hits =  "
              << trkHits.size() << endmsg;
      continue;
    }

    const auto outliers = marlin_trk.getOutliers();

    std::vector<const edm4hep::TrackerHit*> all_hits;
    std::vector<const edm4hep::TrackerHit*> hits_in_fit_ptr;
    all_hits.reserve(hits_in_fit.size() + outliers.size());
    hits_in_fit_ptr.reserve(hits_in_fit.size());

    for (unsigned ihit = 0; ihit < hits_in_fit.size(); ++ihit) {
      all_hits.push_back(hits_in_fit[ihit].first);
      hits_in_fit_ptr.push_back(hits_in_fit[ihit].first);
    }

    for (unsigned ihit = 0; ihit < outliers.size(); ++ihit) {
      all_hits.push_back(outliers[ihit].first);
    }

    std::string cellIDEncodingString = m_geoSvc->constantAsString(m_encodingStringVariable.value());
    dd4hep::DDSegmentation::BitFieldCoder encoder2(cellIDEncodingString);
    std::vector<int32_t> subdetectorHitNumbers;
    trkUtils.addHitNumbersToTrack(subdetectorHitNumbers, all_hits, false, encoder2);
    trkUtils.addHitNumbersToTrack(subdetectorHitNumbers, hits_in_fit_ptr, true, encoder2);
    for (const auto num : subdetectorHitNumbers) {
      lcio_trk.addToSubdetectorHitNumbers(num);
    }

    // debug() << "processEvent: Hit numbers for track " << lcio_trk.id() << ":  " << endmsg;
    int detID = 0;
    for (size_t ip = 0; ip < lcio_trk.getSubdetectorHitNumbers().size(); ip = ip + 2) {
      detID++;
      debug() << "  det id " << detID << " , nhits in track = " << lcio_trk.getSubdetectorHitNumbers()[ip]
              << " , nhits in fit = " << lcio_trk.getSubdetectorHitNumbers()[ip + 1] << endmsg;
      if (lcio_trk.getSubdetectorHitNumbers()[ip] > 0)
        // TODO: is detID - 1 correct?
        lcio_trk.setType(lcio_trk.getType() | (1 << detID));
    }

    // TODO:
    // if (input_rel_col) {
    //   auto mcParticleVec = relation->getRelatedToObjects(track);
    //   auto weightVec     = relation->getRelatedToWeights(track);
    //   for (size_t i = 0; i < mcParticleVec.size(); ++i) {
    //     LCRelationImpl* relationTrack = new LCRelationImpl(lcioTrkPtr, mcParticleVec[i], weightVec[i]);
    //     trackRelationCollection->addElement(relationTrack);
    //   }
    // }

  } // for loop to the tracks

  // TODO:
  // if (input_rel_col) {
  //   evt->addCollection(trackRelationCollection, _output_track_rel_name);
  // }
  return std::make_tuple(std::move(trackVec), std::move(trackRelationCollection));
}

int RefitFinal::FitInit2(const edm4hep::Track& track, GaudiDDKalTestTrack& marlinTrk) const {
  edm4hep::TrackState trackState;

  size_t refPoint = m_refPoint.value() == -1 ? 0 : static_cast<size_t>(m_refPoint.value());
  if (refPoint >= track.getTrackStates().size()) {
    error() << "Cannot find trackstate for " << m_refPoint << endmsg;
    return 1;
  }
  trackState = track.getTrackStates()[refPoint];

  const bool direction = m_extrapolateForward ? true : false;
  marlinTrk.initialise(trackState, direction);

  return 0;
}
