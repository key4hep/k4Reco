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
    : MultiTransformer(
          name, svcLoc,
          {
              KeyValue("InputTrackCollectionName", "TruthTracks"),
              KeyValues("InputRelationCollectionName", {"SiTrackRelations"}), // KeyValues so that it can be empty
          },
          {
              KeyValue("OutputTrackCollectionName", "RefittedTracks"),
              KeyValue("OutputRelationCollectionName", "RefittedRelation"),
          }) {}

StatusCode RefitFinal::initialize() {
  // Setting the streamlog output is necessary to avoid lots of overhead.
  // Otherwise it would be equivalent to running with every debug message
  // being computed
  streamlog::out.init(std::cout, "");
  streamlog::logscope* scope = new streamlog::logscope(streamlog::out);
  setStreamlogOutputLevel(this, scope);

  m_geoSvc = serviceLocator()->service(m_geoSvcName);
  if (!m_geoSvc) {
    error() << "Unable to retrieve GeoSvc" << endmsg;
    return StatusCode::FAILURE;
  }
  std::string cellIDEncodingString = m_geoSvc->constantAsString(m_encodingStringVariable.value());
  m_encoder = dd4hep::DDSegmentation::BitFieldCoder(cellIDEncodingString);

  m_ddkaltest.init(m_MSOn, m_ElossOn);
  m_ddkaltest.setEncoder(m_encoder);

  return StatusCode::SUCCESS;
}

std::tuple<edm4hep::TrackCollection, edm4hep::TrackMCParticleLinkCollection> RefitFinal::operator()(
  const edm4hep::TrackCollection& input_track_col,
  const std::vector<const edm4hep::TrackMCParticleLinkCollection*>& input_rel_col
) const {
  // establish the track collection that will be created
  edm4hep::TrackCollection trackVec;
  edm4hep::TrackMCParticleLinkCollection trackRelationCollection;

  if (input_rel_col.size() == 0) {
    debug() << "No input relation collection, not creating one either" << endmsg;
  }

  std::map<int, edm4hep::MCParticle> trackIndexToMCParticle;
  if (!input_rel_col.empty()) {
    const auto& relationCol = *input_rel_col[0];
    for (const auto& rel : relationCol) {
      edm4hep::Track trk = rel.getFrom();
      trackIndexToMCParticle[trk.getObjectID().index] = rel.getTo();
    }
  }

  const size_t nTracks = input_track_col.size();

  debug() << " Number of Tracks " << nTracks << endmsg;

  // loop over the input tracks and refit
  int counter = 0;
  std::map<int, int> hitInSubDet;
  for (size_t iTrack = 0; iTrack < nTracks; ++iTrack) {
    const auto& track = input_track_col.at(iTrack);
    const auto& trkHits = track.getTrackerHits();
    // TODO: Don't create a copy because finaliseLCIOTrack and marlin_trk need pointers
    std::vector<const edm4hep::TrackerHit*> trkHitsPtr;
    trkHitsPtr.reserve(trkHits.size());
    for (const auto& hit : trkHits) {
      trkHitsPtr.push_back(&hit);
    }

    auto gaudi_trk = GaudiDDKalTestTrack(this, const_cast<GaudiDDKalTest*>(&m_ddkaltest));

    debug() << "---- track n = " << iTrack << "  n hits = " << trkHits.size() << endmsg;
    if (trkHits.size() < 3) {
      debug() << "Track " << iTrack << " has less than 3 hits, skipping" << endmsg;
      continue;
    }

    hitInSubDet.clear();

    for (const auto& ptr : trkHitsPtr) {
      gaudi_trk.addHit(ptr);
      ++hitInSubDet[m_encoder.get(ptr->getCellID(), "system")];
    }


    edm4hep::MutableTrack edm4hep_trk = trackVec.create();

    edm4hep::CovMatrix6f initialCov;
    initialCov[0] = m_initialTrackError_d0;
    initialCov[2] = m_initialTrackError_phi0;
    initialCov[5] = m_initialTrackError_omega;
    initialCov[9] = m_initialTrackError_z0;
    initialCov[14] = m_initialTrackError_tanL;

    const bool fit_direction = m_extrapolateForward;

    GaudiTrkUtils trkUtils(
      static_cast<const Gaudi::Algorithm*>(this), 
      m_ddkaltest, m_geoSvc, 
      m_encodingStringVariable.value()
    );

    int return_code = trkUtils.createFinalisedLCIOTrack(
      gaudi_trk, trkHitsPtr, edm4hep_trk, 
      fit_direction, initialCov, m_bField, 
      m_Max_Chi2_Incr.value()
    );

    if (return_code != 0) {
      debug() << "finaliseLCIOTrack failed" << endmsg;
      continue;
    }

    debug() << " *** created finalized LCIO track - return code " << return_code << endmsg;

    // fit finished - get hits in the fit
    // remember the hits are ordered in the order in which they were fitted

    const auto hits_in_fit = gaudi_trk.getHitsInFit();

    if (int(hits_in_fit.size()) < m_minClustersOnTrackAfterFit) {
      debug() << "Less than " << m_minClustersOnTrackAfterFit
              << " hits in fit: Track "
                 "Discarded. Number of hits =  "
              << trkHits.size() << endmsg;
      continue;
    }

    const auto outliers = gaudi_trk.getOutliers();

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
      edm4hep_trk.addToSubdetectorHitNumbers(num);
    }

    // debug() << "processEvent: Hit numbers for track " << edm4hep_trk.id() << ":  " << endmsg;
    int detID = 0;
    for (size_t ip = 0; ip < edm4hep_trk.getSubdetectorHitNumbers().size(); ip = ip + 2) {
      detID++;
      debug() << "  det id " << detID << " , nhits in track = " << edm4hep_trk.getSubdetectorHitNumbers()[ip]
              << " , nhits in fit = " << edm4hep_trk.getSubdetectorHitNumbers()[ip + 1] << endmsg;
      if (edm4hep_trk.getSubdetectorHitNumbers()[ip] > 0)
        // TODO: is detID - 1 correct?
        edm4hep_trk.setType(edm4hep_trk.getType() | (1 << detID));
    }

    // if required apply the ReducedChi2 cut
    if (m_ReducedChi2Cut > 0.0 && edm4hep_trk.getChi2() / edm4hep_trk.getNdf() > m_ReducedChi2Cut) {
      debug() << "Track Discarded due to ReducedChi2 cut:  Chi2 = " << edm4hep_trk.getChi2()
              << " Ndf = " << edm4hep_trk.getNdf() << " ReducedChi2 = " << edm4hep_trk.getChi2() / edm4hep_trk.getNdf()
              << endmsg;
      ++counter;
      continue;
    }

    //trackVec.push_back(lcio_trk);
    // create the Track-MCParticle relation
    edm4hep::MutableTrackMCParticleLink relation = trackRelationCollection.create();
    relation.setFrom(edm4hep_trk);
    relation.setTo(trackIndexToMCParticle.at(static_cast<int>(iTrack)));
  } // for loop to the tracks

  debug() << "Final number of Tracks after refit = " << trackVec.size()
          << " Number of discarded Tracks = " << counter << endmsg;

  return std::make_tuple(std::move(trackVec), std::move(trackRelationCollection));
}
