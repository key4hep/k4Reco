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
#include "TruthTrackFinder.h"
#include "GaudiDDKalTestTrack.h"
#include "GaudiTrkUtils.h"

#include <streamlog/logstream.h>
#include <streamlog/streamlog.h>

#include <DDSegmentation/BitFieldCoder.h>
#include <edm4hep/TrackCollection.h>
#include <edm4hep/TrackMCParticleLinkCollection.h>
#include <edm4hep/TrackerHitSimTrackerHitLinkCollection.h>
#include <edm4hep/utils/vector_utils.h>

#include <DD4hep/BitFieldCoder.h>
#include <DD4hep/Detector.h>

#include <algorithm>
#include <cfloat>
#include <climits>
#include <cmath>

/*

 This code performs a true pattern recognition by looping over all MC particles and adding all hits
 associated to them onto a track. This is then fitted and output.

*/

inline bool sort_by_radius(const edm4hep::TrackerHitPlane* hit1, const edm4hep::TrackerHitPlane* hit2) {
  return edm4hep::utils::magnitudeTransverse(hit1->getPosition()) <
         edm4hep::utils::magnitudeTransverse(hit2->getPosition());
}

inline bool sort_by_z(const edm4hep::TrackerHitPlane* hit1, const edm4hep::TrackerHitPlane* hit2) {
  const double z1 = fabs(hit1->getPosition()[2]);
  const double z2 = fabs(hit2->getPosition()[2]);
  return z1 < z2;
}

std::vector<const edm4hep::TrackerHitPlane*>
TruthTrackFinder::removeHitsSameLayer(const std::vector<const edm4hep::TrackerHitPlane*>& trackHits) const {
  std::vector<const edm4hep::TrackerHitPlane*> trackFilteredHits;

  trackFilteredHits.push_back(*(trackHits.begin()));

  auto cellID = trackHits[0]->getCellID();
  // originally "subdet" was used, which translates to 0
  auto subdet = m_encoder.get(cellID, 0);
  // originall "layer" was used, which translates to
  auto layer = m_encoder.get(cellID, "layer");

  for (auto it = trackHits.begin() + 1; it != trackHits.end(); ++it) {
    auto currentCellID = (*it)->getCellID();
    auto currentSubdet = m_encoder.get(currentCellID, 0);
    auto currentLayer = m_encoder.get(currentCellID, "layer");
    if (subdet != currentSubdet) {
      trackFilteredHits.push_back(*it);
    } else if (layer != currentLayer) {
      trackFilteredHits.push_back(*it);
    }
    subdet = currentSubdet;
    layer = currentLayer;
  }
  return trackFilteredHits;
}

TruthTrackFinder::TruthTrackFinder(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(name, svcLoc,
                       {
                           KeyValues("TrackerHitCollectionNames", {"VXDTrackerHits"}),
                           KeyValues("SimTrackerHitRelCollectionNames", {"VXDTrackerHitRelations"}),
                           KeyValues("MCParticleCollectionName", {"MCParticle"}),
                       },
                       {
                           KeyValues("SiTrackCollectionName", {"SiTracks"}),
                           KeyValues("SiTrackRelationCollectionName", {"SiTrackRelations"}),
                       }) {}

StatusCode TruthTrackFinder::initialize() {
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

  m_ddkaltest.init();
  m_ddkaltest.setEncoder(m_encoder);

  // // Get the magnetic field
  dd4hep::Detector* lcdd = m_geoSvc->getDetector();
  const double position[3] = {0, 0, 0};      // position to calculate magnetic field at (the origin in this case)
  double magneticFieldVector[3] = {0, 0, 0}; // initialise object to hold magnetic field
  lcdd->field().magneticField(position, magneticFieldVector); // get the magnetic field vector from DD4hep
  m_magneticField = static_cast<float>(magneticFieldVector[2] / dd4hep::tesla); // z component at (0,0,0)

  return StatusCode::SUCCESS;
}

std::tuple<edm4hep::TrackCollection, edm4hep::TrackMCParticleLinkCollection>
TruthTrackFinder::operator()(const std::vector<const edm4hep::TrackerHitPlaneCollection*>& trackerHitCollections,
                             const std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>& relations,
                             const std::vector<const edm4hep::MCParticleCollection*>& particleCollections) const {
  // // set the correct configuration for the tracking system for this event
  // MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::useQMS> mson(trackFactory, true);
  // MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::usedEdx> elosson(trackFactory, true);
  // MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing> smoothon(trackFactory, false);

  edm4hep::TrackCollection trackCollection;
  edm4hep::TrackMCParticleLinkCollection trackRelationCollection;

  if (particleCollections.size() == 0)
    return std::make_tuple(std::move(trackCollection), std::move(trackRelationCollection));

  /*
Now for each MC particle we want the list of hits belonging to it. The most
efficient way is to loop over all hits once, and store the pointers in a
map, with the key a pointer to the MC particle. We can then loop over each
MC particle at the end and get all of the hits, before making a track.
*/

  // Make the container
  std::map<size_t, std::vector<const edm4hep::TrackerHitPlane*>> particleHits;
  std::vector<edm4hep::TrackerHitPlane> hits;
  hits.reserve(std::accumulate(trackerHitCollections.begin(), trackerHitCollections.end(), 0u,
                               [](size_t sum, const auto& collection) { return sum + collection->size(); }));

  // for (size_t collection = 0; collection < trackerHitCollections.size(); collection++) {
  //   const auto& trackerHitCollection = trackerHitCollections[collection];
  //   info() << "TruthTrackFinder: trackerHitCollection.size(): " << trackerHitCollection->size() << endmsg;
  //   info() << "TruthTrackFinder: relations[collection].size(): " << relations[collection]->size() << endmsg;
  //   for (size_t itHit = 0; itHit < trackerHitCollection->size(); itHit++) {
  //     info() << "TruthTrackFinder: hit cellID: " << (*trackerHitCollection)[itHit].getCellID() << endmsg;
  //     info() << "TruthTrackFinder: simHit cellID: " << relations[collection]->at(itHit).getTo().getCellID() <<
  //     endmsg; if (relations[collection]->at(itHit).getFrom().getCellID() !=
  //     (*trackerHitCollection)[itHit].getCellID())
  //       throw std::runtime_error("CellID mismatch between TrackerHit and SimTrackerHit");
  //     hits.push_back({(*trackerHitCollection)[itHit], relations[collection]->at(itHit).getTo()});
  //   }
  // }

  // Loop over all input collections
  for (size_t collection = 0; collection < relations.size(); collection++) {
    const auto& rel = relations[collection];
    for (unsigned int itRel = 0; itRel < rel->size(); itRel++) {
      const auto& hit = rel->at(itRel).getFrom();
      const auto& simHit = rel->at(itRel).getTo();
      // If the hit was produced by a secondary which was not saved to the MCParticle collection
      if (simHit.isProducedBySecondary())
        continue;

      // Get the particle belonging to that hit
      edm4hep::MCParticle particle = simHit.getParticle();

      // Push back the element into the container
      hits.push_back(hit.as<edm4hep::TrackerHitPlane>());
      particleHits[static_cast<size_t>(particle.id().index)].push_back(&hits.back());
    }
  }

  // Now loop over all particles and get the list of hits
  for (size_t itP = 0; itP < particleCollections[0]->size(); itP++) {
    // Get the vector of hits from the container
    // if (std::ranges::find(particleHits, &mcParticle) == particleHits.end())
    //   continue;
    if (!particleHits.contains(itP)) {
      info() << "Particle not found in particleHits" << endmsg;
      continue;
    }
    auto& trackHits = particleHits[itP];

    // Only make tracks with 3 or more hits
    if (trackHits.size() < 3)
      continue;

    // Sort the hits from smaller to larger radius
    std::sort(trackHits.begin(), trackHits.end(), sort_by_radius);

    // Remove the hits on the same layers (removing those with higher R)
    auto trackFilteredByRHits = removeHitsSameLayer(trackHits);
    if (trackFilteredByRHits.size() < 3)
      continue;

    /*
     Fit - this gets complicated.
     Need to pass a series of objects, including some initial states,
     covariance matrix etc. Set these up, then call the fit.
    */

    // Make the track object and relations object
    edm4hep::MutableTrack track;

    // IMarlinTrk used to fit track - IMarlinTrk interface to separete pattern recogition from fit implementation
    auto marlinTrack = GaudiDDKalTestTrack(this, const_cast<GaudiDDKalTest*>(&m_ddkaltest));
    auto marlinTrackZSort = GaudiDDKalTestTrack(this, const_cast<GaudiDDKalTest*>(&m_ddkaltest));

    // Save a vector of the hits to be used (why is this not attached to the track directly?? MarlinTrkUtils to be
    // updated?)
    std::vector<const edm4hep::TrackerHitPlane*> trackfitHits;
    for (unsigned int itTrackHit = 0; itTrackHit < trackFilteredByRHits.size(); itTrackHit++)
      trackfitHits.push_back(trackFilteredByRHits[itTrackHit]);

    // Make an initial covariance matrix with very broad default values
    edm4hep::CovMatrix6f covMatrix{};
    covMatrix[0] = m_initialTrackError_d0;    // sigma_d0^2
    covMatrix[2] = m_initialTrackError_phi0;  // sigma_phi0^2
    covMatrix[5] = m_initialTrackError_omega; // sigma_omega^2
    covMatrix[9] = m_initialTrackError_z0;    // sigma_z0^2
    covMatrix[14] = m_initialTrackError_tanL; // sigma_tanl^2

    bool direction = ((m_fitForward) ? true : false);

    int fitError = 0;

    GaudiTrkUtils trkUtils(static_cast<const Gaudi::Algorithm*>(this), m_ddkaltest, m_geoSvc,
                           m_encodingStringVariable.value());

    // TODO: Implement code for running when m_useTruthInPrefit is true
    if (m_useTruthInPrefit) {
      throw std::runtime_error("UseTruthInPrefit not implemented yet");
      //   HelixTrack helix(mcParticle->getVertex(), mcParticle->getMomentum(), mcParticle->getCharge(),
      //   m_magneticField);

      //   double trueD0        = helix.getD0();
      //   double truePhi       = helix.getPhi0();
      //   double trueOmega     = helix.getOmega();
      //   double trueZ0        = helix.getZ0();
      //   double trueTanLambda = helix.getTanLambda();

      //   // float ref_point[3] = { 0., 0., 0. };
      //   helix.moveRefPoint(trackfitHits.at(0)->getPosition()[0], trackfitHits.at(0)->getPosition()[1],
      //                      trackfitHits.at(0)->getPosition()[2]);
      //   float ref_point[3] = {float(helix.getRefPointX()), float(helix.getRefPointY()), float(helix.getRefPointZ())};
      //   TrackStateImpl* trkState =
      //       new TrackStateImpl(TrackState::AtIP, trueD0, truePhi, trueOmega, trueZ0, trueTanLambda, covMatrix,
      //       ref_point);

      //   // int prefitError =  createFit(trackfitHits, marlinTrack, trkState, m_magneticField,  direction,
      //   // m_maxChi2perHit); streamlog_out(DEBUG2) << "---- createFit - error_fit = " << error_fit << endmsg;

      //   // if (prefitError == 0) {
      //   //   fitError = finaliseLCIOTrack(marlinTrack, track, trackfitHits,  direction );
      //   //   streamlog_out(DEBUG2) << "---- finalisedLCIOTrack - error = " << error << endmsg;
      //   // }

      //   fitError = MarlinTrk::createFinalisedLCIOTrack(marlinTrack, trackfitHits, track, direction, trkState,
      //                                                  m_magneticField, m_maxChi2perHit);

      //   // If first fit attempt fails, try a new fit with hits ordered by z

      //   if (fitError != 0) {
      //     // Sort the hits from smaller to larger z
      //     std::sort(trackfitHits.begin(), trackfitHits.end(), sort_by_z);

      //     // Removing the hits on the same layers (remove those with higher z)
      //     EVENT::TrackerHitVec trackFilteredByZHits;
      //     removeHitsSameLayer(trackfitHits, trackFilteredByZHits);
      //     if (trackFilteredByZHits.size() < 3)
      //       continue;

      //     // If fit with hits ordered by radius has failed, the track is probably a 'spiral' track.
      //     // Fitting 'backward' is very difficult for spiral track, so the default direction here is set as 'forward'
      //     fitError = MarlinTrk::createFinalisedLCIOTrack(marlinTrackZSort, trackFilteredByZHits, track,
      //                                                    MarlinTrk::IMarlinTrack::forward, trkState, m_magneticField,
      //                                                    m_maxChi2perHit);
      //   }

      //   delete trkState;

      //   ////////////////////////

    } // end: helical prefit initialised with info from truth and then track fitted and saved as a lcio track

    else {
      // DEFAULT procedure: Try to fit
      fitError = trkUtils.createFinalisedLCIOTrack(marlinTrack, trackfitHits, track, direction, covMatrix,
                                                   m_magneticField, m_maxChi2perHit);

      // If first fit attempt fails, try a new fit with hits ordered by z

      if (fitError != 0) {
        // we need to clean the track object
        track = edm4hep::MutableTrack();

        // Sort the hits from smaller to larger z
        std::sort(trackfitHits.begin(), trackfitHits.end(), sort_by_z);

        // Removing the hits on the same layers (remove those with higher z)
        const auto trackFilteredByZHits = removeHitsSameLayer(trackfitHits);
        if (trackFilteredByZHits.size() < 3)
          continue;

        // If fit with hits ordered by radius has failed, the track is probably a 'spiral' track.
        // Fitting 'backward' is very difficult for spiral track, so the default directiin here is set as 'forward'
        fitError = trkUtils.createFinalisedLCIOTrack(marlinTrackZSort, trackFilteredByZHits, track, true, covMatrix,
                                                     m_magneticField, m_maxChi2perHit);
      }

    } // end: track fitted (prefit from fit from 3 hits) and saved as a lcio track

    /////////////////////////////////////////////

    debug() << "TruthTrackFinder: fitError " << fitError << endmsg;

    // Check track quality - if fit fails chi2 will be 0
    if (fitError != 0 || track.getChi2() <= 0 || track.getNdf() <= 0) {
      continue;
    }

    const auto hits_in_fit = marlinTrack.getHitsInFit();
    std::vector<const edm4hep::TrackerHitPlane*> hits_in_fit_ptr;

    for (const auto& hit : hits_in_fit) {
      hits_in_fit_ptr.push_back(hit.first);
    }

    /// Fill hits associated to the track by pattern recognition and hits in fit
    std::vector<int32_t> subdetectorHitNumbers;
    trkUtils.addHitNumbersToTrack(subdetectorHitNumbers, trackHits, false, m_encoder);
    trkUtils.addHitNumbersToTrack(subdetectorHitNumbers, hits_in_fit_ptr, true, m_encoder);
    for (const auto num : subdetectorHitNumbers) {
      track.addToSubdetectorHitNumbers(num);
    }

    debug() << "TruthTrackFinder: trackHits.size(): " << trackHits.size()
            << " trackfitHits.size(): " << trackfitHits.size() << " hits_in_fit.size(): " << hits_in_fit.size()
            << endmsg;

    // Push back to the output container
    trackCollection.push_back(track);

    // Make the particle to track link
    auto link = trackRelationCollection.create();
    link.setFrom(track);
    link.setTo((*particleCollections[0])[itP]);
    link.setWeight(1.0);
  }

  return std::make_tuple(std::move(trackCollection), std::move(trackRelationCollection));
}
