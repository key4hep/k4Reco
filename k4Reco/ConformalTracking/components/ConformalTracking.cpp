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
#include "ConformalTracking.h"

#include "GaudiDDKalTest.h"
#include "GaudiDDKalTestTrack.h"
#include "GaudiTrkUtils.h"

#include "Cell.h"
#include "KDCluster.h"
#include "KDTrack.h"
#include "KDTree.h"

#include <streamlog/logstream.h>
#include <streamlog/streamlog.h>

#include <edm4hep/CovMatrix6f.h>
#include <edm4hep/SimTrackerHit.h>
#include <edm4hep/Track.h>
#include <edm4hep/TrackState.h>
#include <edm4hep/TrackerHitPlaneCollection.h>
#include <edm4hep/Vector3d.h>
#include <edm4hep/utils/vector_utils.h>

#include <k4Interface/IGeoSvc.h>

#include <DD4hep/DD4hepUnits.h>
#include <DD4hep/Detector.h>
#include <DDRec/Material.h>
#include <DDRec/SurfaceManager.h>

#include <GaudiKernel/MsgStream.h>

#include <TFile.h>
#include <TLine.h>
#include <TMath.h>
#include <TStopwatch.h>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>

bool sort_by_radius(const edm4hep::TrackerHitPlane& hit1, const edm4hep::TrackerHitPlane& hit2);
bool sort_by_radius(const edm4hep::TrackerHitPlane* hit1, const edm4hep::TrackerHitPlane* hit2);

// Sort kd hits from smaller to larger radius
inline bool sort_by_lower_radiusKD(const SKDCluster& hit1, const SKDCluster& hit2) {
  return hit1->getR() < hit2->getR();
}

ConformalTracking::ConformalTracking(const std::string& name, ISvcLocator* svcLoc)
    : Transformer(name, svcLoc,
                  {
                      KeyValues("TrackerHitCollectionNames", {}),
                      KeyValues("MCParticleCollectionName", {"MCParticle"}),
                      KeyValues("RelationsNames", {}),
                  },
                  {
                      KeyValues("SiTrackCollectionName", {"CATracks"}),
                  }) {}

StatusCode ConformalTracking::initialize() {
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

  const auto& locs = inputLocations(0);

  std::vector<int> m_vertexBarrelHits, m_vertexEndcapHits, m_vertexCombinedHits, m_trackerHits, m_allHits;

  for (const auto& str : m_inputVertexBarrelCollections) {
    auto it = std::find(locs.begin(), locs.end(), str);
    if (it == locs.end()) {
      throw std::runtime_error("Collection " + str + " not found in input collections");
    }
    m_vertexBarrelHits.push_back(std::distance(locs.begin(), it));
  }

  for (const auto& str : m_inputVertexEndcapCollections) {
    auto it = std::find(locs.begin(), locs.end(), str);
    if (it == locs.end()) {
      throw std::runtime_error("Collection " + str + " not found in input collections");
    }
    m_vertexEndcapHits.push_back(std::distance(locs.begin(), it));
  }

  m_vertexCombinedHits = m_vertexBarrelHits;
  m_vertexCombinedHits.insert(m_vertexCombinedHits.end(), m_vertexEndcapHits.begin(), m_vertexEndcapHits.end());

  for (const auto& str : m_inputMainTrackerHitCollections) {
    auto it = std::find(locs.begin(), locs.end(), str);
    if (it == locs.end()) {
      throw std::runtime_error("Collection " + str + " not found in input collections");
    }
    m_trackerHits.push_back(std::distance(locs.begin(), it));
  }

  for (size_t i = 0; i < locs.size(); ++i) {
    m_allHits.push_back(i);
  }

  // Default parsing
  if (m_stepCollections.empty()) {
    int step = 0;
    // Build tracks in the vertex barrel
    m_stepParameters.emplace_back(m_vertexBarrelHits, m_maxCellAngle, m_maxCellAngleRZ, m_chi2cut, m_minClustersOnTrack,
                                  m_maxDistance, m_slopeZRange, m_highPTcut, /*highPT*/ true, /*OnlyZS*/ false,
                                  /*rSearch*/ false,
                                  /*vtt*/ true, /*kalmanFitForward*/ true, step++,
                                  /*combine*/ true, /*build*/ true, /*extend*/ false, /*sort*/ false);
    // Extend through the endcap
    m_stepParameters.emplace_back(m_vertexEndcapHits, m_maxCellAngle, m_maxCellAngleRZ, m_chi2cut, m_minClustersOnTrack,
                                  m_maxDistance, m_slopeZRange, m_highPTcut, /*highPT*/ true, /*OnlyZS*/ false,
                                  /*rSearch*/ false, /*vtt*/ true,
                                  /*kalmanFitForward*/ true, step++,
                                  /*combine*/ true, /*build*/ false, /*extend*/ true, /*sort*/ false);
    // Make combined vertex tracks
    m_stepParameters.emplace_back(m_vertexCombinedHits, m_maxCellAngle, m_maxCellAngleRZ, m_chi2cut,
                                  m_minClustersOnTrack, m_maxDistance, m_slopeZRange, m_highPTcut, /*highPT*/ true,
                                  /*OnlyZS*/ false,
                                  /*rSearch*/ false, /*vtt*/ true,
                                  /*kalmanFitForward*/ true, step++,
                                  /*combine*/ true, /*build*/ true, /*extend*/ false, /*sort*/ false);
    // Make leftover tracks in the vertex with lower requirements
    // 1. open the cell angles
    if (m_enableTCVC) {
      m_stepParameters.emplace_back(m_vertexCombinedHits, m_maxCellAngle * 5.0, m_maxCellAngleRZ * 5.0, m_chi2cut,
                                    m_minClustersOnTrack, m_maxDistance, m_slopeZRange, m_highPTcut,
                                    /*highPT*/ true, /*OnlyZS*/ false,
                                    /*rSearch*/ true, /*vtt*/ true, /*kalmanFitForward*/ true, step++,
                                    /*combine*/ not m_enableTCVC, /*build*/ true, /*extend*/ false, /*sort*/ false);
    }
    // 2. open further the cell angles and increase the chi2cut
    m_stepParameters.emplace_back(std::vector<int>{}, m_maxCellAngle * 10.0, m_maxCellAngleRZ * 10.0, m_chi2cut * 20.0,
                                  m_minClustersOnTrack, m_maxDistance, m_slopeZRange, m_highPTcut,
                                  /*highPT*/ true, /*OnlyZS*/ false,
                                  /*rSearch*/ true, /*vtt*/ true, /*kalmanFitForward*/ true, step++,
                                  /*combine*/ false, /*build*/ true, /*extend*/ false, /*sort*/ false);
    // 3. min number of hits on the track = 4

    m_stepParameters.emplace_back(std::vector<int>{}, m_maxCellAngle * 10.0, m_maxCellAngleRZ * 10.0, m_chi2cut * 20.0,
                                  /*m_minClustersOnTrack*/ 4, m_maxDistance, m_slopeZRange, m_highPTcut,
                                  /*highPT*/ true,
                                  /*OnlyZS*/ false,
                                  /*rSearch*/ true, /*vtt*/ true, /*kalmanFitForward*/ true, step++,
                                  /*combine*/ false, /*build*/ true, /*extend*/ false, /*sort*/ true);
    // Extend through inner and outer trackers
    m_stepParameters.emplace_back(m_trackerHits, m_maxCellAngle * 10.0, m_maxCellAngleRZ * 10.0, m_chi2cut * 20.0,
                                  /*m_minClustersOnTrack*/ 4, m_maxDistance, m_slopeZRange, /*highPTcut*/ 1.0,
                                  /*highPT*/ true,
                                  /*OnlyZS*/ false,
                                  /*rSearch*/ true, /*vtt*/ true, /*kalmanFitForward*/ true, step++,
                                  /*combine*/ true, /*build*/ false, /*extend*/ true, /*sort*/ false);
    // Finally reconstruct displaced tracks
    m_stepParameters.emplace_back(m_allHits, m_maxCellAngle * 10.0, m_maxCellAngleRZ * 10.0, m_chi2cut * 10.0,
                                  /*m_minClustersOnTrack*/ 5, 0.015, m_slopeZRange, m_highPTcut,
                                  /*highPT*/ false, /*OnlyZS*/ true,
                                  /*rSearch*/ true,
                                  /*vtt*/ false, /*kalmanFitForward*/ true, step++,
                                  /*combine*/ true, /*build*/ true, /*extend*/ false, /*sort*/ false);
  } else {
    for (size_t i = 0; i < m_stepCollections.size(); i++) {
      const auto& collections = m_stepCollections[i];
      const auto& parNames = m_stepParametersNames[i];
      const auto& parValues = m_stepParametersValues[i];
      const auto& flags = m_stepParametersFlags[i];
      const auto& functions = m_stepParametersFunctions[i];

      std::vector<int> indexes;
      for (const auto& str : collections) {
        auto it = std::find(locs.begin(), locs.end(), str);
        if (it == locs.end()) {
          throw std::runtime_error("Collection " + str + " not found in input collections");
        }
        indexes.push_back(std::distance(locs.begin(), it));
      }
      auto maxCellAngle = m_maxCellAngle;
      auto maxCellAngleRZ = m_maxCellAngleRZ;
      auto chi2cut = m_chi2cut;
      auto minClustersOnTrack = m_minClustersOnTrack;
      auto maxDistance = m_maxDistance;
      auto slopeZRange = m_slopeZRange;
      auto highPTcut = m_highPTcut;

      auto highPT = false;
      auto OnlyZS = false;
      auto rSearch = false;
      auto vtt = false;
      auto kalmanFitForward = false;

      auto combine = false;
      auto build = false;
      auto extend = false;
      auto sort = false;

      for (size_t j = 0; j < parNames.size(); j++) {
        if (parNames[j] == "MaxCellAngle") {
          maxCellAngle = parValues[j];
        } else if (parNames[j] == "MaxCellAngleRZ") {
          maxCellAngleRZ = parValues[j];
        } else if (parNames[j] == "Chi2Cut") {
          chi2cut = parValues[j];
        } else if (parNames[j] == "MinClustersOnTrack") {
          minClustersOnTrack = parValues[j];
        } else if (parNames[j] == "MaxDistance") {
          maxDistance = parValues[j];
        } else if (parNames[j] == "SlopeZRange") {
          slopeZRange = parValues[j];
        } else if (parNames[j] == "HighPTCut") {
          highPTcut = parValues[j];
        }
      }

      highPT = std::find(flags.begin(), flags.end(), "HighPTFit") != flags.end();
      OnlyZS = std::find(flags.begin(), flags.end(), "OnlyZSchi2cut") != flags.end();
      rSearch = std::find(flags.begin(), flags.end(), "RadialSearch") != flags.end();
      vtt = std::find(flags.begin(), flags.end(), "VertexToTracker") != flags.end();
      kalmanFitForward = std::find(flags.begin(), flags.end(), "KalmanFitForward") != flags.end() ||
                         !(std::find(flags.begin(), flags.end(), "KalmanFitBackward") != flags.end());

      combine = std::find(functions.begin(), functions.end(), "CombineCollections") != functions.end();
      build = std::find(functions.begin(), functions.end(), "BuildNewTracks") != functions.end();
      extend = std::find(functions.begin(), functions.end(), "ExtendTracks") != functions.end();
      sort = std::find(functions.begin(), functions.end(), "SortTracks") != functions.end();

      m_stepParameters.emplace_back(indexes, maxCellAngle, maxCellAngleRZ, chi2cut, minClustersOnTrack, maxDistance,
                                    slopeZRange, highPTcut, highPT, OnlyZS, rSearch, vtt, kalmanFitForward, i, combine,
                                    build, extend, sort);
    }
  }

  m_ddkaltest.init();
  m_ddkaltest.setEncoder(m_encoder);

  // Default values for track fitting
  m_initialTrackError_d0 = 1e6;
  m_initialTrackError_phi0 = 1e2;
  m_initialTrackError_omega = 1e-4;
  m_initialTrackError_z0 = 1e6;
  m_initialTrackError_tanL = 1e2;
  m_maxChi2perHit = 1e2;

  m_magneticField = getBzAtOrigin(); // z component at (0,0,0)

  if (m_debugPlots) {
    //   m_szDistribution  = new TH2F("m_szDistribution", "m_szDistribution", 20000, -100, 100, 200, -10, 10);
    //   m_uvDistribution  = new TH2F("m_uvDistribution", "m_uvDistribution", 1000, -0.05, 0.05, 1000, -0.05, 0.05);
    //   m_xyDistribution  = new TH2F("m_xyDistribution", "m_xyDistribution", 500, -1500, 1500, 500, -1500, 1500);
    //   m_xyzDistribution = new TH3F("m_xyzDistribution", "m_xyzDistribution", 50, 0, 100, 50, 0, 100, 100, 0, 25);

    //   // Histograms for neighbors parameters
    m_neighX = new TH1F("m_neighX", "m_neighX", 500, -1500, 1500);
    m_neighY = new TH1F("m_neighY", "m_neighY", 500, -1500, 1500);
    m_neighZ = new TH1F("m_neighZ", "m_neighZ", 500, -2500, 2500);

    m_slopeZ = new TH1F("m_slopeZ", "m_slopeZ", 1000, -100, 100);
    m_slopeZ_true = new TH1F("m_slopeZ_true", "m_slopeZ_true", 1000, -100, 100);
    m_slopeZ_true_first = new TH1F("m_slopeZ_true_first", "m_slopeZ_true_first", 1000, -100, 100);
    m_slopeZ_vs_pt_true = new TH2F("m_diffZ_pt_true_first", "m_diffZ_pt_true_first", 2000, -500, 500, 400, 0, 100);

    // Histograms for tuning parameters (cell angle cut, cell length cut)
    m_cellAngle = new TH1F("cellAngle", "cellAngle", 1250, 0, 0.05);
    m_cellDOCA = new TH1F("cellDOCA", "cellDOCA", 100., 0, 0.1);
    m_cellAngleRadius = new TH2F("cellAngleRadius", "cellAngleRadius", 400, 0, 0.04, 1000, 0, 0.04);
    m_cellLengthRadius = new TH2F("cellLengthRadius", "cellLengthRadius", 300, 0, 0.03, 1000, 0, 0.04);
    m_cellAngleLength = new TH2F("cellAngleLength", "cellAngleLength", 400, 0, 0.04, 300, 0, 0.03);
    m_conformalChi2 = new TH1F("conformalChi2", "conformalChi2", 100, 0, 100);
    m_conformalChi2real = new TH1F("conformalChi2real", "conformalChi2real", 1000, 0, 1000);
    m_conformalChi2fake = new TH1F("conformalChi2fake", "conformalChi2fake", 1000, 0, 1000);
    m_conformalChi2Purity = new TH2F("conformalChi2Purity", "conformalChi2Purity", 150, 0, 1.5, 1000, 0, 1000);

    m_conformalChi2MC = new TH1F("conformalChi2MC", "conformalChi2MC", 1000, 0, 1000);
    m_conformalChi2PtMC = new TH2F("conformalChi2PtMC", "conformalChi2PtMC", 1000, 0, 1000, 1000, 0, 100);
    m_conformalChi2VertexRMC = new TH2F("conformalChi2VertexRMC", "conformalChi2VertexRMC", 1000, 0, 1000, 100, 0, 100);

    m_conformalChi2SzMC = new TH1F("conformalChi2SzMC", "conformalChi2SzMC", 1000, 0, 1000);
    m_conformalChi2SzPtMC = new TH2F("conformalChi2SzPtMC", "conformalChi2SzPtMC", 1000, 0, 1000, 1000, 0, 100);
    m_conformalChi2SzVertexRMC =
        new TH2F("conformalChi2SzVertexRMC", "conformalChi2SzVertexRMC", 1000, 0, 1000, 100, 0, 100);

    m_cellAngleMC = new TH1F("cellAngleMC", "cellAngleMC", 1250, 0, 0.05);
    m_cellDOCAMC = new TH1F("cellDOCAMC", "cellDOCAMC", 100., 0, 0.1);
    m_cellAngleRadiusMC = new TH2F("cellAngleRadiusMC", "cellAngleRadiusMC", 400, 0, 0.04, 1000, 0, 0.04);
    m_cellLengthRadiusMC = new TH2F("cellLengthRadiusMC", "cellLengthRadiusMC", 300, 0, 0.03, 1000, 0, 0.04);
    m_cellAngleLengthMC = new TH2F("cellAngleLengthMC", "cellAngleLengthMC", 400, 0, 0.04, 300, 0, 0.03);

    m_cellAngleRZMC = new TH1F("cellAngleRZMC", "cellAngleRZMC", 1250, 0, 0.05);

    // Histograms for "event display"
    m_conformalEvents = new TH2F("conformalEvents", "conformalEvents", 1000, -0.05, 0.05, 1000, -0.05, 0.05);
    m_nonconformalEvents = new TH2F("nonconformalEvents", "nonconformalEvents", 500, -1500, 1500, 500, -1500, 1500);
    m_conformalEventsRTheta =
        new TH2F("conformalEventsRTheta", "conformalEventsRTheta", 200, 0, 0.05, 632, -0.02, 6.30);
    m_conformalEventsMC = new TH2F("conformalEventsMC", "conformalEventsMC", 1000, -0.05, 0.05, 1000, -0.05, 0.05);

    m_canvConformalEventDisplay = new TCanvas("canvConformalEventDisplay", "canvConformalEventDisplay");
    m_canvConformalEventDisplayAllCells =
        new TCanvas("canvConformalEventDisplayAllCells", "canvConformalEventDisplayAllCells");
    m_canvConformalEventDisplayAcceptedCells =
        new TCanvas("canvConformalEventDisplayAcceptedCells", "canvConformalEventDisplayAcceptedCells");
    m_canvConformalEventDisplayMC = new TCanvas("canvConformalEventDisplayMC", "canvConformalEventDisplayMC");
    m_canvConformalEventDisplayMCunreconstructed =
        new TCanvas("canvConformalEventDisplayMCunreconstructed", "canvConformalEventDisplayMCunreconstructed");
  }

  return StatusCode::SUCCESS;
}

edm4hep::TrackCollection ConformalTracking::operator()(

    const std::vector<const edm4hep::TrackerHitPlaneCollection*>& trackerHits,
    const std::vector<const edm4hep::MCParticleCollection*>&,
    const std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>&) const {
  //------------------------------------------------------------------------------------------------------------------
  // This pattern recognition algorithm is based on two concepts: conformal mapping and cellular automaton. Broadly
  // speaking, the 2D xy projection of all hits is transformed such that circles (helix projections) become straight
  // lines. The tracking is then considered as a 2D straight line search, using the z information to reduce
  // combinatorics.
  //
  // The hits from each input collection are transformed into conformal hit positions, with some binary trees created
  // to allow fast nearest neighbour calculations. All hits are then considered as seeds (starting from outer radius)
  // and an attempt to make cells leading back to this seed is carried out. If a long enough chain of cells can be
  // produced, this is defined as a track and the hits contained are all removed from further consideration. Once all
  // hits in a collection have been considered, the next collection of hits is added to the unused hits and the search
  // for new tracks begin again (first an attempt to extend existing tracks is performed, followed by a new search
  // using all unused hits as seeding points).
  //
  // Where several paths are possible back to the seed position, the candidate with lowest chi2/ndof is chosen.
  //------------------------------------------------------------------------------------------------------------------

  auto outputTrackCollection = edm4hep::TrackCollection();
  auto debugHitCollection = edm4hep::TrackerHitPlaneCollection();
  debugHitCollection.setSubsetCollection();

  // Make the collection of conformal hits that will be used, with a link back to
  // the corresponding tracker hit.

  // Collections to be stored throughout the tracking
  std::map<size_t, SharedKDClusters> collectionClusters;       // Conformal hits
  std::map<SKDCluster, edm4hep::TrackerHitPlane> kdClusterMap; // Their link to "real" hits

  // Debug collections (not filled if debug off)
  // std::map<MCParticle*, SharedKDClusters> particleHits;   // List of conformal hits on each MC particle
  // std::map<MCParticle*, bool>             reconstructed;  // Check for MC particles
  // SharedKDClusters                        debugHits;      // Debug hits for plotting

  // Create the conformal hit collections for each tracker hit collection (and save the link)
  for (size_t iColl = 0; iColl < trackerHits.size(); iColl++) {
    const auto& collection = trackerHits[iColl];
    info() << "Processing collection " << inputLocations("TrackerHitCollectionNames")[iColl] << "with size "
           << collection->size() << endmsg;
    // Loop over tracker hits and make conformal hit collection
    SharedKDClusters tempClusters;
    for (size_t iHit = 0; iHit < collection->size(); iHit++) {
      const auto& hit = (*collection)[iHit];
      // Get subdetector information and check if the hit is in the barrel or endcaps
      // Hardcoded from the ILD detector
      // int  subdet   = m_encoder[lcio::LCTrackerCellID::subdet()];
      // int  side     = m_encoder[lcio::LCTrackerCellID::side()];
      // int  layer    = m_encoder[lcio::LCTrackerCellID::layer()];
      // int  module   = m_encoder[lcio::LCTrackerCellID::module()];
      // int  sensor   = m_encoder[lcio::LCTrackerCellID::sensor()];
      const auto celId = hit.getCellID();
      const auto subdet = m_encoder.get(celId, 0);
      const auto side = m_encoder.get(celId, 1);
      const auto layer = m_encoder.get(celId, 2);
      const auto module = m_encoder.get(celId, 3);
      const auto sensor = m_encoder.get(celId, 4);
      bool isEndcap = false;
      bool forward = false;

      // TODO: Hardcoded from the ILD detector
      // if (side != ILDDetID::barrel) {
      //   if (side == ILDDetID::fwd)
      if (side != 0) {
        isEndcap = true;
        if (side == 1)
          forward = true;
      }

      // Make a new kd cluster
      auto kdhit = std::make_shared<KDCluster>(hit, isEndcap, forward);

      // Set the subdetector information
      kdhit->setDetectorInfo(subdet, side, layer, module, sensor);

      // Store the link between the two
      kdClusterMap.emplace(kdhit, hit);
      tempClusters.push_back(kdhit);

      // Store the MC link if in debug mode
      // if (m_debugPlots) {
      //   // Get the related simulated hit(s)
      //   const LCObjectVec& simHitVector = relations[collection]->getRelatedToObjects(hit);
      //   // Take the first hit only (TODO: this should be changed? Loop over all related simHits and add an entry for
      //   each mcparticle so that this hit is in each fit?) SimTrackerHit* simHit =
      //   dynamic_cast<SimTrackerHit*>(simHitVector.at(0));
      //   // Get the particle belonging to that hit
      //   MCParticle* particle = simHit->getMCParticle();
      //   // Store the information (not for secondaries)
      //   kdParticles[kdhit] = particle;
      //   kdSimHits[kdhit]   = simHit;
      //   if (!simHit->isProducedBySecondary()) {
      //     particleHits[particle].push_back(kdhit);
      //   }
      //   m_debugger.registerHit(collection, kdhit, hit);
      //   // Draw plots for event 0
      //   if (m_eventNumber == 0) {
      //     m_conformalEvents->Fill(kdhit->getU(), kdhit->getV());
      //     m_nonconformalEvents->Fill(hit->getPosition()[0], hit->getPosition()[1]);
      //     m_conformalEventsRTheta->Fill(kdhit->getR(), kdhit->getTheta());
      //   }
      // }
    }
    collectionClusters[iColl] = tempClusters;
  }

  // WHAT TO DO ABOUT THIS?? POSSIBLY MOVE DEPENDING ON MC RECONSTRUCTION (and in fact, would fit better into the check
  // reconstruction code at present)

  // Now loop over all MC particles and make the cells connecting hits
  // if (m_debugPlots && MCParticles.size() > 0) {
  //   for (const auto& mcParticle : *MCParticles[0]) {
  //     // Get the vector of hits from the container
  //     if (particleHits.count(mcParticle) == 0)
  //       continue;
  //     SharedKDClusters trackHits = m_debugger.getAssociatedHits(mcParticle);  //particleHits[mcParticle];
  //     // Only make tracks with n or more hits
  //     if (trackHits.size() < (unsigned int)m_minClustersOnTrack)
  //       continue;
  //     // Discard low momentum particles
  //     double particlePt = sqrt(mcParticle->getMomentum()[0] * mcParticle->getMomentum()[0] +
  //                              mcParticle->getMomentum()[1] * mcParticle->getMomentum()[1]);
  //     // Cut on stable particles
  //     if (mcParticle->getGeneratorStatus() != 1)
  //       continue;
  //     // Sort the hits from larger to smaller radius
  //     std::sort(trackHits.begin(), trackHits.end(), sort_by_radiusKD);

  //     // Make a track
  //     auto pars    = _stepParameters[0];
  //     auto mcTrack = std::unique_ptr<KDTrack>(new KDTrack(pars));
  //     // Loop over all hits for debugging
  //     for (auto const& cluster : trackHits) {
  //       // Get the conformal clusters
  //       mcTrack->add(cluster);
  //     }

  //     // Fit the track and plot the chi2
  //     mcTrack->linearRegression();
  //     mcTrack->linearRegressionConformal();
  //     m_conformalChi2MC->Fill(mcTrack->chi2ndof());
  //     m_conformalChi2PtMC->Fill(mcTrack->chi2ndof(), particlePt);
  //     m_conformalChi2SzMC->Fill(mcTrack->chi2ndofZS());
  //     m_conformalChi2SzPtMC->Fill(mcTrack->chi2ndofZS(), particlePt);

  //     double mcVertexX = mcParticle->getVertex()[0];
  //     double mcVertexY = mcParticle->getVertex()[1];
  //     double mcVertexR = sqrt(pow(mcVertexX, 2) + pow(mcVertexY, 2));
  //     m_conformalChi2VertexRMC->Fill(mcTrack->chi2ndof(), mcVertexR);
  //     m_conformalChi2SzVertexRMC->Fill(mcTrack->chi2ndofZS(), mcVertexR);

  //     // Now loop over the hits and make cells - filling histograms along the way
  //     int nHits = trackHits.size();
  //     for (int itHit = 0; itHit < (nHits - 2); itHit++) {
  //       // Get the conformal clusters
  //       SKDCluster cluster0 = trackHits[itHit];
  //       SKDCluster cluster1 = trackHits[itHit + 1];
  //       SKDCluster cluster2 = trackHits[itHit + 2];

  //       // Make the two cells connecting these three hits
  //       auto cell = std::make_shared<Cell>(cluster0, cluster1);
  //       cell->setWeight(itHit);
  //       auto cell1 = std::make_shared<Cell>(cluster1, cluster2);
  //       cell1->setWeight(itHit + 1);

  //       if (itHit == 0)
  //         m_cellDOCAMC->Fill(cell->doca());

  //       // Fill the debug/tuning plots
  //       double angleBetweenCells   = cell->getAngle(cell1);
  //       double angleRZBetweenCells = cell->getAngleRZ(cell1);
  //       double cell0Length = sqrt(pow(cluster0->getU() - cluster1->getU(), 2) + pow(cluster0->getV() -
  //       cluster1->getV(), 2)); double cell1Length = sqrt(pow(cluster1->getU() - cluster2->getU(), 2) +
  //       pow(cluster1->getV() - cluster2->getV(), 2));

  //       m_cellAngleMC->Fill(angleBetweenCells);
  //       m_cellAngleRadiusMC->Fill(cluster2->getR(), angleBetweenCells);
  //       m_cellLengthRadiusMC->Fill(cluster0->getR(), cell0Length);
  //       m_cellAngleLengthMC->Fill(cell1Length, angleBetweenCells);
  //       m_cellAngleRZMC->Fill(angleRZBetweenCells);

  //       // Draw cells on the first event
  //       if (m_eventNumber == 0) {
  //         // Fill the event display (hit positions)
  //         m_conformalEventsMC->Fill(cluster0->getU(), cluster0->getV());
  //         m_conformalEventsMC->Fill(cluster1->getU(), cluster1->getV());
  //         m_conformalEventsMC->Fill(cluster2->getU(), cluster2->getV());
  //         // Draw the cell lines on the event display. Use the line style to show
  //         // if the cells would have been cut by some of the search criteria
  //         m_canvConformalEventDisplayMC->cd();
  //         if (itHit == 0) {
  //           drawline(cluster0, cluster1, itHit + 1);
  //         }
  //         // Draw line style differently if the cell angle was too large
  //         if (angleBetweenCells > (m_maxCellAngle)) {
  //           drawline(cluster1, cluster2, itHit + 2, 3);
  //         } else {
  //           drawline(cluster1, cluster2, itHit + 2);
  //         }
  //       }
  //     }
  //   }
  // }

  // Draw the final set of conformal hits (on top of the cell lines)
  // if (m_debugPlots && m_eventNumber) {
  //   m_canvConformalEventDisplayMC->cd();
  //   m_conformalEventsMC->DrawCopy("same");
  //   // Draw the non-MC event display
  //   m_canvConformalEventDisplay->cd();
  //   m_conformalEvents->DrawCopy("");
  //   m_canvConformalEventDisplayAllCells->cd();
  //   m_conformalEvents->DrawCopy("");
  //   m_canvConformalEventDisplayAcceptedCells->cd();
  //   m_conformalEvents->DrawCopy("");
  //   m_canvConformalEventDisplayMC->cd();
  //   m_conformalEvents->DrawCopy("");
  //   m_canvConformalEventDisplayMCunreconstructed->cd();
  //   m_conformalEvents->DrawCopy("");
  // }

  // Now the track reconstruction strategy. Perform a sequential search, with hits
  // removed from the seeding collections once tracks have been built

  // The final vector of conformal tracks
  UniqueKDTracks conformalTracks;
  SharedKDClusters kdClusters;
  UKDTree nearestNeighbours = nullptr;

  debug() << "m_stepParameters.size() " << m_stepParameters.size() << endmsg;
  for (const auto& parameters : m_stepParameters) {
    debug() << "Parameters : ";
    for (const auto& elem : parameters.m_collections) {
      debug() << " " << elem;
    }
    debug() << endmsg;
    debug() << "maxCellAngle " << parameters.m_maxCellAngle << endmsg;
    debug() << "maxCellAngleRZ " << parameters.m_maxCellAngleRZ << endmsg;
    debug() << "chi2cut " << parameters.m_chi2cut << endmsg;
    debug() << "minClustersOnTrack " << parameters.m_minClustersOnTrack << endmsg;
    debug() << "maxDistance " << parameters.m_maxDistance << endmsg;
    debug() << "maxSlopeZ " << parameters.m_maxSlopeZ << endmsg;
    debug() << "highPTcut " << parameters.m_highPTcut << endmsg;
    debug() << "highPTfit " << parameters.m_highPTfit << endmsg;
    debug() << "onlyZSchi2cut " << parameters.m_onlyZSchi2cut << endmsg;
    debug() << "radialSearch " << parameters.m_radialSearch << endmsg;
    debug() << "vertexToTracker " << parameters.m_vertexToTracker << endmsg;
    debug() << "kalmanFitForward " << parameters.m_kalmanFitForward << endmsg;
    debug() << "step " << parameters.m_step << endmsg;
    debug() << "combine " << parameters.m_combine << endmsg;
    debug() << "build " << parameters.m_build << endmsg;
    debug() << "extend " << parameters.m_extend << endmsg;
    debug() << "sortTracks " << parameters.m_sortTracks << endmsg;
    debug() << "tightenStep " << parameters.m_tightenStep << endmsg;

    runStep(kdClusters, nearestNeighbours, conformalTracks, collectionClusters, parameters);

    info() << "conformalTracks.size() " << conformalTracks.size() << endmsg;

    if (msgLevel(MSG::DEBUG)) {
      for (const auto& confTrack : conformalTracks) {
        debug() << "- Track " << &confTrack << " has " << confTrack->m_clusters.size() << " hits" << endmsg;
        for (unsigned int ht = 0; ht < confTrack->m_clusters.size(); ht++) {
          SKDCluster const& kdhit = confTrack->m_clusters.at(ht);
          debug() << "-- Hit " << ht << ": [x,y,z] = [" << kdhit->getX() << ", " << kdhit->getY() << ", "
                  << kdhit->getZ() << "]" << endmsg;
        }
      }
    }
  }

  // Clean up
  nearestNeighbours.reset(nullptr);

  // Now in principle have all conformal tracks, but due to how the check for clones is performed (ish) there is a
  // possibility that clones/fakes are still present. Try to remove them by looking at overlapping hits. Turned off at
  // the moment

  // Now make "real" tracks from all of the conformal tracks
  info() << "*** CA has made " << conformalTracks.size() << (conformalTracks.size() == 1 ? " track ***" : " tracks ***")
         << endmsg;

  // Loop over all track candidates
  for (const auto& conformalTrack : conformalTracks) {
    info() << "- Fitting track " << &conformalTrack << endmsg;

    // Make the LCIO track hit vector
    std::vector<const edm4hep::TrackerHitPlane*> trackHits;
    for (const auto& cluster : conformalTrack->m_clusters) {
      trackHits.push_back(&kdClusterMap[cluster]);
    }

    // Sort the hits from smaller to larger radius
    std::ranges::sort(trackHits,
                      (bool (*)(const edm4hep::TrackerHitPlane*, const edm4hep::TrackerHitPlane*))sort_by_radius);

    // Now we can make the track object and relations object, and fit the track
    edm4hep::MutableTrack track;

    // First, for some reason there are 2 track objects, one which gets saved and one which is used for fitting. Don't
    // ask...
    // TODO: Remove const_cast
    auto marlinTrk = GaudiDDKalTestTrack(this, const_cast<GaudiDDKalTest*>(&m_ddkaltest));

    // Make an initial covariance matrix with very broad default values
    // Track states in EDM4hep are stored in a 6x6 covariance matrix
    edm4hep::CovMatrix6f covMatrix{};
    covMatrix[0] = m_initialTrackError_d0;    // sigma_d0^2
    covMatrix[2] = m_initialTrackError_phi0;  // sigma_phi0^2
    covMatrix[5] = m_initialTrackError_omega; // sigma_omega^2
    covMatrix[9] = m_initialTrackError_z0;    // sigma_z0^2
    covMatrix[14] = m_initialTrackError_tanL; // sigma_tanl^2

    debug() << " Track hits before fit = " << trackHits.size() << endmsg;

    GaudiTrkUtils trkUtils(static_cast<const Gaudi::Algorithm*>(this), m_ddkaltest, m_geoSvc,
                           m_encodingStringVariable.value());

    // Try to fit
    int fitError = trkUtils.createFinalisedLCIOTrack(marlinTrk, trackHits, track, conformalTrack->m_kalmanFitForward,
                                                     covMatrix, m_magneticField, m_maxChi2perHit);

    // debug() << " Fit direction " << ((conformalTrack->m_kalmanFitForward) ? "forward" : "backward")
    //                       << endmsg;
    // debug() << " Track hits after fit = " << track->getTrackerHits().size() << endmsg;

    // // If the track is too short (usually less than 7 hits correspond to a vertex track)
    // // the fit is tried using the inverted direction
    if (int(track.getTrackerHits().size()) < m_maxHitsInvFit || fitError != 0) {
      auto marlinTrack_inv = GaudiDDKalTestTrack(this, const_cast<GaudiDDKalTest*>(&m_ddkaltest));
      edm4hep::MutableTrack track_inv;

      // Try to fit on the other way
      int fitError_inv =
          trkUtils.createFinalisedLCIOTrack(marlinTrack_inv, trackHits, track_inv, !conformalTrack->m_kalmanFitForward,
                                            covMatrix, m_magneticField, m_maxChi2perHit);

      debug() << " Fit direction " << ((!conformalTrack->m_kalmanFitForward) ? "forward" : "backward") << endmsg;
      debug() << " Track hits after inverse fit = " << track_inv.getTrackerHits().size() << endmsg;

      if (track_inv.getTrackerHits().size() > track.getTrackerHits().size()) {
        debug() << " Track is replaced. " << endmsg;
        // TODO
        //  track.swap(track_inv);
        //  marlinTrack.swap(marlinTrack_inv);
        fitError = fitError_inv;
      } else {
        debug() << " Track is not replaced. " << endmsg;
      }
    }

    // Check track quality - if fit fails chi2 will be 0. For the moment add hits by hand to any track that fails the
    // track fit, and store it as if it were ok...
    if (fitError != 0) {
      debug() << "- Fit fail error " << fitError << endmsg;
      continue;
    }

    // Check if track has minimum number of hits
    if (int(track.getTrackerHits().size()) < m_minClustersOnTrackAfterFit) {
      debug() << "- Track has " << track.getTrackerHits().size() << " hits. The minimum required is "
              << m_minClustersOnTrackAfterFit << endmsg;
      continue;
    }

    // Add hit information TODO: this is just a fudge for the moment, since we only use vertex hits. Should do for each
    // subdetector once enabled This is kept for compatibility with the original code in ConformalTracking.cc
    std::vector<int32_t> hitNumbers;
    // hitNumbers.resize(2 * lcio::ILDDetID::ETD);
    hitNumbers.resize(2 * 6);
    // hitNumbers[2 * lcio::ILDDetID::VXD - 2] = trackHits.size();
    hitNumbers[2 * 1 - 2] = trackHits.size();
    for (const auto num : hitNumbers) {
      track.addToSubdetectorHitNumbers(num);
    }

    // calculate purities and check if track has been reconstructed
    // if (m_debugPlots) {
    //   m_conformalChi2->Fill(conformalTrack->chi2ndof());
    //   debug() << "-------------------- New TRACK --------------------" << endmsg;
    //   // debug() << " LCIO track fit chi2 is "<<track->getChi2()<<std::endl;
    //   double purity = checkReal(conformalTrack, reconstructed, particleHits);
    //   if (purity >= m_purity) {
    //     m_conformalChi2real->Fill(conformalTrack->chi2ndof());
    //   }
    //   if (purity < m_purity) {
    //     m_conformalChi2fake->Fill(conformalTrack->chi2ndof());
    //   }
    //   m_conformalChi2Purity->Fill(purity, conformalTrack->chi2ndof());
    // }

    // Push back to the output container
    outputTrackCollection.push_back(track);
  }

  // Draw the cells for all produced tracks
  if (m_debugPlots && m_eventNumber == 0) {
    m_canvConformalEventDisplay->cd();
    for (auto& debugTrack : conformalTracks) {
      SharedKDClusters clusters = debugTrack->m_clusters;
      std::sort(clusters.begin(), clusters.end(), sort_by_lower_radiusKD);
      for (size_t itCluster = 1; itCluster < clusters.size(); itCluster++)
        drawline(clusters[itCluster - 1], clusters[itCluster], clusters.size() - itCluster);
    }
  }

  // Draw the conformal event display hits for debugging
  if (m_debugPlots && m_eventNumber == 0) {
    m_canvConformalEventDisplay->cd();
    m_conformalEvents->DrawCopy("same");
    m_canvConformalEventDisplayAllCells->cd();
    m_conformalEvents->DrawCopy("same");
    m_canvConformalEventDisplayAcceptedCells->cd();
    m_conformalEvents->DrawCopy("same");
    m_canvConformalEventDisplayMCunreconstructed->cd();
    m_conformalEvents->DrawCopy("same");
  }
  // if (m_debugPlots) {
  //   int nReconstructed(0), nUnreconstructed(0);
  //   // Additionally draw all tracks that were not reconstructed
  //   m_canvConformalEventDisplayMCunreconstructed->cd();
  //   int nParticles = particleCollection->getNumberOfElements();

  //   // Make the nearest neighbour tree to debug particle reconstruction issues
  //   SharedKDClusters kdClusters_debug;
  //   for (size_t i = 0; i < trackerHitCollections.size(); i++)
  //     kdClusters_debug.insert(kdClusters_debug.begin(), collectionClusters[i].begin(), collectionClusters[i].end());
  //   auto nearestNeighbours_debug = UKDTree(new KDTree(kdClusters_debug, m_thetaRange, m_sortTreeResults));

  //   for (int itP = 0; itP < nParticles; itP++) {
  //     // Get the particle
  //     MCParticle* mcParticle = dynamic_cast<MCParticle*>(particleCollection->getElementAt(itP));
  //     // Get the conformal hits
  //     if (particleHits.count(mcParticle) == 0)
  //       continue;
  //     SharedKDClusters mcHits = particleHits[mcParticle];
  //     // Cut on the number of hits
  //     int uniqueHits = getUniqueHits(mcHits);
  //     if (uniqueHits < m_minClustersOnTrack)
  //       continue;
  //     // Check if it was stable
  //     if (mcParticle->getGeneratorStatus() != 1)
  //       continue;
  //     // Check if it was reconstructed
  //     debug() << "-------------------- New PARTICLE --------------------" << endmsg;
  //     // List the pt
  //     double particlePt = sqrt(mcParticle->getMomentum()[0] * mcParticle->getMomentum()[0] +
  //                              mcParticle->getMomentum()[1] * mcParticle->getMomentum()[1]);
  //     debug() << "Particle pt: " << particlePt << endmsg;
  //     //checkReconstructionFailure(mcParticle, particleHits, nearestNeighbours_debug, _stepParameters[0]);
  //     if (reconstructed.count(mcParticle)) {
  //       nReconstructed++;
  //       continue;
  //     }
  //     // Draw the cells connecting the hits
  //     std::sort(mcHits.begin(), mcHits.end(), sort_by_radiusKD);
  //     for (size_t iHit = 0; iHit < (mcHits.size() - 1); iHit++) {
  //       drawline(mcHits[iHit], mcHits[iHit + 1], iHit + 1);
  //     }
  //      debug() << "Unreconstructed particle pt: " << particlePt << endmsg;
  //     nUnreconstructed++;

  //     // Check why particles were not reconstructed
  //     //      checkReconstructionFailure(mcParticle, particleHits, used, nearestNeighbours);
  //   }
  //   debug() << "Reconstructed " << nReconstructed << " particles out of " << nReconstructed + nUnreconstructed
  //                         << ". Gives efficiency "
  //                         << 100. * (double)nReconstructed / (double)(nReconstructed + nUnreconstructed) << "%" <<
  //                         endmsg;
  //   nearestNeighbours_debug.reset(nullptr);
  // }

  m_eventNumber++;
  return outputTrackCollection;
}

StatusCode ConformalTracking::finalize() {
  if (m_debugPlots) {
    auto file = TFile::Open("ConformalTracking.root", "RECREATE");
    // for (auto& h : _h) {
    //   h->Write();
    //   delete h;
    // }
    m_neighX->Write();
    delete m_neighX;
    m_neighY->Write();
    delete m_neighY;
    m_neighZ->Write();
    delete m_neighZ;

    m_slopeZ->Write();
    delete m_slopeZ;

    m_slopeZ_true->Write();
    delete m_slopeZ_true;
    m_slopeZ_true_first->Write();
    delete m_slopeZ_true_first;
    m_slopeZ_vs_pt_true->Write();
    delete m_slopeZ_vs_pt_true;

    m_cellAngle->Write();
    delete m_cellAngle;
    m_cellDOCA->Write();
    delete m_cellDOCA;
    m_cellAngleRadius->Write();
    delete m_cellAngleRadius;
    m_cellLengthRadius->Write();
    delete m_cellLengthRadius;
    m_cellAngleLength->Write();
    delete m_cellAngleLength;
    m_conformalChi2->Write();
    delete m_conformalChi2;
    m_conformalChi2real->Write();
    delete m_conformalChi2real;
    m_conformalChi2fake->Write();
    delete m_conformalChi2fake;
    m_conformalChi2Purity->Write();
    delete m_conformalChi2Purity;

    m_conformalChi2MC->Write();
    delete m_conformalChi2MC;
    m_conformalChi2PtMC->Write();
    delete m_conformalChi2PtMC;
    m_conformalChi2VertexRMC->Write();
    delete m_conformalChi2VertexRMC;

    m_conformalChi2SzMC->Write();
    delete m_conformalChi2SzMC;
    m_conformalChi2SzPtMC->Write();
    delete m_conformalChi2SzPtMC;
    m_conformalChi2SzVertexRMC->Write();
    delete m_conformalChi2SzVertexRMC;

    m_cellAngleMC->Write();
    delete m_cellAngleMC;
    m_cellDOCAMC->Write();
    delete m_cellDOCAMC;
    m_cellAngleRadiusMC->Write();
    delete m_cellAngleRadiusMC;
    m_cellLengthRadiusMC->Write();
    delete m_cellLengthRadiusMC;
    m_cellAngleLengthMC->Write();
    delete m_cellAngleLengthMC;

    m_cellAngleRZMC->Write();
    delete m_cellAngleRZMC;

    m_conformalEvents->Write();
    delete m_conformalEvents;
    m_nonconformalEvents->Write();
    delete m_nonconformalEvents;
    m_conformalEventsRTheta->Write();
    delete m_conformalEventsRTheta;
    m_conformalEventsMC->Write();
    delete m_conformalEventsMC;

    m_canvConformalEventDisplay->Write();
    delete m_canvConformalEventDisplay;
    m_canvConformalEventDisplayAllCells->Write();
    delete m_canvConformalEventDisplayAllCells;
    m_canvConformalEventDisplayAcceptedCells->Write();
    delete m_canvConformalEventDisplayAcceptedCells;
    m_canvConformalEventDisplayMC->Write();
    delete m_canvConformalEventDisplayMC;
    m_canvConformalEventDisplayMCunreconstructed->Write();
    delete m_canvConformalEventDisplayMCunreconstructed;

    file->Close();
  }
  return StatusCode::SUCCESS;
}

// Sort tracker hits from smaller to larger radius
inline bool sort_by_radius(const edm4hep::TrackerHitPlane& hit1, const edm4hep::TrackerHitPlane& hit2) {
  return edm4hep::utils::magnitudeTransverse(hit1.getPosition()) <
         edm4hep::utils::magnitudeTransverse(hit2.getPosition());
}

inline bool sort_by_radius(const edm4hep::TrackerHitPlane* hit1, const edm4hep::TrackerHitPlane* hit2) {
  return edm4hep::utils::magnitudeTransverse(hit1->getPosition()) <
         edm4hep::utils::magnitudeTransverse(hit2->getPosition());
}

// Sort kd hits from larger to smaller radius
inline bool sort_by_radiusKD(const SKDCluster& hit1, const SKDCluster& hit2) { return hit1->getR() > hit2->getR(); }

// Sort kdhits by lower to higher layer number
inline bool sort_by_layer(const SKDCluster& hit1, const SKDCluster& hit2) {
  if (hit1->getSubdetector() != hit2->getSubdetector())
    return (hit1->getSubdetector() < hit2->getSubdetector());
  else if (hit1->getLayer() != hit2->getLayer())
    return (hit1->getLayer() < hit2->getLayer());
  else if (hit1->getSide() != hit2->getSide())
    return (hit1->getSide() < hit2->getSide());
  else
    return false;
}

// Sort kdhits by higher to lower layer number
inline bool sort_by_lower_layer(const SKDCluster& hit1, const SKDCluster& hit2) {
  if (hit1->getSubdetector() != hit2->getSubdetector())
    return (hit1->getSubdetector() > hit2->getSubdetector());
  else if (hit1->getLayer() != hit2->getLayer())
    return (hit1->getLayer() > hit2->getLayer());
  else if (hit1->getSide() != hit2->getSide())
    return (hit1->getSide() > hit2->getSide());
  else
    return false;
}

bool ConformalTracking::neighbourIsCompatible(const SKDCluster& neighbourHit, const SKDCluster& seedHit,
                                              const double slopeZRange) const {
  const double distanceX(neighbourHit->getX() - seedHit->getX());
  const double distanceY(neighbourHit->getY() - seedHit->getY());
  const double distance(sqrt(distanceX * distanceX + distanceY * distanceY));
  const double deltaZ(neighbourHit->getZ() - seedHit->getZ());
  if (fabs(deltaZ / distance) > slopeZRange) {
    debug() << "- z condition not met" << endmsg;
    if (debugSeed && seedHit == debugSeed)
      debug() << "- z condition not met" << endmsg;
    return false;
  }

  return true;
}

double ConformalTracking::getBzAtOrigin() const {
  double bfield(0.0);

  const dd4hep::Detector* theDetector = m_geoSvc->getDetector();
  const double position[3] = {0, 0, 0};      // position to calculate magnetic field at (the origin in this case)
  double magneticFieldVector[3] = {0, 0, 0}; // initialise object to hold magnetic field
  theDetector->field().magneticField(position, magneticFieldVector); // get the magnetic field vector from DD4hep
  bfield = magneticFieldVector[2] / dd4hep::tesla;                   // z component at (0,0,0)
  return bfield;
}

// Combine collections
UKDTree ConformalTracking::combineCollections(SharedKDClusters& kdClusters, const std::vector<int>& combination,
                                              const std::map<size_t, SharedKDClusters>& collectionClusters) const {
  // Clear the input objects
  kdClusters.clear();

  // Loop over all given collections
  for (const auto& comb : combination) {
    // Copy the clusters to the output vector
    const SharedKDClusters& clusters = collectionClusters.at(comb);
    kdClusters.insert(kdClusters.end(), clusters.begin(), clusters.end());
  }

  debug() << "*** combineCollections: Collection has " << kdClusters.size() << " hits" << endmsg;

  // Sort the KDClusters from larger to smaller radius
  std::ranges::sort(kdClusters, sort_by_radiusKD);

  // Make the binary search tree. This tree class contains two binary trees - one sorted by u-v and the other by theta
  return std::make_unique<KDTree>(kdClusters, m_thetaRange, m_sortTreeResults);
}

// Given a list of connected cells (so-called cellular tracks), return the candidate(s) with lowest chi2/degrees of
// freedom. Attempt to remove hits to see if there is a significant improvement in the chi2/ndof, to protect against
// noise hits causing a good track to be discarded. If several candidates have low chi2/ndof and are not clones (limited
// sharing of hits) then return all of them. Given that there is no material scattering taken into account this helps
// retain low pt tracks, which may have worse chi2/ndof than ghosts/real tracks with an additional unrelated hit from
// the low pt track.
UniqueKDTracks ConformalTracking::getFittedTracks(UniqueCellularTracks& candidateTracks,
                                                  const Parameters& parameters) const {
  debug() << "***** getFittedTracks" << endmsg;

  UniqueKDTracks finalTracks;

  // Make a container for all tracks being considered, initialise variables
  UniqueKDTracks trackContainer;
  //  std::vector<double> trackChi2ndofs;

  // Loop over all candidate tracks and do an inital fit to get the track angle (needed to calculate the
  // hit errors for the error-weighted fit)
  for (auto& candidateTrack : candidateTracks) {
    // If there are not enough hits on the track, ignore it
    if (int(candidateTrack->size()) < (parameters.m_minClustersOnTrack - 2)) {
      debug() << "- Track " << candidateTrack.get() << ": not enough hits on track." << endmsg;
      candidateTrack.reset();
      continue;
    }

    // Make the fitting object. TGraphErrors used for 2D error-weighted fitting
    UKDTrack track = std::unique_ptr<KDTrack>(new KDTrack(parameters));

    // Loop over all hits and add them to the fitter (and track)
    int npoints = 0;
    track->add((*candidateTrack)[0]->getEnd());
    npoints++;

    for (const auto& trackCell : (*candidateTrack)) {
      track->add(trackCell->getStart());
      npoints++;
    }
    track->linearRegression(parameters.m_highPTfit);
    track->linearRegressionConformal();

    double chi2ndofTOT = track->chi2ndof() + track->chi2ndofZS();

    // We try to see if there are spurious hits causing the chi2 to be very large. This would cause us to throw away
    // good tracks with perhaps just a single bad hit. Try to remove each hit and see if fitting without it causes a
    // significant improvement in chi2/ndof

    // Loop over each hit (starting at the back, since we will use the 'erase' function to get rid of them)
    // and see if removing it improves the chi2/ndof
    int removed = 0;
    if (chi2ndofTOT > parameters.m_chi2cut &&
        chi2ndofTOT < parameters.m_chi2cut) { // CHANGE ME?? Upper limit to get rid of really terrible tracks (temp
                                              // lower changed from 0 to parameters.m_chi2cut)

      debug() << "- Track " << candidateTrack.get() << ": chi2ndofTOT > chi2cut. Try to fit without point" << endmsg;
      for (int point = npoints - 1; point >= 0; point--) {
        // Stop if we would remove too many points on the track to meet the minimum hit requirement (or if the track has
        // more than 2 hits removed)
        if ((npoints - removed - 1) < parameters.m_minClustersOnTrack || removed == 2)
          break;

        // Refit the track without this point
        double newChi2ndofTOT = fitWithoutPoint(*track, point);

        // If the chi2/ndof is significantly better, remove the point permanently CHANGE ME??
        //        if( (chi2ndofTOT - newChi2ndofTOT) > 0 && (chi2ndofTOT - newChi2ndofTOT) > 1. ){
        if ((newChi2ndofTOT - chi2ndofTOT) < chi2ndofTOT) {
          track->remove(point);
          removed++;
          chi2ndofTOT = newChi2ndofTOT;
        }
      }
    }

    debug() << "- Track " << candidateTrack.get() << " has " << track->m_clusters.size() << " hits after fit" << endmsg;
    // for (unsigned int cluster = 0; cluster < track->m_clusters.size() && streamlog_level(DEBUG8); cluster++) {
    //   debug() << "-- Hit " << cluster << ": [x,y] = [" << track->m_clusters.at(cluster)->getX() << ", "
    //                         << track->m_clusters.at(cluster)->getY() << "]" << endmsg;
    // }

    // Store the track information
    trackContainer.push_back(std::move(track));

    candidateTrack.reset();

  } // end for candidateTracks

  // Now have all sets of conformal tracks and their chi2/ndof. Decide which tracks to send back, ie. the one with
  // lowest chi2/ndof, and possibly others if they are not clones and have similar chi2 value
  getLowestChi2(finalTracks, trackContainer);
  debug() << "getFittedTracks *****" << endmsg;

  // Send back the final set of tracks
  return finalTracks;
}

// Pick the lowest chi2/ndof KDTrack from a list of possible tracks, and additionally return other tracks in the
// collection with similar chi2/ndof values that don't share many hits
void ConformalTracking::getLowestChi2(UniqueKDTracks& finalTracks, UniqueKDTracks& trackContainer) const {
  debug() << "***** getLowestChi2" << endmsg;

  // Get the lowest chi2/ndof value from the given tracks
  //  double lowestChi2ndof = *std::min_element(trackChi2ndofs.begin(),trackChi2ndofs.end());
  UKDTrack& lowestChi2ndofTrack = trackContainer[0];
  double lowestChi2ndof = sqrt(lowestChi2ndofTrack->chi2ndof() * lowestChi2ndofTrack->chi2ndof() +
                               lowestChi2ndofTrack->chi2ndofZS() * lowestChi2ndofTrack->chi2ndofZS());

  for (const auto& itTrack : trackContainer) {
    double chi2ndof = sqrt(itTrack->chi2ndof() * itTrack->chi2ndof() + itTrack->chi2ndofZS() * itTrack->chi2ndofZS());

    debug() << "- Track " << itTrack.get() << " has " << chi2ndof << " chi2" << endmsg;
    if (chi2ndof < lowestChi2ndof) {
      lowestChi2ndof = itTrack->chi2ndof();
      // lowestChi2ndofTrack = trackContainer[itTrack];
    }
  }

  // finalTracks.push_back(lowestChi2ndofTrack);
  // return;
  debug() << "Save more than one track if similar chi2" << endmsg;

  // Loop over all other tracks and decide whether or not to save them
  // Look at the difference in chi2/ndof - we want to keep tracks with similar chi2/ndof. If they
  // are clones then take the longest
  copy_if(std::make_move_iterator(trackContainer.begin()), std::make_move_iterator(trackContainer.end()),
          std::back_inserter(finalTracks), [lowestChi2ndof](UKDTrack const& track) {
            return ((sqrt((track->chi2ndof() * track->chi2ndof() + track->chi2ndofZS() * track->chi2ndofZS())) -
                     lowestChi2ndof) < 10.);
          });

  debug() << "Final fitted tracks: " << finalTracks.size() << endmsg;
  // if (streamlog_level(DEBUG8)) {
  //   for (auto& itTrk : finalTracks) {
  //     debug() << "- Track " << itTrk.get() << " has " << itTrk->m_clusters.size() << " hits" << endmsg;
  //     for (unsigned int cluster = 0; cluster < itTrk->m_clusters.size(); cluster++) {
  //       debug() << "-- Hit " << cluster << ": [x,y] = [" << itTrk->m_clusters.at(cluster)->getX() << ", "
  //                             << itTrk->m_clusters.at(cluster)->getY() << "]" << endmsg;
  //     }
  //   }
  // }
  trackContainer.clear();
  debug() << "getLowestChi2 *****" << endmsg;

  return;
}

double ConformalTracking::fitWithoutPoint(KDTrack track, int point) const {
  // Remove the given point from the track
  track.remove(point);

  // Calculate the track chi2 with the final fitted values
  track.linearRegression();
  track.linearRegressionConformal();

  return track.chi2ndof() + track.chi2ndofZS();
}

void ConformalTracking::updateCell(const SCell& cell) const {
  if (cell->getTo().size() == 0)
    return;
  for (size_t i = 0; i < cell->getTo().size(); i++) {
    SCell(cell->getTo()[i])->update(cell);
    updateCell(SCell(cell->getTo()[i]));
  }
}

SKDCluster ConformalTracking::extrapolateCell(const SCell& cell, double distance) const {
  // Fake cluster to be returned
  SKDCluster extrapolatedCluster = std::make_shared<KDCluster>();

  // Linear extrapolation of the cell - TODO: check that cell gradients have correct sign and remove checks here
  double gradient = cell->getGradient();
  double deltaU = sqrt(distance * distance / (1 + gradient * gradient));
  double deltaV = std::abs(gradient) * deltaU;

  if ((cell->getStart()->getU() - cell->getEnd()->getU()) > 0)
    deltaU *= (-1.);
  if ((cell->getStart()->getV() - cell->getEnd()->getV()) > 0)
    deltaV *= (-1.);

  extrapolatedCluster->setU(cell->getEnd()->getU() + deltaU);
  extrapolatedCluster->setV(cell->getEnd()->getV() + deltaV);

  return extrapolatedCluster;
}

// Function to check if two KDtracks contain several hits in common
int ConformalTracking::overlappingHits(const UKDTrack& track1, const UKDTrack& track2) const {
  int nHitsInCommon = 0;
  for (size_t hit = 0; hit < track1->m_clusters.size(); hit++) {
    if (std::find(track2->m_clusters.begin(), track2->m_clusters.end(), track1->m_clusters[hit]) !=
        track2->m_clusters.end())
      nHitsInCommon++;
  }
  return nHitsInCommon;
}

void ConformalTracking::extendSeedCells(SharedCells& cells, UKDTree& nearestNeighbours, bool extendingTrack,
                                        const SharedKDClusters& /*debugHits*/, Parameters const& parameters,
                                        bool vertexToTracker) const {
  debug() << "***** extendSeedCells" << endmsg;

  size_t nCells = 0;
  int depth = 0;
  size_t startPos = 0;

  // Keep track of existing cells in case there are branches in the track
  std::map<SKDCluster, SharedCells> existingCells;

  debug() << "Extending " << cells.size()
          << " cells. Start with seed cells (depth 0) and then move forward (depth > 0)." << endmsg;
  // Try to create all "downstream" cells until no more can be added
  while (cells.size() != nCells) {
    // Extend all cells with depth N. In the next iteration, look at cells with depth N+1
    nCells = cells.size();
    debug() << "Depth = " << depth << endmsg;

    for (size_t itCell = startPos; itCell < nCells; itCell++) {
      debug() << "- Extend cell " << itCell << ": A ([x,y] = [" << cells[itCell]->getStart()->getX() << ", "
              << cells[itCell]->getStart()->getY() << "]) - B ([x,y] = [" << cells[itCell]->getEnd()->getX() << ", "
              << cells[itCell]->getEnd()->getY() << "])" << endmsg;
      // Get the end point of the cell (to search for neighbouring hits to form new cells connected to this one)
      SKDCluster const& hit = cells[itCell]->getEnd();
      double searchDistance = parameters.m_maxDistance; // hit->getR();
      if (searchDistance > hit->getR())
        searchDistance = 1.2 * hit->getR();

      // Extrapolate along the cell and then make a 2D nearest neighbour search at this extrapolated point
      SKDCluster const& fakeHit =
          extrapolateCell(cells[itCell], searchDistance / 2.); // TODO: make this search a function of radius
      SharedKDClusters results;
      nearestNeighbours->allNeighboursInRadius(
          fakeHit, 0.625 * searchDistance, results, [&hit, vertexToTracker](SKDCluster const& nhit) {
            if (nhit->used())
              return true;
            if (hit->sameLayer(nhit))
              return true;
            if (nhit->endcap() && hit->endcap() && (nhit->forward() != hit->forward()))
              return true;
            if ((vertexToTracker && nhit->getR() >= hit->getR()) || (!vertexToTracker && nhit->getR() <= hit->getR()))
              return true;
            return false;
          });

      debug() << "- Found " << results.size() << " neighbours from cell extrapolation " << endmsg;
      if (extendingTrack) {
        debug() << "Extrapolating cell " << itCell << " from u,v: " << hit->getU() << "," << hit->getV() << endmsg;
        debug() << "Extrapolated hit at u,v: " << fakeHit->getU() << "," << fakeHit->getV() << endmsg;
        debug() << "Found " << results.size() << " neighbours from cell extrapolation" << endmsg;
      }

      // Make new cells pointing inwards
      for (unsigned int neighbour = 0; neighbour < results.size(); neighbour++) {
        // Get the neighbouring hit
        SKDCluster const& nhit = results[neighbour];

        debug() << "-- Neighbour " << neighbour << ": [x,y] = [" << nhit->getX() << ", " << nhit->getY() << "]"
                << endmsg;
        if (extendingTrack)
          debug() << "looking at neighbour " << neighbour << " at u,v: " << nhit->getU() << "," << nhit->getV()
                  << endmsg;

        // Check that it is not used, is not on the same detector layer, points inwards and has real z pointing away
        // from IP
        //        if(used.count(nhit)){if(extendingTrack)debug()<<"- used"<<std::endl; continue;}

        if (nhit->used()) {
          debug() << "-- used" << endmsg;
          continue;
        }
        // if (hit->sameLayer(nhit)) {
        //   if (extendingTrack)
        //     debug() << "- same layer" << endmsg;
        //   continue;
        // }
        // if (nhit->endcap() && hit->endcap() && (nhit->forward() != hit->forward()))
        //   continue;
        // if ((vertexToTracker && nhit->getR() >= hit->getR()) || (!vertexToTracker && nhit->getR() <= hit->getR())) {
        //   if (extendingTrack)
        //     debug() << "- " << (vertexToTracker ? "higher radius" : "lower radius") << endmsg;
        //   continue;
        // }

        // Check if this cell already exists (rejoining branch) FIXME - allows rejoining a branch without checking cell
        // angles
        auto const& existingCellsForHit = existingCells.find(hit);
        if (existingCellsForHit != existingCells.end()) {
          bool alreadyExists = false;
          for (const auto& existingCell : existingCellsForHit->second) {
            if (existingCell->getEnd() == nhit) {
              debug() << "-- cell A ([x,y] = [" << hit->getX() << ", " << hit->getY() << "]) - B ([x,y] = ["
                      << nhit->getX() << ", " << nhit->getY() << "]) already exists" << endmsg;
              alreadyExists = true;

              // Check if cell angle is too large to rejoin
              if (cells[itCell]->getAngle(existingCell) > parameters.m_maxCellAngle ||

                  cells[itCell]->getAngleRZ(existingCell) > parameters.m_maxCellAngleRZ) {
                debug() << "-- cell A ([x,y] = [" << hit->getX() << ", " << hit->getY() << "]) - B ([x,y] = ["
                        << nhit->getX() << ", " << nhit->getY() << "]) angle too large" << endmsg;
                continue;
              }
              // Otherwise add the path
              cells[itCell]->setTo(existingCell);
              existingCell->setFrom(cells[itCell]);
              updateCell(existingCell);
            }
          }
          if (alreadyExists)
            continue;
        }

        // Make the new cell
        Cell cell(hit, nhit);

        debug() << "-- made new cell A ([x,y] = [" << hit->getX() << ", " << hit->getY() << "]) - B ([x,y] = ["
                << nhit->getX() << ", " << nhit->getY() << "]) " << endmsg;
        if (extendingTrack)
          debug() << "- make new cell" << endmsg;

        // Check if the new cell is compatible with the previous cell (angle between the two is acceptable)
        //        if( cells[itCell]->getAngle(cell) > (parameters.m_maxCellAngle*exp(-0.001/nhit->getR())) ){
        if (cells[itCell]->getAngle(cell) > parameters.m_maxCellAngle ||
            cells[itCell]->getAngleRZ(cell) > parameters.m_maxCellAngleRZ) {
          // Debug plotting
          //          if(m_debugPlots && m_eventNumber == 0){
          //            m_canvConformalEventDisplayAllCells->cd();
          //            drawline(hit,nhit,cells[itCell]->getWeight()+2,3);
          //          }

          debug() << "-- cell angle too large. Discarded! " << endmsg;

          if (extendingTrack)
            debug() << "-- discarded!" << endmsg;

          continue;
        }

        // Set the information about which cell this new cell is attached to and store it
        cells.emplace_back(std::make_shared<Cell>(std::move(cell)));
        auto const& scell = cells.back();
        existingCells[hit].push_back(scell);
        scell->setFrom(cells[itCell]);
        cells[itCell]->setTo(scell);

        // Debug plotting
        //        if(m_debugPlots && m_eventNumber == 0){
        //          m_canvConformalEventDisplayAllCells->cd();
        //          drawline(hit,nhit,cells[itCell]->getWeight()+2);
        //        }
      }

      // Finished adding new cells to this cell
    }

    // All new cells added at this depth
    startPos = nCells;
    depth++;
  }
  debug() << "extendSeedCells *****" << endmsg;

  // No more downstream cells can be added
}

// Add a point to a track and return the delta chi2
void ConformalTracking::fitWithPoint(KDTrack kdTrack, SKDCluster& hit, double& deltaChi2, double& deltaChi2zs) const {
  double chi2 = kdTrack.chi2ndof();
  double chi2zs = kdTrack.chi2ndofZS();
  /*if(m_debugPlots){

    double xMeasured = hit->getU();
    double yMeasured = hit->getV();
    double dx = hit->getErrorU();
    double dv = hit->getErrorV();

    if(kdTrack.m_rotated){
      double newxMeasured = yMeasured;
      double newyMeasured = -1. * xMeasured;
      double newdx        = dv;
      double newdv        = dx;
      xMeasured           = newxMeasured;
      yMeasured           = newyMeasured;
      dx                  = newdx;
      dv                  = newdv;
    }
    double residualY = (kdTrack.m_gradient * xMeasured + kdTrack.m_quadratic * xMeasured * xMeasured +
  kdTrack.m_intercept) - yMeasured;
    // Get the error on the hit position
    double term = kdTrack.m_gradient + 2. * kdTrack.m_quadratic * xMeasured;
    double dy2  = (dv * dv) + (term * term * dx * dx);

    streamlog_out(DEBUG7)<<"- hit has delta chi2 of "<<(residualY * residualY) / (dy2)<<std::endl;

  }*/
  kdTrack.add(hit);
  kdTrack.linearRegression();
  kdTrack.linearRegressionConformal();
  double newchi2 = kdTrack.chi2ndof();
  double newchi2zs = kdTrack.chi2ndofZS();

  deltaChi2 = newchi2 - chi2;
  deltaChi2zs = newchi2zs - chi2zs;
}

// Take a collection of tracks and try to extend them into the collection of clusters passed.
void ConformalTracking::extendTracks(UniqueKDTracks& conformalTracks, SharedKDClusters& collection,
                                     UKDTree& nearestNeighbours, Parameters const& parameters) const {
  // Loop over all current tracks. At the moment this is a "stupid" algorithm: it will simply try to add every
  // hit in the collection to every track, and keep the ones thta have a good chi2. In fact, it will extrapolate
  // the track and do a nearest neighbours search, but this seemed to fail for some reason, TODO!

  debug() << "EXTENDING " << conformalTracks.size() << " tracks, into " << collection.size() << " hits" << endmsg;
  if (collection.empty())
    return;

  // Sort the hit collection by layer
  std::ranges::sort(collection, sort_by_layer);

  // First extend high pt tracks
  extendTracksPerLayer(conformalTracks, collection, nearestNeighbours, parameters, parameters.m_vertexToTracker);

  // Mark hits from "good" tracks as being used
  for (const auto& conformalTrack : conformalTracks)
    for (const auto& thisCluster : conformalTrack->m_clusters)
      thisCluster->used(true);

  // index just for debug
  [[maybe_unused]] int debug_idxTrack = 0;

  for (const auto& track : conformalTracks) {
    // Make sure that track hits are ordered from largest to smallest radius
    std::ranges::sort(track->m_clusters, sort_by_radiusKD);

    // Get the associated MC particle
    // edm4hep::MCParticle* associatedParticle = nullptr;
    // if (m_debugPlots) {
    //   associatedParticle = m_debugger.getAssociatedParticle(track);
    //   if (associatedParticle == nullptr)
    //     debug() << "- nullptr particle!" << endmsg;

    //   // Check the track pt estimate
    //   TLorentzVector mc_helper;
    //   mc_helper.SetPxPyPzE(associatedParticle->getMomentum()[0], associatedParticle->getMomentum()[1],
    //                        associatedParticle->getMomentum()[2], associatedParticle->getEnergy());

    //   debug() << "- extending track " << debug_idxTrack << " with pt = " << mc_helper.Pt()
    //                         << ". pt estimate: " << track->pt() << " chi2/ndof " << track->chi2ndof() << " and
    //                         chi2/ndof ZS "
    //                         << track->chi2ndofZS() << endmsg;
    // }
    debug_idxTrack++;

    // Create a seed cell (connecting the first two hits in the track vector - those at smallest conformal radius)
    size_t nclusters = track->m_clusters.size();
    Cell seedCell{track->m_clusters[nclusters - 2], track->m_clusters[nclusters - 1]};

    // Loop over all hits found and check if any have a sufficiently good delta Chi2
    SharedKDClusters goodHits;
    SKDCluster bestCluster = nullptr;
    double bestChi2 = 0.;

    for (size_t nKDHit = 0; nKDHit < collection.size(); nKDHit++) {
      // Get the kdHit and check if it has already been used (assigned to a track)
      SKDCluster kdhit = collection[nKDHit];
      bool associated = false;
      // if (m_debugPlots) {
      //   associated = m_debugger.isAssociated(kdhit, associatedParticle);
      //   if (associated)
      //     debug() << "-- hit " << nKDHit << " belongs to this track" << endmsg;
      // }
      // debug()<<"Detector "<<kdhit->getSubdetector()<<", layer "<<kdhit->getLayer()<<", side
      // "<<kdhit->getSide()<<std::endl;

      // If this hit is on a new sensor, then add the hit from the previous sensor and start anew
      if (bestCluster && !(kdhit->sameSensor(bestCluster))) {
        bestCluster->used(true);
        track->add(bestCluster);
        track->linearRegression();
        track->linearRegressionConformal();
        bestCluster = nullptr;
      }

      // Don't reuse hits
      if (kdhit->used()) {
        if (associated)
          debug() << "used" << endmsg;
        continue;
      }

      // Don't pick up hits in the opposite side of the detector
      if ((track->m_clusters[nclusters - 1]->getZ() > 0. && kdhit->getZ() < 0.) ||
          (track->m_clusters[nclusters - 1]->getZ() < 0. && kdhit->getZ() > 0.)) {
        if (associated)
          debug() << "opposite side of detector" << endmsg;
        continue;
      }

      // First check that the hit is not wildly away from the track (make cell and check angle)
      // SCell extensionCell = std::make_shared<Cell>(track->m_clusters[0],results2[newHit]);
      Cell extensionCell(track->m_clusters[nclusters - 1], kdhit);
      double cellAngle = seedCell.getAngle(extensionCell);
      double cellAngleRZ = seedCell.getAngleRZ(extensionCell);
      if (cellAngle > 3. * parameters.m_maxCellAngle || cellAngleRZ > 3. * parameters.m_maxCellAngleRZ) {
        if (associated)
          debug() << "-- killed by cell angle cut" << endmsg;
        continue;
      }

      // Now fit the track with the new hit and check the increase in chi2
      double deltaChi2(0.), deltaChi2zs(0.);

      fitWithPoint(*(track), kdhit, deltaChi2, deltaChi2zs); // track->deltaChi2(results2[newHit]);

      if (associated) {
        debug() << "-- hit was fitted and has a delta chi2 of " << deltaChi2 << " and delta chi2zs of " << deltaChi2zs
                << endmsg;
      }

      // We have an estimate of the pT here, could use it in the chi2 criteria
      double chi2cut = parameters.m_chi2cut;
      // if (track->pt() < 5.)
      //   chi2cut = 1000.;

      if (deltaChi2 > chi2cut || deltaChi2zs > chi2cut)
        continue;

      bool onSameSensor = false;
      for (const auto& clusterOnTrack : track->m_clusters) {
        if (kdhit->sameSensor(clusterOnTrack)) {
          onSameSensor = true;
          break;
        }
      }
      if (onSameSensor)
        continue;

      if (associated)
        debug() << "-- valid candidate!" << endmsg;

      if (!bestCluster || deltaChi2 < bestChi2) {
        bestCluster = kdhit;
        bestChi2 = deltaChi2;
      }
    }

    if (bestCluster) {
      bestCluster->used(true);
      track->add(bestCluster);
      track->linearRegression();
      track->linearRegressionConformal();
      bestCluster = nullptr;
    }

  } // End of loop over tracks
}

void ConformalTracking::extendTracksPerLayer(UniqueKDTracks& conformalTracks, SharedKDClusters& collection,
                                             UKDTree& nearestNeighbours, Parameters const& parameters,
                                             bool vertexToTracker) const {
  debug() << "*** extendTracksPerLayer: extending " << conformalTracks.size() << " tracks layer by layer, into "
          << collection.size() << " hits" << endmsg;

  if (collection.size() == 0)
    return;

  debug() << "Total number of tracks = " << conformalTracks.size() << endmsg;

  // index just for debug
  int debug_idxTrack = 0;

  // Sort the input collection by layer
  std::ranges::sort(collection, vertexToTracker ? sort_by_layer : sort_by_lower_layer);

  // Loop over all tracks
  for (const auto& track : conformalTracks) {
    debug() << "Track " << debug_idxTrack << endmsg;
    debug_idxTrack++;

    // Only look at high pt tracks
    if (track->pt() < parameters.m_highPTcut)
      continue;

    // Make sure that the hits are ordered in KDradius
    std::ranges::sort(track->m_clusters, vertexToTracker ? sort_by_radiusKD : sort_by_lower_radiusKD);

    // Make the cell from the two hits of the track from which to extend
    size_t nclusters = track->m_clusters.size();
    debug() << "Track has " << nclusters << " hits" << endmsg;
    for (size_t i = 0; i < nclusters; i++) {
      debug() << "- Hit " << i << ": [x,y] = [" << track->m_clusters.at(i)->getX() << ", "
              << track->m_clusters.at(i)->getY() << "]" << endmsg;
    }
    Cell seedCell = {track->m_clusters[nclusters - 2], track->m_clusters[nclusters - 1]};
    debug() << "Seed cell A ([x,y] = [" << track->m_clusters[nclusters - 2]->getX() << ","
            << track->m_clusters[nclusters - 2]->getY() << "]) - B ([x,y] = ["
            << track->m_clusters[nclusters - 1]->getX() << "," << track->m_clusters[nclusters - 1]->getY() << "])"
            << endmsg;

    // Initialize variables for nearest neighbours search and subdetector layers
    SharedKDClusters results;
    bool loop = true;
    int extendInSubdet = 0;
    int extendInLayer = 0;
    int final_subdet = 0;
    int final_layer = 0;
    bool firstLayer = true;
    bool skipLayer = false;
    SKDCluster expectedHit = nullptr;

    // Start the extension into input collection, layer by layer
    do {
      debug() << "- Start extension" << endmsg;
      // Find nearest neighbours in theta and sort them by layer
      SKDCluster const& kdhit = seedCell.getEnd();
      debug() << "- Endpoint seed cell: [x,y] = [" << kdhit->getX() << ", " << kdhit->getY() << "]; "
              << "]; r = " << kdhit->getR() << "; [subdet,layer] = [" << kdhit->getSubdetector() << ", "
              << kdhit->getLayer() << "]" << endmsg;
      double theta = kdhit->getTheta();
      nearestNeighbours->allNeighboursInTheta(
          theta, m_thetaRange * 4, results, [&kdhit, vertexToTracker](SKDCluster const& nhit) {
            // if same subdet, take only the hits within two layers
            if (nhit->getSubdetector() == kdhit->getSubdetector()) {
              if ((vertexToTracker && nhit->getLayer() > (kdhit->getLayer() + 2)) ||
                  (!vertexToTracker && nhit->getLayer() < (kdhit->getLayer() - 2))) {
                return true;
              }
            }

            if ((vertexToTracker && nhit->getR() > kdhit->getR()) || (!vertexToTracker && nhit->getR() < kdhit->getR()))
              return true;
            return false;
          });
      // nearestNeighbours->allNeighboursInRadius(kdhit, parameters.m_maxDistance, results);
      debug() << "- Found " << results.size() << " neighbours. " << endmsg;
      // if (streamlog_level(DEBUG9)) {
      //   for (auto const& neighbour : results) {
      //     debug() << "-- Neighbour from allNeighboursInTheta : [x,y] = [" << neighbour->getX() << ", "
      //                           << neighbour->getY() << "]; r = " << neighbour->getR() << "; [subdet,layer] = ["
      //                           << neighbour->getSubdetector() << ", " << neighbour->getLayer() << "]" << endmsg;
      //   }
      // }

      if (results.size() == 0) {
        loop = false;
        continue;
      }

      std::ranges::sort(results, vertexToTracker ? sort_by_layer : sort_by_lower_layer);
      // Get the final values of subdet and layer to stop the loop at the (vertexToTracker? outermost : innermost) layer
      // with neighbours
      final_subdet = results.back()->getSubdetector();
      final_layer = results.back()->getLayer();
      debug() << "- Final layer with neighbours: [subdet,layer] = [" << final_subdet << ", " << final_layer << "]"
              << endmsg;

      // If no hit was found on a layer, we use expectedHit as current position
      // However, the nearest neighbour search is always performed from the last real hit (kdhit)
      const SKDCluster& hitOnCurrentLayer = skipLayer ? expectedHit : kdhit;

      // Set the subdetector layer in which to extend -- for this it is important that the hits are sorted by layer
      // before First layer is based on first neighbour Next layers are based on next wrt current
      for (const auto& neighbour : results) {
        debug() << "-- Neighbour : [x,y] = [" << neighbour->getX() << ", " << neighbour->getY()
                << "]; [subdet,layer] = [" << neighbour->getSubdetector() << ", " << neighbour->getLayer() << "]"
                << endmsg;

        if (firstLayer) { // first neighbour
          extendInSubdet = results.front()->getSubdetector();
          extendInLayer = results.front()->getLayer();
          debug() << "-- firstLayer) [subdet,layer] = [" << extendInSubdet << ", " << extendInLayer << "]" << endmsg;
          firstLayer = false;
          break;
        } else if (neighbour->getSubdetector() ==
                   hitOnCurrentLayer->getSubdetector()) { // next layer is in same subdetector
          if ((vertexToTracker && (neighbour->getLayer() > hitOnCurrentLayer->getLayer())) ||
              (!vertexToTracker && (neighbour->getLayer() < hitOnCurrentLayer->getLayer()))) {
            extendInSubdet = neighbour->getSubdetector();
            extendInLayer = neighbour->getLayer();
            debug() << "-- sameSub,diffLayer) [subdet,layer] = [" << extendInSubdet << ", " << extendInLayer << "]"
                    << endmsg;
            break;
          }
        } else if ((vertexToTracker && (neighbour->getSubdetector() > hitOnCurrentLayer->getSubdetector())) ||
                   (!vertexToTracker && (neighbour->getSubdetector() <
                                         hitOnCurrentLayer->getSubdetector()))) { // next layer is in next subdetector
          extendInSubdet = neighbour->getSubdetector();
          extendInLayer = neighbour->getLayer();
          debug() << "-- diffSub) [subdet,layer] = [" << extendInSubdet << ", " << extendInLayer << "]" << endmsg;
          break;
        } else
          continue;
      }

      // Set the condition to end the loop
      if (extendInSubdet == final_subdet && extendInLayer == final_layer) {
        debug() << "- Ending the loop. Reached final subdet layer" << endmsg;
        loop = false;
      }

      // Initialize variables for choosing the best neighbour in layer
      SKDCluster bestCluster = nullptr;
      std::map<SKDCluster, double> bestClustersWithChi2 = {};
      std::vector<SKDCluster> bestClusters = {};
      double bestChi2 = std::numeric_limits<double>::max();
      double chi2 = track->chi2ndof();
      double chi2zs = track->chi2ndofZS();
      KDTrack tempTrack = *track;

      // Loop over neighbours
      for (const auto& neighbour : results) {
        // only in the layer in which to extend
        if (neighbour->getSubdetector() == extendInSubdet && neighbour->getLayer() == extendInLayer) {
          // Store the hit in case no bestCluster will be found in this layer
          expectedHit = neighbour;

          debug() << "-- Looking at neighbour " << neighbour << ": [x,y] = [" << neighbour->getX() << ", "
                  << neighbour->getY() << "]" << endmsg;

          // Check that the hit has not been used
          if (neighbour->used()) {
            debug() << "-- used" << endmsg;
            continue;
          }

          // Check that the hit is not in the opposite side of the detector (if endcap)
          if (neighbour->endcap() && kdhit->endcap() && (neighbour->forward() != kdhit->forward())) {
            debug() << "-- opposite side of detector" << endmsg;
            continue;
          }

          // Check that radial conditions are met
          if ((vertexToTracker && neighbour->getR() >= kdhit->getR()) ||
              (!vertexToTracker && neighbour->getR() <= kdhit->getR())) {
            debug() << "-- radial conditions not met" << endmsg;
            continue;
          }

          // Check that the hit is not wildly away from the track (make cell and check angle)
          Cell extensionCell(track->m_clusters[nclusters - 1], neighbour);
          double cellAngle = seedCell.getAngle(extensionCell);
          double cellAngleRZ = seedCell.getAngleRZ(extensionCell);
          double maxCellAngle = parameters.m_maxCellAngle;
          double maxCellAngleRZ = parameters.m_maxCellAngleRZ;

          if (cellAngle > maxCellAngle || cellAngleRZ > maxCellAngleRZ) {
            debug() << "-- killed by cell angle cut" << endmsg;
            continue;
          }

          // Now fit the track with the new hit and check the increase in chi2 - use a tempTrack object: add hit, fit,
          // get chi2, remove hit
          double deltaChi2(0.), deltaChi2zs(0.);
          debug() << "-- tempTrack has " << tempTrack.m_clusters.size() << " hits " << endmsg;
          tempTrack.add(neighbour);
          tempTrack.linearRegression();
          tempTrack.linearRegressionConformal();
          double newchi2 = tempTrack.chi2ndof();
          double newchi2zs = tempTrack.chi2ndofZS();
          debug() << "-- tempTrack has now " << tempTrack.m_clusters.size() << " hits " << endmsg;
          deltaChi2 = newchi2 - chi2;
          deltaChi2zs = newchi2zs - chi2zs;
          debug() << "-- hit was fitted and has a delta chi2 of " << deltaChi2 << " and delta chi2zs of " << deltaChi2zs
                  << endmsg;
          tempTrack.remove(tempTrack.m_clusters.size() - 1);
          debug() << "-- tempTrack has now " << tempTrack.m_clusters.size() << " hits " << endmsg;

          double chi2cut = parameters.m_chi2cut;
          if (deltaChi2 > chi2cut || deltaChi2zs > chi2cut) {
            debug() << "-- killed by chi2 cut" << endmsg;
            continue;
          }
          debug() << "-- valid candidate" << endmsg;

          bestClustersWithChi2[neighbour] = deltaChi2;

          // bestCluster still empty - fill it with the first candidate
          // otherwise fill it with the one with best chi2
          if (!bestCluster || deltaChi2 < bestChi2) {
            bestCluster = neighbour;
            bestChi2 = deltaChi2;
          } else {
            continue;
          }

        } // end if on the extension layer

      } // end loop on neighbours

      debug() << "-- this seed cells has " << bestClustersWithChi2.size() << " good candidates." << endmsg;

      // put the best cluster already found
      if (bestCluster) {
        bestClusters.push_back(bestCluster);
      }

      // put all the other cluster with a similar chi2 - probably duplicates
      float chi2window = 100.0;
      if (bestClustersWithChi2.size() > 0) {
        for (const auto& clu : bestClustersWithChi2) {
          debug() << "-- Best cluster candidate: [x,y] = [" << clu.first->getX() << ", " << clu.first->getY()
                  << "]; r = " << clu.first->getR() << "; radius = " << clu.first->getRadius() << endmsg;
          if (clu.second < bestChi2 + chi2window && clu.second > bestChi2 - chi2window) {
            if (!clu.first->sameSensor(bestCluster))
              bestClusters.push_back(clu.first);
            else
              debug() << "-- current candidate and best one are in the same sensor " << endmsg;
          }
        }
      }
      debug() << "-- this seed cells will be updated with " << bestClusters.size() << " candidates." << endmsg;
      std::ranges::sort(bestClusters, vertexToTracker ? sort_by_radiusKD : sort_by_lower_radiusKD);

      // If bestClusters have been found in this layer, add them to the track, update the seed cell and reset
      if (!bestClusters.empty()) {
        skipLayer = false;
        nclusters = track->m_clusters.size();
        debug() << "- nclusters = " << nclusters << endmsg;
        for (size_t i = 0; i < nclusters; i++) {
          debug() << "- Hit " << i << ": [x,y] = [" << track->m_clusters.at(i)->getX() << ", "
                  << track->m_clusters.at(i)->getY() << "]; r = " << track->m_clusters.at(i)->getR()
                  << "; radius = " << track->m_clusters.at(i)->getRadius() << endmsg;
        }
        // create new cell with last track cluster and last bestCluster
        // Best clusters are already ordered depending on vertexToTracker bool
        seedCell = Cell(track->m_clusters[nclusters - 1], bestClusters.at(bestClusters.size() - 1));

        for (const auto& bestClu : bestClusters) {
          debug() << "- Found bestCluster [x,y] = [" << bestClu->getX() << ", " << bestClu->getY()
                  << "]; r = " << bestClu->getR() << endmsg;
          track->add(bestClu);
          track->linearRegression();
          track->linearRegressionConformal();
          bestClu->used(true);
          nclusters = track->m_clusters.size();
          debug() << "- nclusters = " << nclusters << endmsg;
          for (size_t i = 0; i < nclusters; i++) {
            debug() << "- Hit " << i << ": [x,y] = [" << track->m_clusters.at(i)->getX() << ", "
                    << track->m_clusters.at(i)->getY() << "]; r = " << track->m_clusters.at(i)->getR()
                    << "; radius = " << track->m_clusters.at(i)->getRadius() << endmsg;
          }
        }
      }
      // If not bestCluster has been found in this layer, make cell with the expected hit (from extrapolation) and
      // increment the missing hit count
      else {
        debug() << "- Found no bestCluster for [subdet,layer] = [" << extendInSubdet << ", " << extendInLayer << "]"
                << endmsg;
        skipLayer = true;
      }
      debug() << (skipLayer ? "- Still" : "- Updated") << " seed cell A ([x,y] = [" << seedCell.getStart()->getX()
              << ", " << seedCell.getStart()->getY() << "]) - B ([x,y] = [" << seedCell.getEnd()->getX() << ", "
              << seedCell.getEnd()->getY() << ")]" << endmsg;

      // Clear the neighbours tree
      results.clear();
      bestClustersWithChi2.clear();
      bestClusters.clear();

    } while (loop); // end of track extension

    debug() << "Track ends with " << track->m_clusters.size() << " hits" << endmsg;
    // if (streamlog_level(DEBUG9)) {
    //   for (unsigned int i = 0; i < track->m_clusters.size(); i++) {
    //     debug() << "- Hit " << i << ": [x,y] = [" << track->m_clusters.at(i)->getX() << ", "
    //                           << track->m_clusters.at(i)->getY() << "], [subdet,layer] = ["
    //                           << track->m_clusters.at(i)->getSubdetector() << ", " <<
    //                           track->m_clusters.at(i)->getLayer()
    //                           << "]" << endmsg;
    //   }
    // }

  } // end loop on tracks
}

bool ConformalTracking::toBeUpdated(UniqueCellularTracks const& cellularTracks) const {
  for (const auto& aCellularTrack : cellularTracks)
    if (aCellularTrack->back()->getFrom().size() > 0)
      return true;
  return false;
}

// New test at creating cellular tracks. In this variant, don't worry about clones etc, give all possible routes back to
// the seed cell. Then cut on number of clusters on each track, and pass back (good tracks to then be decided based on
// best chi2
UniqueCellularTracks ConformalTracking::createTracksNew(const SCell& seedCell) const {
  debug() << "***** createTracksNew" << endmsg;

  // Final container to be returned
  UniqueCellularTracks cellularTracks;

  // Make the first cellular track using the seed cell
  auto seedTrack = std::make_unique<cellularTrack>();
  seedTrack->push_back(seedCell);
  cellularTracks.push_back(std::move(seedTrack));

  debug() << "Follow all paths from higher weighted cells back to the seed cell" << endmsg;
  // Now start to follow all paths back from this seed cell
  // While there are still tracks that are not finished (last cell weight 0), keep following their path
  while (toBeUpdated(cellularTracks)) {
    //   debug()<<"== Updating "<<cellularTracks.size()<<" tracks"<<std::endl;
    // Loop over all (currently existing) tracks
    size_t nTracks = cellularTracks.size();

    debug() << "Going to create " << nTracks << " tracks" << endmsg;
    if (nTracks > 10000) {
      warning() << "WARNING: Going to create " << nTracks << " tracks " << endmsg;
    }
    if (nTracks > m_tooManyTracks) {
      error() << "Too many tracks (" << nTracks << " > " << m_tooManyTracks
              << " are going to be created, tightening parameters" << endmsg;
      throw std::runtime_error("Too many tracks");
    }
    for (size_t itTrack = 0; itTrack < nTracks; itTrack++) {
      // If the track is finished, do nothing
      //      if(cellularTracks[itTrack].back()->getWeight() == 0) continue;
      if (cellularTracks[itTrack]->back()->getFrom().size() == 0) {
        debug() << "- Cellular track " << itTrack << " is finished " << endmsg;
        //       debug()<<"-- Track "<<itTrack<<" is finished"<<std::endl;
        continue;
      }

      // While there is only one path leading from this cell, follow that path
      SCell cell = cellularTracks[itTrack]->back();
      //     debug()<<"-- Track "<<itTrack<<" has "<<(*(cell->getFrom())).size()<<" cells attached to the end of
      //     it"<<std::endl;
      //      while(cell->getWeight() > 0 && (*(cell->getFrom())).size() == 1){
      while (cell->getFrom().size() == 1) {
        debug() << "- Cellular track " << itTrack << " is a simple extension" << endmsg;
        //       debug()<<"- simple extension"<<std::endl;
        // Get the cell that it attaches to
        auto parentCell = SCell(cell->getFrom()[0]);
        debug() << "- Added parent cell A ([x,y] = [" << parentCell->getStart()->getX() << ", "
                << parentCell->getStart()->getY() << "]) - B ([x,y] = [" << cell->getEnd()->getX() << ", "
                << cell->getEnd()->getY() << "])" << endmsg;
        // Attach it to the track and continue
        cellularTracks[itTrack]->push_back(parentCell);
        cell = parentCell;
      }

      // If the track is finished, do nothing
      //      if(cellularTracks[itTrack].back()->getWeight() == 0) continue;
      if (cellularTracks[itTrack]->back()->getFrom().size() == 0)
        continue;

      // If the weight is != 0 and there is more than one path to follow, branch the track (create a new one for each
      // path)
      debug() << "- Cellular track " << itTrack << " has more than one extension" << endmsg;

      //     debug()<<"- making "<<nBranches<<" branches"<<std::endl;

      // For each additional branch make a new track
      for (size_t itBranch = 1; itBranch < cell->getFrom().size(); itBranch++) {
        debug() << "-- Cellular track " << itTrack << ", extension " << itBranch << endmsg;
        auto branchedTrack = std::unique_ptr<cellularTrack>(new cellularTrack(*(cellularTracks[itTrack].get())));
        auto branchedParentCell = SCell(cell->getFrom()[itBranch]);
        debug() << "-- Added branched parent cell A ([x,y] = [" << branchedParentCell->getStart()->getX() << ", "
                << branchedParentCell->getStart()->getY() << "]) - B ([x,y] = [" << cell->getEnd()->getX() << ", "
                << cell->getEnd()->getY() << "])" << endmsg;
        branchedTrack->push_back(std::move(branchedParentCell));
        cellularTracks.push_back(std::move(branchedTrack));
      }

      // Keep the existing track for the first branch
      cellularTracks[itTrack]->push_back(SCell(cell->getFrom()[0]));
    }
  }

  debug() << "Number of finalcellularTracks = " << cellularTracks.size() << endmsg;
  // if (streamlog_level(DEBUG8)) {
  //   for (auto& cellTrack : finalcellularTracks) {
  //      debug() << "- Finalcelltrack is made of " << cellTrack->size() << " cells " << endmsg;
  //     for (unsigned int finalcell = 0; finalcell < cellTrack->size(); finalcell++) {
  //        debug() << "-- Cell A ([x,y] = [" << (*cellTrack)[finalcell]->getStart()->getX() << ", "
  //                             << (*cellTrack)[finalcell]->getStart()->getY() << "]) - B ([x,y] = ["
  //                             << (*cellTrack)[finalcell]->getEnd()->getX() << ", "
  //                             << (*cellTrack)[finalcell]->getEnd()->getY() << "])" << endmsg;
  //     }
  //   }
  // }
  debug() << "createTracksNew *****" << endmsg;

  return cellularTracks;
}

void ConformalTracking::buildNewTracks(UniqueKDTracks& conformalTracks, SharedKDClusters& collection,
                                       UKDTree& nearestNeighbours, Parameters const& parameters, bool radialSearch,
                                       bool vertexToTracker) const {
  debug() << "*** buildNewTracks" << endmsg;

  // Sort the input collection by radius - higher to lower if starting with the vertex detector (high R in conformal
  // space)
  std::ranges::sort(collection, vertexToTracker ? sort_by_radiusKD : sort_by_lower_radiusKD);

  // Loop over all hits, using each as a seed to produce a new track
  for (unsigned int nKDHit = 0; nKDHit < collection.size(); nKDHit++) {
    auto stopwatch_hit = TStopwatch();
    auto stopwatch_hit_total = TStopwatch();

    // Get the kdHit and check if it has already been used (assigned to a track)
    SKDCluster kdhit = collection[nKDHit];

    debug() << "Seed hit " << nKDHit << ": [x,y,z] = [" << kdhit->getX() << ", " << kdhit->getY() << ", "
            << kdhit->getZ() << "]" << endmsg;
    if (m_debugPlots) {
      ++m_X[kdhit->getX()];
      ++m_Y[kdhit->getY()];
      ++m_Z[kdhit->getZ()];
    }

    if (debugSeed && kdhit == debugSeed)
      debug() << "Starting to seed with debug cluster" << endmsg;
    if (kdhit->used()) {
      debug() << "hit already used" << endmsg;
      continue;
    }

    // Debug: Plot residuals between hit and associated SimTrackerHit
    // TODO
    // debug() << "SimHit : [x,y,z] = [" << kdSimHits[kdhit]->getPosition()[0] << ", "
    //                       << kdSimHits[kdhit]->getPosition()[1] << ", " << kdSimHits[kdhit]->getPosition()[2] << "] "
    //                       << endmsg;
    //      if(kdhit->getR() < 0.003) break; // new cut - once we get to inner radius we will never make tracks. temp?
    //      TODO: make parameter? FCC (0.005 to 0.003)

    // The tracking differentiates between the first and all subsequent hits on a chain.
    // First, take the seed hit and look for sensible hits nearby to make an initial set
    // of cells. Once these are found, extrapolate the cells and look for additional hits
    // to produce a new cell along the chain. This is done mainly for speed: you ignore
    // all combinations which would be excluded when compared with the seed cell. In the
    // end we will try to produce a track starting from this seed, then begin again for
    // the next seed hit.

    // Get the initial seed cells for this hit, by looking at neighbouring hits in a given
    // radial distance. Check that they are at lower radius and further from the IP
    SharedKDClusters results;
    double theta = kdhit->getTheta();
    // Filter already if the neighbour is used, is on the same detector layer,
    // or is in the opposite side of the detector and points inwards
    if (radialSearch)
      nearestNeighbours->allNeighboursInRadius(
          kdhit, parameters.m_maxDistance, results, [&kdhit, vertexToTracker](SKDCluster const& nhit) {
            if (nhit->used())
              return true;
            if (kdhit->sameLayer(nhit))
              return true;
            // not pointing in the same direction
            if (nhit->endcap() && kdhit->endcap() && (nhit->forward() != kdhit->forward()))
              return true;
            // radial conditions not met
            if ((vertexToTracker && nhit->getR() >= kdhit->getR()) ||
                (!vertexToTracker && nhit->getR() <= kdhit->getR()))
              return true;
            return false;
          });
    else
      nearestNeighbours->allNeighboursInTheta(
          theta, m_thetaRange, results, [&kdhit, vertexToTracker](SKDCluster const& nhit) {
            if (nhit->used())
              return true;
            if (kdhit->sameLayer(nhit))
              return true;
            // not pointing in the same direction
            if (nhit->endcap() && kdhit->endcap() && (nhit->forward() != kdhit->forward()))
              return true;
            // radial conditions not met
            if ((vertexToTracker && nhit->getR() >= kdhit->getR()) ||
                (!vertexToTracker && nhit->getR() <= kdhit->getR()))
              return true;

            return false;
          });

    if (m_debugTime)
      debug() << "  Time report: Searching for " << results.size() << " neighbours took "
              << stopwatch_hit.RealTime() * 1000 << std::scientific << " milli-seconds" << endmsg;
    stopwatch_hit.Start(true);

    debug() << "Picked up " << results.size() << " neighbours from " << (radialSearch ? "radial" : "theta") << " search"
            << endmsg;

    // Sort the neighbours by radius
    if (debugSeed && kdhit == debugSeed)
      debug() << "- picked up " << results.size() << " neighbours from " << (radialSearch ? "radial" : "theta")
              << " search" << endmsg;
    if (results.size() == 0)
      continue;
    std::ranges::sort(results, vertexToTracker ? sort_by_radiusKD : sort_by_lower_radiusKD);

    // Objects to hold cells
    SharedCells cells;
    [[maybe_unused]] bool isFirst = true;

    // Make seed cells pointing inwards/outwards (conformal space)
    for (unsigned int neighbour = 0; neighbour < results.size(); neighbour++) {
      // Get the neighbouring hit
      SKDCluster const& nhit = results[neighbour];

      debug() << "- Neighbour " << neighbour << ": [x,y,z] = [" << nhit->getX() << ", " << nhit->getY() << ", "
              << nhit->getZ() << "]" << endmsg;

      if (!neighbourIsCompatible(nhit, kdhit, parameters.m_maxSlopeZ)) {
        continue;
      }

      // Check if the cell would be too long (hit very far away)
      double length2 = ((kdhit->getU() - nhit->getU()) * (kdhit->getU() - nhit->getU()) +
                        (kdhit->getV() - nhit->getV()) * (kdhit->getV() - nhit->getV()));
      if (length2 > parameters.m_maxDistance * parameters.m_maxDistance) {
        debug() << "- cell between A ([x,y] = [" << kdhit->getX() << ", " << kdhit->getY() << "]) and B ([x,y] = ["
                << nhit->getX() << ", " << nhit->getY() << "]) is too long" << endmsg;
        continue;
      }

      if (m_debugPlots) {
        m_neighX->Fill(nhit->getX());
        m_neighY->Fill(nhit->getY());
        m_neighZ->Fill(nhit->getZ());

        double distanceX = nhit->getX() - kdhit->getX();
        double distanceY = nhit->getY() - kdhit->getY();
        double distance = sqrt(distanceX * distanceX + distanceY * distanceY);
        double deltaZ = nhit->getZ() - kdhit->getZ();
        double slopeZ = deltaZ / distance;
        m_slopeZ->Fill(slopeZ);

        // Debug using the seed hit and the associated SimTrackerHit
        // if (m_debugPlots && kdParticles[kdhit] == kdParticles[nhit] && kdhit != nhit) {
        //   double simpt = sqrt(kdParticles[kdhit]->getMomentum()[0] * kdParticles[kdhit]->getMomentum()[0] +
        //                       kdParticles[kdhit]->getMomentum()[1] * kdParticles[kdhit]->getMomentum()[1]);
        //   debug() << "- They were produced by the same MCParticle with pt = " << simpt << endmsg;
        //   debug() << "- SimHit : [x,y,z] = [" << kdSimHits[nhit]->getPosition()[0] << ", "
        //           << kdSimHits[nhit]->getPosition()[1] << ", " << kdSimHits[nhit]->getPosition()[2] << "] " <<
        //           endmsg;
        //   debug() << "- Delta : [x,y,z] = [" << nhit->getX() - kdhit->getX() << ", " << nhit->getY() - kdhit->getY()
        //           << ", " << nhit->getZ() - kdhit->getZ() << "] " << endmsg;

        //   m_slopeZ_true->Fill(slopeZ);
        //   m_slopeZ_vs_pt_true->Fill(slopeZ, simpt);

        //   if (isFirst) {
        //     debug() << "- is first" << endmsg;
        //     m_slopeZ_true_first->Fill(slopeZ);
        //     isFirst = false;
        //   }
        // }
        // Debugging: uncomment in the case you want to create seeds only with hits belonging to the same MCParticle
        // else {
        //  continue;
        //}
      }

      // Create the new seed cell
      cells.emplace_back(std::make_shared<Cell>(kdhit, nhit));

      debug() << "- made cell between A ([x,y] = [" << kdhit->getX() << ", " << kdhit->getY() << "]) and B ([x,y] = ["
              << nhit->getX() << ", " << nhit->getY() << "])" << endmsg;
      if (debugSeed && kdhit == debugSeed)
        debug() << "- made cell with neighbour " << neighbour << " at " << nhit->getU() << "," << nhit->getV()
                << endmsg;

      // Debug plotting
      // if (m_debugPlots) {
      //   const auto& cell = cells.back();
      //   m_cellDOCA->Fill(cell->doca());
      //   if (m_eventNumber == 0) {
      //     m_canvConformalEventDisplayAllCells->cd();
      //     drawline(kdhit, nhit, 1);
      //   }
      // }
    }

    if (m_debugTime)
      debug() << "  Time report: Making " << cells.size() << " seed cells took " << stopwatch_hit.RealTime() * 1000
              << std::scientific << " milli-seconds" << endmsg;
    stopwatch_hit.Start(true);

    debug() << "Produced " << cells.size() << " seed cells from seed hit A ([x,y] = [" << kdhit->getX() << ", "
            << kdhit->getY() << "])" << endmsg;
    if (debugSeed && kdhit == debugSeed)
      debug() << "- produced " << cells.size() << " seed cells" << endmsg;

    // No seed cells produced
    if (cells.size() == 0)
      continue;

    // All seed cells have been created, now try create all "downstream" cells until no more can be added
    SharedKDClusters debugHits;
    extendSeedCells(cells, nearestNeighbours, true, debugHits, parameters, vertexToTracker);

    if (m_debugTime)
      debug() << "  Time report: Extending " << cells.size() << " seed cells took " << stopwatch_hit.RealTime() * 1000
              << std::scientific << " milli-seconds" << endmsg;
    stopwatch_hit.Start(true);

    debug() << "After extension, have " << cells.size() << " cells from seed hit A ([x,y] = [" << kdhit->getX() << ", "
            << kdhit->getY() << "])" << endmsg;
    if (debugSeed && kdhit == debugSeed)
      debug() << "- after extension, have " << cells.size() << " cells" << endmsg;

    // Now have all cells stemming from this seed hit. If it is possible to produce a track (ie. cells with depth X)
    // then we will now...
    //      if(depth < (m_minClustersOnTrack-1)) continue; // TODO: check if this is correct

    // We create all acceptable tracks by looping over all cells with high enough weight to create
    // a track and trace their route back to the seed hit. We then have to choose the best candidate
    // at the end (by minimum chi2 of a linear fit)
    std::map<SCell, bool> usedCells;
    UniqueKDTracks cellTracks;

    // Sort Cells from highest to lowest weight
    std::ranges::sort(cells, [](const SCell& a, const SCell& b) { return a->getWeight() > b->getWeight(); });

    // Create tracks by following a path along cells
    debug() << "Create tracks by following path along cells. Loop over cells, sorted by weight" << endmsg;

    for (size_t itCell = 0; itCell < cells.size(); itCell++) {
      debug() << "- Cell " << itCell << " between A ([x,y] = [" << cells[itCell]->getStart()->getX() << ", "
              << cells[itCell]->getStart()->getY() << "]) and B ([x,y] = [" << cells[itCell]->getEnd()->getX() << ", "
              << cells[itCell]->getEnd()->getY() << "]) has weight " << cells[itCell]->getWeight() << endmsg;
      // Check if this cell has already been used
      if (debugSeed && kdhit == debugSeed)
        debug() << "-- looking at cell " << itCell << endmsg;
      if (usedCells.count(cells[itCell])) {
        debug() << "-- used cell" << endmsg;
        continue;
      }
      // Check if this cell could produce a track (is on a long enough chain)

      if (cells[itCell]->getWeight() < (parameters.m_minClustersOnTrack - 2)) {
        debug() << "-- cell can not produce a track: weight < (minClustersOnTrack - 2)" << endmsg;
        break;
      }
      // Produce all tracks leading back to the seed hit from this cell
      UniqueCellularTracks candidateTracks;
      UniqueCellularTracks candidateTracksTemp =
          createTracksNew(cells[itCell]); // Move back to using used cells here? With low chi2/ndof?

      copy_if(std::make_move_iterator(candidateTracksTemp.begin()), std::make_move_iterator(candidateTracksTemp.end()),
              std::back_inserter(candidateTracks), [&parameters](UcellularTrack const& track) {
                return (int(track->size()) >= (parameters.m_minClustersOnTrack - 1));
              });
      candidateTracksTemp.clear();

      debug() << "- From cell, produced " << candidateTracks.size() << " candidate tracks" << endmsg;
      // if (streamlog_level(DEBUG9)) {
      //   for (auto& candidateTrack : candidateTracks) {
      //      debug() << "--  track is made of " << candidateTrack->size() << " cells " << endmsg;
      //     for (unsigned int candidateCell = 0; candidateCell < candidateTrack->size(); candidateCell++) {
      //        debug() << "--- cell between A ([x,y] = [" << (*candidateTrack)[candidateCell]->getStart()->getX()
      //                             << ", " << (*candidateTrack)[candidateCell]->getStart()->getY() << "]) and B ([x,y]
      //                             = ["
      //                             << (*candidateTrack)[candidateCell]->getEnd()->getX() << ", "
      //                             << (*candidateTrack)[candidateCell]->getEnd()->getY() << "])" << endmsg;
      //     }
      //   }
      // }
      // Debug plotting
      // if (m_debugPlots && m_eventNumber == 2) {
      //   m_canvConformalEventDisplayAcceptedCells->cd();
      //   for (auto& track : candidateTracks) {
      //     for (unsigned int trackCell = 0; trackCell < track->size(); trackCell++) {
      //       drawline((*track)[trackCell]->getStart(), (*track)[trackCell]->getEnd(), track->size() - trackCell);
      //     }
      //   }
      // }

      // Look at the candidate tracks and fit them (+pick the best chi2/ndof)
      //        debug()<<"- produced "<<candidateTracks.size()<<" candidate tracks"<<std::endl;
      if (debugSeed && kdhit == debugSeed)
        debug() << "- produced " << candidateTracks.size() << " candidate tracks" << endmsg;
      if (candidateTracks.size() == 0)
        continue;
      std::vector<double> chi2ndof;

      // Temporary check of how many track candidates should not strictly have been allowed
      // checkUnallowedTracks(candidateTracks, parameters);

      UniqueKDTracks bestTracks =
          getFittedTracks(candidateTracks,
                          parameters); // Returns all tracks at the moment, not lowest chi2 CHANGE ME

      // Store track(s) for later
      cellTracks.insert(cellTracks.end(), std::make_move_iterator(bestTracks.begin()),
                        std::make_move_iterator(bestTracks.end()));
    }

    debug() << "Final number of fitted tracks to this seed hit: " << cellTracks.size() << endmsg;
    if (m_debugTime)
      debug() << "  Time report: " << cellTracks.size() << " tracks reconstructed from cells took "
              << stopwatch_hit.RealTime() * 1000 << std::scientific << " milli-seconds" << endmsg;
    stopwatch_hit.Start(true);

    // All tracks leading back to the seed hit have now been found. Decide which are the feasible candidates (may be
    // more than 1)
    if (debugSeed && kdhit == debugSeed)
      debug() << "== final number of candidate tracks to this seed hit: " << cellTracks.size() << endmsg;
    if (cellTracks.size() == 0) {
      continue;
    }

    UniqueKDTracks bestTracks = std::move(cellTracks); // CHANGE ME - temp to give all tracks
    if (debugSeed && kdhit == debugSeed) {
      debug() << "== final number of stored tracks to this seed hit: " << bestTracks.size() << endmsg;
      for (size_t itBest = 0; itBest < bestTracks.size(); itBest++)
        debug() << "- track " << itBest << " has chi2/ndof " << bestTracks[itBest]->chi2ndof() << endmsg;
    }

    // Could now think to do the full helix fit and apply a chi2 cut. TODO

    // Store the final CA tracks. First sort them by length, so that if we produce a clone with just 1 hit missing, the
    // larger track is taken

    debug() << "Now cut on chi2 and treate the clones. Loop over tracks, sorted by length" << endmsg;

    std::ranges::sort(bestTracks,
                      [](const UKDTrack& a, const UKDTrack& b) { return a->m_clusters.size() > b->m_clusters.size(); });
    for (auto& bestTrack : bestTracks) {
      bool bestTrackUsed = false;

      // Cut on chi2
      double chi2cut = parameters.m_chi2cut;

      if ((parameters.m_onlyZSchi2cut && bestTrack->chi2ndofZS() > chi2cut) ||
          (!parameters.m_onlyZSchi2cut && (bestTrack->chi2ndof() > chi2cut || bestTrack->chi2ndofZS() > chi2cut))) {
        debug() << "- Track has chi2 too large" << endmsg;
        bestTrack.reset();
        continue;
      }

      // Check if the new track is a clone
      bool clone = false;

      for (auto& conformalTrack : conformalTracks) {
        const unsigned int nOverlappingHits = overlappingHits(bestTrack, conformalTrack);
        if (nOverlappingHits >= 2) {
          clone = true;
          debug() << "- Track is a clone" << endmsg;

          // Calculate the new and existing chi2 values
          double newchi2 = bestTrack->chi2ndofZS() + bestTrack->chi2ndof();
          double oldchi2 = conformalTrack->chi2ndofZS() + conformalTrack->chi2ndof();

          // If the new track is an existing track + segment, take the new track
          if (nOverlappingHits == conformalTrack->m_clusters.size()) {
            conformalTrack = std::move(bestTrack);
            bestTrackUsed = true;
            debug() << "- New = existing + segment. Replaced existing with new" << endmsg;
          }

          // If the new track is a subtrack of an existing track, don't consider it further (already try removing bad
          // hits from tracks

          else if (nOverlappingHits == bestTrack->m_clusters.size()) {
            debug() << "- New + segment = existing. New ignored" << endmsg;
            break;
          }
          // Otherwise take the longest
          else if (bestTrack->m_clusters.size() == conformalTrack->m_clusters.size()) { // New track equal in length
            if (newchi2 > oldchi2) {
              debug() << "- New equal length. Worse chi2:" << newchi2 << endmsg;
              break;
            }
            // Take it
            conformalTrack = std::move(bestTrack);
            bestTrackUsed = true;
            debug() << "- New equal. Better chi2 Replaced existing with new" << endmsg;

          } else if (bestTrack->m_clusters.size() > conformalTrack->m_clusters.size()) { // New track longer

            // Take it
            conformalTrack = std::move(bestTrack);
            bestTrackUsed = true;
            debug() << "- New longer. Replaced existing with new" << endmsg;

          } else if (bestTrack->m_clusters.size() < conformalTrack->m_clusters.size()) { // Old track longer
            debug() << "- Old longer. New ignored" << endmsg;
            break;
          }
          break;
        }
      }

      // If not a clone, save the new track
      if (!clone) {
        bestTrackUsed = true;

        debug() << "- Track is not a clone. Pushing back best track with chi2/ndof " << bestTrack->chi2ndof()
                << " and hits: " << endmsg;

        if (debugSeed && kdhit == debugSeed) {
          debug() << "== Pushing back best track with chi2/ndof " << bestTrack->chi2ndof() << endmsg;
        } else {
          debug() << "Pushing back best track with chi2/ndof " << bestTrack->chi2ndof() << endmsg;
        }

        // for (unsigned int cluster = 0; cluster < bestTrack->m_clusters.size() && streamlog_level(DEBUG9); cluster++)
        // {
        //   debug() << "-- Hit " << cluster << ": [x,y] = [" << bestTrack->m_clusters.at(cluster)->getX() << ", "
        //                         << bestTrack->m_clusters.at(cluster)->getY() << "]" << endmsg;
        // }

        conformalTracks.push_back(std::move(bestTrack));
      }

      if (not bestTrackUsed) {
        bestTrack.reset();
      }

    } // end for besttracks
    if (m_debugTime) {
      debug() << "  Time report: Sort best tracks took " << stopwatch_hit.RealTime() * 1000 << std::scientific
              << " milli-seconds" << endmsg;
      debug() << " Time report: Total time for seed hit " << nKDHit << " = " << stopwatch_hit_total.RealTime() * 1000
              << std::scientific << " milli-seconds" << endmsg;
    }
  }
}

void ConformalTracking::runStep(SharedKDClusters& kdClusters, UKDTree& nearestNeighbours,
                                UniqueKDTracks& conformalTracks,
                                std::map<size_t, SharedKDClusters> const& collectionClusters,
                                const Parameters& parameters) const {
  auto stopwatch = TStopwatch();
  stopwatch.Start(false);

  if (parameters.m_combine) {
    nearestNeighbours = combineCollections(kdClusters, parameters.m_collections, collectionClusters);
  }

  if (parameters.m_build) {
    bool caughtException = false;
    Parameters thisParameters(parameters);
    do {
      caughtException = false;
      try {
        buildNewTracks(conformalTracks, kdClusters, nearestNeighbours, thisParameters, parameters.m_radialSearch,
                       parameters.m_vertexToTracker);
      } catch (std::runtime_error&) {
        info() << "caught too many tracks, tightening parameters" << endmsg;
        caughtException = true;
        thisParameters.tighten();
        if (not m_retryTooManyTracks || thisParameters.m_tightenStep > 10) {
          error() << "Skipping event" << endmsg;
          throw;
        }
      }
    } while (caughtException);

    stopwatch.Stop();
    if (m_debugTime)
      debug() << "Step " << parameters.m_step << " buildNewTracks took " << stopwatch.RealTime() << " seconds"
              << endmsg;
    stopwatch.Reset();
  }
  if (parameters.m_extend) {
    extendTracks(conformalTracks, kdClusters, nearestNeighbours, parameters);
    stopwatch.Stop();
    if (m_debugTime)
      debug() << "Step " << parameters.m_step << " extendTracks took " << stopwatch.RealTime() << " seconds" << endmsg;
    stopwatch.Reset();
  }

  if (parameters.m_sortTracks) {
    // FIXME: Still needed?
    // Sort by pt (low to hight)
    std::ranges::sort(conformalTracks,
                      [](const UKDTrack& track1, const UKDTrack& track2) { return (track1->pt() > track2->pt()); });
  }

  // Mark hits from "good" tracks as being used
  for (const auto& conformalTrack : conformalTracks)
    for (const auto& thisCluster : conformalTrack->m_clusters)
      thisCluster->used(true);
}

// Draw a line on the current canvas
void ConformalTracking::drawline(const SKDCluster& hitStart, const SKDCluster& hitEnd, int colour, int style) const {
  TLine* line = new TLine(hitStart->getU(), hitStart->getV(), hitEnd->getU(), hitEnd->getV());
  line->SetLineColor(colour);
  line->SetLineStyle(style);
  line->Draw();
}
