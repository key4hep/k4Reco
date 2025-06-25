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
#ifndef K4RECO_TRUTHTRACKFINDER_H
#define K4RECO_TRUTHTRACKFINDER_H 1

#include "GaudiDDKalTest.h"

#include <DDSegmentation/BitFieldCoder.h>

#include <Gaudi/Property.h>

#include <edm4hep/MCParticleCollection.h>
#include <edm4hep/TrackCollection.h>
#include <edm4hep/TrackMCParticleLinkCollection.h>
#include <edm4hep/TrackerHitPlaneCollection.h>
#include <edm4hep/TrackerHitSimTrackerHitLinkCollection.h>

#include <k4FWCore/Transformer.h>
#include <k4Interface/IGeoSvc.h>

#include <string>
#include <vector>

struct TruthTrackFinder final
    : k4FWCore::MultiTransformer<std::tuple<edm4hep::TrackCollection, edm4hep::TrackMCParticleLinkCollection>(
          const std::vector<const edm4hep::TrackerHitPlaneCollection*>&,
          const std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>&,
          const std::vector<const edm4hep::MCParticleCollection*>&)> {
  TruthTrackFinder(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;

  std::tuple<edm4hep::TrackCollection, edm4hep::TrackMCParticleLinkCollection>
  operator()(const std::vector<const edm4hep::TrackerHitPlaneCollection*>& trackerHitCollections,
             const std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>& relations,
             const std::vector<const edm4hep::MCParticleCollection*>& particleCollections) const override;

  template <typename T>
  std::vector<const T*> removeHitsSameLayer(const std::vector<const T*>& trackHits) const;

  // Track fit parameters
  float m_initialTrackError_d0{1.e6};
  float m_initialTrackError_phi0{1.e2};
  float m_initialTrackError_omega{static_cast<float>(1.e-4)};
  float m_initialTrackError_z0{1.e6};
  float m_initialTrackError_tanL{1.e2};
  double m_maxChi2perHit{1.e2};
  float m_magneticField{};

  Gaudi::Property<bool> m_useTruthInPrefit{
      this, "UseTruthInPrefit", false,
      "If true use the truth information to initialise the helical prefit, otherwise use prefit by fitting 3 hits"};
  Gaudi::Property<bool> m_fitForward{this, "FitForward", false,
                                     "If true fit 'forward' (go forward + smooth back adding last two hits with Kalman "
                                     "FIlter steps), otherwise fit backward "};
  Gaudi::Property<std::string> m_geoSvcName{this, "GeoSvcName", "GeoSvc", "The name of the GeoSvc instance"};
  Gaudi::Property<std::string> m_encodingStringVariable{
      this, "EncodingStringParameterName", "GlobalTrackerReadoutID",
      "The name of the DD4hep constant that contains the Encoding string for tracking detectors"};
  GaudiDDKalTest m_ddkaltest{this};
  SmartIF<IGeoSvc> m_geoSvc;
  dd4hep::DDSegmentation::BitFieldCoder m_encoder;
};

#endif

DECLARE_COMPONENT(TruthTrackFinder)
