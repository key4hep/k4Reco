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
#ifndef K4RECO_REFITFINAL_H
#define K4RECO_REFITFINAL_H 1

#include "GaudiDDKalTest.h"
#include "GaudiDDKalTestTrack.h"

#include <DDSegmentation/BitFieldCoder.h>

#include <edm4hep/TrackCollection.h>
#include <edm4hep/TrackMCParticleLinkCollection.h>

#include <k4FWCore/Transformer.h>
#include <k4Interface/IGeoSvc.h>

#include <Gaudi/Property.h>

#include <limits>
#include <string>

struct RefitFinal final : k4FWCore::MultiTransformer<
    std::tuple<edm4hep::TrackCollection, edm4hep::TrackMCParticleLinkCollection>(
    const edm4hep::TrackCollection&, const std::vector<const edm4hep::TrackMCParticleLinkCollection*>&)> {
  RefitFinal(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;

  std::tuple<edm4hep::TrackCollection, edm4hep::TrackMCParticleLinkCollection> operator()(
    const edm4hep::TrackCollection&,
    const std::vector<const edm4hep::TrackMCParticleLinkCollection*>&
  ) const override;

  Gaudi::Property<bool> m_MSOn{this, "MultipleScatteringOn", true, "Use MultipleScattering in Fit"};
  Gaudi::Property<bool> m_ElossOn{this, "EnergyLossOn", true, "Use Energy Loss in Fit"};
  Gaudi::Property<bool> m_SmoothOn{this, "SmoothOn", false, "Smooth All Mesurement Sites in Fit"};
  
  Gaudi::Property<float> m_initialTrackError_d0{this, "InitialTrackErrorD0", 1.e6f,
        "Value used for the initial d0 variance of the trackfit"};
  Gaudi::Property<float> m_initialTrackError_phi0{this, "InitialTrackErrorPhi0", 1.e2f,
        "Value used for the initial phi0 variance of the trackfit"};
  Gaudi::Property<float> m_initialTrackError_omega{this, "InitialTrackErrorOmega", 1.e-4f,
        "Value used for the initial omega variance of the trackfit"};
  Gaudi::Property<float> m_initialTrackError_z0{this, "InitialTrackErrorZ0", 1.e6f,
        "Value used for the initial z0 variance of the trackfit"};
  Gaudi::Property<float> m_initialTrackError_tanL{this, "InitialTrackErrorTanL", 1.e2f,
        "Value used for the initial tanL variance of the trackfit"};
  Gaudi::Property<double> m_Max_Chi2_Incr{this, "Max_Chi2_Incr", std::numeric_limits<double>::max(),
        "maximum allowable chi2 increment when moving from one site to another"};
  Gaudi::Property<int> m_refPoint{this, "ReferencePoint", -1,
        "Identifier of the reference point to use for the fit initialisation, -1 means at 0 0 0"};
  Gaudi::Property<bool> m_extrapolateForward{this, "extrapolateForward", true,
        "if true extrapolation in the forward direction (in-out), otherwise backward (out-in)"};
  Gaudi::Property<int> m_minClustersOnTrackAfterFit{this, "MinClustersOnTrackAfterFit", 4,
        "Final minimum number of track clusters"};
  Gaudi::Property<int> m_maxOutliersAllowed{this, "MaxOutliersAllowed", 99,
        "Maximum number of outliers allowed on the refitted track"};
  Gaudi::Property<double> m_ReducedChi2Cut{this, "ReducedChi2Cut", -1.0,
        "Cut on maximum allowed reduced chi2"};

  Gaudi::Property<std::string> m_geoSvcName{this, "GeoSvcName", "GeoSvc", "The name of the GeoSvc instance"};
  Gaudi::Property<std::string> m_encodingStringVariable{this, "EncodingStringParameterName", "GlobalTrackerReadoutID",
      "The name of the DD4hep constant that contains the Encoding string for tracking detectors"};

  GaudiDDKalTest m_ddkaltest{this};

  SmartIF<IGeoSvc> m_geoSvc;
  float m_bField = 5.0;
  dd4hep::DDSegmentation::BitFieldCoder m_encoder;
};

#endif

DECLARE_COMPONENT(RefitFinal)
