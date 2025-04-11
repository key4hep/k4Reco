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

struct RefitFinal final
    : k4FWCore::MultiTransformer<std::tuple<edm4hep::TrackCollection, edm4hep::TrackMCParticleLinkCollection>(
          const edm4hep::TrackCollection&, const std::vector<const edm4hep::TrackMCParticleLinkCollection*>&)> {
  RefitFinal(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;

  std::tuple<edm4hep::TrackCollection, edm4hep::TrackMCParticleLinkCollection>
  operator()(const edm4hep::TrackCollection&,
             const std::vector<const edm4hep::TrackMCParticleLinkCollection*>&) const override;

  int FitInit2(const edm4hep::Track& track, GaudiDDKalTestTrack& _marlinTrk) const;

  // /* helper function to get collection using try catch block */
  // lcio::LCCollection* GetCollection(lcio::LCEvent* evt, std::string colName);

  // /** Input track collection name for refitting.
  //  */
  // std::string _input_track_col_name = "TruthTracks";

  // /** output track collection name.
  //  */
  // std::string _output_track_col_name = "RefittedTracks";

  // /** Input track relations name.
  //  */
  // std::string _input_track_rel_name = "SiTrackRelations";

  // /** Output track relations name for refitting.
  //  */
  // std::string _output_track_rel_name = "RefittedRelation";

  // /** pointer to the IMarlinTrkSystem instance
  //  */
  // MarlinTrk::IMarlinTrkSystem* _trksystem = nullptr;

  // int _n_run = -1;

  // bool _MSOn = true;
  // bool _ElossOn = true;
  // bool _SmoothOn = false;
  // double _Max_Chi2_Incr = DBL_MAX;
  // int _refPoint = -1;

  // bool _extrapolateForward = true;
  // int _minClustersOnTrackAfterFit = 0;

  // std::shared_ptr<UTIL::BitField64> _encoder{};

  // registerProcessorParameter("Max_Chi2_Incr", "maximum allowable chi2 increment when moving from one site to
  // another",
  //                            _Max_Chi2_Incr, _Max_Chi2_Incr);

  Gaudi::Property<bool> m_MSOn{this, "MultipleScatteringOn", true, "Use MultipleScattering in Fit"};
  Gaudi::Property<bool> m_ElossOn{this, "EnergyLossOn", true, "Use Energy Loss in Fit"};
  Gaudi::Property<bool> m_SmoothOn{this, "SmoothOn", false, "Smooth All Mesurement Sites in Fit"};
  Gaudi::Property<double> m_Max_Chi2_Incr{this, "Max_Chi2_Incr", std::numeric_limits<double>::max(),
                                          "maximum allowable chi2 increment when moving from one site to another"};
  Gaudi::Property<int> m_refPoint{
      this, "ReferencePoint", -1,
      "Identifier of the reference point to use for the fit initialisation, -1 means at 0 0 0"};
  Gaudi::Property<bool> m_extrapolateForward{
      this, "extrapolateForward", true,
      "if true extrapolation in the forward direction (in-out), otherwise backward (out-in)"};
  Gaudi::Property<int> m_minClustersOnTrackAfterFit{this, "MinClustersOnTrackAfterFit", 4,
                                                    "Final minimum number of track clusters"};

  Gaudi::Property<std::string> m_geoSvcName{this, "GeoSvcName", "GeoSvc", "The name of the GeoSvc instance"};
  Gaudi::Property<std::string> m_encodingStringVariable{
      this, "EncodingStringParameterName", "GlobalTrackerReadoutID",
      "The name of the DD4hep constant that contains the Encoding string for tracking detectors"};
  GaudiDDKalTest m_ddkaltest{this};

  SmartIF<IGeoSvc> m_geoSvc;
  dd4hep::DDSegmentation::BitFieldCoder m_encoder;
};

#endif

DECLARE_COMPONENT(RefitFinal)
