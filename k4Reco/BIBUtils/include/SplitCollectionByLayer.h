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
#ifndef K4RECO_SPLITCOLLECTIONBYLAYER_H
#define K4RECO_SPLITCOLLECTIONBYLAYER_H 1

#include "Gaudi/Property.h"

#include <edm4hep/TrackerHitPlaneCollection.h>

#include <k4FWCore/Transformer.h>
#include <k4Interface/IGeoSvc.h>

#include <GaudiKernel/StatusCode.h>

#include <string>
#include <vector>

/** === SplitCollectionByLayer ===
 *  Splits a single tracker-hit collection into several output collections
 *  according to the layer number decoded from each hit's cellID. Every output
 *  collection i collects the hits whose layer lies in the closed interval
 *  [StartLayers[i], EndLayers[i]]; a hit is copied into every output whose
 *  interval contains its layer, so overlapping intervals duplicate the hit as
 *  in the original processor. The outputs are subset collections referencing
 *  the original hits.
 *
 *  Gaudi-native port of the Marlin SplitCollectionByLayer processor from
 *  MarlinTrkProcessors (MuonColliderSoft). The original processor dispatched on
 *  the LCIO hit type at run time; this port targets edm4hep::TrackerHitPlane
 *  collections, consistent with the other BIBUtils tracker algorithms. The
 *  number of output collections is arbitrary and follows the length of the
 *  OutputCollections list, which must match StartLayers and EndLayers.
 *
 *  @author F. Gaede, DESY (original Marlin processor)
 */
struct SplitCollectionByLayer final : k4FWCore::Transformer<std::vector<edm4hep::TrackerHitPlaneCollection>(
                                          const edm4hep::TrackerHitPlaneCollection&)> {
  SplitCollectionByLayer(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;

  std::vector<edm4hep::TrackerHitPlaneCollection>
  operator()(const edm4hep::TrackerHitPlaneCollection& trackerHits) const override;

private:
  Gaudi::Property<std::vector<int>> m_startLayers{
      this, "StartLayers", {}, "First layer (inclusive) routed to each output collection"};
  Gaudi::Property<std::vector<int>> m_endLayers{
      this, "EndLayers", {}, "Last layer (inclusive) routed to each output collection"};
  Gaudi::Property<std::string> m_encodingStringVariable{
      this, "EncodingStringParameterName", "GlobalTrackerReadoutID",
      "The name of the DD4hep constant that contains the cellID encoding string for tracking detectors"};
  Gaudi::Property<std::string> m_geoSvcName{this, "GeoSvcName", "GeoSvc", "The name of the GeoSvc instance"};

  SmartIF<IGeoSvc> m_geoSvc;
};

#endif
