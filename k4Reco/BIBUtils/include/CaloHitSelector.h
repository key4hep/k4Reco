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
#ifndef K4RECO_CALOHITSELECTOR_H
#define K4RECO_CALOHITSELECTOR_H 1

#include "Gaudi/Property.h"

#include <edm4hep/CaloHitSimCaloHitLinkCollection.h>
#include <edm4hep/CalorimeterHitCollection.h>

#include <k4FWCore/Transformer.h>
#include <k4Interface/IGeoSvc.h>

#include <GaudiKernel/StatusCode.h>

#include <memory>
#include <string>
#include <tuple>

class TH2D;

/** === CaloHitSelector ===
 *  Applies per-cell energy and timing selections to calorimeter hits in order
 *  to suppress beam-induced background (BIB). The energy threshold is read, as
 *  a function of polar angle and calorimeter layer, from two ROOT histograms
 *  (th_2dmode_sym, stddev_sym) stored in the file pointed to by
 *  ThresholdsFilePath: threshold = mode + Nsigma * stddev. A constant
 *  FlatThreshold (GeV) can be used instead. Hits passing the threshold are
 *  further required to fall inside a time window, after correcting the hit time
 *  for the time of flight from the origin. Selected hits are written to a
 *  subset collection together with a freshly built reco-to-sim link collection.
 *
 *  Gaudi-native port of the Marlin CaloHitSelector processor from MyBIBUtils
 *  (https://github.com/madbaron/MyBIBUtils).
 *
 *  @author F. Meloni, DESY (original Marlin processor)
 */
struct CaloHitSelector final
    : k4FWCore::MultiTransformer<
          std::tuple<edm4hep::CalorimeterHitCollection, edm4hep::CaloHitSimCaloHitLinkCollection>(
              const edm4hep::CalorimeterHitCollection&, const edm4hep::CaloHitSimCaloHitLinkCollection&)> {
  CaloHitSelector(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;

  std::tuple<edm4hep::CalorimeterHitCollection, edm4hep::CaloHitSimCaloHitLinkCollection>
  operator()(const edm4hep::CalorimeterHitCollection& caloHits,
             const edm4hep::CaloHitSimCaloHitLinkCollection& caloLinks) const override;

private:
  Gaudi::Property<std::string> m_thFile{this, "ThresholdsFilePath", "",
                                        "Path to the ROOT file holding the threshold maps (th_2dmode_sym, stddev_sym)"};
  Gaudi::Property<int> m_nSigma{this, "Nsigma", 3, "Number of BIB-energy sigmas added on top of the modal threshold"};
  Gaudi::Property<double> m_flatThreshold{this, "FlatThreshold", 0.,
                                          "Constant energy threshold [GeV]; if > 0 it overrides the map-based one"};
  Gaudi::Property<double> m_timeWindowMin{this, "TimeWindowMin", -0.5,
                                          "Lower edge of the time-of-flight corrected acceptance window [ns]"};
  Gaudi::Property<double> m_timeWindowMax{this, "TimeWindowMax", 10.,
                                          "Upper edge of the time-of-flight corrected acceptance window [ns]"};
  Gaudi::Property<bool> m_doBIBsubtraction{this, "DoBIBsubtraction", false,
                                           "Subtract the mean expected BIB energy (modal map) from each cell"};
  Gaudi::Property<std::string> m_encodingStringVariable{
      this, "EncodingStringParameterName", "GlobalCalorimeterReadoutID",
      "Name of the DD4hep constant holding the cellID encoding string for calorimeters"};

  SmartIF<IGeoSvc> m_geoSvc;

  std::unique_ptr<TH2D> m_thresholdMap;
  std::unique_ptr<TH2D> m_stddevMap;
};

#endif
