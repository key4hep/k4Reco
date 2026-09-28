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
#ifndef K4RECO_SPLITCOLLECTIONBYPOLARANGLE_H
#define K4RECO_SPLITCOLLECTIONBYPOLARANGLE_H 1

#include "Gaudi/Accumulators/RootHistogram.h"
#include "Gaudi/Property.h"

#include <edm4hep/SimTrackerHitCollection.h>
#include <edm4hep/TrackerHitPlaneCollection.h>
#include <edm4hep/TrackerHitSimTrackerHitLinkCollection.h>

#include <k4FWCore/Transformer.h>

#include <GaudiKernel/StatusCode.h>

#include <array>
#include <memory>
#include <string>
#include <tuple>

#include "GAUDI_VERSION.h"

#if GAUDI_MAJOR_VERSION < 39
namespace Gaudi::Accumulators {
template <unsigned int ND, atomicity Atomicity = atomicity::full, typename Arithmetic = double>
using StaticRootHistogram =
    Gaudi::Accumulators::RootHistogramingCounterBase<ND, Atomicity, Arithmetic, naming::histogramString>;
}
#endif

/** === SplitCollectionByPolarAngle ===
 *  Selects the tracker hits whose polar angle theta = acos(z/r) lies inside the
 *  configurable window [PolarAngleLowerLimit, PolarAngleUpperLimit] (expressed in
 *  degrees), together with the corresponding simulated hits and reco-to-sim
 *  links. The accepted reco and sim hits are written out as subset collections
 *  and a freshly built link collection connects them.
 *
 *  Gaudi-native port of the Marlin SplitCollectionByPolarAngle processor from
 *  MarlinTrkProcessors (MuonColliderSoft). Each instance handles a single
 *  tracker collection; configure one instance per collection. The simulated
 *  hits are taken from the input reco-to-sim links, so no separate simhit input
 *  collection is needed.
 *
 *  @author M. Casarsa, INFN Trieste (original Marlin processor)
 */
struct SplitCollectionByPolarAngle final
    : k4FWCore::MultiTransformer<std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::SimTrackerHitCollection,
                                            edm4hep::TrackerHitSimTrackerHitLinkCollection>(
          const edm4hep::TrackerHitPlaneCollection&, const edm4hep::TrackerHitSimTrackerHitLinkCollection&)> {
  SplitCollectionByPolarAngle(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;

  std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::SimTrackerHitCollection,
             edm4hep::TrackerHitSimTrackerHitLinkCollection>
  operator()(const edm4hep::TrackerHitPlaneCollection& trackerHits,
             const edm4hep::TrackerHitSimTrackerHitLinkCollection& trackerHitLinks) const override;

private:
  Gaudi::Property<double> m_thetaMin{this, "PolarAngleLowerLimit", 50.,
                                     "Lower limit on the hit polar angle in degrees"};
  Gaudi::Property<double> m_thetaMax{this, "PolarAngleUpperLimit", 130.,
                                     "Upper limit on the hit polar angle in degrees"};
  Gaudi::Property<bool> m_fillHistos{this, "FillHistograms", false, "Flag to fill the diagnostic histograms"};

  enum { hTheta = 0, hSize };
  std::array<std::unique_ptr<Gaudi::Accumulators::StaticRootHistogram<1>>, hSize> m_histograms;
};

#endif
