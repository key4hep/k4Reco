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
#ifndef K4RECO_FILTERCONEHITS_H
#define K4RECO_FILTERCONEHITS_H 1

#include "Gaudi/Accumulators/RootHistogram.h"
#include "Gaudi/Property.h"

#include <edm4hep/MCParticleCollection.h>
#include <edm4hep/SimTrackerHitCollection.h>
#include <edm4hep/TrackerHitPlaneCollection.h>
#include <edm4hep/TrackerHitSimTrackerHitLinkCollection.h>

#include <k4FWCore/Transformer.h>
#include <k4Interface/IGeoSvc.h>

#include <GaudiKernel/StatusCode.h>

#include <array>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "GAUDI_VERSION.h"

#if GAUDI_MAJOR_VERSION < 39
namespace Gaudi::Accumulators {
template <unsigned int ND, atomicity Atomicity = atomicity::full, typename Arithmetic = double>
using StaticRootHistogram =
    Gaudi::Accumulators::RootHistogramingCounterBase<ND, Atomicity, Arithmetic, naming::histogramString>;
}
#endif

/** === FilterConeHits ===
 *  Selects the tracker hits that lie inside a cone opened around the trajectory
 *  of a generator-level MC particle, together with the corresponding simulated
 *  hits and reco-to-sim links. For each MC particle (whose generator status is
 *  in ConeAroundStatus) a helix is built from its production vertex, momentum
 *  and charge in the detector magnetic field. A hit is kept when its angular
 *  distance to the helix is below DeltaRCut and/or its 3D distance to the helix
 *  is below Dist3DCut. The outputs are subset collections of the accepted reco
 *  and sim hits plus a freshly built link collection.
 *
 *  Gaudi-native port of the Marlin FilterConeHits processor from
 *  MarlinTrkProcessors (MuonColliderSoft). Each instance handles a single
 *  tracker subdetector; configure one instance per collection. The simulated
 *  hits are taken from the input reco-to-sim links, so no separate simhit input
 *  collection is needed.
 *
 *  @author M. Casarsa, INFN Trieste (original Marlin processor)
 */
struct FilterConeHits final
    : k4FWCore::MultiTransformer<std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::SimTrackerHitCollection,
                                            edm4hep::TrackerHitSimTrackerHitLinkCollection>(
          const edm4hep::MCParticleCollection&, const edm4hep::TrackerHitPlaneCollection&,
          const edm4hep::TrackerHitSimTrackerHitLinkCollection&)> {
  FilterConeHits(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;

  std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::SimTrackerHitCollection,
             edm4hep::TrackerHitSimTrackerHitLinkCollection>
  operator()(const edm4hep::MCParticleCollection& mcParticles, const edm4hep::TrackerHitPlaneCollection& trackerHits,
             const edm4hep::TrackerHitSimTrackerHitLinkCollection& trackerHitLinks) const override;

private:
  Gaudi::Property<double> m_deltaRCut{this, "DeltaRCut", -1.,
                                      "Maximum angular distance between the hits and the particle trajectory [rad]"};
  Gaudi::Property<double> m_dist3DCut{this, "Dist3DCut", -1.,
                                      "Maximum 3D distance between the hits and the extrapolated helix [mm]"};
  Gaudi::Property<std::vector<int>> m_coneAroundStatus{
      this, "ConeAroundStatus", {1}, "List of MC particle generator statuses to build cones around"};
  Gaudi::Property<bool> m_fillHistos{this, "FillHistograms", false, "Flag to fill the diagnostic histograms"};
  Gaudi::Property<double> m_trackerOuterRadius{this, "TrackerOuterRadius", 1500.,
                                               "Outer radius of the tracker barrel used to clip the helix [mm]"};
  Gaudi::Property<std::string> m_geoSvcName{this, "GeoSvcName", "GeoSvc", "The name of the GeoSvc instance"};

  SmartIF<IGeoSvc> m_geoSvc;
  double m_bField{0.}; ///< Bz at the origin [Tesla]

  enum { hDistXY = 0, hDistZ, hDist3D, hAngle, hTime, hPathLength, hSize };
  std::array<std::unique_ptr<Gaudi::Accumulators::StaticRootHistogram<1>>, hSize> m_histograms;
};

#endif
