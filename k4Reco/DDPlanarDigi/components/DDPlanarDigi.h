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
#ifndef DDPLANARDIGI_H
#define DDPLANARDIGI_H

#include "Gaudi/Accumulators/RootHistogram.h"
#include "Gaudi/Property.h"

#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackerHitPlaneCollection.h"
#include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h"

#include "k4FWCore/Transformer.h"
#include "k4Interface/IGeoSvc.h"
#include "k4Interface/IUniqueIDGenSvc.h"

#include "DDRec/SurfaceManager.h"

#include "TRandom2.h"

#include <string>
#include <vector>

#include "GAUDI_VERSION.h"

#if GAUDI_MAJOR_VERSION < 39
namespace Gaudi::Accumulators {
  template <unsigned int ND, atomicity Atomicity = atomicity::full, typename Arithmetic = double>
  using StaticRootHistogram =
      Gaudi::Accumulators::RootHistogramingCounterBase<ND, Atomicity, Arithmetic, naming::histogramString>;
}
#endif

/** ======= DDPlanarDigi ========== <br>
 * Creates TrackerHits from SimTrackerHits, smearing them according to the input parameters.
 * The positions of "digitized" TrackerHits are obtained by gaussian smearing positions
 * of SimTrackerHits perpendicular and along the ladder according to the specified point resolutions.
 * The geometry of the surface is retreived from DDRec::Surface associated to the hit via cellID.
 *
 *
 * <h4>Input collections and prerequisites</h4>
 * Processor requires a collection of SimTrackerHits <br>
 * <h4>Output</h4>
 * Processor produces collection of smeared TrackerHits<br>
 * @param SimTrackHitCollectionName The name of input collection of SimTrackerHits <br>
 * (default name VXDCollection) <br>
 * @param TrackerHitCollectionName The name of output collection of smeared TrackerHits <br>
 * (default name VTXTrackerHits) <br>
 * @param ResolutionU resolution in direction of u (in mm) <br>
 * (default value 0.004) <br>
 * @param ResolutionV Resolution in direction of v (in mm) <br>
 * (default value 0.004) <br>
 * @param IsStrip whether the hits are 1 dimensional strip measurements <br>
 * (default value false)
 * @param Ladder_Number_encoded_in_cellID ladder number has been encoded in the cellID <br>
 * (default value false) <br>
 * @param Sub_Detector_ID ID of Sub-Detector using UTIL/ILDConf.h from lcio <br>
 * (default value lcio::ILDDetID::VXD) <br>
 * <br>
 *
 * Originally in https://github.com/iLCSoft/MarlinTrkProcessors/blob/master/source/Digitisers/include/DDPlanarDigi.h
 */

enum { hu = 0, hv, hT, hitE, hitsAccepted, diffu, diffv, diffT, hSize };

struct DDPlanarDigi final
    : k4FWCore::MultiTransformer<
          std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::TrackerHitSimTrackerHitLinkCollection>(
              const edm4hep::SimTrackerHitCollection&, const edm4hep::EventHeaderCollection&)> {
  DDPlanarDigi(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;

  std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::TrackerHitSimTrackerHitLinkCollection> operator()(
      const edm4hep::SimTrackerHitCollection& simTrackerHits,
      const edm4hep::EventHeaderCollection&   headers) const override;

private:
  Gaudi::Property<std::string>        m_subDetName{this, "SubDetectorName", "VXD", "Name of the subdetector"};
  Gaudi::Property<bool>               m_isStrip{this, "IsStrip", false, "Whether the hits are 1D strip hits"};
  Gaudi::Property<std::vector<float>> m_resULayer{
      this, "ResolutionU", {0.004}, "Resolution in the direction of u; either one per layer or one for all layers"};
  Gaudi::Property<std::vector<float>> m_resVLayer{
      this, "ResolutionV", {0.004}, "Resolution in the direction of v; either one per layer or one for all layers"};
  Gaudi::Property<std::vector<float>> m_resTLayer{
      this,
      "ResolutionT",
      {-1},
      "Resolution in the direction of t; either one per layer or one for all layers. If the single entry is negative, "
      "disable time smearing. "};
  Gaudi::Property<double> m_efficiency{
      this, "Efficiency", 1.00, "Efficiency of the chamber to detect the charged particle"};    
  Gaudi::Property<bool>   m_forceHitsOntoSurface{this, "ForceHitsOntoSurface", false,
                                               "Project hits onto the surfoce in case they are not yet on the surface"};
  Gaudi::Property<double> m_minEnergy{this, "MinEnergy", 0.0, "Minimum energy (GeV) of SimTrackerHit to be digitized"};

  Gaudi::Property<bool> m_useTimeWindow{
      this, "UseTimeWindow", false,
      "Only accept hits with time (after smearing) within the specified time window (default: false)"};
  Gaudi::Property<bool> m_correctTimesForPropagation{
      this, "CorrectTimesForPropagation", false,
      "Correct hit time for the propagation: radial distance/c (default: false)"};
  Gaudi::Property<std::vector<float>> m_timeWindowMin{
      this, "TimeWindowMin", {-1e9}, "Minimum time (ns) of SimTrackerHit to be digitized"};
  Gaudi::Property<std::vector<float>> m_timeWindowMax{
      this, "TimeWindowMax", {1e9}, "Maximum time (ns) of SimTrackerHit to be digitized"};
  Gaudi::Property<std::string> m_encodingStringVariable{
      this, "EncodingStringParameterName", "GlobalTrackerReadoutID",
      "The name of the DD4hep constant that contains the Encoding string for tracking detectors"};
  Gaudi::Property<std::string> m_geoSvcName{this, "GeoSvcName", "GeoSvc", "The name of the GeoSvc instance"};
  Gaudi::Property<int>    m_maxTries{this, "MaxTries", 10, "Maximum number of tries to find a valid surface for a hit"};
  Gaudi::Property<size_t> m_cellIDBits{this, "CellIDBits", 64, "Number of bits to use for the cellID of the hits"};
  const dd4hep::rec::SurfaceMap*                                                  surfaceMap;
  std::array<std::unique_ptr<Gaudi::Accumulators::StaticRootHistogram<1>>, hSize> m_histograms;
  std::string                                                                     m_collName;

  inline static thread_local TRandom2 m_engine;
  SmartIF<IGeoSvc>                    m_geoSvc;
  SmartIF<IUniqueIDGenSvc>            m_uidSvc;

  std::uint64_t m_mask{static_cast<std::uint64_t>(-1)};
};

DECLARE_COMPONENT(DDPlanarDigi)

#endif
