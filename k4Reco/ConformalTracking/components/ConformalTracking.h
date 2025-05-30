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
#include "Cell.h"
#include "KDCluster.h"
#include "KDTrack.h"
#include "KDTree.h"
#include "Parameters.h"

#include "GaudiDDKalTest.h"

#include <DDSegmentation/BitFieldCoder.h>

#include <edm4hep/MCParticleCollection.h>
#include <edm4hep/TrackCollection.h>
#include <edm4hep/TrackerHitPlaneCollection.h>
#include <edm4hep/TrackerHitSimTrackerHitLinkCollection.h>

#include <k4FWCore/Transformer.h>
#include <k4Interface/IGeoSvc.h>

#include <Gaudi/Accumulators/RootHistogram.h>
#include <Gaudi/Property.h>

#include "GAUDI_VERSION.h"

#if GAUDI_MAJOR_VERSION < 39
namespace Gaudi::Accumulators {
template <unsigned int ND, atomicity Atomicity = atomicity::full, typename Arithmetic = double>
using StaticRootHistogram =
    Gaudi::Accumulators::RootHistogramingCounterBase<ND, Atomicity, Arithmetic, naming::histogramString>;
}
#endif

#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>

#include <map>
#include <string>
#include <vector>

class IGeoSvc;

struct ConformalTracking final : k4FWCore::Transformer<edm4hep::TrackCollection(
                                     const std::vector<const edm4hep::TrackerHitPlaneCollection*>&,
                                     const std::vector<const edm4hep::MCParticleCollection*>&,
                                     const std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>&)> {
  ConformalTracking(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;
  StatusCode finalize() override;

  edm4hep::TrackCollection
  operator()(const std::vector<const edm4hep::TrackerHitPlaneCollection*>&,
             const std::vector<const edm4hep::MCParticleCollection*>&,
             const std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>&) const override;

private:
  // Cell creation
  SKDCluster extrapolateCell(Cell::SCell const&, double) const;
  void extendSeedCells(SharedCells&, UKDTree&, bool, const SharedKDClusters&, Parameters const&,
                       bool vertexToTracker = true) const;

  void runStep(SharedKDClusters&, UKDTree&, UniqueKDTracks&, std::map<size_t, SharedKDClusters> const&,
               Parameters const&) const;

  // Track finding
  void buildNewTracks(UniqueKDTracks&, SharedKDClusters&, UKDTree&, Parameters const&, bool radialSearch = false,
                      bool vertexToTracker = true) const;
  bool neighbourIsCompatible(const SKDCluster& neighbourHit, const SKDCluster& seedHit, const double slopeZRange) const;
  void extendTracks(UniqueKDTracks&, SharedKDClusters&, UKDTree&, Parameters const&) const;
  UKDTree combineCollections(SharedKDClusters&, std::vector<int> const&,
                             std::map<size_t, SharedKDClusters> const&) const;

  void extendHighPT(UniqueKDTracks&, SharedKDClusters&, UKDTree&, Parameters const&, bool radialSearch = false);

  void extendTracksPerLayer(UniqueKDTracks&, SharedKDClusters&, UKDTree&, Parameters const&,
                            bool vertexToTracker = true) const;

  UniqueCellularTracks createTracksNew(const Cell::SCell&) const;
  bool toBeUpdated(UniqueCellularTracks const&) const;
  void updateCell(Cell::SCell const&) const;

  // Track fitting
  UniqueKDTracks getFittedTracks(UniqueCellularTracks&, const Parameters&) const;
  void getLowestChi2(UniqueKDTracks&, UniqueKDTracks&) const;

  double fitWithoutPoint(KDTrack, int) const;
  int overlappingHits(const UKDTrack&, const UKDTrack&) const;

  void fitWithPoint(KDTrack kdTrack, SKDCluster& hit, double& deltaChi2, double& deltaChi2zs) const;

  // Used for debugging
  void drawline(const SKDCluster&, const SKDCluster&, int, int style = 1) const;

  double getBzAtOrigin() const;

  Gaudi::Property<double> m_maxCellAngle{this, "MaxCellAngle", 0.035,
                                         "Cut on angle between two cells for cell to be valid"};
  Gaudi::Property<double> m_maxCellAngleRZ{this, "MaxCellAngleRZ", 0.035,
                                           "Cut on angle between two cells in RZ for cell to be valid"};
  Gaudi::Property<double> m_maxDistance{this, "MaxDistance", 0.015,
                                        "Maximum length of a cell (max. distance between two hits)"};
  Gaudi::Property<double> m_highPTcut{this, "HighPTCut", 10.0,
                                      "pT threshold (in GeV) for enabling extendHighPT in extendTracks"};
  Gaudi::Property<double> m_chi2cut{this, "MaxChi2", 300., "Maximum chi2/ndof for linear conformal tracks"};
  Gaudi::Property<int> m_minClustersOnTrack{this, "MinClustersOnTrack", 6,
                                            "Minimum number of clusters to create a track in pattern recognition"};
  Gaudi::Property<bool> m_enableTCVC{this, "EnableTightCutsVertexCombined", true,
                                     "Enabled tight cuts as first step of reconstruction in vertex b+e [TMP!!]"};

  Gaudi::Property<bool> m_debugPlots{this, "DebugPlots", false, "Plots for debugging the tracking"};
  Gaudi::Property<bool> m_debugTime{this, "DebugTiming", false, "Print out time profile"};

  Gaudi::Property<bool> m_retryTooManyTracks{this, "RetryTooManyTracks", true,
                                             "retry with tightened parameters, when too many tracks are being created"};
  Gaudi::Property<size_t> m_tooManyTracks{this, "TooManyTracks", 1000, "Number of tracks that is considered too many"};
  Gaudi::Property<bool> m_sortTreeResults{this, "SortTreeResults", true, "sort results from kdtree search"};
  Gaudi::Property<double> m_purity{this, "trackPurity", 0.75,
                                   "Purity value used for checking if tracks are real or not"};
  Gaudi::Property<double> m_thetaRange{this, "ThetaRange", 0.1, "Angular range for initial cell seeding"};
  Gaudi::Property<int> m_minClustersOnTrackAfterFit{this, "MinClustersOnTrackAfterFit", 4,
                                                    "Final minimum number of track clusters"};
  Gaudi::Property<int> m_maxHitsInvFit{this, "MaxHitInvertedFit", 0,
                                       "Maximum number of track hits to try the inverted fit"};

  Gaudi::Property<std::string> m_encodingStringVariable{
      this, "EncodingStringParameterName", "GlobalTrackerReadoutID",
      "The name of the DD4hep constant that contains the Encoding string for tracking detectors"};

  Gaudi::Property<std::string> m_geoSvcName{this, "GeoSvcName", "GeoSvc", "The name of the GeoSvc instance"};

  Gaudi::Property<std::vector<std::vector<std::string>>> m_stepCollections{
      this, "stepCollections", {}, "Integer to add to the dummy values written to the edm"};
  Gaudi::Property<std::vector<std::vector<std::string>>> m_stepParametersNames{
      this, "stepParametersNames", {}, "Integer to add to the dummy values written to the edm"};
  Gaudi::Property<std::vector<std::vector<double>>> m_stepParametersValues{
      this, "stepParametersValues", {}, "Integer to add to the dummy values written to the edm"};
  Gaudi::Property<std::vector<std::vector<std::string>>> m_stepParametersFlags{
      this, "stepParametersFlags", {}, "Integer to add to the dummy values written to the edm"};
  Gaudi::Property<std::vector<std::vector<std::string>>> m_stepParametersFunctions{
      this, "stepParametersFunctions", {}, "Integer to add to the dummy values written to the edm"};

  Gaudi::Property<std::vector<std::string>> m_inputMainTrackerHitCollections{
      this, "MainTrackerHitCollectionNames", {}, "Name of the TrackerHit input collections from the Main Tracker"};
  Gaudi::Property<std::vector<std::string>> m_inputVertexBarrelCollections{
      this, "VertexBarrelHitCollectionNames", {}, "Name of the TrackerHit input collections from the Vertex Barrel"};
  Gaudi::Property<std::vector<std::string>> m_inputVertexEndcapCollections{
      this, "VertexEndcapHitCollectionNames", {}, "Name of the TrackerHit input collections from the Vertex Endcap"};

  SmartIF<IGeoSvc> m_geoSvc;

  // Track fit parameters
  mutable double m_initialTrackError_d0 = 0.0;
  mutable double m_initialTrackError_phi0 = 0.0;
  mutable double m_initialTrackError_omega = 0.0;
  mutable double m_initialTrackError_z0 = 0.0;
  mutable double m_initialTrackError_tanL = 0.0;
  mutable double m_maxChi2perHit = 0.0;
  mutable double m_magneticField = 0.0;

  // Histograms
  mutable Gaudi::Accumulators::StaticRootHistogram<1> m_X{this, "m_X", "m_X", {500, -1500, 1500}};
  mutable Gaudi::Accumulators::StaticRootHistogram<1> m_Y{this, "m_Y", "m_Y", {500, -1500, 1500}};
  mutable Gaudi::Accumulators::StaticRootHistogram<1> m_Z{this, "m_Z", "m_Z", {500, -2500, 2500}};

  TH1F* m_neighX = nullptr;
  TH1F* m_neighY = nullptr;
  TH1F* m_neighZ = nullptr;
  TH1F* m_slopeZ = nullptr;
  TH1F* m_slopeZ_true = nullptr;
  TH1F* m_slopeZ_true_first = nullptr;
  TH2F* m_slopeZ_vs_pt_true = nullptr;

  TH1F* m_cellAngle = nullptr;
  TH1F* m_cellDOCA = nullptr;
  TH2F* m_cellAngleRadius = nullptr;
  TH2F* m_cellLengthRadius = nullptr;
  TH2F* m_cellAngleLength = nullptr;
  TH1F* m_conformalChi2 = nullptr;
  TH1F* m_conformalChi2real = nullptr;
  TH1F* m_conformalChi2fake = nullptr;
  TH2F* m_conformalChi2Purity = nullptr;

  TH1F* m_absZ = nullptr;

  TH1F* m_conformalChi2MC = nullptr;
  TH2F* m_conformalChi2PtMC = nullptr;
  TH2F* m_conformalChi2VertexRMC = nullptr;

  TH1F* m_conformalChi2SzMC = nullptr;
  TH2F* m_conformalChi2SzPtMC = nullptr;
  TH2F* m_conformalChi2SzVertexRMC = nullptr;

  TH1F* m_cellAngleMC = nullptr;
  TH1F* m_cellDOCAMC = nullptr;
  TH1F* m_cellAngleRZMC = nullptr;
  TH2F* m_cellAngleRadiusMC = nullptr;
  TH2F* m_cellLengthRadiusMC = nullptr;
  TH2F* m_cellAngleLengthMC = nullptr;

  TH2F* m_conformalEvents = nullptr;
  TH2F* m_nonconformalEvents = nullptr;
  TH2F* m_conformalEventsRTheta = nullptr;
  TH2F* m_conformalEventsMC = nullptr;

  TCanvas* m_canvConformalEventDisplay = nullptr;
  TCanvas* m_canvConformalEventDisplayAllCells = nullptr;
  TCanvas* m_canvConformalEventDisplayAcceptedCells = nullptr;
  TCanvas* m_canvConformalEventDisplayMC = nullptr;
  TCanvas* m_canvConformalEventDisplayMCunreconstructed = nullptr;

  TH2F* m_szDistribution = nullptr;
  TH2F* m_uvDistribution = nullptr;
  TH2F* m_xyDistribution = nullptr;
  TH3F* m_xyzDistribution = nullptr;

  // Other constants
  SKDCluster debugSeed = nullptr;
  double m_slopeZRange = 1000.0;
  // ConformalDebugger                    m_debugger{};
  std::map<SKDCluster, edm4hep::MCParticle*> kdParticles{};  // Link from conformal hit to MC particle
  std::map<SKDCluster, edm4hep::SimTrackerHit*> kdSimHits{}; // Link from conformal hit to SimHit

  std::vector<Parameters> m_stepParameters;

  GaudiDDKalTest m_ddkaltest{this};
  mutable int m_eventNumber = 0;

  dd4hep::DDSegmentation::BitFieldCoder m_encoder;
};

DECLARE_COMPONENT(ConformalTracking)
