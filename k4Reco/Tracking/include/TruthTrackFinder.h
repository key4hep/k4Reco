#ifndef K4RECO_TRUTHTRACKFINDER_H
#define K4RECO_TRUTHTRACKFINDER_H 1

#include "GaudiDDKalTest.h"

#include <DD4hep/BitFieldCoder.h>

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

  std::tuple<edm4hep::TrackCollection, edm4hep::TrackMCParticleLinkCollection> operator()(
      const std::vector<const edm4hep::TrackerHitPlaneCollection*>&             trackerHitCollections,
      const std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>& relations,
      const std::vector<const edm4hep::MCParticleCollection*>&                  particleCollections) const override;

  std::vector<const edm4hep::TrackerHitPlane*> removeHitsSameLayer(
      const std::vector<const edm4hep::TrackerHitPlane*>& trackHits) const;

  // Track fit parameters
  float  m_initialTrackError_d0{1.e6};
  float  m_initialTrackError_phi0{1.e2};
  float  m_initialTrackError_omega{static_cast<float>(1.e-4)};
  float  m_initialTrackError_z0{1.e6};
  float  m_initialTrackError_tanL{1.e2};
  double m_maxChi2perHit{1.e2};
  float  m_magneticField{};

  Gaudi::Property<bool> m_useTruthInPrefit{
      this, "UseTruthInPrefit", false,
      "If true use the truth information to initialise the helical prefit, otherwise use prefit by fitting 3 hits"};
  Gaudi::Property<bool>        m_fitForward{this, "FitForward", false,
                                     "If true fit 'forward' (go forward + smooth back adding last two hits with Kalman "
                                            "FIlter steps), otherwise fit backward "};
  Gaudi::Property<std::string> m_geoSvcName{this, "GeoSvcName", "GeoSvc", "The name of the GeoSvc instance"};
  Gaudi::Property<std::string> m_encodingStringVariable{
      this, "EncodingStringParameterName", "GlobalTrackerReadoutID",
      "The name of the DD4hep constant that contains the Encoding string for tracking detectors"};
  GaudiDDKalTest                                         m_ddkaltest{this};
  SmartIF<IGeoSvc>                                       m_geoSvc;
  std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder> m_encoder;
};

#endif

DECLARE_COMPONENT(TruthTrackFinder)
