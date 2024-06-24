#include "podio/Frame.h"
#include "podio/ROOTReader.h"

#include "edm4hep/CaloHitContributionCollection.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"

#include "k4FWCore/Transformer.h"
#include "k4Interface/IUniqueIDGenSvc.h"

#include "Gaudi/Property.h"
#include "Gaudi/Parsers/Factory.h"

#include <random>

struct EventHolder {
  std::vector<std::vector<std::string>>           m_fileNames;
  std::vector<std::unique_ptr<podio::ROOTReader>> m_rootFileReaders;
  std::vector<int>                                m_totalNumberOfEvents;
  std::map<int, podio::Frame>                     m_events;

  int m_nextEntry = 0;

  EventHolder(const std::vector<std::vector<std::string>>& fileNames) : m_fileNames(fileNames) {
    for (auto& names : m_fileNames) {
      m_rootFileReaders.emplace_back(std::make_unique<podio::ROOTReader>());
      m_rootFileReaders.back()->openFiles(names);
      m_totalNumberOfEvents.push_back(m_rootFileReaders.back()->getEntries("events"));
    }
  }
  EventHolder() = default;

  size_t size() const { return m_fileNames.size(); }
};

using retType =
    std::tuple<edm4hep::MCParticleCollection, std::map<std::string, edm4hep::SimTrackerHitCollection>,
               std::map<std::string, edm4hep::SimCalorimeterHitCollection>, edm4hep::CaloHitContributionCollection>;

struct OverlayTiming : public k4FWCore::MultiTransformer<retType(
                           const edm4hep::EventHeaderCollection& headers, const edm4hep::MCParticleCollection&,
                           const std::map<std::string, const edm4hep::SimTrackerHitCollection&>&,
                           const std::map<std::string, const edm4hep::SimCalorimeterHitCollection&>&,
                           const edm4hep::CaloHitContributionCollection&)> {
  OverlayTiming(const std::string& name, ISvcLocator* svcLoc)
      : MultiTransformer(
            name, svcLoc,
            {KeyValues("EventHeader", {"EventHeader"}), KeyValues("MCParticles", {"DefaultMCParticles"}),
             KeyValues("SimTrackerHits", {"DefaultSimTrackerHits"}),
             KeyValues("SimCalorimeterHits", {"DefaultSimCalorimeterHits"}),
             KeyValues("CaloHitContributions", {"CaloHitContributions"})},
            {KeyValues("OutputMCParticles", {"NewMCParticles"}), KeyValues("OutputSimTrackerHits", {"MCParticles1"}),
             KeyValues("OutputSimCalorimeterHits", {"MCParticles2"}),
             KeyValues("OutputCaloHitContributions", {"OverlayCaloHitContributions"})}) {}

  template <typename T> void overlayCollection(std::string collName, const podio::CollectionBase& inColl);

  virtual StatusCode initialize() final;
  /**  Execute.
   *   @return status code
   */
  virtual StatusCode finalize() final;

  retType virtual operator()(
      const edm4hep::EventHeaderCollection& headers, const edm4hep::MCParticleCollection& mcParticles,
      const std::map<std::string, const edm4hep::SimTrackerHitCollection&>&     simTrackerHits,
      const std::map<std::string, const edm4hep::SimCalorimeterHitCollection&>& simCalorimeterHits,
      const edm4hep::CaloHitContributionCollection&                             caloHitContribs) const final;

  std::pair<float, float> define_time_windows(const std::string& Collection_name) const;

private:

  Gaudi::Property<std::vector<std::string>> m_SimTrackerHitNames{
      this, "SimTrackerHitNames", {}, "List of names of the SimTrackerHit outputs"};
  Gaudi::Property<std::vector<std::string>> m_SimCalorimeterHitNames{
      this, "SimCalorimeterHitNames", {}, "List of names of the SimCalorimeterHit outputs"};

  Gaudi::Property<bool>        _randomBX{this, "RandomBx", false,
                                  "Place the physics event at an random position in the train: overrides PhysicsBX"};
  mutable Gaudi::Property<int> _BX_phys{this, "PhysicsBX", 1, "Number of the Bunch crossing of the physics event"};
  Gaudi::Property<float>       _tpcVdrift_mm_ns{this, "TPCDriftvelocity", float(5.0e-2),
                                          "Drift velocity in TPC [mm/ns] - default 5.0e-2 (5cm/us)"};
  Gaudi::Property<int>         _nBunchTrain{this, "NBunchtrain", 1, "Number of bunches in a bunch train"};
  Gaudi::Property<int> m_startWithBackgroundFile{this, "StartBackgroundFileIndex", -1,
                                                         "Which background file to startWith"};
  Gaudi::Property<int>         m_startWithBackgroundEvent{this, "StartBackgroundEventIndex", -1,
                                                  "Which background event to startWith"};

  Gaudi::Property<std::vector<std::vector<std::string>>> m_inputFileNames{
      this, "BackgroundFileNames", {}, "Name of the edm4hep input file(s) with background."};

  Gaudi::Property<std::vector<double>> m_Noverlay{
      this, "NumberBackground", {}, "Number of Background events to overlay - either fixed or Poisson mean"};

  Gaudi::Property<std::vector<bool>> m_Poisson{
      this,
      "Poisson_random_NOverlay",
      {},
      "Draw random number of Events to overlay from Poisson distribution with mean value NumberBackground"};

  Gaudi::Property<std::string> _mcParticleCollectionName{
      this, "MCParticleCollectionName", "MCParticle", "The name of the MCParticle collection in the background files"};

  Gaudi::Property<float> _DefaultStart_int{
      this, "Start_Integration_Time", float(-0.25),
      "Starting integration time.  Should be shortly before the BX, but may need to be shifted "
      "earlier if the vertex is smeared in time."};

  Gaudi::Property<float> _T_diff{this, "Delta_t", float(0.5), "Time difference between BXs in the BXtrain"};

  mutable std::unique_ptr<EventHolder> m_bkgEvents{};

  Gaudi::Property<std::map<std::string, std::vector<float>>> m_timeWindows{
    this, "TimeWindows", std::map<std::string, std::vector<float>>(), "Time windows for the different collections"};
private:
  // Collections with Integration Times

  Gaudi::Property<float> _BeamCal_int{
      this, "BeamCalCollection_IntegrationTime", {float(10)}, "Integration time for the BeamCalCollection"};
  Gaudi::Property<float> _LumiCal_int{
      this, "LumiCalCollection_IntegrationTime", {float(10)}, "Integration time for the LumiCalCollection"};
  Gaudi::Property<float> _EcalBarrel_int{
      this, "EcalBarrelCollection_IntegrationTime", {float(10)}, "Integration time for the EcalBarrelCollection"};
  Gaudi::Property<float> _EcalEndcap_int{
      this, "EcalEndcapCollection_IntegrationTime", {float(10)}, "Integration time for the EcalEndcapCollection"};
  Gaudi::Property<float> _HcalBarrelReg_int{this,
                                            "HcalBarrelRegCollection_IntegrationTime",
                                            {float(10)},
                                            "Integration time for the HcalBarrelRegCollection / HCalBarrelCollection"};
  Gaudi::Property<float> _HcalBarrel_int{
      this, "HcalBarrelCollection_IntegrationTime", {float(10)}, "Integration time for the HcalBarrelCollection"};
  Gaudi::Property<float> _HcalEndCaps_int{
      this, "HcalEndCapsCollection_IntegrationTime", {float(10)}, "Integration time for the HcalEndCapsCollection"};
  Gaudi::Property<float> _HcalEndCapRings_int{this,
                                              "HcalEndCapRingsCollection_IntegrationTime",
                                              {float(10)},
                                              "Integration time for the HcalEndCapRingsCollection"};
  Gaudi::Property<float> _MuonBarrel_int{
      this, "MuonBarrelCollection_Integration_Time", {float(10)}, "Integration time for the MuonBarrelCollection"};
  Gaudi::Property<float> _MuonEndCap_int{this,
                                         "MuonEndCapCollection_Integration_Time",
                                         {float(10)},
                                         "Integration time for the MuonEndCapCollection / YokeEndcapCollection"};

  // ILD specific
  Gaudi::Property<float> _LHcal_int{
      this, "LHcalCollection_Integration_Time", {float(10)}, "Integration time for the LHcalCollection"};

  Gaudi::Property<float> _EcalEndcapRing_int{this,
                                             "EcalEndcapRingCollection_Integration_Time",
                                             {float(10)},
                                             "Integration time for the EcalEndcapRingCollection"};

  Gaudi::Property<float> _EcalBarrelPreShower_int{this,
                                                  "EcalBarrelPreShowerCollection_Integration_Time",
                                                  {float(10)},
                                                  "Integration time for the EcalBarrelPreShowerCollection"};

  Gaudi::Property<float> _EcalEndcapPreShower_int{this,
                                                  "EcalEndcapPreShowerCollection_Integration_Time",
                                                  {float(10)},
                                                  "Integration time for the EcalEndcapPreShowerCollection"};

  Gaudi::Property<float> _EcalEndcapRingPreShower_int{this,
                                                      "EcalEndcapRingPreShowerCollection_Integration_Time",
                                                      {float(10)},
                                                      "Integration time for the EcalEndcapRingPreShowerCollection"};

  Gaudi::Property<float> _ETD_int{
      this, "ETDCollection_Integration_Time", {float(10)}, "Integration time for the ETDCollection"};
  Gaudi::Property<float> _FTD_int{
      this, "FTDCollection_Integration_Time", {float(10)}, "Integration time for the FTDCollection"};
  Gaudi::Property<float> _SET_int{
      this, "SETCollection_Integration_Time", {float(10)}, "Integration time for the SETCollection"};
  Gaudi::Property<float> _SIT_int{
      this, "SITCollection_Integration_Time", {float(10)}, "Integration time for the SITCollection"};
  Gaudi::Property<float> _VXD_int{
      this, "VXDCollection_Integration_Time", {float(10)}, "Integration time for the VXDCollection"};
  Gaudi::Property<float> _TPC_int{
      this, "TPCCollection_Integration_Time", {float(10)}, "Integration time for the TPCCollection"};
  Gaudi::Property<float> _TPCSpacePoint_int{this,
                                            "TPCSpacePointCollection_Integration_Time",
                                            {float(10)},
                                            "Integration time for the TPCSpacePointCollection"};

  // CLIC specific

  Gaudi::Property<float> _EcalPlug_int{
      this, "EcalPlugCollection_Integration_Time", {float(10)}, "Integration time for the ECalPlugCollection"};
  Gaudi::Property<float> _VXDB_int{
      this, "VertexBarrelCollection_Integration_Time", {float(10)}, "Integration time for the VertexBarrelCollection"};
  Gaudi::Property<float> _VXDE_int{
      this, "VertexEndcapCollection_Integration_Time", {float(10)}, "Integration time for the VertexEndcapCollection"};
  Gaudi::Property<float> _ITB_int{this,
                                  "InnerTrackerBarrelCollection_Integration_Time",
                                  {float(10)},
                                  "Integration time for the InnerTrackerBarrelCollection"};
  Gaudi::Property<float> _ITE_int{this,
                                  "InnerTrackerEndcapCollection_Integration_Time",
                                  {float(10)},
                                  "Integration time for the InnerTrackerEndcapCollection"};
  Gaudi::Property<float> _OTB_int{this,
                                  "OuterTrackerBarrelCollection_Integration_Time",
                                  {float(10)},
                                  "Integration time for the OuterTrackerBarrelCollection"};
  Gaudi::Property<float> _OTE_int{this,
                                  "OuterTrackerEndcapCollection_Integration_Time",
                                  {float(10)},
                                  "Integration time for the OuterTrackerEndcapCollection"};
  Gaudi::Property<bool>  m_allowReusingBackgroundFiles{
      this, "AllowReusingBackgroundFiles", false, "If true the same background file can be used for the same event"};

private:
  mutable std::mt19937     m_engine;
  SmartIF<IUniqueIDGenSvc> m_uidSvc;

};
