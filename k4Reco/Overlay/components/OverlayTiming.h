#include "podio/Frame.h"
#include "podio/ROOTReader.h"
#include "podio/Reader.h"

#include "edm4hep/CaloHitContributionCollection.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"

#include "k4FWCore/Transformer.h"
#include "k4Interface/IUniqueIDGenSvc.h"

// Needed for some of the more complex properties
#include "Gaudi/Parsers/Factory.h"
#include "Gaudi/Property.h"

#include <random>

struct EventHolder {
  std::vector<std::vector<std::string>>           m_fileNames;
  std::vector<podio::Reader>                      m_rootFileReaders;
  std::vector<int>                                m_totalNumberOfEvents;
  std::map<int, podio::Frame>                     m_events;

  int m_nextEntry = 0;

  EventHolder(const std::vector<std::vector<std::string>>& fileNames) : m_fileNames(fileNames) {
    for (auto& names : m_fileNames) {
      m_rootFileReaders.emplace_back(podio::makeReader(names));
      m_totalNumberOfEvents.push_back(m_rootFileReaders.back().getEntries("events"));
    }
  }
  EventHolder() = default;

  size_t size() const { return m_fileNames.size(); }
};

using retType = std::tuple<edm4hep::MCParticleCollection, std::vector<edm4hep::SimTrackerHitCollection>,
                           std::vector<edm4hep::SimCalorimeterHitCollection>,
                           std::vector<edm4hep::CaloHitContributionCollection>>;

struct OverlayTiming : public k4FWCore::MultiTransformer<retType(
                           const edm4hep::EventHeaderCollection& headers, const edm4hep::MCParticleCollection&,
                           const std::vector<const edm4hep::SimTrackerHitCollection*>&,
                           const std::vector<const edm4hep::SimCalorimeterHitCollection*>&
                           // const std::map<std::string, const edm4hep::CaloHitContributionCollection&>&
                                                                 )> {
  OverlayTiming(const std::string& name, ISvcLocator* svcLoc)
      : MultiTransformer(
            name, svcLoc,
            {KeyValues("EventHeader", {"EventHeader"}), KeyValues("MCParticles", {"DefaultMCParticles"}),
             KeyValues("SimTrackerHits", {"DefaultSimTrackerHits"}),
             KeyValues("SimCalorimeterHits", {"DefaultSimCalorimeterHits"})
            },
            {KeyValues("OutputMCParticles", {"NewMCParticles"}), KeyValues("OutputSimTrackerHits", {"MCParticles1"}),
             KeyValues("OutputSimCalorimeterHits", {"MCParticles2"}),
             KeyValues("OutputCaloHitContributions", {"OverlayCaloHitContributions"})}) {}

  template <typename T> void overlayCollection(std::string collName, const podio::CollectionBase& inColl);

  virtual StatusCode initialize() final;

  retType virtual operator()(
      const edm4hep::EventHeaderCollection& headers, const edm4hep::MCParticleCollection& mcParticles,
      const std::vector<const edm4hep::SimTrackerHitCollection*>&       simTrackerHits,
      const std::vector<const edm4hep::SimCalorimeterHitCollection*>&   simCalorimeterHits
                             ) const final;

  std::pair<float, float> define_time_windows(const std::string& Collection_name) const;

private:

  // These correspond to the index position in the argument list
  constexpr static int TRACKERHIT_INDEX_POSITION = 2;
  constexpr static int SIMCALOHIT_INDEX_POSITION = 3;

  Gaudi::Property<bool>        _randomBX{this, "RandomBx", false,
                                  "Place the physics event at an random position in the train: overrides PhysicsBX"};
  mutable Gaudi::Property<int> _BX_phys{this, "PhysicsBX", 1, "Number of the Bunch crossing of the physics event"};
  Gaudi::Property<float>       _tpcVdrift_mm_ns{this, "TPCDriftvelocity", float(5.0e-2),
                                          "Drift velocity in TPC [mm/ns] - default 5.0e-2 (5cm/us)"};
  Gaudi::Property<int>         _nBunchTrain{this, "NBunchtrain", 1, "Number of bunches in a bunch train"};
  Gaudi::Property<int>         m_startWithBackgroundFile{this, "StartBackgroundFileIndex", -1,
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
      this, "BackgroundMCParticleCollectionName", "MCParticle",
      "The name of the MCParticle collection in the background files"};

  Gaudi::Property<float> _DefaultStart_int{
      this, "Start_Integration_Time", float(-0.25),
      "Starting integration time.  Should be shortly before the BX, but may need to be shifted "
      "earlier if the vertex is smeared in time."};

  Gaudi::Property<float> _T_diff{this, "Delta_t", float(0.5), "Time difference between BXs in the BXtrain"};

  mutable std::unique_ptr<EventHolder> m_bkgEvents{};

  Gaudi::Property<std::map<std::string, std::vector<float>>> m_timeWindows{
      this, "TimeWindows", std::map<std::string, std::vector<float>>(), "Time windows for the different collections"};
  Gaudi::Property<bool> m_allowReusingBackgroundFiles{
      this, "AllowReusingBackgroundFiles", false, "If true the same background file can be used for the same event"};

private:
  mutable std::mt19937     m_engine;
  SmartIF<IUniqueIDGenSvc> m_uidSvc;
};
