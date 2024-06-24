#include "OverlayTiming.h"
#include <GaudiKernel/MsgStream.h>

#include "podio/Frame.h"

#include "edm4hep/MutableCaloHitContribution.h"

#include <limits>
#include <random>

DECLARE_COMPONENT(OverlayTiming)

template <typename T> inline float time_of_flight(const T& pos) {
  // Returns the time of flight to the radius in ns
  // mm/m/s = 10^{-3}s = 10^6 ns d.h. 299 mm/ns
  return std::sqrt((pos[0] * pos[0]) + (pos[1] * pos[1]) + (pos[2] * pos[2])) / 299.792458;
}

std::tuple<float, float, bool> OverlayTiming::define_time_windows(const std::string& Collection_name) const {
  // Default values for collections not named below
  // if the collection is found, they are overwritten
  auto this_start = _DefaultStart_int;
  auto this_stop  = std::numeric_limits<float>::max();
  auto TPC_hits   = false;

  const std::map<std::string, std::tuple<float, float, bool>> parameters = {
      // CLIC / ILD common name collections
      {"BeamCalCollection", {0., _BeamCal_int, false}},
      {"LumiCalCollection", {0., _LumiCal_int, false}},
      // ILD
      // calo
      {"EcalBarrelCollection", {0., _EcalBarrel_int, false}},
      {"EcalBarrelPreShowerCollection", {0., _EcalBarrelPreShower_int, false}},
      {"EcalEndcapCollection", {0., _EcalEndcap_int, false}},
      {"EcalEndcapPreShowerCollection", {0., _EcalEndcapPreShower_int, false}},
      {"EcalEndcapRingCollection", {0., _EcalEndcapRing_int, false}},
      {"EcalEndcapRingPreShowerCollection", {0., _EcalEndcapRingPreShower_int, false}},
      {"HcalBarrelRegCollection", {0., _HcalBarrelReg_int, false}},
      {"HcalEndCapRingsCollection", {0., _HcalEndCapRings_int, false}},
      {"HcalEndCapsCollection", {0., _HcalEndCaps_int, false}},
      {"LHcalCollection", {0., _LHcal_int, false}},
      // muon system
      {"MuonBarrelCollection", {0., _MuonBarrel_int, false}},
      {"MuonEndCapCollection", {0., _MuonEndCap_int, false}},
      // tracker
      {"ETDCollection", {0., _ETD_int, false}},
      {"FTDCollection", {0., _FTD_int, false}},
      {"SETCollection", {0., _SET_int, false}},
      {"SITCollection", {0., _SIT_int, false}},
      {"VXDCollection", {0., _VXD_int, false}},
      {"TPCCollection", {-_TPC_int / 2, _TPC_int / 2, true}},
      {"TPCSpacePointCollection", {-_TPCSpacePoint_int / 2, _TPCSpacePoint_int / 2, true}},
      // CLIC
      // calo
      {"ECalBarrelCollection", {0., _EcalBarrel_int, false}},
      {"ECalEndcapCollection", {0., _EcalEndcap_int, false}},
      {"ECalPlugCollection", {0., _EcalPlug_int, false}},
      {"HCalBarrelCollection", {0., _HcalBarrelReg_int, false}},
      {"HCalEndcapCollection", {0., _HcalEndCaps_int, false}},
      {"HCalRingCollection", {0., _HcalEndCapRings_int, false}},
      // muon system
      {"YokeBarrelCollection", {0., _MuonBarrel_int, false}},
      {"YokeEndcapCollection", {0., _MuonEndCap_int, false}},
      // tracker
      {"VertexBarrelCollection", {0., _VXDB_int, false}},
      {"VertexEndcapCollection", {0., _VXDE_int, false}},
      {"InnerTrackerBarrelCollection", {0., _ITB_int, false}},
      {"InnerTrackerEndcapCollection", {0., _ITE_int, false}},
      {"OuterTrackerBarrelCollection", {0., _OTB_int, false}},
      {"OuterTrackerEndcapCollection", {0., _OTE_int, false}},
  };

  if (parameters.find(Collection_name) != parameters.end()) {
    this_start = std::get<0>(parameters.at(Collection_name));
    this_stop  = std::get<1>(parameters.at(Collection_name));
    TPC_hits   = std::get<2>(parameters.at(Collection_name));
  }
  return std::make_tuple(this_start, this_stop, TPC_hits);
}

StatusCode OverlayTiming::initialize() {
  m_uidSvc = service<IUniqueIDGenSvc>("UniqueIDGenSvc", true);
  if (!m_uidSvc) {
    error() << "Unable to get UniqueIDGenSvc" << endmsg;
  }

  std::vector<std::vector<std::string>> inputFiles;
  inputFiles = m_inputFileNames.value();
  // if (m_startWithBackgroundFile >= 0) {
  //   inputFiles = std::vector<std::string>(m_inputFileNames.begin() + m_startWithBackgroundFile, m_inputFileNames.end());
  // } else {
  //   inputFiles = m_inputFileNames;
  // }
  // TODO:: shuffle input files
  // std::shuffle(inputFiles.begin(), inputFiles.end(), m_engine);

  m_bkgEvents = make_unique<EventHolder>(inputFiles);
  for (auto& val : m_bkgEvents->m_totalNumberOfEvents) {
    info() << "Total number of events in group: " << val << endmsg;
  }

  return StatusCode::SUCCESS;
}

retType OverlayTiming::operator()(const edm4hep::EventHeaderCollection&                                 headers,
                                  const edm4hep::MCParticleCollection&                                  particles,
                                  const std::map<std::string, const edm4hep::SimTrackerHitCollection&>& simTrackerHits,
                                  const std::map<std::string, const edm4hep::SimCalorimeterHitCollection&>& simCaloHits,
                                  const edm4hep::CaloHitContributionCollection& caloHitContribs) const {
  const auto seed = m_uidSvc->getUniqueID(headers[0].getEventNumber(), headers[0].getRunNumber(), this->name());
  m_engine.seed(seed);

  // Output collections
  auto oparticles       = edm4hep::MCParticleCollection();
  auto osimTrackerHits  = std::map<std::string, edm4hep::SimTrackerHitCollection>();
  auto osimCaloHits     = std::map<std::string, edm4hep::SimCalorimeterHitCollection>();
  auto ocaloHitContribs = edm4hep::CaloHitContributionCollection();

  // copy MCParticles for physics event into a new collection
  for (auto&& part : particles) {
    oparticles->push_back(part.clone());
  }

  // Crop the SimTrackerHits
  for (auto& [name, coll] : simTrackerHits) {
    auto [this_start, this_stop, TPCHits] = define_time_windows(name);
    auto ocoll                            = edm4hep::SimTrackerHitCollection();
    for (auto elem : coll) {
      const float tof = time_of_flight(elem.getPosition());
      if ((elem.getTime() > this_start + tof) && (elem.getTime() < this_stop + tof)) {
        auto nhit = elem.clone();
        ocoll->push_back(nhit);
      }
    }
    osimTrackerHits["New" + name] = std::move(ocoll);
  }

  std::map<std::string, std::map<uint64_t, edm4hep::MutableSimCalorimeterHit>> cellIDsMap;
  for (auto& [name, coll] : simCaloHits) {
    auto [this_start, this_stop, TPCHits] = define_time_windows(name);
    auto& mapp                            = cellIDsMap[name];
    for (auto&& elem : coll) {
      const float tof                = time_of_flight(elem.getPosition());
      int         within_time_window = 0;
      for (auto contrib : elem.getContributions()) {
        if (!((contrib.getTime() > this_start + tof) && (contrib.getTime() < this_stop + tof)))
          continue;
        within_time_window++;
        // TODO: Make sure a contribution is not added twice
        auto newcontrib = contrib.clone();
        // ocaloHitContribs.push_back(newcontrib);
      }
      if (within_time_window) {
        auto newhit            = elem.clone();
        mapp.emplace(elem.getCellID(), newhit);
      }
    }
  }

  // Iterate over each group of files and parameters
  info() << "m_bkgEvents->size() = " << m_bkgEvents->size() << endmsg;
  for (size_t groupIndex = 0; groupIndex < m_bkgEvents->size(); groupIndex++) {
    if (_randomBX) {
      _BX_phys = std::uniform_int_distribution<int>(0, _nBunchTrain - 1)(m_engine);
      debug() << "Physics Event was placed in the " << _BX_phys << " bunch crossing!" << endmsg;
    }

    // define a permutation for the events to overlay -- the physics event is per definition at position 0
    std::vector<int> permutation;

    // Permutation has negative values and the last one is 0
    // if (!_randomBX) then _BX_phys (default = 1)
    for (int i = -(_BX_phys - 1); i < _nBunchTrain - (_BX_phys - 1); ++i) {
      permutation.push_back(i);
    }
    std::shuffle(permutation.begin(), permutation.end(), m_engine);

    // TODO: Check that there is anything to overlay

    info() << "Starting at event: " << m_bkgEvents->m_nextEntry << endmsg;

    if (m_startWithBackgroundEvent >= 0) {
      info() << "Skipping to event: " << m_startWithBackgroundEvent << endmsg;
      m_bkgEvents->m_nextEntry = m_startWithBackgroundEvent;
    }

    // Now overlay the background evnts to each bunchcrossing in the bunch train
    for (int bxInTrain = 0; bxInTrain < _nBunchTrain; ++bxInTrain) {
      const int BX_number_in_train = permutation.at(bxInTrain);

      int NOverlay_to_this_BX = 0;

      if (_Poisson) {
        NOverlay_to_this_BX = std::poisson_distribution<>(_NOverlay)(m_engine);
      } else {
        NOverlay_to_this_BX = _NOverlay;
      }

      debug() << "Will overlay " << NOverlay_to_this_BX << " events to BX number " << BX_number_in_train + _BX_phys
              << endmsg;

      for (int k = 0; k < NOverlay_to_this_BX; ++k) {
        info() << "Overlaying event " << m_bkgEvents->m_nextEntry << " to BX number " << BX_number_in_train + _BX_phys
               << endmsg;
        auto ptr =
            m_bkgEvents->m_rootFileReaders[groupIndex]->readEntry(podio::Category::Event, m_bkgEvents->m_nextEntry);
        if (!ptr) {
          Gaudi::Algorithm::warning() << "No more events in background file " << endmsg;
          break;
        }
        auto backgroundEvent      = podio::Frame(std::move(ptr));
        auto availableCollections = backgroundEvent.getAvailableCollections();

        // Either 0 or negative
        auto timeOffset = BX_number_in_train * _T_diff;

        if (std::find(availableCollections.begin(), availableCollections.end(), _mcParticleCollectionName) !=
            availableCollections.end()) {
          for (auto part : backgroundEvent.get<edm4hep::MCParticleCollection>(_mcParticleCollectionName)) {
            auto npart = part.clone();
            npart.setTime(part.getTime() + timeOffset);
            npart.setOverlay(true);
            oparticles->push_back(npart);
          }
        } else {
          Gaudi::Algorithm::warning() << "Collection " << _mcParticleCollectionName << " not found in background event"
                                      << endmsg;
        }

        for (auto& [name, coll] : simTrackerHits) {
          info() << "Processing collection " << name << endmsg;
          if (std::find(availableCollections.begin(), availableCollections.end(), name) == availableCollections.end()) {
            Gaudi::Algorithm::warning() << "Collection " << name << " not found in background event" << endmsg;
            continue;
          }
          auto [this_start, this_stop, TPCHits] = define_time_windows(name);
          // There are only contributions to the readout if the hits are in the integration window
          if (this_stop <= (BX_number_in_train - _BX_phys) * _T_diff) {
            info() << "Skipping collection " << name << " as it is not in the integration window" << endmsg;
            continue;
          }
          auto& ocoll = osimTrackerHits["New" + name];
          if (std::abs(timeOffset) < std::numeric_limits<float>::epsilon() || !TPCHits) {
            for (auto elem : backgroundEvent.get<edm4hep::SimTrackerHitCollection>(name)) {
              const float tof = time_of_flight(elem.getPosition());

              if (!((elem.getTime() + timeOffset > this_start + tof) &&
                    (elem.getTime() + timeOffset < this_stop + tof)))
                continue;
              auto nhit = elem.clone();
              nhit.setOverlay(true);
              ocoll->push_back(nhit);
            }
          } else if (TPCHits) {
            for (auto elem : backgroundEvent.get<edm4hep::SimTrackerHitCollection>(name)) {
              const float tof = time_of_flight(elem.getPosition());
              if (!((elem.getTime() + timeOffset > this_start + tof) &&
                    (elem.getTime() + timeOffset < this_stop + tof)))
                continue;
              auto nhit     = elem.clone();
              auto position = elem.getPosition();
              if (position.z <= 0.) {
                position.z -= timeOffset * _tpcVdrift_mm_ns;
              } else {
                position.z += timeOffset * _tpcVdrift_mm_ns;
              }
              nhit.setPosition(position);
              nhit.setTime(nhit.getTime() + timeOffset);
              nhit.setOverlay(true);
              ocoll->push_back(nhit);
            }
          }
        }

        for (auto& [name, coll] : simCaloHits) {
          info() << "Processing collection " << name << endmsg;
          if (std::find(availableCollections.begin(), availableCollections.end(), name) == availableCollections.end()) {
            Gaudi::Algorithm::warning() << "Collection " << name << " not found in background event" << endmsg;
            continue;
          }
          auto [this_start, this_stop, TPCHits] = define_time_windows(name);
          // There are only contributions to the readout if the hits are in the integration window
          if (this_stop <= (BX_number_in_train - _BX_phys) * _T_diff) {
            info() << "Skipping collection " << name << " as it is not in the integration window" << endmsg;
            continue;
          }

          auto ocoll = edm4hep::SimCalorimeterHitCollection();
          auto mapp  = cellIDsMap[name];
          for (auto elem : backgroundEvent.get<edm4hep::SimCalorimeterHitCollection>(name)) {
            if (mapp.find(elem.getCellID()) == mapp.end()) {
              // There is no Hit at this position -- the new hit can be added, if it is not outside the window
              auto calhit = edm4hep::MutableSimCalorimeterHit();
              bool add    = false;
              for (auto& contrib : elem.getContributions()) {
                if ((contrib.getTime() + timeOffset > this_start) && (contrib.getTime() + timeOffset < this_stop)) {
                  add = true;
                  // TODO: Make sure a contribution is not added twice
                  auto newcontrib = contrib.clone();
                  calhit.addToContributions(newcontrib);
                  ocaloHitContribs.push_back(newcontrib);
                }
              }
              if (add) {
                calhit.setCellID(elem.getCellID());
                calhit.setEnergy(elem.getEnergy());
                calhit.setPosition(elem.getPosition());
                mapp[calhit.getCellID()] = calhit;
                info() << "Adding new hit with cellID: " << elem.getCellID() << endmsg;
                info() << calhit.getContributions().size() << endmsg;
              }
            } else {
              // there is already a hit at this position....
              // TODO: get the hit
              auto& calhit = mapp[elem.getCellID()];
              for (auto& contrib : elem.getContributions()) {
                if ((contrib.getTime() + timeOffset > this_start) && (contrib.getTime() + timeOffset < this_stop)) {
                  // TODO: Make sure a contribution is not added twice
                  auto newcontrib = contrib.clone();
                  calhit.addToContributions(newcontrib);
                  ocaloHitContribs.push_back(newcontrib);
                  // info() << "Hit found, adding contribution with cellID: " << elem.getCellID() << endmsg;
                  // info() << calhit.getContributions().size() << endmsg;
                }
              }
            }
          }
          for (auto& [cellID, hit] : mapp) {
            // why clone?
            info() << "CellID: " << cellID << ", contributions: " << hit.getContributions().size() << endmsg;
            // info() << "begin: " << hit.getContributions().begin() << endmsg;
            // info() << "end: " << hit.getContributions().end() << endmsg;
            ocoll->push_back(hit.clone());
          }
          osimCaloHits["New" + name] = std::move(ocoll);
        }
      }
    }
  }
  return std::make_tuple(std::move(oparticles), std::move(osimTrackerHits), std::move(osimCaloHits),
                         std::move(ocaloHitContribs));
}

StatusCode OverlayTiming::finalize() { return Gaudi::Algorithm::finalize(); }
