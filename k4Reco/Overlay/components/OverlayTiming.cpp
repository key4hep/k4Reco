#include "OverlayTiming.h"
#include <GaudiKernel/MsgStream.h>

#include "podio/Frame.h"

#include "edm4hep/MutableCaloHitContribution.h"

#include <TMath.h>

#include <limits>
#include <random>

DECLARE_COMPONENT(OverlayTiming)

template <typename T> inline float time_of_flight(const T& pos) {
  // Returns the time of flight to the radius in ns
  // Assumming positions in mm, then mm/m/s = 10^-3 s = 10^6* 10^-9 s = 10^6 ns
  return std::sqrt((pos[0] * pos[0]) + (pos[1] * pos[1]) + (pos[2] * pos[2])) / TMath::C() * 1e6;
}

std::pair<float, float> OverlayTiming::define_time_windows(const std::string& collection_name) const {
  return {m_timeWindows.value().at(collection_name)[0], m_timeWindows.value().at(collection_name)[1]};
  // // Default values for collections not named below
  // // if the collection is found, they are overwritten
  // auto this_start = _DefaultStart_int;
  // auto this_stop  = std::numeric_limits<float>::max();
  // auto TPC_hits   = false;

  // const std::map<std::string, std::tuple<float, float, bool>> parameters = {
  //     // CLIC / ILD common name collections
  //     {"BeamCalCollection", {0., _BeamCal_int, false}},
  //     {"LumiCalCollection", {0., _LumiCal_int, false}},
  //     // ILD
  //     // calo
  //     {"EcalBarrelCollection", {0., _EcalBarrel_int, false}},
  //     {"EcalBarrelPreShowerCollection", {0., _EcalBarrelPreShower_int, false}},
  //     {"EcalEndcapCollection", {0., _EcalEndcap_int, false}},
  //     {"EcalEndcapPreShowerCollection", {0., _EcalEndcapPreShower_int, false}},
  //     {"EcalEndcapRingCollection", {0., _EcalEndcapRing_int, false}},
  //     {"EcalEndcapRingPreShowerCollection", {0., _EcalEndcapRingPreShower_int, false}},
  //     {"HcalBarrelRegCollection", {0., _HcalBarrelReg_int, false}},
  //     {"HcalEndCapRingsCollection", {0., _HcalEndCapRings_int, false}},
  //     {"HcalEndCapsCollection", {0., _HcalEndCaps_int, false}},
  //     {"LHcalCollection", {0., _LHcal_int, false}},
  //     // muon system
  //     {"MuonBarrelCollection", {0., _MuonBarrel_int, false}},
  //     {"MuonEndCapCollection", {0., _MuonEndCap_int, false}},
  //     // tracker
  //     {"ETDCollection", {0., _ETD_int, false}},
  //     {"FTDCollection", {0., _FTD_int, false}},
  //     {"SETCollection", {0., _SET_int, false}},
  //     {"SITCollection", {0., _SIT_int, false}},
  //     {"VXDCollection", {0., _VXD_int, false}},
  //     {"TPCCollection", {-_TPC_int / 2, _TPC_int / 2, true}},
  //     {"TPCSpacePointCollection", {-_TPCSpacePoint_int / 2, _TPCSpacePoint_int / 2, true}},
  //     // CLIC
  //     // calo
  //     {"ECalBarrelCollection", {0., _EcalBarrel_int, false}},
  //     {"ECalEndcapCollection", {0., _EcalEndcap_int, false}},
  //     {"ECalPlugCollection", {0., _EcalPlug_int, false}},
  //     {"HCalBarrelCollection", {0., _HcalBarrelReg_int, false}},
  //     {"HCalEndcapCollection", {0., _HcalEndCaps_int, false}},
  //     {"HCalRingCollection", {0., _HcalEndCapRings_int, false}},
  //     // muon system
  //     {"YokeBarrelCollection", {0., _MuonBarrel_int, false}},
  //     {"YokeEndcapCollection", {0., _MuonEndCap_int, false}},
  //     // tracker
  //     {"VertexBarrelCollection", {0., _VXDB_int, false}},
  //     {"VertexEndcapCollection", {0., _VXDE_int, false}},
  //     {"InnerTrackerBarrelCollection", {0., _ITB_int, false}},
  //     {"InnerTrackerEndcapCollection", {0., _ITE_int, false}},
  //     {"OuterTrackerBarrelCollection", {0., _OTB_int, false}},
  //     {"OuterTrackerEndcapCollection", {0., _OTE_int, false}},
  // };

  // if (parameters.find(Collection_name) != parameters.end()) {
  //   this_start = std::get<0>(parameters.at(Collection_name));
  //   this_stop  = std::get<1>(parameters.at(Collection_name));
  //   TPC_hits   = std::get<2>(parameters.at(Collection_name));
  // }
  // return std::make_tuple(this_start, this_stop, TPC_hits);
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
    if (val == 0) {
      std::string err = "No events found in the background files";
      for (auto& file : m_inputFileNames.value()) {
        err += " " + file[0];
      }
      return StatusCode::FAILURE;
    }
  }

  if (m_Noverlay.empty()) {
    info() << "Using the default number of overlay events (1) for each group, since none was specified with "
              "NumberBackground "
           << endmsg;
    m_Noverlay.value() = std::vector<double>(m_bkgEvents->size(), 1);
  }

  std::vector<bool> actualPoisson;
  if (m_Poisson.empty()) {
    info() << "Using the default overlay mode (no Poission distribution) for each group, since none was specified with "
              "Poisson_random_Noverlay"
           << endmsg;
    m_Poisson.value() = std::vector<bool>(m_bkgEvents->size(), false);
  }

  return StatusCode::SUCCESS;
}

retType OverlayTiming::operator()(const edm4hep::EventHeaderCollection&                                 headers,
                                  const edm4hep::MCParticleCollection&                                  particles,
                                  const std::map<std::string, const edm4hep::SimTrackerHitCollection&>& simTrackerHits,
                                  const std::map<std::string, const edm4hep::SimCalorimeterHitCollection&>& simCaloHits,
                                  const edm4hep::CaloHitContributionCollection& caloHitContribs) const {

  if (m_SimTrackerHitNames.size() != simTrackerHits.size()) {
    throw std::runtime_error(
        "The number of names for the output SimTrackerHit collections is not the expected one, expected: " +
        std::to_string(m_SimTrackerHitNames.size()) + " but the input has " + std::to_string(simTrackerHits.size()) +
        " collections");
  }
  if (m_SimCalorimeterHitNames.size() != simCaloHits.size()) {
    throw std::runtime_error(
        "The number of names for the output SimCalorimeterHit collections is not the expected one, expected: " +
        std::to_string(m_SimCalorimeterHitNames.size()) + " but the input has " + std::to_string(simCaloHits.size()) +
        " collections");
  }

  const auto seed = m_uidSvc->getUniqueID(headers[0].getEventNumber(), headers[0].getRunNumber(), this->name());
  m_engine.seed(seed);

  // Output collections
  auto oparticles       = edm4hep::MCParticleCollection();
  auto osimTrackerHits  = std::map<std::string, edm4hep::SimTrackerHitCollection>();
  auto osimCaloHits     = std::map<std::string, edm4hep::SimCalorimeterHitCollection>();
  auto ocaloHitContribs = edm4hep::CaloHitContributionCollection();

  // Copy MCParticles for physics event into a new collection
  for (auto&& part : particles) {
    oparticles->push_back(part.clone());
  }

  // Crop the SimTrackerHits and copy them into a new collection
  for (auto it = simTrackerHits.begin(); it != simTrackerHits.end(); ++it) {
    auto& [name, coll]                    = *it;
    auto [this_start, this_stop] = define_time_windows(name);
    auto ocoll                            = edm4hep::SimTrackerHitCollection();
    for (auto elem : coll) {
      const float tof = time_of_flight(elem.getPosition());
      if ((elem.getTime() > this_start + tof) && (elem.getTime() < this_stop + tof)) {
        auto nhit = elem.clone();
        ocoll->push_back(nhit);
      }
    }
    osimTrackerHits[m_SimTrackerHitNames[std::distance(simTrackerHits.begin(), it)]] = std::move(ocoll);
  }

  // Crop the SimCalorimeterHits and copy them into a new collection,
  // together with the contributions
  std::map<std::string, std::map<uint64_t, edm4hep::MutableSimCalorimeterHit>> cellIDsMap;
  for (auto& [name, coll] : simCaloHits) {
    auto [this_start, this_stop] = define_time_windows(name);
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
        ocaloHitContribs.push_back(newcontrib);
      }
      if (within_time_window) {
        auto newhit = elem.clone();
        mapp.emplace(elem.getCellID(), newhit);
      }
    }
  }



  // Iterate over each group of files and parameters
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

      if (m_Poisson[groupIndex]) {
        NOverlay_to_this_BX = std::poisson_distribution<>(m_Noverlay[groupIndex])(m_engine);
      } else {
        NOverlay_to_this_BX = m_Noverlay[groupIndex];
      }

      debug() << "Will overlay " << NOverlay_to_this_BX << " events to BX number " << BX_number_in_train + _BX_phys
              << endmsg;

      for (int k = 0; k < NOverlay_to_this_BX; ++k) {
        info() << "Overlaying event " << m_bkgEvents->m_nextEntry << " to BX number " << BX_number_in_train + _BX_phys
               << endmsg;
        auto ptr =
            m_bkgEvents->m_rootFileReaders[groupIndex]->readEntry(podio::Category::Event, m_bkgEvents->m_nextEntry);
        if (!ptr) {
          warning() << "No more events in background file " << endmsg;
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
          warning() << "Collection " << _mcParticleCollectionName << " not found in background event"
                                      << endmsg;
        }

        for (auto it = simTrackerHits.begin(); it != simTrackerHits.end(); ++it) {
          auto& [name, coll] = *it;
          info() << "Processing collection " << name << endmsg;
          if (std::find(availableCollections.begin(), availableCollections.end(), name) == availableCollections.end()) {
            warning() << "Collection " << name << " not found in background event" << endmsg;
            continue;
          }
          auto [this_start, this_stop] = define_time_windows(name);
          // There are only contributions to the readout if the hits are in the integration window
          if (this_stop <= (BX_number_in_train - _BX_phys) * _T_diff) {
            info() << "Skipping collection " << name << " as it is not in the integration window" << endmsg;
            continue;
          }
          auto& ocoll = osimTrackerHits[m_SimTrackerHitNames[std::distance(simTrackerHits.begin(), it)]];
          if (std::abs(timeOffset) < std::numeric_limits<float>::epsilon()) {
            for (auto elem : backgroundEvent.get<edm4hep::SimTrackerHitCollection>(name)) {
              const float tof = time_of_flight(elem.getPosition());

              if (!((elem.getTime() + timeOffset > this_start + tof) &&
                    (elem.getTime() + timeOffset < this_stop + tof)))
                continue;
              auto nhit = elem.clone();
              nhit.setOverlay(true);
              ocoll->push_back(nhit);
            }
          }
          // else if (TPCHits) {
          //   for (auto elem : backgroundEvent.get<edm4hep::SimTrackerHitCollection>(name)) {
          //     const float tof = time_of_flight(elem.getPosition());
          //     if (!((elem.getTime() + timeOffset > this_start + tof) &&
          //           (elem.getTime() + timeOffset < this_stop + tof)))
          //       continue;
          //     auto nhit     = elem.clone();
          //     auto position = elem.getPosition();
          //     if (position.z <= 0.) {
          //       position.z -= timeOffset * _tpcVdrift_mm_ns;
          //     } else {
          //       position.z += timeOffset * _tpcVdrift_mm_ns;
          //     }
          //     nhit.setPosition(position);
          //     nhit.setTime(nhit.getTime() + timeOffset);
          //     nhit.setOverlay(true);
          //     ocoll->push_back(nhit);
          //   }
          // }
        }

        for (auto it = simCaloHits.begin(); it != simCaloHits.end(); ++it) {
          auto& [name, coll] = *it;
          info() << "Processing collection " << name << endmsg;
          if (std::find(availableCollections.begin(), availableCollections.end(), name) == availableCollections.end()) {
            warning() << "Collection " << name << " not found in background event" << endmsg;
            continue;
          }
          auto [this_start, this_stop] = define_time_windows(name);
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
                }
              }
            }
          }
          for (auto& [cellID, hit] : mapp) {
            ocoll->push_back(hit.clone());
          }
          osimCaloHits[m_SimCalorimeterHitNames[std::distance(simCaloHits.begin(), it)]] = std::move(ocoll);
        }
      }
    }
  }
  return std::make_tuple(std::move(oparticles), std::move(osimTrackerHits), std::move(osimCaloHits),
                         std::move(ocaloHitContribs));
}

StatusCode OverlayTiming::finalize() { return Gaudi::Algorithm::finalize(); }
