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
  // Assumming positions in mm, then mm/m/s = 10^-3 s = 10^6 ns
  return std::sqrt((pos[0] * pos[0]) + (pos[1] * pos[1]) + (pos[2] * pos[2])) / TMath::C() * 1e6;
}

std::pair<float, float> OverlayTiming::define_time_windows(const std::string& collection_name) const {
  try {
    return {m_timeWindows.value().at(collection_name)[0], m_timeWindows.value().at(collection_name)[1]};
  } catch (const std::out_of_range& e) {
    error() << "No time window defined for collection " << collection_name << endmsg;
    throw e;
  }
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

retType OverlayTiming::operator()(
    const edm4hep::EventHeaderCollection& headers, const edm4hep::MCParticleCollection& particles,
    const std::map<std::string, const edm4hep::SimTrackerHitCollection&>&       simTrackerHits,
    const std::map<std::string, const edm4hep::SimCalorimeterHitCollection&>&   simCaloHits
                                  ) const {
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
  // if (m_CaloHitContributionNames.size() != caloHitContribs.size()) {
  //   throw std::runtime_error(
  //       "The number of names for the output CaloHitContribution collections is not the expected one, expected: " +
  //       std::to_string(m_CaloHitContributionNames.size()) + " but the input has " +
  //       std::to_string(caloHitContribs.size()) + " collections");
  // }
  std::map<std::string, std::string> simCaloToContribution;
  for (size_t i = 0; i < simCaloHits.size(); ++i) {
    simCaloToContribution[std::next(simCaloHits.begin(), i)->first] = m_CaloHitContributionNames[i];
  }

  const auto seed = m_uidSvc->getUniqueID(headers[0].getEventNumber(), headers[0].getRunNumber(), this->name());
  m_engine.seed(seed);

  // Output collections
  auto oparticles       = edm4hep::MCParticleCollection();
  auto osimTrackerHits  = std::map<std::string, edm4hep::SimTrackerHitCollection>();
  auto osimCaloHits     = std::map<std::string, edm4hep::SimCalorimeterHitCollection>();
  auto ocaloHitContribs = std::map<std::string, edm4hep::CaloHitContributionCollection>();

  // Copy MCParticles for physics event into a new collection
  for (auto&& part : particles) {
    oparticles->push_back(part.clone(false));
  }
  // Fix relations to point to the new particles
  for (size_t i = 0; i < particles.size(); ++i) {
    for (auto& parent : particles[i].getParents()) {
      oparticles[i].addToParents(oparticles[parent.getObjectID().index]);
    }
    for (auto& daughter : particles[i].getDaughters()) {
      oparticles[i].addToDaughters(oparticles[daughter.getObjectID().index]);
    }
  }

  // Copy the SimTrackerHits and crop them
  for (auto it = simTrackerHits.begin(); it != simTrackerHits.end(); ++it) {
    auto& [name, coll]           = *it;
    auto [this_start, this_stop] = define_time_windows(name);
    auto ocoll                   = edm4hep::SimTrackerHitCollection();
    for (auto&& elem : coll) {
      const float tof = time_of_flight(elem.getPosition());
      if ((elem.getTime() > this_start + tof) && (elem.getTime() < this_stop + tof)) {
        auto nhit = elem.clone(false);
        nhit.setParticle(oparticles[elem.getParticle().getObjectID().index]);
        ocoll->push_back(nhit);
      }
    }
    osimTrackerHits[m_SimTrackerHitNames[std::distance(simTrackerHits.begin(), it)]] = std::move(ocoll);
  }

  // Copy the SimCalorimeterHits and crop them together with the contributions
  std::map<std::string, std::map<uint64_t, edm4hep::MutableSimCalorimeterHit>> cellIDsMap;
  for (auto& [name, coll] : simCaloHits) {
    auto [this_start, this_stop] = define_time_windows(name);
    auto& calHitMap              = cellIDsMap[name];
    auto& caloHitContribs        = ocaloHitContribs[simCaloToContribution[name]];
    for (auto&& elem : coll) {
      const float      tof                = time_of_flight(elem.getPosition());
      bool             within_time_window = false;
      std::vector<int> thisContribs;
      for (auto&& contrib : elem.getContributions()) {
        if (!((contrib.getTime() > this_start + tof) && (contrib.getTime() < this_stop + tof)))
          continue;
        within_time_window = true;
        // TODO: Make sure a contribution is not added twice
        auto newcontrib = contrib.clone(false);
        newcontrib.setParticle(oparticles[contrib.getParticle().getObjectID().index]);
        thisContribs.push_back(caloHitContribs.size());
        caloHitContribs.push_back(std::move(newcontrib));
      }
      if (within_time_window) {
        auto newhit = elem.clone(false);
        for (auto& contrib : thisContribs) {
          newhit.addToContributions(caloHitContribs[contrib]);
        }
        // if (elem.getCellID() == 18444492282485805327ULL) {
        //   info() << "Adding (crop) hit with cellID " << elem.getCellID() << " to collection " << name << " with index "
        //          << newhit.id().index << endmsg;
        // }
        calHitMap.emplace(elem.getCellID(), std::move(newhit));
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

        if (std::find(availableCollections.begin(), availableCollections.end(), _mcParticleCollectionName) ==
            availableCollections.end()) {
          warning() << "Collection " << _mcParticleCollectionName << " not found in background event" << endmsg;
        }

        // To fix the relations we will need to have a map from old to new particle index
        std::map<int, int>                                           oldToNewMap;
        std::map<int, std::pair<std::vector<int>, std::vector<int>>> parentDaughterMap;

        auto& particles = backgroundEvent.get<edm4hep::MCParticleCollection>(_mcParticleCollectionName);
        int   j         = oparticles.size();
        for (size_t i = 0; i < particles.size(); ++i) {
          auto npart = particles[i].clone(false);

          npart.setTime(particles[i].getTime() + timeOffset);
          npart.setOverlay(true);
          oparticles->push_back(npart);
          for (auto& parent : particles[i].getParents()) {
            parentDaughterMap[j].first.push_back(parent.getObjectID().index);
          }
          for (auto& daughter : particles[i].getDaughters()) {
            parentDaughterMap[j].second.push_back(daughter.getObjectID().index);
          }
          oldToNewMap[i] = j;
          j++;
        }
        for (auto& [index, parentsDaughters] : parentDaughterMap) {
          auto& [parents, daughters] = parentsDaughters;
          for (auto& parent : parents) {
            if (parentDaughterMap.find(oldToNewMap[parent]) == parentDaughterMap.end()) {
              // warning() << "Parent " << parent << " not found in background event" << endmsg;
              continue;
            }
            oparticles[index].addToParents(oparticles[oldToNewMap[parent]]);
          }
          for (auto& daughter : daughters) {
            if (parentDaughterMap.find(oldToNewMap[daughter]) == parentDaughterMap.end()) {
              // warning() << "Parent " << daughter << " not found in background event" << endmsg;
              continue;
            }
            // info() << "Adding (daughter) " << daughter << " to " << index << endmsg;
            oparticles[index].addToDaughters(oparticles[oldToNewMap[daughter]]);
          }
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
            for (auto&& elem : backgroundEvent.get<edm4hep::SimTrackerHitCollection>(name)) {
              const float tof = time_of_flight(elem.getPosition());

              if (!((elem.getTime() + timeOffset > this_start + tof) &&
                    (elem.getTime() + timeOffset < this_stop + tof)))
                continue;
              auto nhit = elem.clone(false);
              nhit.setOverlay(true);
              nhit.setParticle(oparticles[oldToNewMap[elem.getParticle().getObjectID().index]]);
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

          auto& calHitMap = cellIDsMap[name];
          auto& calHitContribs = ocaloHitContribs[simCaloToContribution[name]];
          for (auto&& elem : backgroundEvent.get<edm4hep::SimCalorimeterHitCollection>(name)) {
            if (calHitMap.find(elem.getCellID()) == calHitMap.end()) {
              // There is no hit at this position. The new hit can be added, if it is not outside the window
              auto calhit = edm4hep::MutableSimCalorimeterHit();
              bool add    = false;
              for (auto& contrib : elem.getContributions()) {
                if ((contrib.getTime() + timeOffset > this_start) && (contrib.getTime() + timeOffset < this_stop)) {
                  add = true;
                  // TODO: Make sure a contribution is not added twice
                  auto newcontrib = contrib.clone(false);
                  newcontrib.setParticle(oparticles[oldToNewMap[contrib.getParticle().getObjectID().index]]);
                  calhit.addToContributions(newcontrib);
                  calHitContribs.push_back(newcontrib);
                }
              }
              if (add) {
                calhit.setCellID(elem.getCellID());
                calhit.setEnergy(elem.getEnergy());
                calhit.setPosition(elem.getPosition());
                // if (calhit.getCellID() == 18444492282485805327ULL) {
                //   info() << "Adding hit with cellID " << calhit.getCellID() << " to collection " << name
                //          << " with index " << calhit.id().index << endmsg;
                // }
                calHitMap[calhit.getCellID()] = calhit;
              }
            } else {
              // there is already a hit at this position
              auto& calhit = calHitMap[elem.getCellID()];
              for (auto& contrib : elem.getContributions()) {
                if ((contrib.getTime() + timeOffset > this_start) && (contrib.getTime() + timeOffset < this_stop)) {
                  // TODO: Make sure a contribution is not added twice
                  auto newcontrib = contrib.clone(false);
                  newcontrib.setParticle(oparticles[oldToNewMap[contrib.getParticle().getObjectID().index]]);
                  calhit.addToContributions(newcontrib);
                  calHitContribs.push_back(newcontrib);
                }
              }
            }
          }
        }
      }
    }
  }
  for (auto& [name, calHitMap] : cellIDsMap) {
    info() << "There are " << calHitMap.size() << " hits in collection " << name << endmsg;
    auto ocoll = edm4hep::SimCalorimeterHitCollection();
    for (auto& [cellID, hit] : calHitMap) {
      if (hit.getCellID() == 18444492282485805327ULL) {
        info() << "Adding hit with cellID " << hit.getCellID() << " to collection " << name << " with index "
               << hit.id().index << endmsg;
      }
      ocoll->push_back(std::move(hit));
    }
    osimCaloHits[m_SimCalorimeterHitNames[std::distance(simCaloHits.begin(), simCaloHits.find(name))]] =
        std::move(ocoll);
  }

  info() << "oparticles.size() = " << oparticles.size() << endmsg;

  return std::make_tuple(std::move(oparticles), std::move(osimTrackerHits), std::move(osimCaloHits),
                         std::move(ocaloHitContribs));
}

StatusCode OverlayTiming::finalize() { return Gaudi::Algorithm::finalize(); }
