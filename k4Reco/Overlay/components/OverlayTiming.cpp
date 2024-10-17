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
#include "OverlayTiming.h"
#include <GaudiKernel/MsgStream.h>

#include "podio/Frame.h"

#include "edm4hep/Constants.h"
#include "edm4hep/MutableCaloHitContribution.h"

#include "k4FWCore/MetadataUtils.h"

#include <TMath.h>

#include <limits>
#include <random>

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
      error() << err << endmsg;
      return StatusCode::FAILURE;
    }
  }

  if (std::any_of(m_bkgEvents->m_totalNumberOfEvents.begin(), m_bkgEvents->m_totalNumberOfEvents.end(),
                  [this](const int& val) { return this->m_startWithBackgroundEvent >= val; })) {
    throw GaudiException("StartWithBackgroundEvent is larger than the number of events in the background files", name(),
                         StatusCode::FAILURE);
  }

  if (m_Noverlay.empty()) {
    info() << "Using the default number of overlay events (1) for each group, since none was specified with "
              "NumberBackground "
           << endmsg;
    m_Noverlay = std::vector<double>(m_bkgEvents->size(), 1);
  }

  if (m_Poisson.empty()) {
    info() << "Using the default overlay mode (no Poission distribution) for each group, since none was specified with "
              "Poisson_random_NOverlay"
           << endmsg;
    m_Poisson = std::vector<bool>(m_bkgEvents->size(), false);
  }

  return StatusCode::SUCCESS;
}

retType OverlayTiming::operator()(const edm4hep::EventHeaderCollection&                           headers,
                                  const edm4hep::MCParticleCollection&                            particles,
                                  const std::vector<const edm4hep::SimTrackerHitCollection*>&     simTrackerHits,
                                  const std::vector<const edm4hep::SimCalorimeterHitCollection*>& simCaloHits) const {
  const auto seed = m_uidSvc->getUniqueID(headers[0].getEventNumber(), headers[0].getRunNumber(), this->name());
  m_engine.seed(seed);

  // Output collections
  auto oparticles       = edm4hep::MCParticleCollection();
  auto osimTrackerHits  = std::vector<edm4hep::SimTrackerHitCollection>();
  auto osimCaloHits     = std::vector<edm4hep::SimCalorimeterHitCollection>();
  auto ocaloHitContribs = std::vector<edm4hep::CaloHitContributionCollection>();
  for (size_t i = 0; i < simCaloHits.size(); ++i) {
    ocaloHitContribs.emplace_back(edm4hep::CaloHitContributionCollection());
  }

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
  for (size_t i = 0; i < simTrackerHits.size(); ++i) {
    auto& coll                   = simTrackerHits[i];
    auto  name                   = inputLocations(SIMTRACKERHIT_INDEX_POSITION)[i];
    auto [this_start, this_stop] = define_time_windows(name);
    auto ocoll                   = edm4hep::SimTrackerHitCollection();
    for (auto&& elem : *coll) {
      const float tof = time_of_flight(elem.getPosition());
      if ((elem.getTime() > this_start + tof) && (elem.getTime() < this_stop + tof)) {
        auto nhit = elem.clone(false);
        if (elem.getParticle().getObjectID().index != -1)
          nhit.setParticle(oparticles[elem.getParticle().getObjectID().index]);
        ocoll->push_back(nhit);
      }
    }
    osimTrackerHits.emplace_back(std::move(ocoll));
  }

  // Copy the SimCalorimeterHits and crop them together with the contributions
  std::map<int, std::map<uint64_t, edm4hep::MutableSimCalorimeterHit>> cellIDsMap;
  for (size_t i = 0; i < simCaloHits.size(); ++i) {
    auto& coll                   = simCaloHits[i];
    auto  name                   = inputLocations(SIMCALOHIT_INDEX_POSITION)[i];
    auto [this_start, this_stop] = define_time_windows(name);
    auto& calHitMap              = cellIDsMap[i];
    auto& caloHitContribs        = ocaloHitContribs[i];
    for (auto&& elem : *coll) {
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
        calHitMap.emplace(elem.getCellID(), std::move(newhit));
      }
    }
  }

  // Iterate over each group of files and parameters
  for (size_t groupIndex = 0; groupIndex < m_bkgEvents->size(); groupIndex++) {
    if (m_randomBX) {
      m_physBX = std::uniform_int_distribution<int>(0, _nBunchTrain - 1)(m_engine);
      debug() << "Physics Event was placed in the " << m_physBX << " bunch crossing!" << endmsg;
    }

    // define a permutation for the events to overlay -- the physics event is per definition at position 0
    std::vector<int> permutation;

    // Permutation has negative values and the last one is 0
    // if (!m_randomBX) then m_physBX (default = 1)
    for (int i = -(m_physBX - 1); i < _nBunchTrain - (m_physBX - 1); ++i) {
      permutation.push_back(i);
    }
    std::shuffle(permutation.begin(), permutation.end(), m_engine);

    // TODO: Check that there is anything to overlay

    debug() << "Starting overlay at event: " << m_bkgEvents->m_nextEntry[groupIndex] << " for the background group "
            << groupIndex << endmsg;

    if (m_startWithBackgroundEvent >= 0) {
      info() << "Skipping to event: " << m_startWithBackgroundEvent << endmsg;
      for (auto& entry : m_bkgEvents->m_nextEntry) {
        entry = m_startWithBackgroundEvent;
      }
    }

    // Overlay the background events to each bunchcrossing in the bunch train
    for (int bxInTrain = 0; bxInTrain < _nBunchTrain; ++bxInTrain) {
      const int BX_number_in_train = permutation.at(bxInTrain);

      int NOverlay_to_this_BX = 0;

      if (m_Poisson[groupIndex]) {
        NOverlay_to_this_BX = std::poisson_distribution<>(m_Noverlay[groupIndex])(m_engine);
      } else {
        NOverlay_to_this_BX = m_Noverlay[groupIndex];
      }

      debug() << "Will overlay " << NOverlay_to_this_BX << " events to BX number " << BX_number_in_train + m_physBX
              << endmsg;

      for (int k = 0; k < NOverlay_to_this_BX; ++k) {
        info() << "Overlaying background event " << m_bkgEvents->m_nextEntry[groupIndex] << " from group " << groupIndex
               << " to BX " << bxInTrain << endmsg;
        if (m_bkgEvents->m_nextEntry[groupIndex] >= m_bkgEvents->m_totalNumberOfEvents[groupIndex] &&
            !m_allowReusingBackgroundFiles) {
          throw GaudiException("No more events in background file", name(), StatusCode::FAILURE);
        }
        auto backgroundEvent =
            m_bkgEvents->m_rootFileReaders[groupIndex].readEvent(m_bkgEvents->m_nextEntry[groupIndex]);
        m_bkgEvents->m_nextEntry[groupIndex]++;
        m_bkgEvents->m_nextEntry[groupIndex] %= m_bkgEvents->m_totalNumberOfEvents[groupIndex];
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

        auto& bgParticles = backgroundEvent.get<edm4hep::MCParticleCollection>(_mcParticleCollectionName);
        int   j           = oparticles.size();
        for (size_t i = 0; i < bgParticles.size(); ++i) {
          auto npart = bgParticles[i].clone(false);

          npart.setTime(bgParticles[i].getTime() + timeOffset);
          npart.setOverlay(true);
          oparticles->push_back(npart);
          for (auto& parent : bgParticles[i].getParents()) {
            parentDaughterMap[j].first.push_back(parent.getObjectID().index);
          }
          for (auto& daughter : bgParticles[i].getDaughters()) {
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

        for (size_t i = 0; i < simTrackerHits.size(); ++i) {
          auto name = inputLocations(SIMTRACKERHIT_INDEX_POSITION)[i];
          debug() << "Processing collection " << name << endmsg;
          if (std::find(availableCollections.begin(), availableCollections.end(), name) == availableCollections.end()) {
            warning() << "Collection " << name << " not found in background event" << endmsg;
            continue;
          }
          auto [this_start, this_stop] = define_time_windows(name);
          // There are only contributions to the readout if the hits are in the integration window
          if (this_stop <= (BX_number_in_train - m_physBX) * _T_diff) {
            info() << "Skipping collection " << name << " as it is not in the integration window" << endmsg;
            continue;
          }
          auto& ocoll = osimTrackerHits[i];
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
        }

        for (size_t i = 0; i < simCaloHits.size(); ++i) {
          auto name = inputLocations(SIMCALOHIT_INDEX_POSITION)[i];
          debug() << "Processing collection " << name << endmsg;
          if (std::find(availableCollections.begin(), availableCollections.end(), name) == availableCollections.end()) {
            warning() << "Collection " << name << " not found in background event" << endmsg;
            continue;
          }
          auto [this_start, this_stop] = define_time_windows(name);
          // There are only contributions to the readout if the hits are in the integration window
          if (this_stop <= (BX_number_in_train - m_physBX) * _T_diff) {
            info() << "Skipping collection " << name << " as it is not in the integration window" << endmsg;
            continue;
          }

          auto& calHitMap      = cellIDsMap[i];
          auto& calHitContribs = ocaloHitContribs[i];
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
  // Move the SimCalorimeterHitCollections to the output vector
  // So far they are stored in a map with the cellID as key
  // but they don't belong to any collection yet
  for (auto& [index, calHitMap] : cellIDsMap) {
    auto ocoll = edm4hep::SimCalorimeterHitCollection();
    for (auto& [cellID, hit] : calHitMap) {
      ocoll->push_back(std::move(hit));
    }
    osimCaloHits.emplace_back(std::move(ocoll));
  }

  return std::make_tuple(std::move(oparticles), std::move(osimTrackerHits), std::move(osimCaloHits),
                         std::move(ocaloHitContribs));
}

StatusCode OverlayTiming::finalize() {
  if (m_copyCellIDMetadata) {
    for (auto& [input, output] :
         {std::make_pair(inputLocations("SimTrackerHits"), outputLocations("OutputSimTrackerHits")),
          std::make_pair(inputLocations("SimCalorimeterHits"), outputLocations("OutputSimCalorimeterHits"))}) {
      for (size_t i = 0; i < input.size(); ++i) {
        auto name  = input[i];
        auto value = k4FWCore::getParameter<std::string>(name + "__" + edm4hep::labels::CellIDEncoding, this);
        if (value.has_value()) {
          k4FWCore::putParameter(output[i] + "__" + edm4hep::labels::CellIDEncoding, value.value(), this);
        } else {
          warning() << "No metadata found for " << name << " when copying CellID metadata was requested" << endmsg;
        }
      }
    }
  }

  return Gaudi::Algorithm::finalize();
}

DECLARE_COMPONENT(OverlayTiming)
