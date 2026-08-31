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

/*
 * FastJetAlg.cpp
 *
 *  Created on: 25.05.2010
 *      Author: Lars Weuste (MPP Munich) - weuste@mpp.mpg.de
 *			iterative inclusive algorithm based on design by Marco Battaglia (CERN) -
 * Marco.Battaglia@cern.ch Converted to Gaudi on: 25.08.2025 Conversion: Samuel Ferraro - samuel.rowles.ferraro@cern.ch
 */

#include "FastJetAlg.hxx"

#include <edm4hep/MutableReconstructedParticle.h>
#include <edm4hep/ReconstructedParticle.h>

#include <fmt/format.h>
#include <fmt/ranges.h>

#include <ranges>

#include <k4FWCore/GaudiChecks.h>

using namespace k4Reco::FastJet;

// The registry maps each algorithm name onto the JetDefinition constructor it needs
JetDefinitionFactory::JetDefinitionFactory() {
  registry["kt_algorithm"] = useOneParams;
  registry["cambridge_algorithm"] = useOneParams;
  registry["antikt_algorithm"] = useOneParams;
  registry["genkt_algorithm"] = useTwoParams;
  registry["cambridge_for_passive_algorithm"] = useOneParams;
  registry["genkt_for_passive_algorithm"] = useOneParams;
  registry["ee_kt_algorithm"] = useZeroParams;
  registry["ee_genkt_algorithm"] = useTwoParams;
  registry["SISConePlugin"] = [](fastjet::JetAlgorithm, const std::vector<float>& params, fastjet::RecombinationScheme,
                                 fastjet::Strategy) {
    fastjet::SISConePlugin* pl;
    pl = new fastjet::SISConePlugin(params.at(0), params.at(1));
    auto jetAlgo = std::make_unique<fastjet::JetDefinition>(pl);
    jetAlgo->delete_plugin_when_unused();
    return jetAlgo;
  };
  registry["SISConeSphericalPlugin"] = [](fastjet::JetAlgorithm, const std::vector<float>& params,
                                          fastjet::RecombinationScheme, fastjet::Strategy) {
    fastjet::SISConeSphericalPlugin* pl;
    pl = new fastjet::SISConeSphericalPlugin(params.at(0), params.at(1));
    auto jetAlgo = std::make_unique<fastjet::JetDefinition>(pl);
    jetAlgo->delete_plugin_when_unused();
    return jetAlgo;
  };
}

std::unique_ptr<fastjet::JetDefinition> JetDefinitionFactory::create(const std::string& type,
                                                                     fastjet::JetAlgorithm m_jetAlgoType,
                                                                     const std::vector<float>& params,
                                                                     fastjet::RecombinationScheme m_jetRecoScheme,
                                                                     fastjet::Strategy m_strategy) {
  if (const auto it = registry.find(type); it != registry.end()) {
    return it->second(m_jetAlgoType, params, m_jetRecoScheme, m_strategy);
  }
  throw GaudiException("Jet algorithm type not in factory registry", "JetDefinitionFactory", StatusCode::FAILURE);
}

FastJetAlg::FastJetAlg(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(name, svcLoc, {KeyValue("recParticleIn", "MCParticle")},
                       {KeyValue("jetOut", "JetOut"), KeyValue("recParticleOut", "Constituents")}),
      m_jetAlgo(nullptr), m_jetAlgoType(), m_clusterMode(NONE), m_jetRecoScheme(), m_strategyName(""), m_strategy(),
      m_requestedNumberOfJets(0), m_yCut(0.0), m_minPt(0.0), m_minE(0.0) {}

StatusCode FastJetAlg::initialize() {
  m_strategy = fastjet::Best;
  m_strategyName = "Best";
  info() << "Strategy: " << m_strategyName << endmsg;

  if (m_jetRecoSchemeName.value().compare("E_scheme") == 0)
    m_jetRecoScheme = fastjet::E_scheme;
  else if (m_jetRecoSchemeName.value().compare("pt_scheme") == 0)
    m_jetRecoScheme = fastjet::pt_scheme;
  else if (m_jetRecoSchemeName.value().compare("pt2_scheme") == 0)
    m_jetRecoScheme = fastjet::pt2_scheme;
  else if (m_jetRecoSchemeName.value().compare("Et_scheme") == 0)
    m_jetRecoScheme = fastjet::Et_scheme;
  else if (m_jetRecoSchemeName.value().compare("Et2_scheme") == 0)
    m_jetRecoScheme = fastjet::Et2_scheme;
  else if (m_jetRecoSchemeName.value().compare("BIpt_scheme") == 0)
    m_jetRecoScheme = fastjet::BIpt_scheme;
  else if (m_jetRecoSchemeName.value().compare("BIpt2_scheme") == 0)
    m_jetRecoScheme = fastjet::BIpt2_scheme;
  else {
    error() << "Unknown recombination scheme: " << m_jetRecoSchemeName << endmsg;
    return StatusCode::FAILURE;
  }
  info() << "recombination scheme: " << m_jetRecoSchemeName << endmsg;

  K4_GAUDI_CHECK(validateClusterModeParams());

  // The parameters have been validated above, so they can simply be assigned here
  m_clusterMode = NAME_TO_CLUSTER_MODE_PARAMS_MAP.at(m_clusterModeName).mode;
  switch (m_clusterMode) {
  case FJ_inclusive:
    m_minPt = m_clusterModeParams[0];
    break;
  case OWN_inclusiveIteration:
    m_requestedNumberOfJets = (int)m_clusterModeParams[0];
    m_minE = m_clusterModeParams[1];
    break;
  case FJ_exclusive_nJets:
    m_requestedNumberOfJets = (int)m_clusterModeParams[0];
    break;
  case FJ_exclusive_yCut:
    m_yCut = m_clusterModeParams[0];
    break;
  case NONE:
    break;
  }
  info() << "Cluster mode: " << m_clusterMode << endmsg;

  // The SISCone plugins are not built from a fastjet::JetAlgorithm, so they have no type to look up
  if (m_jetAlgoName != "SISConePlugin" and m_jetAlgoName != "SISConeSphericalPlugin")
    m_jetAlgoType = getAlgoType();

  K4_GAUDI_CHECK(validateParams());
  K4_GAUDI_CHECK(validateClusterModes());

  m_jetAlgo =
      theJetDefinitionFactory->create(m_jetAlgoName, m_jetAlgoType, m_jetAlgoParams, m_jetRecoScheme, m_strategy);

  return StatusCode::SUCCESS;
}

std::tuple<edm4hep::ReconstructedParticleCollection, edm4hep::ReconstructedParticleCollection>
FastJetAlg::operator()(const edm4hep::ReconstructedParticleCollection& inputCollection) const {
  edm4hep::ReconstructedParticleCollection outputCollection;
  outputCollection.setSubsetCollection(true);
  edm4hep::ReconstructedParticleCollection jetCollection;

  PseudoJetList jets;
  PseudoJetList pjList;
  for (size_t i = 0; i < inputCollection.size(); ++i) {
    edm4hep::ReconstructedParticle par = inputCollection.at(i);
    pjList.push_back(
        fastjet::PseudoJet(par.getMomentum().x, par.getMomentum().y, par.getMomentum().z, par.getEnergy()));
    pjList.back().set_user_index(i); // save the id of this recParticle
  }

  fastjet::ClusterSequence cs = fastjet::ClusterSequence(pjList, *m_jetAlgo);

  if (m_clusterMode == FJ_inclusive) {
    jets = cs.inclusive_jets(m_minPt);
  } else if (m_clusterMode == FJ_exclusive_yCut) {
    jets = cs.exclusive_jets_ycut(m_yCut);
  } else if (m_clusterMode == FJ_exclusive_nJets) {
    // sanity check: if we have not enough particles, FJ will cause an assert
    if (inputCollection.size() < (size_t)m_requestedNumberOfJets) {
      warning() << "Not enough elements in the input collection to create " << m_requestedNumberOfJets << " jets."
                << endmsg;
    } else {
      jets = cs.exclusive_jets((int)(m_requestedNumberOfJets));
    }
  } else if (m_clusterMode == OWN_inclusiveIteration) {
    // sanity check: if we have not enough particles, FJ will cause an assert
    if (inputCollection.size() < (size_t)m_requestedNumberOfJets) {
      warning() << "Not enough elements in the input collection to create " << m_requestedNumberOfJets << " jets."
                << endmsg;

    } else {
      // lets do a iterative procedure until we found the correct number of jets
      // for that we will do inclusive clustering, modifying the R parameter in some kind of minimization
      // this is based on Marco Battaglia's FastJetClustering
      double R = M_PI_4;    // maximum of R is Pi/2, minimum is 0. So we start hat Pi/4
      double RDiff = R / 2; // the step size we modify the R parameter at each iteration. Its size for the n-th step is
                            // R/(2n), i.e. starts with R/2
      PseudoJetList jets_it;
      unsigned nJets;
      int iIter = 0; // nr of current iteration
      // these variables are only used if the SisCone(Spherical)Plugin is selected
      // This is necessary, as these are plugins and hence use a different constructor than
      // the built in fastjet algorithms
      // here we save pointer to the plugins, so that if they are created, we can delete them again after usage
      std::unique_ptr<fastjet::SISConePlugin> pluginSisCone;
      std::unique_ptr<fastjet::SISConeSphericalPlugin> pluginSisConeSph;
      // check if we use the siscones
      bool useSisCone = m_jetAlgoName.value().compare("SISConePlugin") == 0;
      bool useSisConeSph = m_jetAlgoName.value().compare("SISConeSphericalPlugin") == 0;
      // save the 2nd parameter of the SisCones
      double sisConeOverlapThreshold = 0;
      if (useSisCone || useSisConeSph)
        sisConeOverlapThreshold = m_jetAlgoParams[1];
      // do a maximum of N iterations
      for (iIter = 0; iIter < ITERATIVE_INCLUSIVE_MAX_ITERATIONS; iIter++) {
        // do the clustering for this value of R. For this we need to re-initialize the JetDefinition, as it takes the R
        // parameter
        std::unique_ptr<fastjet::JetDefinition> jetDefinition;
        // unfortunately SisCone(spherical) are being initialized differently, so we have to check for this
        if (useSisCone) {
          pluginSisCone = std::make_unique<fastjet::SISConePlugin>(R, sisConeOverlapThreshold);
          jetDefinition = std::make_unique<fastjet::JetDefinition>(pluginSisCone.get());
        } else if (useSisConeSph) {
          pluginSisConeSph = std::make_unique<fastjet::SISConeSphericalPlugin>(R, sisConeOverlapThreshold);
          jetDefinition = std::make_unique<fastjet::JetDefinition>(pluginSisConeSph.get());
        } else {
          jetDefinition = std::make_unique<fastjet::JetDefinition>(m_jetAlgoType, R, m_jetRecoScheme, m_strategy);
        }
        // now we can finally create the cluster sequence
        fastjet::ClusterSequence cs_it(pjList, *jetDefinition);
        jets_it = cs_it.inclusive_jets(0); // no pt cut, we will do an energy cut
        jets.clear();
        // keep only the jets above the energy threshold
        for (unsigned j = 0; j < jets_it.size(); j++) {
          if (jets_it[j].E() > m_minE) {
            jets.push_back(jets_it[j]);
          }
        }

        nJets = jets.size();
        debug() << iIter << " " << R << " " << jets_it.size() << " " << nJets << endmsg;
        if (nJets == m_requestedNumberOfJets) { // if the number of jets is correct: success!
          break;
        } else if (nJets < m_requestedNumberOfJets) {
          // if number of jets is too small: we need a smaller Radius per jet (so
          // that we get more jets)
          R -= RDiff;
        } else if (nJets > m_requestedNumberOfJets) { // if the number of jets is too
          // high: increase the Radius
          R += RDiff;
        }
        RDiff /= 2;
      }
    }
  }

  for (const auto& jet : jets) {
    edm4hep::MutableReconstructedParticle rec = jetCollection.create();
    rec.setEnergy(jet.E());
    rec.setMass(jet.m());
    rec.setMomentum(edm4hep::Vector3f(jet.px(), jet.py(), jet.pz()));

    for (const auto& constituent : cs.constituents(jet)) {
      const auto particle = inputCollection.at(constituent.user_index());
      rec.addToParticles(particle);
      outputCollection.push_back(particle);
    }
  }

  return std::make_tuple(std::move(jetCollection), std::move(outputCollection));
}

fastjet::JetAlgorithm FastJetAlg::getAlgoType() const {
  if (!k4Reco::FastJet::NAME_TO_ALGORITHM_MAP.contains(m_jetAlgoName)) {
    throw GaudiException("Could not find jet algorithm name in internal map", name(), StatusCode::FAILURE);
  }
  return k4Reco::FastJet::NAME_TO_ALGORITHM_MAP.at(m_jetAlgoName);
}

// Checks the number of algorithm parameters against the expected number
bool FastJetAlg::validateParams() {
  const int nParams = static_cast<int>(m_jetAlgoParams.size());

  // ee_genkt_algorithm may be given R only, in which case the exponent p defaults to 1
  if (m_jetAlgoName == "ee_genkt_algorithm") {
    if (nParams == 1) {
      m_jetAlgoParams.value().push_back(1.);
      info() << "When only 1 parameter is provided for ee_genkt_algorithm it is assumed to be R, and the exponent p is "
                "assumed to be equal to 1"
             << endmsg;
      return true;
    }
    return nParams == 2;
  }

  const auto it = k4Reco::FastJet::NAME_TO_NR_PARAMS_MAP.find(m_jetAlgoName);
  return it != k4Reco::FastJet::NAME_TO_NR_PARAMS_MAP.end() && it->second == nParams;
}

// Checks the clustering mode and the number of parameters it expects. The names of the
// expected parameters are part of that mapping and are used for the error messages
bool FastJetAlg::validateClusterModeParams() const {
  const auto it = k4Reco::FastJet::NAME_TO_CLUSTER_MODE_PARAMS_MAP.find(m_clusterModeName);
  if (it == k4Reco::FastJet::NAME_TO_CLUSTER_MODE_PARAMS_MAP.end()) {
    error() << fmt::format("Unknown cluster mode: '{}'. Available modes are: {}", m_clusterModeName.value(),
                           fmt::join(k4Reco::FastJet::NAME_TO_CLUSTER_MODE_PARAMS_MAP | std::views::keys, ", "))
            << endmsg;
    return false;
  }

  const auto& paramNames = it->second.paramNames;
  if (m_clusterModeParams.size() != paramNames.size()) {
    error() << fmt::format("Wrong number of parameters for clustering mode '{}'. Expected {} ({}) but got {}",
                           m_clusterModeName.value(), paramNames.size(), fmt::join(paramNames, ", "),
                           m_clusterModeParams.size())
            << endmsg;
    return false;
  }

  return true;
}

// Checks that the requested clustering mode is supported by the requested algorithm
bool FastJetAlg::validateClusterModes() const {
  const auto it = k4Reco::FastJet::NAME_TO_CLUSTER_MODE_MAP.find(m_jetAlgoName);
  return it != k4Reco::FastJet::NAME_TO_CLUSTER_MODE_MAP.end() && (it->second & m_clusterMode) == m_clusterMode;
}
