/*
 * FastJetAlg.h
 *
 *  Created on: 25.05.2010
 *      Author: Lars Weuste (MPP Munich) - weuste@mpp.mpg.de
 *			iterative inclusive algorithm based on design by Marco Battaglia (CERN) - Marco.Battaglia@cern.ch
 *  Converted to Gaudi on: 25.08.2025
 *      Conversion: Samuel Ferraro - samuel.rowles.ferraro@cern.ch
 */

#ifndef FASTJETALG_H_
#define FASTJETALG_H_

#include "EClusterMode.h"

#include <k4FWCore/Transformer.h>
#include <edm4hep/ReconstructedParticleCollection.h>

//FastJet
#include <fastjet/PseudoJet.hh>
#include <fastjet/JadePlugin.hh>
#include <fastjet/SISConePlugin.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/CDFJetCluPlugin.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/NestedDefsPlugin.hh>
#include <fastjet/CDFMidPointPlugin.hh>
#include <fastjet/EECambridgePlugin.hh>
#include <fastjet/SISConeSphericalPlugin.hh>
// #include <fastjet/contrib/ValenciaPlugin.hh>

#include <vector>
#include <string>
#include <stdexcept>
#include <map>

// some namespacing for organization
namespace k4Reco::FastJet {
  constexpr static int ITERATIVE_INCLUSIVE_MAX_ITERATIONS = 20;
  
  static const std::map<std::string, fastjet::JetAlgorithm> NAME_TO_ALGORITHM_MAP = {
    {"kt_algorithm", fastjet::kt_algorithm},
    {"cambridge_algorithm", fastjet::cambridge_algorithm},
    {"antikt_algorithm", fastjet::antikt_algorithm},
    {"genkt_algorithm", fastjet::genkt_algorithm},
    {"cambridge_for_passive_algorithm", fastjet::cambridge_for_passive_algorithm},
    {"genkt_for_passive_algorithm", fastjet::genkt_for_passive_algorithm},
    {"ee_kt_algorithm", fastjet::ee_kt_algorithm},
    {"ee_genkt_algorithm", fastjet::ee_genkt_algorithm},
  };
  static const std::map<std::string, int> NAME_TO_NR_PARAMS_MAP = {
    {"kt_algorithm", 1},
    {"cambridge_algorithm", 1},
    {"antikt_algorithm", 1},
    {"genkt_algorithm", 2},
    {"cambridge_for_passive_algorithm", 1},
    {"genkt_for_passive_algorithm", 1},
    {"ee_kt_algorithm", 0},
    {"ee_genkt_algorithm", 2},
    {"SISConePlugin", 2},
    {"SISConeSphericalPlugin", 2},
  };
  static const std::map<std::string, int> NAME_TO_CLUSTER_MODE_MAP = {
    {"kt_algorithm", FJ_inclusive | FJ_exclusive_nJets | FJ_exclusive_yCut | OWN_inclusiveIteration},
    {"cambridge_algorithm", FJ_inclusive | FJ_exclusive_nJets | FJ_exclusive_yCut | OWN_inclusiveIteration},
    {"antikt_algorithm", FJ_inclusive | OWN_inclusiveIteration},
    {"genkt_algorithm", FJ_inclusive | OWN_inclusiveIteration | FJ_exclusive_nJets | FJ_exclusive_yCut},
    {"cambridge_for_passive_algorithm", FJ_inclusive | OWN_inclusiveIteration | FJ_exclusive_nJets | FJ_exclusive_yCut},
    {"genkt_for_passive_algorithm", FJ_inclusive | OWN_inclusiveIteration},
    {"ee_kt_algorithm", FJ_exclusive_nJets | FJ_exclusive_yCut},
    {"ee_genkt_algorithm", FJ_inclusive | FJ_exclusive_nJets | FJ_exclusive_yCut},
    {"SISConePlugin", FJ_inclusive | OWN_inclusiveIteration},
    {"SISConeSphericalPlugin", FJ_inclusive | OWN_inclusiveIteration},
  };

  class JetDefinitionFactory {
    using creatorFunc = std::function<std::unique_ptr<fastjet::JetDefinition>(
						       fastjet::JetAlgorithm,
						       const std::vector<float>&,
						       fastjet::RecombinationScheme,
						       fastjet::Strategy
						       )>;
    std::map<std::string, creatorFunc> registry;
  public:
    JetDefinitionFactory();
    std::unique_ptr<fastjet::JetDefinition> create(const std::string& type,
						   fastjet::JetAlgorithm m_jetAlgoType,
						   const std::vector<float>& params,
						   fastjet::RecombinationScheme m_jetRecoScheme,
						   fastjet::Strategy m_strategy);
  private:
    static std::unique_ptr<fastjet::JetDefinition> useZeroParams(fastjet::JetAlgorithm m_jetAlgoType,
							  const std::vector<float>& params,
							  fastjet::RecombinationScheme m_jetRecoScheme,
							  fastjet::Strategy m_strategy) {
      return std::make_unique<fastjet::JetDefinition>(m_jetAlgoType, m_jetRecoScheme, m_strategy);
    }
    static std::unique_ptr<fastjet::JetDefinition> useOneParams(
							 fastjet::JetAlgorithm m_jetAlgoType,
							 const std::vector<float>& params,
							 fastjet::RecombinationScheme m_jetRecoScheme,
							 fastjet::Strategy m_strategy
							 ){
      return std::make_unique<fastjet::JetDefinition>(m_jetAlgoType, params.at(0), m_jetRecoScheme, m_strategy);
    }
    static std::unique_ptr<fastjet::JetDefinition> useTwoParams(
							 fastjet::JetAlgorithm m_jetAlgoType,
							 const std::vector<float>& params,
							 fastjet::RecombinationScheme m_jetRecoScheme,
							 fastjet::Strategy m_strategy
							 ) {
      return std::make_unique<fastjet::JetDefinition>(m_jetAlgoType, params.at(0), params.at(1), m_jetRecoScheme, m_strategy);
    }
  };

  //Forward declaration
typedef std::vector< fastjet::PseudoJet > PseudoJetList;

class SkippedFixedNrJetException: public std::runtime_error {
public:
    SkippedFixedNrJetException():std::runtime_error("") {}
};

class SkippedMaxIterationException: public std::runtime_error {
public:
    SkippedMaxIterationException(PseudoJetList& jets) :std::runtime_error(""), m_jets(jets) {}
    PseudoJetList m_jets;
};
}

struct FastJetAlg : k4FWCore::MultiTransformer<
  std::tuple<edm4hep::ReconstructedParticleCollection,
	     edm4hep::ReconstructedParticleCollection>(const edm4hep::ReconstructedParticleCollection&)> {
public:
  FastJetAlg(const std::string& name, ISvcLocator* svcLoc);

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  StatusCode initialize();

  /** Called for every run.
   */
  std::tuple<edm4hep::ReconstructedParticleCollection, edm4hep::ReconstructedParticleCollection> operator()(
													    const edm4hep::ReconstructedParticleCollection& inputCollection) const;

private:
  std::vector<std::string> defaultJetAlgoNameAndParams{"kt_algorithm", "0.7"};
  std::vector<std::string> defaultClusterMode{"Inclusive", "0.0"};

  Gaudi::Property<std::string> m_jetAlgoName{this, "algorithm", "kt_algorithm", "Name of the algorithm to use for making jets. E.g. kt_algorithm. Full list of algorithms can be seen in the FastJet code."};
  Gaudi::Property<std::vector<float>> m_jetAlgoParams{this, "algorithmParameters", {0.7}, "Parameters required by each specific algorithm. The amount of parameters and interpretation varies per algorithm."};
  Gaudi::Property<std::string> m_clusterModeName{this, "clusteringMode", "Inclusive", "Clustering mode for the algorithm. Not all clustering modes are available for all algorithms. See Fast Jet code for availabilities."};
  Gaudi::Property<std::vector<float>> m_clusterModeParams{this, "clusteringParams", {0.0}, "Cluster mode parameters. One of 'Inclusive <minPt>', 'InclusiveIterativeNJets <nrJets> <minE>', 'ExclusiveNJets <nrJets>', 'ExclusiveYCut <yCut>'. Note: not all modes are available for all algorithms. Some parameters are input as floats, despite being eventual conversion to integers."};
  
  Gaudi::Property<std::string> m_jetRecoSchemeName{this, "recombinationScheme", std::string("E_scheme"), "The recombination scheme used when merging 2 particles. Usually there is no need to use anything else than 4-Vector addition: E_scheme."};

  // jet algorithm
  //std::string m_jetAlgoName;
  std::unique_ptr<fastjet::JetDefinition> m_jetAlgo;
  fastjet::JetAlgorithm m_jetAlgoType;

  // clustering mode
  //std::string m_clusterModeName;
  k4Reco::FastJet::EClusterMode m_clusterMode;

  // jet reco scheme
  fastjet::RecombinationScheme m_jetRecoScheme;

  // jet strategy
  std::string m_strategyName;
  fastjet::Strategy m_strategy;

  // parameters
  unsigned m_requestedNumberOfJets;
  double m_yCut;
  double m_minPt;
  double m_minE;

private:
  fastjet::JetAlgorithm getAlgoType() const;
  bool validateParams();
  bool validateClusterModes() const;
  std::unique_ptr<k4Reco::FastJet::JetDefinitionFactory> theJetDefinitionFactory = std::make_unique<k4Reco::FastJet::JetDefinitionFactory>();
  
};

DECLARE_COMPONENT(FastJetAlg)

#endif /* FASTJETALG_H_ */
