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

#define ITERATIVE_INCLUSIVE_MAX_ITERATIONS 20

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

    Gaudi::Property<std::string> m_jetRecoSchemeName{this, "recombinationScheme", std::string("E_scheme"), "The recombination scheme used when merging 2 particles. Usually there is no need to use anything else than 4-Vector addition: E_scheme."};
    Gaudi::Property<std::vector<std::string>> m_jetAlgoNameAndParams{this, "algorithm", defaultJetAlgoNameAndParams, "Selects the algorithm and its parameters. E.g. 'kt_algorithm 0.7' or 'ee_kt_algorithm'. For a full list of supported algorithms, see the logfile after execution."};
    Gaudi::Property<std::vector<std::string>> m_clusterModeNameAndParam{this, "clusteringMode", defaultClusterMode, "One of 'Inclusive <minPt>', 'InclusiveIterativeNJets <nrJets> <minE>', 'ExclusiveNJets <nrJets>', 'ExclusiveYCut <yCut>'. Note: not all modes are available for all algorithms."};

    // jet algorithm
    std::string m_jetAlgoName;
    fastjet::JetDefinition* m_jetAlgo;
    fastjet::JetAlgorithm m_jetAlgoType;

    // clustering mode
    std::string m_clusterModeName;
    EClusterMode m_clusterMode;

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
    bool isJetAlgo(std::string algo, int nrParams, int supportedModes) const;
};

DECLARE_COMPONENT(FastJetAlg)

bool FastJetAlg::isJetAlgo(std::string algo, int nrParams, int supportedModes) const {
    info() << " " << algo;

    // check if the chosen algorithm is the same as it was passed
    if (m_jetAlgoName.compare(algo) != 0) {
    return false;
    }
    info() << "*"; // mark the current algorithm as the selected
        // one, even before we did our checks on nr of
        // parameters etc.

    // check if we have enough number of parameters
    if ((int)m_jetAlgoNameAndParams.size() - 1 != nrParams) {
        error() << std::endl
        << "Wrong numbers of parameters for algorithm: " << algo << std::endl
        << "We need " << nrParams << " params, but we got " << m_jetAlgoNameAndParams.size() - 1 << std::endl;
        throw GaudiException("You have insufficient number of parameters for this algorithm! See log for more details.", name(), StatusCode::FAILURE);
    }

    // check if the mode is supported via a binary AND
    if ((supportedModes & m_clusterMode) != m_clusterMode) {
        error() << std::endl
        << "This algorithm is not capable of running in this clustering mode ("
        << m_clusterMode << "). Sorry!" << std::endl;
        throw GaudiException("This algorithm is not capable of running in this mode", name(), StatusCode::FAILURE);
    }

    return true;
}


#endif /* FASTJETALG_H_ */