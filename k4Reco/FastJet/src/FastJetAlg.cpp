/*
 * FastJetAlg.cpp
 *
 *  Created on: 25.05.2010
 *      Author: Lars Weuste (MPP Munich) - weuste@mpp.mpg.de
 *			iterative inclusive algorithm based on design by Marco Battaglia (CERN) - Marco.Battaglia@cern.ch
 *  Converted to Gaudi on: 25.08.2025
 *      Conversion: Samuel Ferraro - samuel.rowles.ferraro@cern.ch
 */

#include "FastJetAlg.hxx"

#include <edm4hep/MCParticle.h>
#include <edm4hep/MutableVertex.h>
#include <edm4hep/MutableParticleID.h>
#include <edm4hep/ReconstructedParticle.h>
#include <edm4hep/MutableReconstructedParticle.h>

#include <sstream>

FastJetAlg::FastJetAlg(const std::string& name, ISvcLocator* svcLoc) : MultiTransformer(name, svcLoc, 
    { KeyValues("recParticleIn", {"MCParticle"}) },
    { KeyValues("jetOut", {"JetOut"}),
      KeyValues("recParticleOut", {"Constituents"}) }),
            m_jetAlgoName(""),
            m_jetAlgo(NULL),
            m_jetAlgoType(),
            m_clusterModeName(""),
            m_clusterMode( NONE ),
            m_jetRecoScheme(),
            m_strategyName(""),
            m_strategy(),
            m_requestedNumberOfJets(0),
            m_yCut(0.0),
            m_minPt(0.0),
            m_minE(0.0) {}

StatusCode FastJetAlg::initialize() {
    // ------------------ Init Strategy ------------------
    m_strategy = fastjet::Best;
    m_strategyName = "Best";
    info() << "Strategy: " << m_strategyName << endmsg;

    // ------------------ Init Reco Scheme ------------------
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
        throw GaudiException("Unknown FastJet recombination scheme! See log for more details.", name(), StatusCode::FAILURE);
    }
    info() << "recombination scheme: " << m_jetRecoSchemeName << endmsg;

    // ------------------ Init Cluster Mode ------------------
    // at least a name has to be given
    if (m_clusterModeNameAndParam.size() == 0)
        throw GaudiException("Cluster mode not specified", name(), StatusCode::FAILURE);
    // save the name of the cluster mode
    m_clusterModeName = m_clusterModeNameAndParam[0];
    m_clusterMode = NONE;
    // check the different cluster mode possibilities, and check if the number of parameters are correct
    if (m_clusterModeName.compare("Inclusive") == 0) {
        if (m_clusterModeNameAndParam.size() != 2) {
            error()  << "Wrong number of values for parameter clusteringMode 'Inclusive': missing minPt" << endmsg;
            throw GaudiException("Wrong Parameter(s) for Clustering Mode. Expected:\n <parameter name=\"clusteringMode\" type=\"StringVec\"> Inclusive <minPt> </parameter>", name(), StatusCode::FAILURE);
        }
        m_minPt = atof(m_clusterModeNameAndParam[1].c_str());
        m_clusterMode = FJ_inclusive;
    } else if (m_clusterModeName.compare("InclusiveIterativeNJets") == 0) {
        if (m_clusterModeNameAndParam.size() != 3) {
            throw GaudiException("Wrong Parameter(s) for Clustering Mode. Expected:\n <parameter name=\"clusteringMode\" type=\"StringVec\"> InclusiveIterativeNJets <NJets> <minE> </parameter>", name(), StatusCode::FAILURE);
        }
        m_requestedNumberOfJets = atoi(m_clusterModeNameAndParam[1].c_str());
        m_minE = atoi(m_clusterModeNameAndParam[2].c_str());
        m_clusterMode = OWN_inclusiveIteration;
    }  else if (m_clusterModeName.compare("ExclusiveNJets") == 0) {
        if (m_clusterModeNameAndParam.size() != 2) {
            throw GaudiException("Wrong Parameter(s) for Clustering Mode. Expected:\n <parameter name=\"clusteringMode\" type=\"StringVec\"> ExclusiveNJets <NJets> </parameter>", name(), StatusCode::FAILURE);
        }
        m_requestedNumberOfJets = atoi(m_clusterModeNameAndParam[1].c_str());
        m_clusterMode = FJ_exclusive_nJets;
    } else if (m_clusterModeName.compare("ExclusiveYCut") == 0) {
        if (m_clusterModeNameAndParam.size() != 2) {
            throw GaudiException("Wrong Parameter(s) for Clustering Mode. Expected:\n <parameter name=\"clusteringMode\" type=\"StringVec\"> ExclusiveYCut <yCut> </parameter>", name(), StatusCode::FAILURE);
        }
        m_yCut = atof(m_clusterModeNameAndParam[1].c_str());
        m_clusterMode = FJ_exclusive_yCut;
    } else {
        throw GaudiException("Unknown cluster mode.", name(), StatusCode::FAILURE);
    }
    info() << "Cluster mode: " << m_clusterMode << endmsg;

    // ------------------ Init Jet Algorithm ------------------
    // sanity check
    if (m_jetAlgoNameAndParams.size() < 1)
        throw GaudiException("No Jet algorithm provided!", name(), StatusCode::FAILURE);
    // save the name
    m_jetAlgoName = m_jetAlgoNameAndParams[0];
    // check all supported algorithms and create the appropriate FJ instance
    m_jetAlgo = NULL;
    info() << "Algorithms: ";	// the isJetAlgo function will write to streamlog_out(MESSAGE), so that we get a list of available algorithms in the log
    // example: kt_algorithm, needs 1 parameter, supports inclusive, inclusiveIterative, exlusiveNJets and exlusiveYCut clustering
    if (isJetAlgo("kt_algorithm", 1, FJ_inclusive | FJ_exclusive_nJets | FJ_exclusive_yCut | OWN_inclusiveIteration)) {
        m_jetAlgoType = fastjet::kt_algorithm;
        m_jetAlgo = new fastjet::JetDefinition(
        m_jetAlgoType, atof(m_jetAlgoNameAndParams[1].c_str()), m_jetRecoScheme, m_strategy);
    }
    if (isJetAlgo("cambridge_algorithm", 1, FJ_inclusive | FJ_exclusive_nJets | FJ_exclusive_yCut | OWN_inclusiveIteration)) {
        m_jetAlgoType = fastjet::cambridge_algorithm;
        m_jetAlgo = new fastjet::JetDefinition(
        m_jetAlgoType, atof(m_jetAlgoNameAndParams[1].c_str()), m_jetRecoScheme, m_strategy);
    }
    if (isJetAlgo("antikt_algorithm", 1, FJ_inclusive | OWN_inclusiveIteration)) {
        m_jetAlgoType = fastjet::antikt_algorithm;
        m_jetAlgo = new fastjet::JetDefinition(
        m_jetAlgoType, atof(m_jetAlgoNameAndParams[1].c_str()), m_jetRecoScheme, m_strategy);
    }
    if (isJetAlgo("genkt_algorithm", 2, FJ_inclusive | OWN_inclusiveIteration | FJ_exclusive_nJets | FJ_exclusive_yCut)) {
        m_jetAlgoType = fastjet::genkt_algorithm;
        m_jetAlgo = new fastjet::JetDefinition(
        m_jetAlgoType, atof(m_jetAlgoNameAndParams[1].c_str()), atof(m_jetAlgoNameAndParams[2].c_str()), m_jetRecoScheme, m_strategy);
    }
    if (isJetAlgo("cambridge_for_passive_algorithm", 1, FJ_inclusive | OWN_inclusiveIteration | FJ_exclusive_nJets | FJ_exclusive_yCut)) {
        m_jetAlgoType = fastjet::cambridge_for_passive_algorithm;
        m_jetAlgo = new fastjet::JetDefinition(
        m_jetAlgoType, atof(m_jetAlgoNameAndParams[1].c_str()), m_jetRecoScheme, m_strategy);
    }
    if (isJetAlgo("genkt_for_passive_algorithm", 1, FJ_inclusive | OWN_inclusiveIteration)) {
        m_jetAlgoType = fastjet::genkt_for_passive_algorithm;
        m_jetAlgo = new fastjet::JetDefinition(
        m_jetAlgoType, atof(m_jetAlgoNameAndParams[1].c_str()), m_jetRecoScheme, m_strategy);
    }
    if (isJetAlgo("ee_kt_algorithm", 0, FJ_exclusive_nJets | FJ_exclusive_yCut)) {
        m_jetAlgoType = fastjet::ee_kt_algorithm;
        m_jetAlgo = new fastjet::JetDefinition(
        m_jetAlgoType, m_jetRecoScheme, m_strategy);
    }
    // backwards compatibility for using 1 parameter only assuming exponent to be 1.
    bool commentOnAlgo = false;
    if ((m_jetAlgoNameAndParams[0]=="ee_genkt_algorithm") && ((int)m_jetAlgoNameAndParams.size() == 2)){
        m_jetAlgoNameAndParams.value().push_back("1.");
        commentOnAlgo = true;
    }
    if (isJetAlgo("ee_genkt_algorithm", 2, FJ_inclusive | FJ_exclusive_nJets | FJ_exclusive_yCut)) {
        m_jetAlgoType = fastjet::ee_genkt_algorithm;
        m_jetAlgo = new fastjet::JetDefinition(
        m_jetAlgoType, atof(m_jetAlgoNameAndParams[1].c_str()), atof(m_jetAlgoNameAndParams[2].c_str()), m_jetRecoScheme, m_strategy);
    }
    if (isJetAlgo("SISConePlugin", 2, FJ_inclusive | OWN_inclusiveIteration)) {
        fastjet::SISConePlugin* pl;
        pl = new fastjet::SISConePlugin(
        atof(m_jetAlgoNameAndParams[1].c_str()),
        atof(m_jetAlgoNameAndParams[2].c_str())
        );
        m_jetAlgo = new fastjet::JetDefinition(pl);
        m_jetAlgo->delete_plugin_when_unused();
    }
    if (isJetAlgo("SISConeSphericalPlugin", 2, FJ_inclusive | OWN_inclusiveIteration)) {
        fastjet::SISConeSphericalPlugin* pl;
        pl = new fastjet::SISConeSphericalPlugin(
            atof(m_jetAlgoNameAndParams[1].c_str()),
            atof(m_jetAlgoNameAndParams[2].c_str())
        );
        m_jetAlgo = new fastjet::JetDefinition(pl);
        m_jetAlgo->delete_plugin_when_unused();
    }
    if (isJetAlgo("ValenciaPlugin", 3, FJ_exclusive_nJets | FJ_exclusive_yCut)) {
        fastjet::contrib::ValenciaPlugin* pl;
        pl = new fastjet::contrib::ValenciaPlugin(
            atof(m_jetAlgoNameAndParams[1].c_str()),  // R value
            atof(m_jetAlgoNameAndParams[2].c_str()),  // beta value
            atof(m_jetAlgoNameAndParams[3].c_str())   // gamma value
        );
        m_jetAlgo = new fastjet::JetDefinition(pl);
        m_jetAlgo->delete_plugin_when_unused();
    }
    info() << endmsg; // end of list of available algorithms
    //ee_genkt_algorithm
    if (commentOnAlgo) {info() << "When only 1 parameter is provided for ee_genkt_algorithm it is assumed to be R, and the exponent p is assumed to be equal to 1" << endmsg;}
    if (!m_jetAlgo) {
        error() << "The given algorithm \"" << m_jetAlgoName << "\" is unknown to me!" << endmsg;
        throw GaudiException("Unknown FastJet algorithm.", name(), StatusCode::FAILURE);
    }
    info() << "jet algorithm: " << m_jetAlgo->description() << endmsg;

    return StatusCode::SUCCESS;
}

std::tuple<edm4hep::ReconstructedParticleCollection, edm4hep::ReconstructedParticleCollection> FastJetAlg::operator()(
        const edm4hep::ReconstructedParticleCollection& inputCollection) const {
    edm4hep::ReconstructedParticleCollection outputCollection;
    outputCollection->setSubsetCollection(true);
    edm4hep::ReconstructedParticleCollection jetCollection;

    PseudoJetList jets;
    PseudoJetList pjList;
    for (int i = 0; i < inputCollection.size(); ++i) {
        edm4hep::ReconstructedParticle par = inputCollection.at(i);
        pjList.push_back( fastjet::PseudoJet(
            par.getMomentum().x,
            par.getMomentum().y,
            par.getMomentum().z,
            par.getEnergy() ) );
        pjList.back().set_user_index(i);	// save the id of this recParticle
    }
    
    fastjet::ClusterSequence cs = fastjet::ClusterSequence(pjList, *m_jetAlgo);
    try {
        if (m_clusterMode == FJ_inclusive) {
            jets = cs.inclusive_jets(m_minPt);
        } else if (m_clusterMode == FJ_exclusive_yCut) {
            jets = cs.exclusive_jets_ycut(m_yCut);
        } else if (m_clusterMode == FJ_exclusive_nJets) {
            // sanity check: if we have not enough particles, FJ will cause an assert
            if (inputCollection.size() < (int)m_requestedNumberOfJets) {
                warning() << "Not enough elements in the input collection to create " << m_requestedNumberOfJets << " jets." << endmsg;
                throw SkippedFixedNrJetException();
            } else {
                jets = cs.exclusive_jets((int)(m_requestedNumberOfJets));
            }
        } else if (m_clusterMode == OWN_inclusiveIteration) {
            // sanity check: if we have not enough particles, FJ will cause an assert
            if (inputCollection.size() < (int)m_requestedNumberOfJets) {
                warning() << "Not enough elements in the input collection to create " << m_requestedNumberOfJets << " jets." << endmsg;
                throw SkippedFixedNrJetException();
            } else {
                // lets do a iterative procedure until we found the correct number of jets
                // for that we will do inclusive clustering, modifying the R parameter in some kind of minimization
                // this is based on Marco Battaglia's FastJetClustering
                double R = M_PI_4;	// maximum of R is Pi/2, minimum is 0. So we start hat Pi/4
                double RDiff = R / 2;	// the step size we modify the R parameter at each iteration. Its size for the n-th step is R/(2n), i.e. starts with R/2
                PseudoJetList jets_it;
                unsigned nJets;
                int iIter = 0;	// nr of current iteration
                // these variables are only used if the SisCone(Spherical)Plugin is selected
                // This is necessary, as these are plugins and hence use a different constructor than
                // the built in fastjet algorithms
                // here we save pointer to the plugins, so that if they are created, we can delete them again after usage
                fastjet::SISConePlugin* pluginSisCone = NULL;
                fastjet::SISConeSphericalPlugin* pluginSisConeSph = NULL;
                // check if we use the siscones
                bool useSisCone = m_jetAlgoName.compare("SISConePlugin") == 0;
                bool useSisConeSph = m_jetAlgoName.compare("SISConeSphericalPlugin") == 0;
                // save the 2nd parameter of the SisCones
                double sisConeOverlapThreshold = 0;
                if (useSisCone || useSisConeSph)
                sisConeOverlapThreshold = atof(m_jetAlgoNameAndParams[2].c_str());
                // do a maximum of N iterations
                for (iIter=0; iIter<ITERATIVE_INCLUSIVE_MAX_ITERATIONS; iIter++) {
                    // do the clustering for this value of R. For this we need to re-initialize the JetDefinition, as it takes the R parameter
                    fastjet::JetDefinition* jetDefinition = NULL;
                    // unfortunately SisCone(spherical) are being initialized differently, so we have to check for this
                    if (useSisCone) {
                        pluginSisCone = new fastjet::SISConePlugin(R, sisConeOverlapThreshold);
                        jetDefinition = new fastjet::JetDefinition(pluginSisCone);
                    } else if (useSisConeSph) {
                        pluginSisConeSph = new fastjet::SISConeSphericalPlugin(R, sisConeOverlapThreshold);
                        jetDefinition = new fastjet::JetDefinition(pluginSisConeSph);
                    } else {
                        jetDefinition = new fastjet::JetDefinition(m_jetAlgoType, R, m_jetRecoScheme, m_strategy);
                    }
                    // now we can finally create the cluster sequence
                    fastjet::ClusterSequence cs_it(pjList, *jetDefinition);
                    jets_it = cs_it.inclusive_jets(0);	// no pt cut, we will do an energy cut
                    jets.clear();
                    // count the number of jets above threshold
                    nJets = 0;
                    for (unsigned j=0; j<jets_it.size(); j++)
                        if (jets_it[j].E() > m_minE)
                            jets.push_back(jets_it[j]);
                    nJets = jets.size();
                    debug() << iIter << " " << R << " " << jets_it.size() << " " << nJets << endmsg;
                    if (nJets == m_requestedNumberOfJets) { // if the number of jets is correct: success!
                        delete pluginSisCone; pluginSisCone = NULL;
                        delete pluginSisConeSph; pluginSisConeSph = NULL;
                        delete jetDefinition;
                        break;
                    } else if (nJets < m_requestedNumberOfJets) {
                        // if number of jets is too small: we need a smaller Radius per jet (so
                        // that we get more jets)
                        R -= RDiff;
                    } else if (nJets > m_requestedNumberOfJets) {	// if the number of jets is too
                        // high: increase the Radius
                        R += RDiff;
                    }
                    RDiff /= 2;
                    // clean up
                    delete pluginSisCone; pluginSisCone = NULL;
                    delete pluginSisConeSph; pluginSisConeSph = NULL;
                    delete jetDefinition;
                }
                if (iIter == ITERATIVE_INCLUSIVE_MAX_ITERATIONS) {
                    warning() << "Maximum number of iterations reached. Canceling" << endmsg;
                    throw SkippedMaxIterationException( jets );
                    // Currently we will return the latest results, independent if the number is actually matched
                    // jets.clear();
                }
            }
        }
    } catch(const SkippedFixedNrJetException& e ) {
    } catch(const SkippedMaxIterationException& e ) {
        jets = e.m_jets;
    }

    PseudoJetList::iterator it;
    for (it=jets.begin(); it != jets.end(); it++) {
        // create a reconstructed particle for this jet, and add all the containing particles to it
        edm4hep::MutableReconstructedParticle rec = jetCollection->create();
        rec.setEnergy( (*it).E() );
        rec.setMass( (*it).m() );
        edm4hep::Vector3f mom((*it).px(), (*it).py(), (*it).pz());
        rec.setMomentum(mom);
        for (unsigned int n = 0; n < cs.constituents(*it).size(); ++n) {
            rec.addToParticles(inputCollection.at(cs.constituents(*it)[n].user_index()));
        }

        // add jet constituents to output collection
        for (unsigned int n = 0; n < cs.constituents(*it).size(); ++n) {
            edm4hep::ReconstructedParticle p = inputCollection.at((cs.constituents(*it))[n].user_index());
            outputCollection.push_back( p );
        }
    }

    return std::make_tuple(std::move(outputCollection), std::move(jetCollection));
}