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
#include "LcfiplusAlgorithm.h"

#include <edm4hep/MCParticleCollection.h>
#include <edm4hep/RecoMCParticleLinkCollection.h>
#include <edm4hep/ReconstructedParticleCollection.h>

#include <string>

// geometry utils from MarlinUtil
// #include <GeometryUtil.h>

// static object initialization
// LCIOStorer* LcfiplusAlgorithm::_lcio = 0;

LcfiplusAlgorithm::LcfiplusAlgorithm(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(name, svcLoc,
                       {
                           KeyValues("PFOCollection", {""}),
                           KeyValues("MCPCollection", {""}),
                           KeyValues("MCPFORelation", {""}),
                       },
                       {
                           KeyValues("TrackerHitCollectionName", {"VTXTrackerHits"}),
                       }) {

  _inInit = false;
  _printPeriod = 0;
}

StatusCode LcfiplusAlgorithm::initialize() {
  // try {

  //   debug() << "   init called  " << endmsg;

  //   // usually a good idea to
  //   printParameters();

  //   std::shared_ptr<StringParameters> parameter = parameters();
  //   /*
  //               // obtain algorithm name
  //               if(!parameter->isParameterSet("algorithm")){
  //                       info() << "Lcfiplus algorithm not set. run nothing." << endl;
  //                       return;
  //               }
  //               vector<string> algos;
  //               parameter->getStringVals("algorithm", algos);

  //   */

  //   // algorithm check
  //   if (_algonames.size() == 0) {
  //     info() << "No algorithms given to run. run nothing." << endl;
  //     return;
  //   }

  //   // set globals
  //   if (_magneticField != 0) {

  //     Globals::Instance()->setBField(_magneticField);

  //     streamlog_out(WARNING) << " overriding the B-field at the origin (" << MarlinUtil::getBzAtOrigin()
  //                            << ") with parameter MagneticField : " << _magneticField << endmsg;
  //   } else {

  //     Globals::Instance()->setBField(MarlinUtil::getBzAtOrigin());
  //   }

  //   Globals::Instance()->setBeamSizeX(_beamSizeX);
  //   Globals::Instance()->setBeamSizeY(_beamSizeY);
  //   Globals::Instance()->setBeamSizeZ(_beamSizeZ);

  //   // register observer
  //   Event::Instance()->RegisterObserver(this);

  //   // conversion StringParameters -> Parameters
  //   _param = new Parameters;

  //   StringVec keys;
  //   parameter->getStringKeys(keys);

  //   for (unsigned int i = 0; i < keys.size(); i++) {
  //     vector<string> vals;
  //     parameter->getStringVals(keys[i], vals);
  //     _param->add(keys[i].c_str(), vals);
  //   }

  //   // initialize LCIOStorer
  //   if (!_lcio) {
  //     _lcio = new LCIOStorer(0, 0, true, false, 0); // no file
  //     _lcio->setReadSubdetectorEnergies(_readSubdetectorEnergies);
  //     _lcio->setTrackHitOrdering(_trackHitOrdering);
  //     _lcio->setUpdateVertexRPDaughters(_updateVertexRPDaughters);
  //     _lcio->setIgnoreLackOfVertexRP(_ignoreLackOfVertexRP);
  //     _lcio->setParticleIDAlgorithmName(_pidAlgoName.c_str());

  //     _lcioowner = true;
  //   } else {
  //     _lcioowner = false;

  //     if ((_lcio->getReadSubdetectorEnergies() != _readSubdetectorEnergies) ||
  //         (_lcio->getTrackHitOrdering() != _trackHitOrdering) ||
  //         (_lcio->getUpdateVertexRPDaughters() != _updateVertexRPDaughters) ||
  //         (_lcio->getIgnoreLackOfVertexRP() != _ignoreLackOfVertexRP)) {
  //       throw(lcfiplus::Exception(
  //           "Global parameters do not match to previous processors: specify the same for all LcfiplusProcessors."));
  //     }
  //   }

  //   // load basic collection
  //   if (_useMcp)
  //     _lcio->InitMCPPFOCollections(_pfoCollectionName.c_str(), _mcpCollectionName.c_str(),
  //     _mcpfoRelationName.c_str());
  //   else
  //     _lcio->InitPFOCollections(_pfoCollectionName.c_str());

  //   _inInit = true;

  //   // make algorithm classes and pass param to init
  //   for (unsigned int i = 0; i < _algonames.size(); i++) {
  //     string s = "lcfiplus::";
  //     s += _algonames[i];
  //     TClass* cl = TClass::GetClass(s.c_str());
  //     if (!cl || !cl->InheritsFrom("lcfiplus::Algorithm")) {
  //       info() << "Algorithm " << _algonames[i] << " is not valid. skip." << endl;
  //       continue;
  //     }
  //     Algorithm* newalgo = (Algorithm*)cl->New();
  //     if (!newalgo) {
  //       info() << "Initialization failed!." << endl;
  //       break;
  //     }
  //     _algos.push_back(newalgo);
  //     info() << "Algorithm " << _algonames[i] << " is being initialized." << endl;
  //     newalgo->init(_param);

  //     info() << "Algorithm " << _algonames[i] << " successfully initialized." << endl;
  //   }

  //   _inInit = false;

  //   // printing EventStore collections
  //   Event::Instance()->Print();

  //   _nRun = 0;
  //   _nEvt = 0;

  // } catch (lcfiplus::Exception& e) {
  //   error() << e.what() << endl;
  //   throw(marlin::StopProcessingException(this));
  // } catch (lcfiplus::Exception*& e) {
  //   error() << e->what() << endl;
  //   throw(marlin::StopProcessingException(this));
  // }
  return StatusCode::SUCCESS;
}

std::tuple<edm4hep::ReconstructedParticleCollection>
LcfiplusAlgorithm::operator()(const edm4hep::ReconstructedParticleCollection&, const edm4hep::MCParticleCollection&,
                              const edm4hep::RecoMCParticleLinkCollection&) const {
  try {

    // set LCEvent
    if (_lcioowner)
      _lcio->SetEvent(evt);

    // set deafult mcp/track/netural
    Event::Instance()->setDefaultTracks(_pfoCollectionName.c_str());
    Event::Instance()->setDefaultNeutrals(_pfoCollectionName.c_str());
    Event::Instance()->setDefaultMCParticles(_mcpCollectionName.c_str());

    std::ranges::for_each(_algos, [](auto& algo) { algo->process(); });

    // convert persistent collections
    //		_lcio->AutoConvert();
    for (unsigned int i = 0; i < _vertexColNamesToWrite.size(); i++) {
      bool persist = _vertexColNamesToWriteFlags[i] & EventStore::PERSIST;
      if (persist)
        _lcio->WriteVertices(_vertexColNamesToWrite[i].c_str());
    }
    for (unsigned int i = 0; i < _jetColNamesToWrite.size(); i++) {
      bool persist = _jetColNamesToWriteFlags[i] & EventStore::PERSIST;
      bool writeVertex = _jetColNamesToWriteFlags[i] & EventStore::JET_WRITE_VERTEX;
      if (persist)
        _lcio->WriteJets(_jetColNamesToWrite[i].c_str(), 0, writeVertex);
    }

  } catch (lcfiplus::Exception& e) {
    error() << e.what() << endl;
    throw(marlin::StopProcessingException(this));
  }
}

StatusCode LcfiplusAlgorithm::finalize() {
  // for (unsigned int i = 0; i < _algos.size(); i++) {
  //   _algos[i]->end();
  //   delete _algos[i];
  // }
  // _algos.clear();

  // Event::Instance()->ClearObjects();

  // std::cout << "LcfiplusAlgorithm::end()  " << name() << " processed " << _nEvt << " events in " << _nRun << " runs "
  //           << endmsg;

  // delete _param;
  return StatusCode::SUCCESS;
}

// TODO: What to do with this?
// void LcfiplusAlgorithm::RegisterCallback(const char* name, const char* classname, int flags) {
//   if (!_inInit)
//     return; // no-op other than init

//   if (string("vector<lcfiplus::Vertex*>") == classname) {
//     cout << "Registering output LCIO vertex collection: " << name << endl;
//     _vertexColNamesToWrite.push_back(name);
//     _vertexColNamesToWriteFlags.push_back(flags);
//   } else if (string("vector<lcfiplus::Jet*>") == classname) {
//     cout << "Registering output LCIO jet collection: " << name << endl;
//     _jetColNamesToWrite.push_back(name);
//     _jetColNamesToWriteFlags.push_back(flags);
//   }
// }
