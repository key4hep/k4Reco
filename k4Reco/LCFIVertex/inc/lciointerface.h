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
#ifndef LCFILCIOINTERFACE_h
#define LCFILCIOINTERFACE_h 1

#include "lcio.h"
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Track.h>
#include <EVENT/Vertex.h>
#include <IMPL/AccessChecked.h>
#include <IMPL/ReconstructedParticleImpl.h>

#include <inc/event.h>
#include <inc/jet.h>
#include <inc/track.h>
#include <inc/vertex.h>

namespace vertex_lcfi {
using namespace lcio;

vertex_lcfi::Track* trackFromLCIORP(Event* MyEvent, lcio::ReconstructedParticle* RP);
vertex_lcfi::Jet* jetFromLCIORP(Event* MyEvent, lcio::ReconstructedParticle* RP);
ReconstructedParticle* addDecayChainToLCIOEvent(LCEvent* MyLCIOEvent, DecayChain* MyDecayChain,
                                                std::string VertexCollectionName, std::string TrackRPCollectionName,
                                                bool StoreTrackChiSquareds = false);
DecayChain* decayChainFromLCIORP(Jet* MyJet, ReconstructedParticle* DecayChainRP);
lcio::Vertex* vertexFromLCFIVertex(vertex_lcfi::Vertex* MyLCFIVertex);
vertex_lcfi::Vertex* vertexFromLCIOVertex(lcio::Vertex* LCIOVertex, Event* MyEvent);

class ReconstructedParticleLCFI : private IMPL::ReconstructedParticleImpl {
public:
  void removeParticle(EVENT::ReconstructedParticle* particle);
  void wipePIDs();
  void copyPIDsFrom(ReconstructedParticle* InputRP);
  void makeWritable();
};

} // namespace vertex_lcfi

#endif
