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
#include "../inc/event.h"
#include "../inc/vertex.h"
#include "../util/inc/memorymanager.h"
#include "../util/inc/vector3.h"

namespace vertex_lcfi {
using namespace util;

Event::Event() {
  SymMatrix3x3 IPError;
  IPError(0, 0) = 10.0 / 1000.0;
  IPError(1, 1) = 10.0 / 1000.0;
  IPError(2, 2) = 10.0 / 1000.0;
  _IPVertex = new Vertex((this), std::vector<Track*>(), Vector3(0, 0, 0), IPError, true, 0, 1);
  MemoryManager<Vertex>::Event()->registerObject(_IPVertex);
}

Event::Event(const Vector3& Position, const SymMatrix3x3& Error) {
  _IPVertex = new Vertex(const_cast<Event*>(this), std::vector<Track*>(), Position, Error, true, 0, 1);
  MemoryManager<Vertex>::Event()->registerObject(_IPVertex);
}

Event::Event(Vertex* ipVertex) : _IPVertex(ipVertex) {}

void Event::addJet(Jet* AJet) { _Jets.push_back(AJet); }

const std::vector<Jet*>& Event::jets() const { return _Jets; }

const std::vector<Track*>& Event::tracks() const { return _Tracks; }

void Event::addTrack(Track* AddTrack) { _Tracks.push_back(AddTrack); }

const Vector3& Event::interactionPoint() const { return _IPVertex->position(); }

const SymMatrix3x3& Event::interactionPointError() const { return _IPVertex->positionError(); }

Vertex* Event::ipVertex() const { return _IPVertex; }
} // namespace vertex_lcfi
