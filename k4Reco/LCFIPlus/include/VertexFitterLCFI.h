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

#ifndef K4RECO_VERTEXFITTERLCFI_H
#define K4RECO_VERTEXFITTERLCFI_H

#include "lcfiplus.h"

#include <inc/event.h>
#include <inc/track.h>
#include <util/inc/helixrep.h>
#include <util/inc/memorymanager.h>
#include <util/inc/util.h>
#include <zvtop/include/VertexFitterKalman.h>
#include <zvtop/include/interactionpoint.h>

#include <list>
#include <map>
#include <vector>

namespace lcfiplus {

namespace detail {

  class LCFIVertexEvent final {
  public:
    LCFIVertexEvent() : m_event(new vertex_lcfi::Event()) {
      vertex_lcfi::MemoryManager<vertex_lcfi::Event>::Event()->registerObject(m_event);
    }

    ~LCFIVertexEvent() { vertex_lcfi::MetaMemoryManager::Event()->delAllObjects(); }

    vertex_lcfi::Track* makeTrack(const Track* track) const {
      vertex_lcfi::util::HelixRep helix;
      helix.d0() = track->getD0();
      helix.z0() = track->getZ0();
      helix.phi() = track->getPhi();
      helix.invR() = track->getOmega();
      helix.tanLambda() = track->getTanLambda();

      vertex_lcfi::util::Vector3 momentum;
      momentum.x() = track->Px();
      momentum.y() = track->Py();
      momentum.z() = track->Pz();

      vertex_lcfi::util::SymMatrix5x5 covariance;
      const double* cov = track->getCovMatrix();
      // LCFIPlus's classic covariance ordering is the ordering expected by the
      // original LCFIPlus LcfiInterface adapter.
      covariance(0, 0) = cov[0];
      covariance(3, 0) = cov[6];
      covariance(3, 3) = cov[9];
      covariance(1, 0) = cov[1];
      covariance(1, 1) = cov[2];
      covariance(0, 2) = cov[3];
      covariance(1, 2) = cov[4];
      covariance(2, 2) = cov[5];
      covariance(1, 3) = cov[7];
      covariance(2, 3) = cov[8];
      covariance(0, 4) = cov[10];
      covariance(1, 4) = cov[11];
      covariance(2, 4) = cov[12];
      covariance(3, 4) = cov[13];
      covariance(4, 4) = cov[14];

      const double charge = helix.invR() > 0.0 ? 1.0 : -1.0;
      auto* result =
          new vertex_lcfi::Track(m_event, helix, momentum, charge, covariance, {}, const_cast<Track*>(track));
      vertex_lcfi::MemoryManager<vertex_lcfi::Track>::Event()->registerObject(result);
      return result;
    }

  private:
    vertex_lcfi::Event* m_event;
  };

} // namespace detail

template <class Iterator>
class VertexFitterLCFI {
public:
  Vertex* operator()(Iterator tracksBegin, Iterator tracksEnd, Vertex* pointConstraint = nullptr) const {
    detail::LCFIVertexEvent event;
    std::vector<vertex_lcfi::TrackState*> trackStates;
    std::map<vertex_lcfi::TrackState*, const Track*> trackMap;

    for (auto track = tracksBegin; track != tracksEnd; ++track) {
      auto* lcfiTrack = event.makeTrack(*track);
      auto* state = lcfiTrack->makeState();
      trackStates.push_back(state);
      trackMap.emplace(state, *track);
    }

    vertex_lcfi::ZVTOP::InteractionPoint* constraint = nullptr;
    if (pointConstraint != nullptr) {
      vertex_lcfi::util::Vector3 position;
      position(0) = pointConstraint->getX();
      position(1) = pointConstraint->getY();
      position(2) = pointConstraint->getZ();

      vertex_lcfi::util::SymMatrix3x3 covariance;
      const double* cov = pointConstraint->getCov();
      covariance(0, 0) = cov[Vertex::xx];
      covariance(0, 1) = cov[Vertex::xy];
      covariance(0, 2) = cov[Vertex::xz];
      covariance(1, 1) = cov[Vertex::yy];
      covariance(1, 2) = cov[Vertex::yz];
      covariance(2, 2) = cov[Vertex::zz];
      constraint = new vertex_lcfi::ZVTOP::InteractionPoint(position, covariance);
    }

    vertex_lcfi::util::Vector3 position;
    vertex_lcfi::util::Matrix3x3 covariance;
    double chi2 = 0.0;
    double constraintChi2 = 0.0;
    std::map<vertex_lcfi::TrackState*, double> trackChi2;
    vertex_lcfi::ZVTOP::VertexFitterKalman{}.fitVertex(trackStates, constraint, position, covariance, chi2, trackChi2,
                                                       constraintChi2);
    if (constraint != nullptr) {
      chi2 -= constraint->chi2(position);
      delete constraint;
    }

    double vertexCovariance[6] = {covariance(0, 0), covariance(0, 1), covariance(1, 1),
                                  covariance(0, 2), covariance(1, 2), covariance(2, 2)};
    auto* vertex = new Vertex(chi2, vertex_lcfi::util::prob(chi2, trackStates.size() * 2 - 3), position(0), position(1),
                              position(2), vertexCovariance, false);
    for (auto* state : trackStates)
      vertex->add(trackMap.at(state), trackChi2.at(state));
    return vertex;
  }
};

using VertexFitterLCFI_V = VertexFitterLCFI<std::vector<const Track*>::const_iterator>;
using VertexFitterLCFI_L = VertexFitterLCFI<std::list<const Track*>::const_iterator>;

} // namespace lcfiplus

#endif
