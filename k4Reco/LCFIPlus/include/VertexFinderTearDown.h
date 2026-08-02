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
#ifndef K4RECO_VERTEXFINDERTEARDOWN_H
#define K4RECO_VERTEXFINDERTEARDOWN_H 1

#include "VertexFitterLCFI.h"
#include "VertexFitterSimple.h"
#include "lcfiplus.h"

#include <iostream>
#include <list>
#include <map>
#include <utility>
#include <vector>

namespace lcfiplus {

class SortTracksByChi2 { // decending order
public:
  bool operator()(const std::pair<const Track*, double>& p1, const std::pair<const Track*, double>& p2) {
    return p1.second > p2.second;
  }
};

// Function for recursive search of vertices using TearDown method
// std::vector<lcfiplus::Vertex*>* findTearDownVertices(const Event& evt, const Jet& jet);

// Primary Vertex finder with TearDown method
lcfiplus::Vertex* findPrimaryVertex(TrackVec& tracks, double chi2 = 9.0, bool beamspotConstraint = true,
                                    bool smearBeamspot = true);
//	lcfiplus::Vertex * findPrimaryVertex(const   std::vector<Track *> &tracks, const   std::vector<Track *>
//&beamTracks, double chi2 = 9.0);

// implementation of TearDown method
template <template <class T, class Allocator = std::allocator<T>> class Container = std::vector,
          template <class Iterator> class VertexFitter = VertexFitterLCFI>
class VertexFinderTearDown {
public:
  Vertex* operator()(const Container<const Track*>& tracks, const Container<const Track*>* fixedTracks = 0,
                     double chiSquareThreshold = 9.0, Container<const Track*>* residual = 0,
                     Vertex* pointConstraint = nullptr) {
    // copy tracks into a list
    std::list<const Track*> trackList;
    trackList.resize(tracks.size() + (fixedTracks ? fixedTracks->size() : 0));
    std::list<const Track*>::iterator listIt = copy(tracks.begin(), tracks.end(), trackList.begin());
    if (fixedTracks) {
      copy(fixedTracks->begin(), fixedTracks->end(), listIt);
    }

    Vertex* resultVertex = nullptr;
    double worstChi2;
    while (trackList.size() >= 2) {
      resultVertex =
          VertexFitter<std::list<const Track*>::iterator>()(trackList.begin(), trackList.end(), pointConstraint);
      const Track* worstTrack = resultVertex->getWorstTrack();
      if (fixedTracks && find(fixedTracks->begin(), fixedTracks->end(), worstTrack) != fixedTracks->end()) {
        // sort Chi2Tracks
        std::vector<std::pair<const Track*, double>> vpair;
        const std::map<const Track*, double>& mpair = resultVertex->getTracksChi2Map();
        vpair.resize(mpair.size());
        partial_sort_copy(mpair.begin(), mpair.end(), vpair.begin(), vpair.end(), SortTracksByChi2());

        unsigned int nworst = 1;
        do {
          worstTrack = vpair[nworst++].first;
        } while (nworst < vpair.size() &&
                 find(fixedTracks->begin(), fixedTracks->end(), worstTrack) != fixedTracks->end());

        std::cout << "The worst track is fixed, " << nworst << "th track will be removed." << std::endl;
      }

      worstChi2 = resultVertex->getChi2Track(worstTrack);
      if (worstChi2 > chiSquareThreshold) {
        trackList.remove(worstTrack);
        if (residual)
          residual->push_back(worstTrack);
        /*
        cout << "Track removed: chi2 = " << worstChi2 << ", vpos = ("
                                << resultVertex->getX() << ", " << resultVertex->getY() << ", " << resultVertex->getZ()
        << ")" <<  endl;
         //*/
        delete resultVertex;
        resultVertex = nullptr;
      } else
        break; // all tracks below chi2 threshold
    }

    return resultVertex;
  }
};
} // namespace lcfiplus

#endif
