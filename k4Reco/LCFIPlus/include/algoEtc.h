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
#ifndef AlgoEtc_h
#define AlgoEtc_h 1

#include "lcfiplus.h"

namespace lcfiplus {
namespace algoEtc {

  // beam pseudo-tracks for primary vertex fitter
  // void makeBeamTracks(Track*& t1, Track*& t2, bool smear=true);
  void makeBeamVertex(Vertex*& vtx, bool smear = true);
  void connectVerticesToJets(const JetVec& jets, const std::vector<Vertex*>& vertices,
                             std::vector<std::vector<Vertex*>>& jetVertices,
                             std::vector<std::vector<const Track*>>& jetResidualTracks, const Vertex* ip = nullptr);
  std::vector<const Track*> extractTracks(const VertexVec& vertices);
  // double calcThrust(   std::vector<TVector3>& list, TVector3& taxis );

  bool SimpleSecMuonFinder(const Track* tr, double d0sigth, double z0sigth, double maxpos, double mudepmin,
                           double ecaldepmin, double ecaldepmax, double hcaldepmin, double hcaldepmax,
                           double maxclusterpertrackenergy = 10., const Vertex* ip = nullptr);
  // bool SimpleSecElectronFinder(const Track* tr, double d0sigth, double z0sigth, double maxpos, double emin,
  //                              double minfracecal, double minecalpertrackenergy, double maxecalpertrackenergy,
  //                              const Vertex* ip = 0);

} // namespace algoEtc
} // namespace lcfiplus

#endif
