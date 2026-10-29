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
#include "VertexFinderTearDown.h"

#include "algoEtc.h"
#include "lcfiplus.h"

// std::vector<lcfiplus::Vertex*>* lcfiplus::findTearDownVertices(const Event& evt, const Jet& jet) {
//   double chi2 = 9.0;
//   bool verbose(false);
//   std::vector<lcfiplus::Vertex*>* pvertices;
//   pvertices = new std::vector<lcfiplus::Vertex*>;

//   // copy tracks in the jet into a list
//   // const std::vector<Track *> &v = jets[nj]->getTracks();
//   const std::vector<const Track*>& v = jet.getTracks();
//   std::list<const Track*> tracksInJet;
//   tracksInJet.resize(v.size());
//   copy(v.begin(), v.end(), tracksInJet.begin());

//   while (tracksInJet.size() >= 2) {
//     std::list<const Track*> residuals;
//     // TODO: This is using the default of VertexFinderTearDown, which is VertexFitterLCFI originally
//     lcfiplus::Vertex* secvtx = lcfiplus::VertexFinderTearDown<std::list>()(tracksInJet, nullptr, chi2, &residuals);
//     if (!secvtx)
//       break;

//     pvertices->push_back(secvtx);

//     if (verbose)
//       std::cout << "    Secondary vertex found! pos = (" << secvtx->getX() << "," << secvtx->getY() << ","
//                 << secvtx->getZ() << "), chi2 = " << secvtx->getChi2() << std::endl;
//     for (unsigned int i = 0; i < secvtx->getTracks().size(); i++) {
//       const Track* tr = secvtx->getTracks()[i];
//       if (verbose) {
//         std::cout << "        Track #" << i << ": p = (" << tr->Px() << "," << tr->Py() << "," << tr->Pz()
//                   << "), chi2 = " << secvtx->getChi2Track(tr)
//                   << ", mcFlavorType = " << evt.getMCParticle(tr)->getFlavorTagCategory() << "\n";
//       }
//     }
//     tracksInJet = residuals;
//   }

//   return pvertices;
// }

lcfiplus::Vertex* lcfiplus::findPrimaryVertex(TrackVec& tracks, double chi2, bool beamspotConstraint,
                                              bool smearBeamspot) {
  Vertex* ip = nullptr;
  if (beamspotConstraint)
    algoEtc::makeBeamVertex(ip, smearBeamspot);

  Vertex* ret = VertexFinderTearDown<std::vector, VertexFitterSimple>()(tracks, nullptr, chi2, nullptr, ip);
  if (ret)
    ret->setPrimary(true);
  else
    return ip; // FIXME: this is safety procedure in case primary vertex is not found; need to confirm this is good
               // behavior

  delete ip;
  return ret;
}

// lcfiplus::Vertex * lcfiplus::findPrimaryVertex(const std::vector<Track *> &tracks, const std::vector<Track *>
// &beamTracks, double chi2)
// {
// 	// make point constraint with beam tracks for initial condition
// 	Vertex * ipv = VertexFitterSimple_V()(beamTracks.begin(), beamTracks.end());
//
// 	// fit with ipv, without beam tracks
// 	Vertex * ret =  VertexFinderTearDown<vector, VertexFitterSimple>()(tracks, 0, chi2, 0, ipv);
// 	delete ipv;
// 	return ret;
// }
