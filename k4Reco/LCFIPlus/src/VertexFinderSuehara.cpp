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

// #include "LcfiInterface.h"
#include "VertexFinderSuehara.h"
#include "VertexFinderTearDown.h"
#include "VertexFitterSimple.h"
#include "lcfiplus.h"

#include "algoEtc.h"
// #include "algoSigProb.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <ranges>
#include <vector>

using namespace lcfiplus;
using namespace lcfiplus::VertexFinderSuehara;

// // vertex compare functions
bool VertexFinderSuehara::VertexNearer(const Vertex* vtx1, const Vertex* vtx2) {
  return vtx1->getPos().Mag() < vtx2->getPos().Mag();
}
bool VertexFinderSuehara::VertexProbLarger(const Vertex* vtx1, const Vertex* vtx2) {
  return vtx1->getProb() > vtx2->getProb();
}

void VertexFinderSuehara::GetVertexList(std::list<const Track*>& tracks, const Vertex* ip, std::vector<Vertex*>& vtx,
                                        std::vector<Vertex*>& v0vtx, VertexFinderSueharaConfig& cfg) {
  // lists of found vertices
  std::list<Vertex*> tr3list;
  std::list<Vertex*> tr2list;

  std::list<const Track*>::iterator trkit1, trkit2;
  std::list<Vertex*>::iterator vit;

  bool verbose = false;

  std::list<const Track*> v0tracks;

  // find all vertices

  int nv = 0;
  int ntr = tracks.size();
  int ntrmax = (ntr - 1) * ntr / 2;

  for (trkit1 = tracks.begin(); trkit1 != tracks.end(); trkit1++) {
    for (trkit2 = trkit1, trkit2++; trkit2 != tracks.end(); trkit2++) {
      nv++;
      if (verbose) {
        std::cerr << "Vertex producing: " << nv << "/" << ntrmax << std::endl;
      }

      // v0tracks rejection
      if (find(v0tracks.begin(), v0tracks.end(), *trkit1) != v0tracks.end())
        break;
      if (find(v0tracks.begin(), v0tracks.end(), *trkit2) != v0tracks.end())
        continue;

      // mass threshold
      TLorentzVector v1 = **trkit1;
      TLorentzVector v2 = **trkit2;

      double mass = (v1 + v2).M();

      // 110216 do not accept tracks with opposite direction
      //			if(v1.Vect().Dot(v2.Vect()) < 0.)continue;
      if (!cfg.avf && mass > std::min(v1.E(), v2.E()))
        continue;
      if (cfg.avf && 0.6 * mass > std::min(v1.E(), v2.E()))
        continue; // 0.7
      // if(cfg.avf && mass > std::max(v1.E(), v2.E())) continue;

      if (mass > cfg.massth)
        continue;

      // obtain vertex
      std::vector<const Track*> vttmp;
      vttmp.push_back(*trkit1);
      vttmp.push_back(*trkit2);

      Vertex* candidate = VertexFitterSimple_V()(vttmp.begin(), vttmp.end(), nullptr);

      double chi2 = std::max(candidate->getChi2Track(*trkit1), candidate->getChi2Track(*trkit2));
      TVector3 vpos = candidate->getPos();

      if (chi2 < cfg.chi2thV0SelTrack && !VertexSelector().passesCut(candidate, cfg.v0selTrack, ip)) {
        //   std::cout << "V0 found at ( " << vtx->getPos().x() << " " << vtx->getPos().y() << " " << vtx->getPos().z()
        //   << ") :
        // 2 tracks removed." << std::endl;
        v0tracks.push_back(*trkit1);
        v0tracks.push_back(*trkit2);
        v0vtx.push_back(candidate);
        break;
      }

      if (!VertexSelector().passesCut(candidate, cfg.v0selVertex, ip)) {
        delete candidate;
        continue;
      }

      // direction & chi2 condition
      if (vpos.Dot((v1 + v2).Vect()) > 0 && chi2 < cfg.chi2th) { // trying 3+ vertex
        if (verbose)
          std::cout << "Vertex accepted." << std::endl;
        Vertex* vtx2 = associateTracks(candidate, constVector(v0vtx), tracks, cfg);
        if (vtx2 != candidate) { // 3+ tracks
          delete candidate;
          tr3list.push_back(vtx2);
        } else { // 2 tracks
          tr2list.push_back(candidate);
        }
      } else { // bad vertex
        delete candidate;
      }
    }
  }

  // v0 rejection again
  for (vit = tr3list.begin(); vit != tr3list.end();) {
    bool deleted = false;
    for (unsigned int n = 0; n < (*vit)->getTracks().size(); n++) {
      const Track* tr = (*vit)->getTracks()[n];
      if (find(v0tracks.begin(), v0tracks.end(), tr) != v0tracks.end()) {
        // new vertex with smaller # of tracks
        std::vector<const Track*> vttmp = (*vit)->getTracks();
        std::vector<const Track*>::iterator itt = find(vttmp.begin(), vttmp.end(), tr);
        vttmp.erase(itt);

        Vertex* refittedVertex = VertexFitterSimple_V()(vttmp.begin(), vttmp.end(), 0);
        if (vttmp.size() > 2)
          tr3list.push_back(refittedVertex);
        else
          tr2list.push_back(refittedVertex);

        delete *vit;
        vit = tr3list.erase(vit);
        deleted = true;
        //   std::cout << "Vertex in tr3list removed due to v0 rejection." << std::endl;
        break;
      }
    }
    if (!deleted)
      vit++;
  }
  for (vit = tr2list.begin(); vit != tr2list.end();) {
    bool deleted = false;
    for (unsigned int n = 0; n < (*vit)->getTracks().size(); n++) {
      if (find(v0tracks.begin(), v0tracks.end(), (*vit)->getTracks()[n]) != v0tracks.end()) {
        vit = tr2list.erase(vit);
        deleted = true;
        std::cout << "Vertex in tr2distlist removed due to v0 rejection." << std::endl;
        break;
      }
    }
    if (!deleted)
      vit++;
  }

  if (verbose) {
    std::cerr << "Vertex produced. 3tr: " << tr3list.size() << ", 2tr: " << tr2list.size() << std::endl;
  }

  while (true) {
    // sort found vertices
    Vertex* curvtx = 0;
    // get tracks to save
    std::vector<const Track*> savetrk;
    if (tr3list.size() > 0) {
      //			tr3list.sort(VertexNearer);
      tr3list.sort(VertexProbLarger);
      curvtx = tr3list.front();
      tr3list.pop_front();
    } else if (tr2list.size() > 0) {
      tr2list.sort(VertexProbLarger);
      curvtx = tr2list.front();
      tr2list.pop_front();
    } else {
      // try to add savetrk to final vtx
      for (unsigned int l = 0; l < savetrk.size(); l++) {
        float okchi2 = 1.0e+6;
        std::vector<Vertex*>::iterator tmpit, okit;

        // first check whether this trk is unused
        bool uflg = true;
        for (tmpit = vtx.begin(); tmpit != vtx.end(); tmpit++) {
          if (find((*tmpit)->getTracks().begin(), (*tmpit)->getTracks().end(), savetrk[l]) !=
              (*tmpit)->getTracks().end()) {
            uflg = false;
            break;
          }
        }
        if (uflg == false)
          continue;

        // start to save omitted track
        for (tmpit = vtx.begin(); tmpit != vtx.end(); tmpit++) {
          // make track vector
          std::vector<const Track*> tmptrk;
          for (unsigned int ll = 0; ll < (*tmpit)->getTracks().size(); ll++)
            tmptrk.push_back((*tmpit)->getTracks()[ll]);
          tmptrk.push_back(savetrk[l]);
          Vertex* tmpvtx = VertexFitterSimple_V()(tmptrk.begin(), tmptrk.end(), (*tmpit), true);
          // check std::max chi2
          float maxchi2 = 0.0; // tmpvtx->getChi2Track(savetrk[ll]);   //0.0;
          for (unsigned int ll = 0; ll < tmptrk.size(); ll++)
            maxchi2 = std::max(maxchi2, (float)tmpvtx->getChi2Track(tmptrk[ll]));

          if (okchi2 > maxchi2) {
            okchi2 = maxchi2;
            okit = tmpit;
          }
          delete tmpvtx;
        }

        if (okchi2 < cfg.chi2th) { // can add this track!
          std::vector<const Track*> oktrk;
          for (unsigned int ll = 0; ll < (*okit)->getTracks().size(); ll++)
            oktrk.push_back((*okit)->getTracks()[ll]);
          oktrk.push_back(savetrk[l]);
          Vertex* okvtx = VertexFitterSimple_V()(oktrk.begin(), oktrk.end(), (*okit), true);
          vtx.erase(okit);
          vtx.push_back(okvtx);
        }
      }

      break; // all vertices gone
    }

    if (verbose)
      std::cerr << "Sort finished." << std::endl;

    // vertex selected
    // set vertexing type
    if (cfg.avf)
      curvtx->setVertexingName("AVF");
    else
      curvtx->setVertexingName("chi2");
    vtx.push_back(curvtx);
    if (verbose) {
      std::cout << "Vertex found with " << curvtx->getTracks().size() << " tracks, prob = " << curvtx->getProb()
                << ", mass = " << curvtx->getVertexMass();
      std::cout << ", pos = ( " << curvtx->getX() << ", " << curvtx->getY() << ", " << curvtx->getZ() << ")"
                << std::endl;
      std::cout << "err = ( xx: " << curvtx->getCov()[Vertex::xx] << ", yy: " << curvtx->getCov()[Vertex::yy]
                << ", zz:" << curvtx->getCov()[Vertex::zz];
      std::cout << ", xy: " << curvtx->getCov()[Vertex::xy] << ", xz: " << curvtx->getCov()[Vertex::xz]
                << ", yz:" << curvtx->getCov()[Vertex::yz] << ")" << std::endl;
    }

    // remove vertices/tracks
    std::list<Vertex*>::iterator itv;

    // tr3list
    if (tr3list.size() > 0) {
      for (itv = tr3list.begin(); itv != tr3list.end();) {
        bool deleted = false;
        Vertex* v = (*itv);
        for (unsigned int itr = 0; itr < curvtx->getTracks().size(); itr++) {
          const Track* tr = curvtx->getTracks()[itr];
          if (find(v->getTracks().begin(), v->getTracks().end(), tr) != v->getTracks().end()) {
            // new vertex with smaller # of tracks
            std::vector<const Track*> vttmp = v->getTracks();
            std::vector<const Track*>::iterator itt = find(vttmp.begin(), vttmp.end(), tr);
            vttmp.erase(itt);

            Vertex* refittedVertex = VertexFitterSimple_V()(vttmp.begin(), vttmp.end(), nullptr);
            if (vttmp.size() > 2)
              tr3list.push_back(refittedVertex);
            else
              tr2list.push_back(refittedVertex);

            // vertex removed
            delete v;
            itv = tr3list.erase(itv);
            deleted = true;
            break;
          }
        }
        if (tr3list.size() == 0)
          break;
        if (!deleted)
          itv++;
      }
    }

    // tr2list
    if (tr2list.size() > 0) {
      for (itv = tr2list.begin(); itv != tr2list.end();) {
        bool deleted = false;
        Vertex* v = (*itv);
        for (unsigned int itr = 0; itr < curvtx->getTracks().size(); itr++) {
          const Track* tr = curvtx->getTracks()[itr];
          if (find(v->getTracks().begin(), v->getTracks().end(), tr) != v->getTracks().end()) {

            // look for candidates to be saved
            if (v->getTracks()[0] == tr) {
              // const_cast<Track*>(v->getTracks()[1])->setChi2(-150.0);
              if (find(curvtx->getTracks().begin(), curvtx->getTracks().end(), v->getTracks()[1]) ==
                  curvtx->getTracks().end()) {
                bool bflg = true;
                for (unsigned int j = 0; j < savetrk.size(); j++) {
                  if (savetrk[j] == v->getTracks()[1]) {
                    bflg = false;
                    break;
                  }
                }
                if (bflg == true)
                  savetrk.push_back(v->getTracks()[1]);
              }
            } else {
              if (find(curvtx->getTracks().begin(), curvtx->getTracks().end(), v->getTracks()[0]) ==
                  curvtx->getTracks().end()) {
                bool bflg = true;
                for (unsigned int j = 0; j < savetrk.size(); j++) {
                  if (savetrk[j] == v->getTracks()[0]) {
                    bflg = false;
                    break;
                  }
                }
                if (bflg == true)
                  savetrk.push_back(v->getTracks()[0]);
              }
            }
            // end

            // vertex removed
            delete v;
            itv = tr2list.erase(itv);
            deleted = true;
            break;
          }
        }
        if (tr2list.size() == 0)
          break;
        if (!deleted)
          itv++;
      }
    }
  }
}

Vertex* VertexFinderSuehara::associateTracks(Vertex* vertex, const VertexVec& v0vtx, std::list<const Track*>& tracks,
                                             VertexFinderSueharaConfig& cfg, std::list<const Track*>* residualTracks) {
  std::vector<const Track*> vt;
  TLorentzVector vsum;

  if (residualTracks) {
    residualTracks->resize(tracks.size());
    copy(tracks.begin(), tracks.end(), residualTracks->begin());
  }

  std::vector<const Track*> v0tracks = algoEtc::extractTracks(v0vtx);

  for (unsigned int i = 0; i < vertex->getTracks().size(); i++) {
    vt.push_back(vertex->getTracks()[i]);
    vsum += *(TLorentzVector*)(vertex->getTracks()[i]);
    if (residualTracks)
      residualTracks->remove(vertex->getTracks()[i]);
  }

  Vertex* curvtx = vertex;

  std::list<const Track*>::iterator trkit1;
  do {
    const Track* tra = nullptr;
    double minchi2 = cfg.chi2th;
    for (trkit1 = tracks.begin(); trkit1 != tracks.end(); trkit1++) {
      if (find(vt.begin(), vt.end(), *trkit1) != vt.end())
        continue;

      // rejecting v0tracks
      if (find(v0tracks.begin(), v0tracks.end(), *trkit1) != v0tracks.end())
        continue;

      TLorentzVector va = **trkit1;
      double mass = (vsum + va).M();

      // 110216 do not accept tracks with opposite direction to current tracks
      //			if(va.Vect().Dot(vsum.Vect()) < 0.)continue;
      // if (!cfg.avf && mass - vsum.M() > std::min(vsum.E(), va.E()))continue;
      if (mass - vsum.M() > std::min(vsum.E(), va.E()))
        continue;

      if (cfg.avf) {
        // cal. distance between vertex and helix track
        Helix hel(*trkit1, PointBase::SECVTX);

        VertexLine line(TVector3(0., 0., 0.), vertex->getPos(), PointBase::SECVTX);
        double linedist = 0;
        TVector3 pos = hel.ClosePoint(line, &linedist);
        TVector3 dist = pos - vertex->getPos();

        // errors
        TVector3 y0 = vertex->getPos();
        TVector3 err(sqrt(vertex->getCov()[Vertex::xx]), sqrt(vertex->getCov()[Vertex::yy]),
                     sqrt(vertex->getCov()[Vertex::zz]));
        double vtxerr =
            sqrt(pow(err.X() * y0.X(), 2.0) + pow(err.Y() * y0.Y(), 2.0) + pow(err.Z() * y0.Z(), 2.0)) / y0.Mag();
        double trkerr = sqrt((*trkit1)->getCovMatrix()[tpar::d0d0] + (*trkit1)->getCovMatrix()[tpar::z0z0]);

        if (dist.Mag() > 20.0 * sqrt(vtxerr * vtxerr + trkerr * trkerr))
          continue; // kill over 20 sigma
      }

      // 110216 do not accept tracks with opposite direction to vertex
      if (va.Vect().Dot(curvtx->getPos()) < 0.)
        continue;

      if (mass > cfg.massth)
        continue;

      vt.push_back(*trkit1);
      Vertex* vtx = VertexFitterSimple_V()(vt.begin(), vt.end(), curvtx, true);
      size_t ntr = vtx->getTracks().size();

      double maxchi2 = 0;
      for (size_t i = 0; i < ntr; i++) {
        double chi2 = vtx->getChi2Track(vtx->getTracks()[i]);
        if (maxchi2 < chi2)
          maxchi2 = chi2;
      }
      if (minchi2 > maxchi2) {
        minchi2 = maxchi2;
        tra = *trkit1;
        if (curvtx != vertex)
          delete curvtx;
        curvtx = vtx;
      } else {
        delete vtx;
      }
      vt.pop_back();
    }
    if (!tra)
      return curvtx;

    vsum += *(TLorentzVector*)tra;

    vt.push_back(tra);
    if (residualTracks)
      residualTracks->remove(tra);
  } while (tracks.size() > 0);

  return curvtx;
}

void VertexFinderSuehara::associateIPTracks(std::vector<Vertex*>& vertices, Vertex* ip,
                                            VertexFinderSueharaConfig& cfg) {
  // Keep the EDM4hep-derived primary vertex owned by the caller. The legacy
  // implementation replaced its local primary-vertex pointer while iterating.
  auto currentIP = std::make_unique<Vertex>(*ip);
  Vertex* beamVertex;
  algoEtc::makeBeamVertex(beamVertex, cfg.beamspotSmearing);

  for (auto*& vertex : vertices) {
    if (vertex->getPos().Mag() < cfg.minimumdistIP)
      continue;

    TLorentzVector momentum;
    for (const auto* track : vertex->getTracks())
      momentum += *track;

    std::vector<const Track*> retainedIPTracks;
    for (const auto* track : currentIP->getTracks()) {
      const TLorentzVector trackMomentum = *track;
      if ((momentum + trackMomentum).M() - momentum.M() > std::min(momentum.E(), trackMomentum.E()) ||
          trackMomentum.Vect().Dot(vertex->getPos()) < 0.) {
        retainedIPTracks.push_back(track);
        continue;
      }

      const double chi2IP = currentIP->getChi2Track(track);
      const double chi2Secondary = VertexFitterSimple_V{}.getChi2(vertex, track, 0, PointBase::SECVTX);
      if (chi2IP <= chi2Secondary * cfg.chi2ratioIP) {
        retainedIPTracks.push_back(track);
        continue;
      }

      auto tracks = vertex->getTracks();
      tracks.push_back(track);
      auto* refittedVertex = VertexFitterSimple_V{}(tracks.begin(), tracks.end(), vertex, true);
      delete vertex;
      vertex = refittedVertex;
    }

    if (retainedIPTracks.size() < currentIP->getTracks().size()) {
      currentIP.reset(VertexFitterSimple_V{}(retainedIPTracks.begin(), retainedIPTracks.end(), beamVertex));
    }
  }

  delete beamVertex;
}

void VertexFinderSuehara::associateIPTracksAVF(std::vector<Vertex*>& vertices, Vertex* ip,
                                               VertexFinderSueharaConfig& cfg) {
  // This is the adaptive variant of the legacy IP-track association. As above,
  // work on a private primary-vertex copy because EDM4hep owns the input one.
  auto currentIP = std::make_unique<Vertex>(*ip);
  Vertex* beamVertex;
  algoEtc::makeBeamVertex(beamVertex, cfg.beamspotSmearing);

  bool changed = true;
  while (changed) {
    changed = false;
    std::vector<const Track*> orderedTracks = currentIP->getTracks();
    std::ranges::sort(orderedTracks, [currentIP = currentIP.get()](const auto* lhs, const auto* rhs) {
      return currentIP->getChi2Track(lhs) > currentIP->getChi2Track(rhs);
    });

    std::vector<const Track*> retainedIPTracks;
    for (const auto* track : orderedTracks) {
      const double chi2IP = currentIP->getChi2Track(track);
      std::vector<double> chi2Secondary(vertices.size(), 1.e50);

      for (std::size_t i = 0; i < vertices.size(); ++i) {
        auto* vertex = vertices[i];
        if (vertex->getPos().Mag() < cfg.minimumdistIP)
          continue;

        TLorentzVector vertexMomentum;
        for (const auto* vertexTrack : vertex->getTracks())
          vertexMomentum += *vertexTrack;
        const TLorentzVector trackMomentum = *track;
        if (trackMomentum.Vect().Dot(vertex->getPos()) < 0. || (vertexMomentum + trackMomentum).M() > 10.)
          continue;

        Helix helix(track, PointBase::SECVTX);
        VertexLine flightLine(currentIP->getPos(), vertex->getPos(), PointBase::SECVTX);
        double lineDistance = 0.;
        const TVector3 closestPoint = helix.ClosePoint(flightLine, &lineDistance);
        const TVector3 separation = closestPoint - vertex->getPos();
        const TVector3 position = vertex->getPos();
        const double positionMagnitude = position.Mag();
        if (positionMagnitude == 0.)
          continue;
        const double vertexError = std::sqrt(std::pow(std::sqrt(vertex->getCov()[Vertex::xx]) * position.X(), 2.) +
                                             std::pow(std::sqrt(vertex->getCov()[Vertex::yy]) * position.Y(), 2.) +
                                             std::pow(std::sqrt(vertex->getCov()[Vertex::zz]) * position.Z(), 2.)) /
                                   positionMagnitude;
        const double trackError = std::sqrt(track->getCovMatrix()[tpar::d0d0] + track->getCovMatrix()[tpar::z0z0]);
        if (separation.Mag() > 20. * std::sqrt(trackError * trackError + vertexError * vertexError))
          continue;

        auto candidateTracks = vertex->getTracks();
        candidateTracks.push_back(track);
        auto* candidate = VertexFitterSimple_V{}(candidateTracks.begin(), candidateTracks.end(), vertex, true);
        chi2Secondary[i] = candidate->getChi2Track(track);
        delete candidate;
      }

      const double temperature = cfg.temperature;
      double normalisation = TMath::Exp(-chi2IP / (2. * temperature));
      for (const auto chi2 : chi2Secondary)
        normalisation += TMath::Exp(-chi2 / (2. * temperature));
      const double primaryWeight =
          TMath::Exp(-chi2IP / (2. * temperature)) / (TMath::Exp(-cfg.chi2th / (2. * temperature)) + normalisation);

      const auto best = std::ranges::max_element(chi2Secondary, {}, [temperature, normalisation, &cfg](double chi2) {
        return TMath::Exp(-chi2 / (2. * temperature)) / (TMath::Exp(-cfg.chi2th / (2. * temperature)) + normalisation);
      });
      if (best == chi2Secondary.end()) {
        retainedIPTracks.push_back(track);
        continue;
      }
      const auto vertexIndex = static_cast<std::size_t>(std::distance(chi2Secondary.begin(), best));
      const double secondaryWeight =
          TMath::Exp(-*best / (2. * temperature)) / (TMath::Exp(-cfg.chi2th / (2. * temperature)) + normalisation);
      if (primaryWeight >= 0.5 || secondaryWeight <= 0.5) {
        retainedIPTracks.push_back(track);
        continue;
      }

      auto candidateTracks = vertices[vertexIndex]->getTracks();
      candidateTracks.push_back(track);
      auto* refittedVertex =
          VertexFitterSimple_V{}(candidateTracks.begin(), candidateTracks.end(), vertices[vertexIndex], true);
      double maximumChi2 = 0.;
      for (const auto* candidateTrack : candidateTracks)
        maximumChi2 = std::max(maximumChi2, refittedVertex->getChi2Track(candidateTrack));
      if (maximumChi2 > cfg.chi2th) {
        delete refittedVertex;
        retainedIPTracks.push_back(track);
        continue;
      }
      delete vertices[vertexIndex];
      vertices[vertexIndex] = refittedVertex;
      changed = true;
    }

    if (retainedIPTracks.size() < currentIP->getTracks().size())
      currentIP.reset(VertexFitterSimple_V{}(retainedIPTracks.begin(), retainedIPTracks.end(), beamVertex));
  }

  delete beamVertex;
}

// // AVF
// void VertexFinderSuehara::associateIPTracksAVF(std::vector<Vertex*>& vertices, Vertex* ip,
//                                                VertexFinderSueharaConfig& cfg) {

//   bool verbose = false;

//   std::vector<const Track*>::const_iterator it;

//   Vertex* vbeam;
//   algoEtc::makeBeamVertex(vbeam, cfg.beamspotSmearing);

//   bool lflg = false;
//   Vertex* ip2 = ip;       // reserve ip
//   while (lflg == false) { // looping until fitting becomes stable(it is "adaptive")
//     lflg = true;
//     std::vector<const Track*> loopiptracks; // for looping tracks(original ip tracks)
//     // arrange ip tracks in chi2 descending order
//     for (it = ip->getTracks().begin(); it != ip->getTracks().end(); it++) {
//       if (loopiptracks.size() == 0)
//         loopiptracks.push_back(*it);
//       else {
//         std::vector<const Track*>::iterator it3 = loopiptracks.begin();
//         while (it3 != loopiptracks.end() && ip->getChi2Track(*it) < ip->getChi2Track(*it3)) {
//           it3++;
//         }
//         if (it3 == loopiptracks.end())
//           loopiptracks.push_back(*it);
//         else
//           loopiptracks.insert(it3, *it);
//       }
//     }

//     std::vector<const Track*> iptracks; // for updating tracks
//     for (it = loopiptracks.begin(); it != loopiptracks.end(); it++) {
//       // get ip chisquare of this track
//       double chi2ip = ip->getChi2Track(*it);

//       // initialize
//       std::vector<double> chi2s(vertices.size(), 1.0e+10);
//       std::vector<double> maxchi2s(vertices.size(), 0.0);

//       for (unsigned int i = 0; i < vertices.size(); i++) { // loop over vertices
//         // vertex quality cut
//         if (vertices[i]->getPos().Mag() < cfg.minimumdistIP) {
//           chi2s[i] = 1.0e+50;
//           continue;
//         }

//         TLorentzVector v;
//         for (unsigned int j = 0; j < vertices[i]->getTracks().size(); j++) {
//           v += *(vertices[i]->getTracks()[j]);
//         }

//         TLorentzVector v2 = **it;

//         // cal. distance between vertex and helix track
//         Helix hel(*it, PointBase::SECVTX);

//         VertexLine line(ip->getPos(), vertices[i]->getPos(), PointBase::SECVTX);
//         double linedist = 0;
//         TVector3 pos = hel.ClosePoint(line, &linedist);
//         TVector3 dist = pos - vertices[i]->getPos();

//         // errors
//         TVector3 y0 = vertices[i]->getPos();
//         TVector3 err(sqrt(vertices[i]->getCov()[Vertex::xx]), sqrt(vertices[i]->getCov()[Vertex::yy]),
//                      sqrt(vertices[i]->getCov()[Vertex::zz]));
//         double vtxerr =
//             sqrt(pow(err.X() * y0.X(), 2.0) + pow(err.Y() * y0.Y(), 2.0) + pow(err.Z() * y0.Z(), 2.0)) / y0.Mag();
//         double trkerr = sqrt((*it)->getCovMatrix()[tpar::d0d0] + (*it)->getCovMatrix()[tpar::z0z0]);

//         // trk selection   give large chi2 value for this vertex
//         if (dist.Mag() > 20.0 * sqrt(trkerr * trkerr + vtxerr * vtxerr)) {
//           chi2s[i] = 1.0e+50;
//           continue;
//         } // kill over 20 sigma
//         if (v2.Vect().Dot(vertices[i]->getPos()) < 0.) {
//           chi2s[i] = 1.0e+50;
//           continue;
//         }
//         if ((v + v2).M() > 10.0) {
//           chi2s[i] = 1.0e+50;
//           continue;
//         } // this track looks like very hard track!!

//         // vertex fitter
//         std::vector<const Track*> trks = vertices[i]->getTracks();
//         // for(int ll=0;ll<trks.size();ll++) const_cast<Track*> (trks[ll])->setChi2()
//         trks.push_back(*it);
//         Vertex* vtx2 = VertexFitterSimple_V()(trks.begin(), trks.end(), vertices[i], true);
//         // double chi2new = vtx2->getChi2Track(*it);   //old
//         // check std::max chi2
//         for (unsigned int j = 0; j < trks.size(); j++) {
//           if (maxchi2s[i] < vtx2->getChi2Track(trks[j]))
//             maxchi2s[i] = vtx2->getChi2Track(trks[j]);
//         }

//         chi2s[i] = vtx2->getChi2Track(*it);

//         delete vtx2;
//       }

//       // cal.weight function of this track!
//       double T = cfg.temperature; // temperature 5.0 to suppress fakes to same extent
//       double core = TMath::Exp(-chi2ip / (2.0 * T));
//       for (unsigned int i = 0; i < vertices.size(); i++) // loop over vertices
//         core += TMath::Exp(-chi2s[i] / (2.0 * T));

//       // start to check weight function!
//       // primary weight
//       double pip = TMath::Exp(-chi2ip / (2.0 * T)) / (TMath::Exp(-cfg.chi2th / (2.0 * T)) + core);
//       double pmax1 = 0.0, pmax2 = 0.0;
//       int maxi1 = 0;
//       // int maxi2=0;
//       // double vangle=3.14;
//       for (unsigned int i = 0; i < vertices.size(); i++) {
//         double tmpp = TMath::Exp(-chi2s[i] / (2.0 * T)) / (TMath::Exp(-cfg.chi2th / (2.0 * T)) + core);
//         if (tmpp > pmax1) {
//           pmax2 = pmax1;
//           // maxi2=maxi1;

//           pmax1 = tmpp;
//           maxi1 = i;
//         } else if (tmpp < pmax1 && tmpp > pmax2) {
//           pmax2 = tmpp;
//           // maxi2=i;
//         }
//       }

//       // calculate angle between vertices
//       // if(maxi1!=maxi2) vangle=vertices[maxi1]->getPos().Angle(vertices[maxi2]->getPos());

//       // if((((*it)->getChi2() > chi2s[maxi1] && chi2s[maxi1]>0.0) || (*it)->getChi2() < 0.0) && chi2s[maxi1]
//       < 1.0e+20)
//       // 	const_cast<Track*> (*it)->setChi2(chi2s[maxi1]);

//       // get good track to be attached!
//       if (!((pip < 0.5 &&
//              pmax1 > 0.5))) // good one
//                             //  || (pmax1>pip && pmax2>pip && pmax1+pmax2>0.8 && vangle<0.2)))  //allow this
//                             situation
//       {
//         iptracks.push_back(*it);
//         continue;
//       }

//       // invoke vertex fitter
//       std::vector<const Track*> tracks;
//       for (unsigned int i = 0; i < vertices[maxi1]->getTracks().size(); i++)
//         tracks.push_back(vertices[maxi1]->getTracks()[i]);
//       tracks.push_back(*it);
//       Vertex* vtx = VertexFitterSimple_V()(tracks.begin(), tracks.end(), vertices[maxi1], true);

//       // check chi2 threshold
//       double maxchi2 = 0.0;
//       for (unsigned int i = 0; i < tracks.size(); i++) {
//         double tmpchi2 = vtx->getChi2Track(tracks[i]);
//         if (maxchi2 < tmpchi2)
//           maxchi2 = tmpchi2;
//       }
//       if (maxchi2 > cfg.chi2th) {
//         delete vtx;
//         iptracks.push_back(*it);
//         continue;
//       }

//       // move track into new vertex
//       delete vertices[maxi1];
//       vertices[maxi1] = vtx;

//       if (verbose)
//         std::cout << "Track # " << (*it)->getId() << " moved to vertex " << maxi1 << ", chi2ip = " << chi2ip
//                   << ", chi2new = " << chi2s[maxi1] << std::endl;

//       lflg = false;
//     }

//     if (iptracks.size() < ip->getTracks().size()) {
//       // tracks removed
//       if (ip != ip2)
//         delete ip;
//       ip = VertexFitterSimple_V()(iptracks.begin(), iptracks.end(), vbeam);
//     }
//   }

//   delete vbeam;

//   if (ip2 != ip) {
//     delete ip;
//     ip = ip2;
//   }

//   return;
// }

void VertexFinderSuehara::buildUp(TrackVec& tracks, std::vector<Vertex*>& vtx, std::vector<Vertex*>& v0vtx,
                                  double chi2thpri, VertexFinderSueharaConfig& cfg, Vertex** pIP) {
  // obtain primary vertex
  std::unique_ptr<Vertex> nip;
  Vertex*& ip = *pIP;
  if (!ip) {
    ip = findPrimaryVertex(tracks, chi2thpri, cfg.beamspotConstraint, cfg.beamspotSmearing);
    //		vtx.push_back(ip);
  } else {
    nip = std::unique_ptr<Vertex>(
        new Vertex(ip->getChi2(), ip->getProb(), ip->getX(), ip->getY(), ip->getZ(), ip->getCov(), true));
    //		vtx.push_back(nip);
  }

  // pickup residuals
  std::list<const Track*> residualTracks;
  for (unsigned int i = 0; i < tracks.size(); i++) {
    if (find(ip->getTracks().begin(), ip->getTracks().end(), tracks[i]) == ip->getTracks().end()) {
      residualTracks.push_back(tracks[i]);
    } else if (nip) {
      nip->add(tracks[i], ip->getChi2Track(tracks[i]));
    }
  }

  //   std::cout << "buildUp: secondary tracks " << residualTracks.size() << "/" << tracks.size() << std::endl;

  // secondary vertex
  GetVertexList(residualTracks, ip, vtx, v0vtx, cfg);
  /*
        Vertex *v;
        do{
  //	lcfiplus::Vertex* VertexFinderSuehara::findOne2( std::list<Track *> &tracks, double chi2th, double massth, bool
  removeTracks){ v = findOne2(residualTracks, 9, massth, true); if(v){ vtx.push_back(v);
                }
        }while(v);
  */
}

// void VertexFinderSuehara::buildUpForJetClustering(TrackVec& tracks, std::vector<Vertex*>& vtx) {
//   VertexFinderSueharaConfig cfg;
//   cfg.chi2th = 25.;

//   std::vector<Vertex*> v0vtx;
//   buildUp(tracks, vtx, v0vtx, cfg.chi2th, cfg);
// }

// std::vector<Vertex*> VertexFinderSuehara::makeSingleTrackVertices(Jet* jet, TrackVec& tracks, VertexVec& v0vtx,
//                                                                   const Vertex* ip, VertexFinderSueharaConfig& cfg) {
//   return makeSingleTrackVertices(jet->getVertices(), tracks, v0vtx, ip, cfg);
// }

// std::vector<Vertex*> VertexFinderSuehara::makeSingleTrackVertices(VertexVec& vtcs, TrackVec& tracks, VertexVec&
// v0vtx,
//                                                                   const Vertex* ip, VertexFinderSueharaConfig& cfg) {
//   std::vector<Vertex*> singlevtcs;

//   std::vector<const Track*> v0tracks = algoEtc::extractTracks(v0vtx);

//   for (unsigned int ntr = 0; ntr < tracks.size(); ntr++) {
//     const Track* track = tracks[ntr];
//     if (find(v0tracks.begin(), v0tracks.end(), track) != v0tracks.end())
//       continue;
//     if (track->E() < cfg.minEnergySingle)
//       continue;
//     // if (!cfg.useBNess && track->E() < cfg.minEnergySingle)continue;

//     // d0/z0 cut
//     double d0 = track->getD0();
//     double d0err = sqrt(track->getCovMatrix()[tpar::d0d0]);
//     double z0 = (ip ? track->getZ0() - ip->getZ() : track->getZ0());
//     double z0err = sqrt(track->getCovMatrix()[tpar::z0z0]);

//     if (!cfg.useBNess && fabs(d0 / d0err) < cfg.mind0SigSingle && fabs(z0 / z0err) < cfg.minz0SigSingle)
//       continue;
//     if (cfg.useBNess && fabs(d0 / d0err) < 2.0 && fabs(z0 / z0err) < 2.0)
//       continue;

//     Helix hel(track, PointBase::SECVTX);

//     for (unsigned int nvtcs = 0; nvtcs < vtcs.size(); nvtcs++) {
//       const Vertex* vtx = vtcs[nvtcs];

//       // rejecting opposite direction
//       if (vtx->getPos().Dot(track->Vect()) < 0)
//         continue;

//       // angular preselection
//       double angle = vtx->getPos().Angle(track->Vect());
//       if (angle > cfg.maxAngleSingle)
//         continue;

//       // calculate closest point
//       VertexLine line((ip ? ip->getPos() : TVector3(0, 0, 0)), vtx->getPos(), PointBase::SECVTX);
//       double linedist = 0;
//       TVector3 pos = hel.ClosePoint(line, &linedist);

//       // selection cuts
//       if (pos.Mag() < cfg.minPosSingle || pos.Mag() > cfg.maxPosSingle)
//         continue;
//       // rejecting opposite vtx position
//       if (pos.Dot(vtx->getPos()) < 0.)
//         continue;

//       if (linedist / pos.Mag() > cfg.maxSeparationPerPosSingle)
//         continue;

//       // BNess cut
//       // if (cfg.useBNess && track->getBNess() < cfg.cutBNess) continue;

//       // all selection passed: make single track vertex
//       double cov[6] = {0., 0., 0., 0., 0., 0.};
//       Vertex* newvtx = new Vertex(0, 1, pos.x(), pos.y(), pos.z(), cov, false);
//       newvtx->add(track);

//       singlevtcs.push_back(newvtx);
//       break; // end searching for this track
//     }
//   }

//   //   std::cout << "makeSingleTrackVertices: " << singlevtcs.size() << " vertices found." << std::endl;
//   return singlevtcs;
// }

// void VertexFinderSuehara::optimizeTwoVertices(Vertex*& v1, Vertex*& v2, int nvr) {
//   bool _verbose = false;
//   double chi2sum = 0.;

//   std::vector<const Track*> n1tracks = v1->getTracks();
//   std::vector<const Track*> n2tracks = v2->getTracks();

//   if (n1tracks.size() > 1)
//     chi2sum += v1->getChi2();
//   if (n2tracks.size() > 1)
//     chi2sum += v2->getChi2();

//   bool move = true;
//   Vertex *oldv1 = v1, *oldv2 = v2;

//   // vertex recombination
//   for (int n = 0; n < nvr && move; n++) {
//     move = false;
//     for (unsigned int ntr = 0; ntr < v1->getTracks().size(); ntr++) {
//       if (n1tracks.size() <= 1)
//         break;

//       const Track* curt = v1->getTracks()[ntr];

//       n2tracks.push_back(curt);
//       Vertex* v2d = VertexFitterSimple_V()(n2tracks.begin(), n2tracks.end());

//       if (v2d->getChi2Track(curt) < v1->getChi2Track(curt)) {
//         if (v2 != oldv2)
//           delete v2;
//         v2 = v2d;
//         n1tracks.erase(find(n1tracks.begin(), n1tracks.end(), curt));
//         if (v1 != oldv1)
//           delete v1;
//         v1 = VertexFitterSimple_V()(n1tracks.begin(), n1tracks.end());
//         move = true;
//         if (_verbose)
//           std::cout << "Track " << ntr << " moved from v1 to v2. v1 chi2: " << v1->getChi2Track(curt)
//                     << ", v2 chi2: " << v2->getChi2Track(curt) << std::endl;

//         ntr--;
//       } else {
//         n2tracks.pop_back();
//         delete v2d;
//       }
//     }
//     for (unsigned int ntr = 0; ntr < v2->getTracks().size(); ntr++) {
//       if (n2tracks.size() <= 1)
//         break;

//       const Track* curt = v2->getTracks()[ntr];

//       n1tracks.push_back(curt);
//       Vertex* v1d = VertexFitterSimple_V()(n1tracks.begin(), n1tracks.end());

//       if (v1d->getChi2Track(curt) < v2->getChi2Track(curt)) {
//         if (v1 != oldv1)
//           delete v1;
//         v1 = v1d;
//         n2tracks.erase(find(n2tracks.begin(), n2tracks.end(), curt));
//         if (v2 != oldv2)
//           delete v2;
//         v2 = VertexFitterSimple_V()(n2tracks.begin(), n2tracks.end());
//         move = true;
//         if (_verbose)
//           std::cout << "Track " << ntr << " moved from v2 to v1. v1 chi2: " << v1->getChi2Track(curt)
//                     << ", v2 chi2: " << v2->getChi2Track(curt) << std::endl;
//       } else {
//         n1tracks.pop_back();
//         delete v1d;
//       }
//     }

//     if (_verbose) {
//       std::cout << "[VR" << n << "]";
//       std::cout << "Tracks = " << n1tracks.size() << ", " << n2tracks.size();
//       std::cout << ", prob = " << v1->getProb() << ", " << v2->getProb();
//       std::cout << ", chi2 = " << v1->getChi2() << ", " << v2->getChi2() << std::endl;
//     }

//     if (!move)
//       break;

//     double chi2sumold = chi2sum;

//     chi2sum = 0;
//     if (n1tracks.size() > 1)
//       chi2sum += v1->getChi2();
//     if (n2tracks.size() > 1)
//       chi2sum += v2->getChi2();

//     if (chi2sumold < chi2sum) {
//       if (_verbose)
//         std::cout << "Chi2 by current VR is bigger than old one; rollback to old one." << std::endl;
//       if (v1 != oldv1) {
//         delete v1;
//         v1 = oldv1;
//       }
//       if (v2 != oldv2) {
//         delete v2;
//         v2 = oldv2;
//       }
//       break;
//     } else {
//       delete oldv1;
//       delete oldv2;
//       oldv1 = v1;
//       oldv2 = v2;
//     }
//   }
// }

// void VertexFinderSuehara::recombineVertices(std::vector<Vertex*>& vertices, std::vector<Vertex*>& singleVertices) {
//   bool _verbose = false;
//   if (vertices.size() + singleVertices.size() == 0) {
//     if (_verbose)
//       std::cout << "VertexFinderSuehara::recombineVertices: no vertices & singleVertices found. no-op." << std::endl;
//     return;
//   }
//   if (vertices.size() + singleVertices.size() == 1) {
//     if (_verbose)
//       std::cout << "VertexFinderSuehara::recombineVertices: just one vertices + singleVertices found. Currently
//       no-op. "
//                    "(should be divided...)"
//                 << std::endl;

//     if (vertices.size() == 0)
//       vertices.push_back(singleVertices[0]);
//     return;
//   }

//   // firstly prepare combined vertex collection
//   std::vector<Vertex*> v = vertices;
//   v.insert(v.end(), singleVertices.begin(), singleVertices.end());

//   if (v.size() > 2) {
//     if (_verbose)
//       std::cout << "VertexFinderSuehara::recombineVertices: " << v.size()
//                 << " vertices + singleVertices found. Trying to combine..." << std::endl;

//     int nvtx = v.size();

//     double chi2min = 1e+300;
//     Vertex* v1opt = 0;
//     Vertex* v2opt = 0;

//     for (int n1 = 0; n1 < nvtx - 1; n1++) {
//       for (int n2 = n1 + 1; n2 < nvtx; n2++) {
//         std::vector<const Track*> n1tracks = v[n1]->getTracks();
//         std::vector<const Track*> n2tracks = v[n2]->getTracks();
//         // n1: first reserved vertex; n2: second reserved vertex

//         for (int n = 0; n < nvtx; n++) {
//           if (n == n1 || n == n2)
//             continue;

//           const std::vector<const Track*>& ntracks = v[n]->getTracks();
//           for (unsigned int ntr = 0; ntr < ntracks.size(); ntr++) {
//             n1tracks.push_back(ntracks[ntr]);
//             Vertex* n1v = VertexFitterSimple_V()(n1tracks.begin(), n1tracks.end());
//             double n1chi2 = n1v->getChi2Track(ntracks[ntr]);
//             delete n1v;

//             n2tracks.push_back(ntracks[ntr]);
//             Vertex* n2v = VertexFitterSimple_V()(n2tracks.begin(), n2tracks.end());
//             double n2chi2 = n2v->getChi2Track(ntracks[ntr]);
//             delete n2v;

//             if (n1chi2 < n2chi2) { // select n1chi2
//               n2tracks.pop_back(); // n1tracks remain added
//             } else {
//               n1tracks.pop_back(); // n2tracks remain added
//             }
//           }
//         }
//         // n1tracks / n2tracks now obtained
//         Vertex* v1 = VertexFitterSimple_V()(n1tracks.begin(), n1tracks.end());
//         Vertex* v2 = VertexFitterSimple_V()(n2tracks.begin(), n2tracks.end());

//         // optimizeTwoVertices here is removed: no better effect, comparing to do last

//         double prob = 1.;
//         if (n1tracks.size() > 1)
//           prob *= v1->getProb();
//         if (n2tracks.size() > 1)
//           prob *= v2->getProb();

//         if (_verbose) {
//           std::cout << "n1 = " << n1 << ", n2 = " << n2 << ", tracks = " << n1tracks.size() << ", " <<
//           n2tracks.size(); std::cout << ", prob = " << v1->getProb() << ", " << v2->getProb(); std::cout << ", chi2 =
//           " << v1->getChi2() << ", " << v2->getChi2() << std::endl;
//         }

//         // save current vertices
//         double chi2 = v1->getChi2() + v2->getChi2();
//         if (chi2 < chi2min) {
//           chi2min = chi2;
//           if (v1opt)
//             delete v1opt;
//           if (v2opt)
//             delete v2opt;
//           v1opt = v1;
//           v2opt = v2;
//         } else {
//           delete v1;
//           delete v2;
//         }
//       }
//     }

//     for (unsigned int n = 0; n < vertices.size(); n++) {
//       delete vertices[n];
//     }
//     vertices.clear();
//     for (unsigned int n = 0; n < singleVertices.size(); n++) {
//       delete singleVertices[n];
//     }
//     singleVertices.clear();

//     optimizeTwoVertices(v1opt, v2opt, 3);

//     if (v1opt->getPos().Mag() < v2opt->getPos().Mag()) {
//       vertices.push_back(v1opt);
//       vertices.push_back(v2opt);
//     } else {
//       vertices.push_back(v2opt);
//       vertices.push_back(v1opt);
//     }
//   } else { // vertex number = 2
//     optimizeTwoVertices(v[0], v[1], 3);

//     vertices.clear();
//     singleVertices.clear();

//     if (v[0]->getPos().Mag() < v[1]->getPos().Mag()) {
//       vertices.push_back(v[0]);
//       vertices.push_back(v[1]);
//     } else {
//       vertices.push_back(v[1]);
//       vertices.push_back(v[0]);
//     }
//   }
// }

// // new function for BNess tagger
// void VertexFinderSuehara::recombineVertices(std::vector<Vertex*>& vertices, std::vector<Vertex*>& singleVertices,
//                                             VertexFinderSueharaConfig& cfg) {
//   bool _verbose = false;
//   if (vertices.size() + singleVertices.size() == 0) {
//     if (_verbose)
//       std::cout << "VertexFinderSuehara::recombineVertices: no vertices & singleVertices found. no-op." << std::endl;
//     return;
//   }
//   if (vertices.size() + singleVertices.size() == 1) {
//     if (_verbose)
//       std::cout << "VertexFinderSuehara::recombineVertices: just one vertices + singleVertices found. Currently
//       no-op. "
//                    "(should be divided...)"
//                 << std::endl;

//     if (vertices.size() == 0)
//       vertices.push_back(singleVertices[0]);
//     return;
//   }

//   // firstly prepare combined vertex collection
//   std::vector<Vertex*> v = vertices;
//   v.insert(v.end(), singleVertices.begin(), singleVertices.end());

//   if (v.size() > 2) {
//     if (_verbose)
//       std::cout << "VertexFinderSuehara::recombineVertices: " << v.size()
//                 << " vertices + singleVertices found. Trying to combine..." << std::endl;

//     int nvtx = v.size();

//     double chi2min = 1e+300;
//     Vertex* v1opt = 0;
//     Vertex* v2opt = 0;

//     for (int n1 = 0; n1 < nvtx - 1; n1++) {
//       for (int n2 = n1 + 1; n2 < nvtx; n2++) {
//         std::vector<const Track*> n1tracks; // = v[n1]->getTracks();
//         std::vector<const Track*> n2tracks; // = v[n2]->getTracks();
//                                             // n1: first reserved vertex; n2: second reserved vertex

//         if (cfg.useBNess) { // fake track rejection here
//           // making vertex
//           double maxbness = -1.0;
//           int hbtr = -1;
//           for (unsigned int ntr = 0; ntr < v[n1]->getTracks().size(); ntr++) {
//             if (maxbness < v[n1]->getTracks()[ntr]->getBNess()) {
//               maxbness = v[n1]->getTracks()[ntr]->getBNess();
//               hbtr = ntr;
//             }

//             // BNess fake rejection
//             if (v[n1]->getTracks()[ntr]->E() >= 1.0 && v[n1]->getTracks()[ntr]->getBNess() < cfg.cutBNess)
//               continue;
//             if (v[n1]->getTracks()[ntr]->E() < 1.0 && v[n1]->getTracks()[ntr]->getBNess() < cfg.cutBNessE1)
//               continue;
//             n1tracks.push_back(v[n1]->getTracks()[ntr]);
//           }

//           if (n1tracks.size() == 0)
//             n1tracks.push_back(v[n1]->getTracks()[hbtr]);
//           if (n1tracks.size() <= 1 && v[n1]->getTracks().size() >= 2) { // make 2 track vtx if possible
//             std::vector<const Track*> tmptracks;
//             double maxchi2 = 1.0e+300;
//             int tid = -1;
//             for (unsigned int i = 0; i < v[n1]->getTracks().size(); i++) {
//               if ((int)i == hbtr)
//                 continue;
//               tmptracks.clear();
//               tmptracks.push_back(n1tracks[0]);
//               tmptracks.push_back(v[n1]->getTracks()[i]);
//               Vertex* tmpvtx = VertexFitterSimple_V()(tmptracks.begin(), tmptracks.end());
//               double tmpmaxchi2 = std::max(tmpvtx->getChi2Track(tmptracks[0]), tmpvtx->getChi2Track(tmptracks[1]));
//               if (tmpmaxchi2 < maxchi2) {
//                 maxchi2 = tmpmaxchi2;
//                 tid = i;
//               }
//               delete tmpvtx;
//             }
//             n1tracks.push_back(v[n1]->getTracks()[tid]);
//           }

//           maxbness = -1.0;
//           hbtr = -1;
//           for (unsigned int ntr = 0; ntr < v[n2]->getTracks().size(); ntr++) {
//             if (maxbness < v[n2]->getTracks()[ntr]->getBNess()) {
//               maxbness = v[n2]->getTracks()[ntr]->getBNess();
//               hbtr = ntr;
//             }

//             // BNess fake rejection
//             if (v[n2]->getTracks()[ntr]->E() >= 1.0 && v[n2]->getTracks()[ntr]->getBNess() < cfg.cutBNess)
//               continue;
//             if (v[n2]->getTracks()[ntr]->E() < 1.0 && v[n2]->getTracks()[ntr]->getBNess() < cfg.cutBNessE1)
//               continue;
//             n2tracks.push_back(v[n2]->getTracks()[ntr]);
//           }

//           if (n2tracks.size() == 0)
//             n2tracks.push_back(v[n2]->getTracks()[hbtr]);

//           if (n2tracks.size() <= 1 && v[n2]->getTracks().size() >= 2) { // make 2 track vtx if possible
//             std::vector<const Track*> tmptracks;
//             double maxchi2 = 1.0e+300;
//             int tid = -1;
//             for (unsigned int i = 0; i < v[n2]->getTracks().size(); i++) {
//               if ((int)i == hbtr)
//                 continue;
//               tmptracks.clear();
//               tmptracks.push_back(n2tracks[0]);
//               tmptracks.push_back(v[n2]->getTracks()[i]);
//               Vertex* tmpvtx = VertexFitterSimple_V()(tmptracks.begin(), tmptracks.end());
//               double tmpmaxchi2 = std::max(tmpvtx->getChi2Track(tmptracks[0]), tmpvtx->getChi2Track(tmptracks[1]));
//               if (tmpmaxchi2 < maxchi2) {
//                 maxchi2 = tmpmaxchi2;
//                 tid = i;
//               }
//               delete tmpvtx;
//             }
//             n2tracks.push_back(v[n2]->getTracks()[tid]);
//           }
//         }

//         for (int n = 0; n < nvtx; n++) {
//           if (n == n1 || n == n2)
//             continue;

//           const std::vector<const Track*>& ntracks = v[n]->getTracks();
//           for (unsigned int ntr = 0; ntr < ntracks.size(); ntr++) {

//             // fake track rejection here
//             if (ntracks[ntr]->E() >= 1.0 && ntracks[ntr]->getBNess() < cfg.cutBNess)
//               continue;
//             if (ntracks[ntr]->E() < 1.0 && ntracks[ntr]->getBNess() < cfg.cutBNessE1)
//               continue;

//             n1tracks.push_back(ntracks[ntr]);
//             Vertex* n1v = VertexFitterSimple_V()(n1tracks.begin(), n1tracks.end());
//             double n1chi2 = n1v->getChi2Track(ntracks[ntr]);
//             delete n1v;

//             n2tracks.push_back(ntracks[ntr]);
//             Vertex* n2v = VertexFitterSimple_V()(n2tracks.begin(), n2tracks.end());
//             double n2chi2 = n2v->getChi2Track(ntracks[ntr]);
//             delete n2v;

//             if (n1chi2 < n2chi2) { // select n1chi2
//               n2tracks.pop_back(); // n1tracks remain added
//             } else {
//               n1tracks.pop_back(); // n2tracks remain added
//             }
//           }
//         }
//         // n1tracks / n2tracks now obtained
//         Vertex* v1 = VertexFitterSimple_V()(n1tracks.begin(), n1tracks.end());
//         Vertex* v2 = VertexFitterSimple_V()(n2tracks.begin(), n2tracks.end());

//         // optimizeTwoVertices here is removed: no better effect, comparing to do last

//         double prob = 1.;
//         if (n1tracks.size() > 1)
//           prob *= v1->getProb();
//         if (n2tracks.size() > 1)
//           prob *= v2->getProb();

//         if (_verbose) {
//           std::cout << "n1 = " << n1 << ", n2 = " << n2 << ", tracks = " << n1tracks.size() << ", " <<
//           n2tracks.size(); std::cout << ", prob = " << v1->getProb() << ", " << v2->getProb(); std::cout << ", chi2 =
//           " << v1->getChi2() << ", " << v2->getChi2() << std::endl;
//         }

//         // save current vertices
//         double chi2 = v1->getChi2() + v2->getChi2();
//         if (chi2 < chi2min) {
//           chi2min = chi2;
//           if (v1opt)
//             delete v1opt;
//           if (v2opt)
//             delete v2opt;
//           v1opt = v1;
//           v2opt = v2;
//         } else {
//           delete v1;
//           delete v2;
//         }
//       }
//     }

//     for (unsigned int n = 0; n < vertices.size(); n++) {
//       delete vertices[n];
//     }
//     vertices.clear();
//     for (unsigned int n = 0; n < singleVertices.size(); n++) {
//       delete singleVertices[n];
//     }
//     singleVertices.clear();

//     optimizeTwoVertices(v1opt, v2opt, 3);

//     if (v1opt->getPos().Mag() < v2opt->getPos().Mag()) {
//       vertices.push_back(v1opt);
//       vertices.push_back(v2opt);
//     } else {
//       vertices.push_back(v2opt);
//       vertices.push_back(v1opt);
//     }
//   } else { // vertex number = 2

//     if (cfg.useBNess) { // check BNess of the two vertices
//       // for updated vertices
//       Vertex *okvtx1 = NULL, *okvtx2 = NULL;
//       for (unsigned int n = 0; n < 2; n++) {
//         const std::vector<const Track*>& ntracks = v[n]->getTracks();
//         if (ntracks.size() >= 2) {
//           std::vector<const Track*> n1tracks;
//           unsigned int hbtr = 0;
//           double maxbness = -1.0;
//           for (unsigned int ntr = 0; ntr < ntracks.size(); ntr++) {
//             if (maxbness < ntracks[ntr]->getBNess()) { // save highest bness score track
//               maxbness = ntracks[ntr]->getBNess();
//               hbtr = ntr;
//             }

//             if (ntracks[ntr]->E() >= 1.0 && ntracks[ntr]->getBNess() < cfg.cutBNess)
//               continue;
//             if (ntracks[ntr]->E() < 1.0 && ntracks[ntr]->getBNess() < cfg.cutBNessE1)
//               continue;
//             n1tracks.push_back(ntracks[ntr]);
//           }

//           // reconstruct vertex
//           if (n1tracks.size() >= 2) {
//             Vertex* n1v = VertexFitterSimple_V()(n1tracks.begin(), n1tracks.end());
//             if (n == 0) {
//               if (okvtx1 != NULL)
//                 delete okvtx1;
//               okvtx1 = n1v;
//             } else if (n == 1) {
//               if (okvtx2 != NULL)
//                 delete okvtx2;
//               okvtx2 = n1v;
//             }
//           } else {
//             // make 2 track vertex using highest score bness track
//             double maxchi2 = 1.0e+300;
//             for (unsigned int i = 0; i < ntracks.size(); i++) {
//               n1tracks.clear();
//               if (i == hbtr)
//                 continue;
//               n1tracks.push_back(ntracks[hbtr]);
//               n1tracks.push_back(ntracks[i]);
//               Vertex* n1v = VertexFitterSimple_V()(n1tracks.begin(), n1tracks.end());
//               double tmpmaxchi2 = std::max(n1v->getChi2Track(n1tracks[0]), n1v->getChi2Track(n1tracks[1]));
//               if (tmpmaxchi2 < maxchi2) {
//                 maxchi2 = n1v->getChi2();
//                 if (n == 0) {
//                   if (okvtx1 != NULL)
//                     delete okvtx1;
//                   okvtx1 = n1v;
//                 } else if (n == 1) {
//                   if (okvtx2 != NULL)
//                     delete okvtx2;
//                   okvtx2 = n1v;
//                 }
//               } else
//                 delete n1v;
//             }
//           }
//         }
//       }

//       // update vertex
//       if (okvtx1 != NULL) {
//         delete v[0];
//         v[0] = okvtx1;
//       }
//       if (okvtx2 != NULL) {
//         delete v[1];
//         v[1] = okvtx2;
//       }
//     }

//     optimizeTwoVertices(v[0], v[1], 3);

//     vertices.clear();
//     singleVertices.clear();

//     if (v[0]->getPos().Mag() < v[1]->getPos().Mag()) {
//       vertices.push_back(v[0]);
//       vertices.push_back(v[1]);
//     } else {
//       vertices.push_back(v[1]);
//       vertices.push_back(v[0]);
//     }
//   }
// }
