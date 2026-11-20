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
#include "VertexFinderSuehara.h"
#include "VertexFitterSimple.h"
#include "algoEtc.h"
#include <algorithm>
#include <cmath>
#include <vector>
using namespace std;
namespace lcfiplus {
vector<Vertex*> VertexFinderSuehara::makeSingleTrackVertices(VertexVec& vtcs, TrackVec& tracks, VertexVec& v0vtx,
                                                             const Vertex* ip, VertexFinderSueharaConfig& cfg) {
  vector<Vertex*> singlevtcs;

  vector<const Track*> v0tracks = algoEtc::extractTracks(v0vtx);

  for (unsigned int ntr = 0; ntr < tracks.size(); ntr++) {
    const Track* track = tracks[ntr];
    if (find(v0tracks.begin(), v0tracks.end(), track) != v0tracks.end())
      continue;
    if (track->E() < cfg.minEnergySingle)
      continue;
    // if (!cfg.useBNess && track->E() < cfg.minEnergySingle)continue;

    // d0/z0 cut
    double d0 = track->getD0();
    double d0err = sqrt(track->getCovMatrix()[tpar::d0d0]);
    double z0 = (ip ? track->getZ0() - ip->getZ() : track->getZ0());
    double z0err = sqrt(track->getCovMatrix()[tpar::z0z0]);

    if (!cfg.useBNess && fabs(d0 / d0err) < cfg.mind0SigSingle && fabs(z0 / z0err) < cfg.minz0SigSingle)
      continue;
    if (cfg.useBNess && fabs(d0 / d0err) < 2.0 && fabs(z0 / z0err) < 2.0)
      continue;

    Helix hel(track, PointBase::SECVTX);

    for (unsigned int nvtcs = 0; nvtcs < vtcs.size(); nvtcs++) {
      const Vertex* vtx = vtcs[nvtcs];

      // rejecting opposite direction
      if (vtx->getPos().Dot(track->Vect()) < 0)
        continue;

      // angular preselection
      double angle = vtx->getPos().Angle(track->Vect());
      if (angle > cfg.maxAngleSingle)
        continue;

      // calculate closest point
      VertexLine line((ip ? ip->getPos() : TVector3(0, 0, 0)), vtx->getPos(), PointBase::SECVTX);
      double linedist = 0;
      TVector3 pos = hel.ClosePoint(line, &linedist);

      // selection cuts
      if (pos.Mag() < cfg.minPosSingle || pos.Mag() > cfg.maxPosSingle)
        continue;
      // rejecting opposite vtx position
      if (pos.Dot(vtx->getPos()) < 0.)
        continue;

      if (linedist / pos.Mag() > cfg.maxSeparationPerPosSingle)
        continue;

      // BNess cut
      // if (cfg.useBNess && track->getBNess() < cfg.cutBNess) continue;

      // all selection passed: make single track vertex
      double cov[6] = {0., 0., 0., 0., 0., 0.};
      Vertex* newvtx = new Vertex(0, 1, pos.x(), pos.y(), pos.z(), cov, false);
      newvtx->add(track);

      singlevtcs.push_back(newvtx);
      break; // end searching for this track
    }
  }

  // cout << "makeSingleTrackVertices: " << singlevtcs.size() << " vertices found." << endl;
  return singlevtcs;
}

void VertexFinderSuehara::optimizeTwoVertices(Vertex*& v1, Vertex*& v2, int nvr) {
  bool _verbose = false;
  double chi2sum = 0.;

  vector<const Track*> n1tracks = v1->getTracks();
  vector<const Track*> n2tracks = v2->getTracks();

  if (n1tracks.size() > 1)
    chi2sum += v1->getChi2();
  if (n2tracks.size() > 1)
    chi2sum += v2->getChi2();

  bool move = true;
  Vertex *oldv1 = v1, *oldv2 = v2;

  // vertex recombination
  for (int n = 0; n < nvr && move; n++) {
    move = false;
    for (unsigned int ntr = 0; ntr < v1->getTracks().size(); ntr++) {
      if (n1tracks.size() <= 1)
        break;

      const Track* curt = v1->getTracks()[ntr];

      n2tracks.push_back(curt);
      Vertex* v2d = VertexFitterSimple_V()(n2tracks.begin(), n2tracks.end());

      if (v2d->getChi2Track(curt) < v1->getChi2Track(curt)) {
        if (v2 != oldv2)
          delete v2;
        v2 = v2d;
        n1tracks.erase(find(n1tracks.begin(), n1tracks.end(), curt));
        if (v1 != oldv1)
          delete v1;
        v1 = VertexFitterSimple_V()(n1tracks.begin(), n1tracks.end());
        move = true;
        if (_verbose)
          cout << "Track " << ntr << " moved from v1 to v2. v1 chi2: " << v1->getChi2Track(curt)
               << ", v2 chi2: " << v2->getChi2Track(curt) << endl;

        ntr--;
      } else {
        n2tracks.pop_back();
        delete v2d;
      }
    }
    for (unsigned int ntr = 0; ntr < v2->getTracks().size(); ntr++) {
      if (n2tracks.size() <= 1)
        break;

      const Track* curt = v2->getTracks()[ntr];

      n1tracks.push_back(curt);
      Vertex* v1d = VertexFitterSimple_V()(n1tracks.begin(), n1tracks.end());

      if (v1d->getChi2Track(curt) < v2->getChi2Track(curt)) {
        if (v1 != oldv1)
          delete v1;
        v1 = v1d;
        n2tracks.erase(find(n2tracks.begin(), n2tracks.end(), curt));
        if (v2 != oldv2)
          delete v2;
        v2 = VertexFitterSimple_V()(n2tracks.begin(), n2tracks.end());
        move = true;
        if (_verbose)
          cout << "Track " << ntr << " moved from v2 to v1. v1 chi2: " << v1->getChi2Track(curt)
               << ", v2 chi2: " << v2->getChi2Track(curt) << endl;
      } else {
        n1tracks.pop_back();
        delete v1d;
      }
    }

    if (_verbose) {
      cout << "[VR" << n << "]";
      cout << "Tracks = " << n1tracks.size() << ", " << n2tracks.size();
      cout << ", prob = " << v1->getProb() << ", " << v2->getProb();
      cout << ", chi2 = " << v1->getChi2() << ", " << v2->getChi2() << endl;
    }

    if (!move)
      break;

    double chi2sumold = chi2sum;

    chi2sum = 0;
    if (n1tracks.size() > 1)
      chi2sum += v1->getChi2();
    if (n2tracks.size() > 1)
      chi2sum += v2->getChi2();

    if (chi2sumold < chi2sum) {
      if (_verbose)
        cout << "Chi2 by current VR is bigger than old one; rollback to old one." << endl;
      if (v1 != oldv1) {
        delete v1;
        v1 = oldv1;
      }
      if (v2 != oldv2) {
        delete v2;
        v2 = oldv2;
      }
      break;
    } else {
      delete oldv1;
      delete oldv2;
      oldv1 = v1;
      oldv2 = v2;
    }
  }
}

void VertexFinderSuehara::recombineVertices(vector<Vertex*>& vertices, vector<Vertex*>& singleVertices) {
  bool _verbose = false;
  if (vertices.size() + singleVertices.size() == 0) {
    if (_verbose)
      cout << "VertexFinderSuehara::recombineVertices: no vertices & singleVertices found. no-op." << endl;
    return;
  }
  if (vertices.size() + singleVertices.size() == 1) {
    if (_verbose)
      cout << "VertexFinderSuehara::recombineVertices: just one vertices + singleVertices found. Currently no-op. "
              "(should be divided...)"
           << endl;

    if (vertices.size() == 0)
      vertices.push_back(singleVertices[0]);
    return;
  }

  // firstly prepare combined vertex collection
  vector<Vertex*> v = vertices;
  v.insert(v.end(), singleVertices.begin(), singleVertices.end());

  if (v.size() > 2) {
    if (_verbose)
      cout << "VertexFinderSuehara::recombineVertices: " << v.size()
           << " vertices + singleVertices found. Trying to combine..." << endl;

    int nvtx = v.size();

    double chi2min = 1e+300;
    Vertex* v1opt = 0;
    Vertex* v2opt = 0;

    for (int n1 = 0; n1 < nvtx - 1; n1++) {
      for (int n2 = n1 + 1; n2 < nvtx; n2++) {
        vector<const Track*> n1tracks = v[n1]->getTracks();
        vector<const Track*> n2tracks = v[n2]->getTracks();
        // n1: first reserved vertex; n2: second reserved vertex

        for (int n = 0; n < nvtx; n++) {
          if (n == n1 || n == n2)
            continue;

          const vector<const Track*>& ntracks = v[n]->getTracks();
          for (unsigned int ntr = 0; ntr < ntracks.size(); ntr++) {
            n1tracks.push_back(ntracks[ntr]);
            Vertex* n1v = VertexFitterSimple_V()(n1tracks.begin(), n1tracks.end());
            double n1chi2 = n1v->getChi2Track(ntracks[ntr]);
            delete n1v;

            n2tracks.push_back(ntracks[ntr]);
            Vertex* n2v = VertexFitterSimple_V()(n2tracks.begin(), n2tracks.end());
            double n2chi2 = n2v->getChi2Track(ntracks[ntr]);
            delete n2v;

            if (n1chi2 < n2chi2) { // select n1chi2
              n2tracks.pop_back(); // n1tracks remain added
            } else {
              n1tracks.pop_back(); // n2tracks remain added
            }
          }
        }
        // n1tracks / n2tracks now obtained
        Vertex* v1 = VertexFitterSimple_V()(n1tracks.begin(), n1tracks.end());
        Vertex* v2 = VertexFitterSimple_V()(n2tracks.begin(), n2tracks.end());

        // optimizeTwoVertices here is removed: no better effect, comparing to do last

        if (_verbose) {
          cout << "n1 = " << n1 << ", n2 = " << n2 << ", tracks = " << n1tracks.size() << ", " << n2tracks.size();
          cout << ", prob = " << v1->getProb() << ", " << v2->getProb();
          cout << ", chi2 = " << v1->getChi2() << ", " << v2->getChi2() << endl;
        }

        // save current vertices
        double chi2 = v1->getChi2() + v2->getChi2();
        if (chi2 < chi2min) {
          chi2min = chi2;
          if (v1opt)
            delete v1opt;
          if (v2opt)
            delete v2opt;
          v1opt = v1;
          v2opt = v2;
        } else {
          delete v1;
          delete v2;
        }
      }
    }

    for (unsigned int n = 0; n < vertices.size(); n++) {
      delete vertices[n];
    }
    vertices.clear();
    for (unsigned int n = 0; n < singleVertices.size(); n++) {
      delete singleVertices[n];
    }
    singleVertices.clear();

    optimizeTwoVertices(v1opt, v2opt, 3);

    if (v1opt->getPos().Mag() < v2opt->getPos().Mag()) {
      vertices.push_back(v1opt);
      vertices.push_back(v2opt);
    } else {
      vertices.push_back(v2opt);
      vertices.push_back(v1opt);
    }
  } else { // vertex number = 2
    optimizeTwoVertices(v[0], v[1], 3);

    vertices.clear();
    singleVertices.clear();

    if (v[0]->getPos().Mag() < v[1]->getPos().Mag()) {
      vertices.push_back(v[0]);
      vertices.push_back(v[1]);
    } else {
      vertices.push_back(v[1]);
      vertices.push_back(v[0]);
    }
  }
}

} // namespace lcfiplus
