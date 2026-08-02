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
#ifndef K4RECO_VERTEXFITTERSIMPLE_H
#define K4RECO_VERTEXFITTERSIMPLE_H

#include "geometry.h"
#include "lcfiplus.h"

#include <iostream>
#include <list>
#include <vector>

namespace lcfiplus {

template <class Iterator>
class VertexFitterSimple {
public:
  Vertex* operator()(Iterator tracksBegin, Iterator tracksEnd, Vertex* pointConstraint = nullptr,
                     bool pointInitialOnly = false) {
    bool verbose = false;

    GeometryHandler* gh = GeometryHandler::Instance();
    if (pointConstraint) {
      std::vector<PointBase*> tracks;
      if (!pointInitialOnly) {
        Point* ip = new Point(pointConstraint, PointBase::PRIVTX);
        tracks.push_back(ip);
      }
      int ntracks = 0;
      for (Iterator it = tracksBegin; it != tracksEnd; it++, ntracks++) {
        tracks.push_back(new Helix(*it, PointBase::PRIVTX));
      }
      if (verbose)
        std::cout << "VertexFitterSimple: number of tracks is " << ntracks << std::endl;

      TVector3 initial = pointConstraint->getPos();
      Point* result = new Point(PointBase::NOTUSED);
      double chi2 = -gh->PointFit(tracks, initial, result);

      TVector3 vresult = result->GetPos();
      double cov[6];
      cov[Vertex::xx] = result->GetErr(0, 0);
      cov[Vertex::xy] = result->GetErr(0, 1);
      cov[Vertex::xz] = result->GetErr(0, 2);
      cov[Vertex::yy] = result->GetErr(1, 1);
      cov[Vertex::yz] = result->GetErr(1, 2);
      cov[Vertex::zz] = result->GetErr(2, 2);

      if (verbose) {
        std::cout << "Vertex cov matrix:" << std::endl;
        std::cout << std::scientific << cov[Vertex::xx] << "  ";
        std::cout << cov[Vertex::yy] << "  ";
        std::cout << cov[Vertex::zz] << "  ";
        std::cout << cov[Vertex::xy] << "  ";
        std::cout << cov[Vertex::yz] << "  ";
        std::cout << cov[Vertex::xz] << std::endl << std::fixed;
      }

      if (verbose) {
        std::cout << "VertexFitterSimple: vertex position is " << std::endl;
        vresult.Print();
      }

      Vertex* vtx =
          new Vertex(chi2, TMath::Prob(chi2, ntracks * 2 - 3), vresult.x(), vresult.y(), vresult.z(), cov, false);
      for (Iterator it = tracksBegin; it != tracksEnd; it++, ntracks++) {
        Helix hel(*it, PointBase::PRIVTX);
        double ll = hel.LogLikelihood(vresult); // need to incorporate vertex error??

        if (verbose)
          std::cout << "VertexFitterSimple: track loglikelihood is " << ll << std::endl;

        vtx->add(*it, -ll);
      }
      for (std::vector<PointBase*>::iterator it = tracks.begin(); it != tracks.end(); it++) {
        delete *it;
      }

      delete result;
      return vtx;
    }

    // without point constraint
    std::vector<Helix*> tracks;
    int ntracks = 0;
    for (Iterator it = tracksBegin; it != tracksEnd; it++, ntracks++) {
      tracks.push_back(new Helix(*it, PointBase::SECVTX));
    }
    if (verbose)
      std::cout << "VertexFitterSimple: number of tracks is " << ntracks << std::endl;

    Point* result = new Point(PointBase::NOTUSED);
    double chi2 = -gh->HelixPointFit(tracks, result);

    TVector3 vresult = result->GetPos();
    double cov[6];
    cov[Vertex::xx] = result->GetErr(0, 0);
    cov[Vertex::xy] = result->GetErr(0, 1);
    cov[Vertex::xz] = result->GetErr(0, 2);
    cov[Vertex::yy] = result->GetErr(1, 1);
    cov[Vertex::yz] = result->GetErr(1, 2);
    cov[Vertex::zz] = result->GetErr(2, 2);

    if (verbose) {
      std::cout << "VertexFitterSimple: vertex position is " << std::endl;
      vresult.Print();
    }

    Vertex* vtx = new Vertex(chi2, (ntracks > 1 ? TMath::Prob(chi2, ntracks * 2 - 3) : 1), vresult.x(), vresult.y(),
                             vresult.z(), cov, false);
    for (Iterator it = tracksBegin; it != tracksEnd; it++, ntracks++) {
      Helix hel(*it, PointBase::SECVTX);
      double ll = hel.LogLikelihood(vresult); // need to incorporate vertex error??
      if (verbose)
        std::cout << "VertexFitterSimple: track loglikelihood is " << ll << std::endl;

      vtx->add(*it, -ll);
    }
    for (std::vector<Helix*>::iterator it = tracks.begin(); it != tracks.end(); it++) {
      delete *it;
    }

    delete result;
    return vtx;
  }

  double getChi2(const Vertex* vtx, const Track* trk, int mode = 1, PointBase::FITFLAG flag = PointBase::NOTUSED) {
    // 110510 suehara for IPassoc study
    if (mode == 0) {
      // mode 0: no fit at all
      Helix hel(trk, flag);
      TVector3 v = vtx->getPos();
      return -hel.LogLikelihood(v);
    } else if (mode == 1) {
      // mode 1: vertex treated as errored point

      Point ptVtx(vtx, flag);
      Helix hel(trk, flag);

      std::vector<PointBase*> vpt;
      vpt.push_back(&ptVtx);
      vpt.push_back(&hel);

      GeometryHandler* gh = GeometryHandler::Instance();
      return -gh->PointFit(vpt, vtx->getPos(), nullptr);
    } else {
      // mode 2: vertex treated as track list
      std::vector<PointBase*> vpt;
      for (unsigned int n = 0; n < vtx->getTracks().size(); n++) {
        vpt.push_back(new Helix(vtx->getTracks()[n], flag));
      }
      vpt.push_back(new Helix(trk, flag));

      GeometryHandler* gh = GeometryHandler::Instance();
      Point ret(PointBase::NOTUSED);
      double ll = -gh->PointFit(vpt, vtx->getPos(), &ret);

      for (unsigned int n = 0; n < vpt.size(); n++) {
        delete vpt[n];
      }
      if (mode == 2) {
        return ll;
      } else {
        Helix hel(trk, flag);
        return -hel.LogLikelihood(ret.GetPos());
      }
    }
  }
};

typedef VertexFitterSimple<std::vector<const Track*>::const_iterator> VertexFitterSimple_V;
typedef VertexFitterSimple<std::list<const Track*>::const_iterator> VertexFitterSimple_L;
} // namespace lcfiplus

#endif
