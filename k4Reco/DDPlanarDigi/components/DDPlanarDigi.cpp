/*
 * Copyright (c) 2014-2023 Key4hep-Project.
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
#include "DDPlanarDigi.h"

#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/MCRecoTrackerHitPlaneAssociationCollection.h"
#include "edm4hep/SimTrackerHit.h"
#include "edm4hep/TrackerHitPlaneCollection.h"

#include "Gaudi/Accumulators/Histogram.h"
#include "Gaudi/Histograming/Sink/Utils.h"

#include "DD4hep/BitFieldCoder.h"
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"

#include <TFile.h>

#include <cmath>
#include <iostream>
#include <random>
#include <sstream>

enum { hu = 0, hv, hT, hitE, hitsAccepted, diffu, diffv, diffT, hSize };

DDPlanarDigi::DDPlanarDigi(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(name, svcLoc,
                       {
                           KeyValue("SimTrackerHitCollectionName", "SimTrackerHits"),
                           KeyValue("HeaderName", "EventHeader"),
                       },
                       {KeyValue("TrackerHitCollectionName", "VTXTrackerHits"),
                        KeyValue("SimTrkHitRelCollection", "VTXTrackerHitRelations")}) {
  m_uidSvc = service<IUniqueIDGenSvc>("UniqueIDGenSvc", true);
  if (!m_uidSvc) {
    error() << "Unable to get UniqueIDGenSvc" << endmsg;
  }

  if (m_resULayer.size() != m_resVLayer.size()) {
    error() << "DDPlanarDigi - Inconsistent number of resolutions given for U and V coordinate: "
            << "ResolutionU  :" << m_resULayer.size() << " != ResolutionV : " << m_resVLayer.size();

    throw std::runtime_error("DDPlanarDigi: Inconsistent number of resolutions given for U and V coordinate");
  }

  m_histograms.resize(hSize);

  m_histograms[hu].reset(new Gaudi::Accumulators::Histogram<1>{this, "hu", "smearing u", {50, -5., +5.}});
  m_histograms[hv].reset(new Gaudi::Accumulators::Histogram<1>{this, "hv", "smearing v", {50, -5., +5.}});
  m_histograms[hT].reset(new Gaudi::Accumulators::Histogram<1>{this, "hT", "smearing time", {50, -5., +5.}});

  m_histograms[diffu].reset(new Gaudi::Accumulators::Histogram<1>{this, "diffu", "diff u", {1000, -.1, +.1}});
  m_histograms[diffv].reset(new Gaudi::Accumulators::Histogram<1>{this, "diffv", "diff v", {1000, -.1, +.1}});
  m_histograms[diffT].reset(new Gaudi::Accumulators::Histogram<1>{this, "diffT", "diff time", {1000, -5., +5.}});

  m_histograms[hitE].reset(new Gaudi::Accumulators::Histogram<1>{this, "hitE", "hitEnergy in keV", {1000, 0, 200}});
  m_histograms[hitsAccepted].reset(
      new Gaudi::Accumulators::Histogram<1>{this, "hitsAccepted", "Fraction of accepted hits [%]", {201, 0, 100.5}});

  m_geoSvc = serviceLocator()->service(m_geoSvcName);
}

StatusCode DDPlanarDigi::initialize() {
  const auto detector = m_geoSvc->getDetector();

  const auto         surfMan = detector->extension<dd4hep::rec::SurfaceManager>();
  dd4hep::DetElement det     = detector->detector(m_subDetName.value());
  surfaceMap                 = surfMan->map(m_subDetName.value());

  if (!surfaceMap) {
    std::stringstream err;
    err << " Could not find surface map for detector: " << det.name() << " in SurfaceManager ";
    throw std::runtime_error(err.str());
  }

  // Get and store the name for a debug message later
  (void)this->getProperty("SimTrackerHitCollectionName", m_collName);

  return StatusCode::SUCCESS;
}

std::tuple<TrackerHitPlaneColl, Association> DDPlanarDigi::operator()(
    const SimTrackerHitCollection& simTrackerHits, const Header& headers) const {
  auto seed = m_uidSvc->getUniqueID(headers[0].getEventNumber(), headers[0].getRunNumber(), this->name());
  info() << "Using seed " << seed << " for event " << headers[0].getEventNumber() << " and run "
         << headers[0].getRunNumber() << endmsg;
  m_engine.seed(seed);
  // For rng calls, use the fact that drawing from a
  // Gaussian with mean mu and variance sigma^2 is the same as drawing from
  // a normal distribution and then multiplying by sigma and adding mu
  auto dist = std::normal_distribution<double>(0, 1);

  int nCreatedHits   = 0;
  int nDismissedHits = 0;

  auto trkhitVec = edm4hep::TrackerHitPlaneCollection();
  auto thsthcol  = edm4hep::MCRecoTrackerHitPlaneAssociationCollection();

  std::string cellIDEncodingString = m_geoSvc->constantAsString(m_encodingStringVariable.value());
  dd4hep::DDSegmentation::BitFieldCoder bitFieldCoder(cellIDEncodingString);

  int nSimHits = simTrackerHits.size();
  debug() << "Processing collection " << m_collName << " with " << simTrackerHits.size() << " hits ... " << endmsg;

  for (const auto& hit : simTrackerHits) {
    ++(*m_histograms[hitE])[hit.getEDep() * (dd4hep::GeV / dd4hep::keV)];

    if (hit.getEDep() < m_minEnergy) {
      debug() << "Hit with insufficient energy " << hit.getEDep() * (dd4hep::GeV / dd4hep::keV) << " keV" << endmsg;
      continue;
    }

    const int cellID0 = hit.getCellID();

    // get the measurement surface for this hit using the CellID
    dd4hep::rec::SurfaceMap::const_iterator sI = surfaceMap->find(cellID0);

    if (sI == surfaceMap->end()) {
      std::stringstream err;
      err << " DDPlanarDigi::processEvent(): no surface found for cellID : " << std::endl;
      // <<   cellid_decoder( hit ).valueString() ;
      throw std::runtime_error(err.str());
    }

    const dd4hep::rec::ISurface* surf  = sI->second;
    int                          layer = bitFieldCoder.get(cellID0, "layer");

    dd4hep::rec::Vector3D oldPos(hit.getPosition()[0], hit.getPosition()[1], hit.getPosition()[2]);
    dd4hep::rec::Vector3D newPos;

    //  Check if Hit is inside sensitive
    if (!surf->insideBounds(dd4hep::mm * oldPos)) {
      // debug() << "  hit at " << oldPos
      //         << " " << cellid_decoder( hit).valueString()
      //         << " is not on surface "
      //         << *surf
      //         << " distance: " << surf->distance(  dd4hep::mm * oldPos )
      //         << endmsg;

      if (m_forceHitsOntoSurface) {
        dd4hep::rec::Vector2D lv           = surf->globalToLocal(dd4hep::mm * oldPos);
        dd4hep::rec::Vector3D oldPosOnSurf = (1. / dd4hep::mm) * surf->localToGlobal(lv);

        debug() << " moved to " << oldPosOnSurf << " distance " << (oldPosOnSurf - oldPos).r() << endmsg;

        oldPos = oldPosOnSurf;

      } else {
        ++nDismissedHits;
        continue;
      }
    }

    // Smear time of the hit and apply the time window cut if needed
    double hitT = hit.getTime();

    if (m_resTLayer.size() and m_resTLayer[0] > 0) {
      float resT = m_resTLayer.size() > 1 ? m_resTLayer[layer] : m_resTLayer[0];

      double tSmear = resT > 0 ? dist(m_engine) * resT : 0;
      ++(*m_histograms[hT])[resT > 0 ? tSmear / resT : 0];
      ++(*m_histograms[diffT])[tSmear];

      hitT += tSmear;
      debug() << "smeared hit at T: " << hit.getTime() << " ns to T: " << hitT
              << " ns according to resolution: " << resT << " ns" << endmsg;
    }

    // Correcting for the propagation time
    if (m_correctTimesForPropagation) {
      double dt = oldPos.r() / (TMath::C() / 1e6);
      hitT -= dt;
      debug() << "corrected hit at R: " << oldPos.r() << " mm by propagation time: " << dt << " ns to T: " << hitT
              << " ns" << endmsg;
    }

    // Skip the hit if its time is outside the acceptance time window
    if (m_useTimeWindow) {
      // TODO: Check that the length of the time window is OK
      float timeWindow_min = m_timeWindowMin.size() > 1 ? m_timeWindowMin[layer] : m_timeWindowMin[0];
      float timeWindow_max = m_timeWindowMax.size() > 1 ? m_timeWindowMax[layer] : m_timeWindowMax[0];
      if (hitT < timeWindow_min || hitT > timeWindow_max) {
        debug() << "hit at T: " << hit.getTime() << " smeared to: " << hitT
                << " is outside the time window: hit dropped" << endmsg;
        ++nDismissedHits;
        continue;
      }
    }

    // Try to smear the hit position but ensure the hit is inside the sensitive region
    dd4hep::rec::Vector3D u = surf->u();
    dd4hep::rec::Vector3D v = surf->v();

    // Get local coordinates on surface
    dd4hep::rec::Vector2D lv = surf->globalToLocal(dd4hep::mm * oldPos);
    double                uL = lv[0] / dd4hep::mm;
    double                vL = lv[1] / dd4hep::mm;

    bool acceptHit = false;
    int  tries     = 0;

    // TODO: check lengths
    float resU = m_resULayer.size() > 1 ? m_resULayer[layer] : m_resULayer[0];
    float resV = m_resVLayer.size() > 1 ? m_resVLayer[layer] : m_resVLayer[0];


    while (tries < m_maxTries) {
      // if( tries > 0 ) debug() << "retry smearing for " <<  cellid_decoder( hit ).valueString() << " : retries " << tries << endmsg;

      double uSmear = dist(m_engine) * resU;
      double vSmear = dist(m_engine) * resV;

      dd4hep::rec::Vector3D newPosTmp;
      if (m_isStrip) {
        if (m_subDetName == "SET") {
          double xStripPos, yStripPos, zStripPos;
          //Find intersection of the strip with the z=centerOfSensor plane to set it as the center of the SET strip
          dd4hep::rec::Vector3D simHitPosSmeared =
              (1. / dd4hep::mm) * (surf->localToGlobal(dd4hep::rec::Vector2D((uL + uSmear) * dd4hep::mm, 0.)));
          zStripPos        = surf->origin()[2] / dd4hep::mm;
          double lineParam = (zStripPos - simHitPosSmeared[2]) / v[2];
          xStripPos        = simHitPosSmeared[0] + lineParam * v[0];
          yStripPos        = simHitPosSmeared[1] + lineParam * v[1];
          newPosTmp        = dd4hep::rec::Vector3D(xStripPos, yStripPos, zStripPos);
        } else {
          newPosTmp = (1. / dd4hep::mm) * (surf->localToGlobal(dd4hep::rec::Vector2D((uL + uSmear) * dd4hep::mm, 0.)));
        }
      } else {
        newPosTmp =
            (1. / dd4hep::mm) *
            (surf->localToGlobal(dd4hep::rec::Vector2D((uL + uSmear) * dd4hep::mm, (vL + vSmear) * dd4hep::mm)));
      }

      debug() << " hit at    : " << oldPos << " smeared to: " << newPosTmp << " uL: " << uL << " vL: " << vL
              << " uSmear: " << uSmear << " vSmear: " << vSmear << endmsg;

      if (surf->insideBounds(dd4hep::mm * newPosTmp)) {
        acceptHit = true;
        newPos    = newPosTmp;

        ++(*m_histograms[hu])[uSmear / resU];
        ++(*m_histograms[hv])[vSmear / resV];

        ++(*m_histograms[diffu])[uSmear];
        ++(*m_histograms[diffv])[vSmear];

        break;
      }

      // debug() << "  hit at " << newPosTmp
      //         << " " << cellid_decoder( hit).valueString()
      //         << " is not on surface "
      //         << " distance: " << surf->distance( dd4hep::mm * newPosTmp )
      //         << endmsg;

      ++tries;
    }

    if (!acceptHit) {
      debug() << "hit could not be smeared within ladder after " << m_maxTries << "  tries: hit dropped" << endmsg;
      ++nDismissedHits;
      continue;
    }

    auto trkHit = trkhitVec.create();

    trkHit.setCellID(cellID0);

    trkHit.setPosition(newPos.const_array());
    trkHit.setTime(hitT);
    trkHit.setEDep(hit.getEDep());

    float u_direction[2];
    u_direction[0] = u.theta();
    u_direction[1] = u.phi();

    float v_direction[2];
    v_direction[0] = v.theta();
    v_direction[1] = v.phi();

    debug() << " U[0] = " << u_direction[0] << " U[1] = " << u_direction[1] << " V[0] = " << v_direction[0]
            << " V[1] = " << v_direction[1] << endmsg;

    trkHit.setU(u_direction);
    trkHit.setV(v_direction);
    trkHit.setDu(resU);

    if (m_isStrip) {
      // store the resolution from the length of the wafer - in case a fitter might want to treat this as 2d hit ....
      double stripRes = surf->length_along_v() / dd4hep::mm / std::sqrt(12);
      trkHit.setDv(stripRes);
      // TODO: Set type?
      // trkHit.setType( UTIL::set_bit( trkHit.getType() ,  UTIL::ILDTrkHitTypeBit::ONE_DIMENSIONAL ) );

    } else {
      trkHit.setDv(resV);
    }

    auto association = thsthcol.create();
    association.setSim(hit);
    association.setRec(trkHit);

    ++nCreatedHits;
  }

  // Filling the fraction of accepted hits in the event
  float accFraction = nSimHits > 0 ? float(nCreatedHits) / float(nSimHits) * 100.0 : 0.0;
  ++(*m_histograms[hitsAccepted])[accFraction];

  debug() << "Created " << nCreatedHits << " hits, " << nDismissedHits << " hits  dismissed" << endmsg;

  return std::make_tuple(std::move(trkhitVec), std::move(thsthcol));
}

StatusCode DDPlanarDigi::finalize() {
  auto file = TFile::Open(m_outputFileName.value().c_str(), "RECREATE");
  auto names = {"hu", "hv", "hT", "hitE", "hitsAccepted", "diffu", "diffv", "diffT", "hSize"};
  auto it = names.begin();
  for (auto& h : m_histograms) {
    std::string name = "";
    // Name that will appear in the stats table
    std::string    histName = *it;
    nlohmann::json json = *h;
    auto [histo, dir] =
        Gaudi::Histograming::Sink::jsonToRootHistogram<Gaudi::Histograming::Sink::Traits<false, TH1D, 1>>(
            name, histName, json);
    histo.Write(*it);
    ++it;
  }
  file->Close();
  return StatusCode::SUCCESS;
}