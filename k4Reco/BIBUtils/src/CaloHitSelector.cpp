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
#include "CaloHitSelector.h"

#include "BIBUtilsHelpers.h"

#include <edm4hep/SimCalorimeterHit.h>

#include <DDSegmentation/BitFieldCoder.h>

#include <TFile.h>
#include <TH2D.h>
#include <TMath.h>

#include <cmath>
#include <unordered_map>

using k4reco::bibutils::objectKey;

CaloHitSelector::CaloHitSelector(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(name, svcLoc,
                       {KeyValues("CaloHitCollectionName", {"EcalBarrelCollectionRec"}),
                        KeyValues("CaloRelationCollectionName", {"EcalBarrelRelationsSimRec"})},
                       {KeyValues("GoodHitCollection", {"EcalBarrelCollectionSel"}),
                        KeyValues("GoodRelationCollection", {"EcalBarrelRelationsSimSel"})}) {}

StatusCode CaloHitSelector::initialize() {
  m_geoSvc = serviceLocator()->service("GeoSvc");
  if (!m_geoSvc) {
    error() << "Unable to retrieve the GeoSvc" << endmsg;
    return StatusCode::FAILURE;
  }

  // Load the threshold maps if a file was provided. They are detached from the
  // file (SetDirectory(nullptr)) so they survive after it is closed.
  if (!m_thFile.value().empty()) {
    std::unique_ptr<TFile> thFile(TFile::Open(m_thFile.value().c_str(), "READ"));
    if (!thFile || thFile->IsZombie()) {
      error() << "Could not open the thresholds file: " << m_thFile.value() << endmsg;
      return StatusCode::FAILURE;
    }
    m_thresholdMap.reset(dynamic_cast<TH2D*>(thFile->Get("th_2dmode_sym")));
    m_stddevMap.reset(dynamic_cast<TH2D*>(thFile->Get("stddev_sym")));
    if (!m_thresholdMap || !m_stddevMap) {
      error() << "Could not find the histograms th_2dmode_sym / stddev_sym in " << m_thFile.value() << endmsg;
      return StatusCode::FAILURE;
    }
    m_thresholdMap->SetDirectory(nullptr);
    m_stddevMap->SetDirectory(nullptr);
    thFile->Close();
  }

  if (!m_thresholdMap && m_flatThreshold <= 0.) {
    error() << "Neither a thresholds file nor a positive FlatThreshold was provided" << endmsg;
    return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS;
}

std::tuple<edm4hep::CalorimeterHitCollection, edm4hep::CaloHitSimCaloHitLinkCollection>
CaloHitSelector::operator()(const edm4hep::CalorimeterHitCollection& caloHits,
                            const edm4hep::CaloHitSimCaloHitLinkCollection& caloLinks) const {
  edm4hep::CalorimeterHitCollection outHits;
  outHits.setSubsetCollection();
  edm4hep::CaloHitSimCaloHitLinkCollection outLinks;

  const std::string encoderString = m_geoSvc->constantAsString(m_encodingStringVariable.value());
  const dd4hep::DDSegmentation::BitFieldCoder bitFieldCoder(encoderString);

  // Map each reconstructed hit to its simulated hit through the input links.
  std::unordered_map<std::uint64_t, edm4hep::SimCalorimeterHit> hitToSim;
  hitToSim.reserve(caloLinks.size());
  for (const auto& link : caloLinks) {
    hitToSim.emplace(objectKey(link.getFrom().getObjectID()), link.getTo());
  }

  std::size_t nAccepted = 0;
  for (const auto& hit : caloHits) {
    const unsigned int layer = bitFieldCoder.get(hit.getCellID(), "layer");

    const auto& pos = hit.getPosition();

    // Polar angle, symmetrized around pi/2 to match the threshold maps.
    double hitTheta = std::atan2(std::sqrt(pos.x * pos.x + pos.y * pos.y), static_cast<double>(pos.z));
    if (hitTheta > TMath::Pi() / 2.) {
      hitTheta = TMath::Pi() - hitTheta;
    }

    double modeThreshold = 0.;
    double stddev = 0.;
    if (m_thresholdMap) {
      const int binx = m_thresholdMap->GetXaxis()->FindBin(hitTheta);
      const int biny = m_thresholdMap->GetYaxis()->FindBin(layer);
      modeThreshold = m_thresholdMap->GetBinContent(binx, biny);
      stddev = m_stddevMap->GetBinContent(binx, biny);
    }

    double threshold = modeThreshold + m_nSigma * stddev;
    if (m_flatThreshold > 0.) {
      threshold = m_flatThreshold;
    }

    double hitEnergy = hit.getEnergy();
    if (m_doBIBsubtraction) {
      hitEnergy -= modeThreshold;
    }

    if (hitEnergy <= threshold) {
      continue;
    }

    // The stored hit time is already corrected for the time of flight from
    // the IP by the digitizer (RealisticCaloDigi with
    // timingCorrectForPropagation, which keeps the earliest accepted
    // contribution's corrected time), so it is used directly here.
    // NB: earlier versions subtracted r/TMath::C() again; with r in mm and
    // TMath::C() in m/s that term was ~1e-6 ns, i.e. numerically a no-op.
    // Subtracting a correctly-computed TOF here would double-count it.
    const double relativeTime = hit.getTime();
    if (relativeTime <= m_timeWindowMin || relativeTime >= m_timeWindowMax) {
      continue;
    }

    outHits.push_back(hit);
    const auto simIt = hitToSim.find(objectKey(hit.getObjectID()));
    if (simIt != hitToSim.end()) {
      auto link = outLinks.create();
      link.setFrom(hit);
      link.setTo(simIt->second);
      link.setWeight(1.0);
    }
    ++nAccepted;
  }

  debug() << nAccepted << " / " << caloHits.size() << " calo hits passed the BIB selection" << endmsg;

  return std::make_tuple(std::move(outHits), std::move(outLinks));
}

DECLARE_COMPONENT(CaloHitSelector)
