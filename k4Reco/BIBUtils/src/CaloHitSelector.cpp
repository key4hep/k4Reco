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

#include <edm4hep/SimCalorimeterHit.h>
#include <edm4hep/utils/vector_utils.h>

#include <podio/ObjectID.h>

#include <k4Interface/IGeoSvc.h>

#include <TFile.h>
#include <TH2D.h>
#include <TMath.h>

#include <exception>
#include <initializer_list>
#include <unordered_map>
#include <utility>

CaloHitSelector::CaloHitSelector(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(name, svcLoc,
                       {KeyValue("CaloHitCollectionName", "EcalBarrelCollectionRec"),
                        KeyValue("CaloRelationCollectionName", "EcalBarrelRelationsSimRec")},
                       {KeyValue("GoodHitCollection", "EcalBarrelCollectionSel"),
                        KeyValue("GoodRelationCollection", "EcalBarrelRelationsSimSel")}) {}

StatusCode CaloHitSelector::initialize() {
  const auto geoSvc = serviceLocator()->service<IGeoSvc>("GeoSvc");
  if (!geoSvc) {
    error() << "Unable to retrieve the GeoSvc" << endmsg;
    return StatusCode::FAILURE;
  }

  // The cellID encoding is fixed by the geometry, so the decoder and the index
  // of the layer field are set up once here rather than for every event.
  try {
    const std::string encoderString = geoSvc->constantAsString(m_encodingStringVariable.value());
    m_bitFieldCoder = std::make_unique<dd4hep::DDSegmentation::BitFieldCoder>(encoderString);
    m_layerIndex = m_bitFieldCoder->index("layer");
  } catch (const std::exception& e) {
    error() << "Could not set up the cellID decoder from " << m_encodingStringVariable.value() << ": " << e.what()
            << endmsg;
    return StatusCode::FAILURE;
  }

  // Load the threshold maps if a file was provided. Each map is detached from
  // the file (SetDirectory(nullptr)) as soon as it is read, so that it is owned
  // only by its unique_ptr and survives the file being closed. The maps are then
  // stored as const: operator() only reads them, which keeps it thread-safe.
  if (!m_thFile.value().empty()) {
    std::unique_ptr<TFile> thFile(TFile::Open(m_thFile.value().c_str(), "READ"));
    if (!thFile || thFile->IsZombie()) {
      error() << "Could not open the thresholds file: " << m_thFile.value() << endmsg;
      return StatusCode::FAILURE;
    }
    std::unique_ptr<TH2D> thresholdMap(dynamic_cast<TH2D*>(thFile->Get("th_2dmode_sym")));
    std::unique_ptr<TH2D> stddevMap(dynamic_cast<TH2D*>(thFile->Get("stddev_sym")));
    for (auto* map : {thresholdMap.get(), stddevMap.get()}) {
      if (map) {
        map->SetDirectory(nullptr);
      }
    }
    if (!thresholdMap || !stddevMap) {
      error() << "Could not find the histograms th_2dmode_sym / stddev_sym in " << m_thFile.value() << endmsg;
      return StatusCode::FAILURE;
    }
    thFile->Close();
    m_thresholdMap = std::move(thresholdMap);
    m_stddevMap = std::move(stddevMap);
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

  // Map each reconstructed hit to its simulated hit through the input links.
  std::unordered_map<podio::ObjectID, edm4hep::SimCalorimeterHit> hitToSim;
  hitToSim.reserve(caloLinks.size());
  for (const auto& link : caloLinks) {
    hitToSim.emplace(link.getFrom().getObjectID(), link.getTo());
  }

  std::size_t nAccepted = 0;
  for (const auto& hit : caloHits) {
    const unsigned int layer = m_bitFieldCoder->get(hit.getCellID(), m_layerIndex);

    // Polar angle, symmetrized around pi/2 to match the threshold maps.
    double hitTheta = edm4hep::utils::anglePolar(hit.getPosition());
    if (hitTheta > TMath::Pi() / 2.) {
      hitTheta = TMath::Pi() - hitTheta;
    }

    double modeThreshold = 0.;
    double stddev = 0.;
    if (m_thresholdMap) {
      // FindFixBin (unlike FindBin) never extends the axis, so the lookup is const.
      const int binx = m_thresholdMap->GetXaxis()->FindFixBin(hitTheta);
      const int biny = m_thresholdMap->GetYaxis()->FindFixBin(layer);
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
    const auto simIt = hitToSim.find(hit.getObjectID());
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
