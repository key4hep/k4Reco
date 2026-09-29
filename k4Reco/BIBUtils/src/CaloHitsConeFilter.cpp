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
#include "CaloHitsConeFilter.h"

#include <edm4hep/SimCalorimeterHit.h>
#include <edm4hep/utils/vector_utils.h>

#include <podio/ObjectID.h>

#include <unordered_map>

CaloHitsConeFilter::CaloHitsConeFilter(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(name, svcLoc,
                       {KeyValue("MCParticleCollectionName", "MCParticle"),
                        KeyValue("CaloHitCollectionName", "EcalBarrelCollectionRec"),
                        KeyValue("CaloRelationCollectionName", "EcalBarrelRelationsSimRec")},
                       {KeyValue("GoodHitCollection", "EcalBarrelCollectionConed"),
                        KeyValue("GoodRelationCollection", "EcalBarrelRelationsSimConed")}) {}

std::tuple<edm4hep::CalorimeterHitCollection, edm4hep::CaloHitSimCaloHitLinkCollection>
CaloHitsConeFilter::operator()(const edm4hep::MCParticleCollection& mcParticles,
                               const edm4hep::CalorimeterHitCollection& caloHits,
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
    const auto& hitPos = hit.getPosition();
    const edm4hep::Vector3d pos{hitPos.x, hitPos.y, hitPos.z};

    bool save = false;
    for (const auto& part : mcParticles) {
      // Keep only the generator-level particles.
      if (part.getGeneratorStatus() != 1) {
        continue;
      }
      const auto& mom = part.getMomentum();
      const double deltaR = edm4hep::utils::angleBetween(mom, pos);
      if (deltaR < m_coneSize) {
        save = true;
        break;
      }
    }

    if (!save) {
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

  debug() << nAccepted << " / " << caloHits.size() << " calo hits kept inside the cone" << endmsg;

  return std::make_tuple(std::move(outHits), std::move(outLinks));
}

DECLARE_COMPONENT(CaloHitsConeFilter)
