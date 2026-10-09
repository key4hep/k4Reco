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
#include "CaloHitConeFilter.h"

#include "LinkUtils.h"

#include <edm4hep/utils/vector_utils.h>

CaloHitConeFilter::CaloHitConeFilter(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(name, svcLoc,
                       {KeyValue("MCParticleCollectionName", "MCParticle"),
                        KeyValue("CaloHitCollectionName", "EcalBarrelCollectionRec"),
                        KeyValue("CaloRelationCollectionName", "EcalBarrelRelationsSimRec")},
                       {KeyValue("GoodHitCollection", "EcalBarrelCollectionConed"),
                        KeyValue("GoodRelationCollection", "EcalBarrelRelationsSimConed")}) {}

std::tuple<edm4hep::CalorimeterHitCollection, edm4hep::CaloHitSimCaloHitLinkCollection>
CaloHitConeFilter::operator()(const edm4hep::MCParticleCollection& mcParticles,
                              const edm4hep::CalorimeterHitCollection& caloHits,
                              const edm4hep::CaloHitSimCaloHitLinkCollection& caloLinks) const {
  edm4hep::CalorimeterHitCollection outHits;
  outHits.setSubsetCollection();
  edm4hep::CaloHitSimCaloHitLinkCollection outLinks;

  const auto linksByHit = k4reco::bibutils::linksByFrom(caloLinks);

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
    if (const auto linksIt = linksByHit.find(hit.getObjectID()); linksIt != linksByHit.end()) {
      k4reco::bibutils::copyLinks(linksIt->second, outLinks);
    }
    ++nAccepted;
  }

  debug() << nAccepted << " / " << caloHits.size() << " calo hits kept inside the cone" << endmsg;

  return std::make_tuple(std::move(outHits), std::move(outLinks));
}

DECLARE_COMPONENT(CaloHitConeFilter)
