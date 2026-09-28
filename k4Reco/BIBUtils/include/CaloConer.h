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
#ifndef K4RECO_CALOCONER_H
#define K4RECO_CALOCONER_H 1

#include "Gaudi/Property.h"

#include <edm4hep/CaloHitSimCaloHitLinkCollection.h>
#include <edm4hep/CalorimeterHitCollection.h>
#include <edm4hep/MCParticleCollection.h>

#include <k4FWCore/Transformer.h>

#include <string>
#include <tuple>

/** === CaloConer ===
 *  Keeps only the calorimeter hits that fall within a fixed angular cone
 *  (ConeWidth, in radians) around the direction of any generator-level
 *  (generatorStatus == 1) MC particle. The selected hits are written to a
 *  subset collection together with a freshly built reco-to-sim link collection.
 *
 *  Gaudi-native port of the Marlin CaloConer processor from MyBIBUtils
 *  (https://github.com/madbaron/MyBIBUtils), used to clean beam-induced
 *  background (BIB) out of the calorimeters before particle flow.
 *
 *  @author F. Meloni, DESY (original Marlin processor)
 */
struct CaloConer final : k4FWCore::MultiTransformer<
                             std::tuple<edm4hep::CalorimeterHitCollection, edm4hep::CaloHitSimCaloHitLinkCollection>(
                                 const edm4hep::MCParticleCollection&, const edm4hep::CalorimeterHitCollection&,
                                 const edm4hep::CaloHitSimCaloHitLinkCollection&)> {
  CaloConer(const std::string& name, ISvcLocator* svcLoc);

  std::tuple<edm4hep::CalorimeterHitCollection, edm4hep::CaloHitSimCaloHitLinkCollection>
  operator()(const edm4hep::MCParticleCollection& mcParticles, const edm4hep::CalorimeterHitCollection& caloHits,
             const edm4hep::CaloHitSimCaloHitLinkCollection& caloLinks) const override;

private:
  Gaudi::Property<double> m_coneSize{this, "ConeWidth", 0.2,
                                     "Half-opening angle of the cone around MC particles [rad]"};
};

#endif
