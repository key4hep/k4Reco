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

// Local
#include "Global.h"
#include "LumiCalClusterer.h"
#include "LumiCalHit.h"

#include <edm4hep/CalorimeterHitCollection.h>
#include <edm4hep/SimCalorimeterHitCollection.h>

// Stdlib
#include <memory>
#include <string>

/* --------------------------------------------------------------------------
   Loop over al hits in the LCCollection and write the hits into vectors
   of CalorimeterHitImpl. Hits are split in two vectors, one for each arm
   of LumiCal.
   -------------------------------------------------------------------------- */
// TODO: CalorimeterHitCollection or SimCalorimeterHitCollection?
std::pair<int, edm4hep::CalorimeterHitCollection>
LumiCalClustererClass::getCalHits(const edm4hep::SimCalorimeterHitCollection& col, MapIntMapIntVCalHit& calHits) {

  if (col.size() < static_cast<size_t>(m_clusterMinNumHits))
    return {0, edm4hep::CalorimeterHitCollection()};

  auto calohits = createCaloHitCollection(col);

  for (size_t i = 0; i < calohits.size(); ++i) {
    const auto& calHitIn = calohits[i];

    int arm(0), layer(0);
    int rCell(0), phiCell(0);

    const double engyHit = static_cast<double>(calHitIn.getEnergy());

    if (engyHit < m_hitMinEnergy)
      continue;

    arm = m_mydecoder->get(calHitIn.getCellID(), "barrel"); // from 1 and 2
    if (arm == 2)
      arm = -1;

    phiCell = m_mydecoder->get(calHitIn.getCellID(), "phi"); // goes from -phiMax/2 to +phiMax/2-1
    if (arm < 0) {
      // for rotation around the Y-axis, so that the phiCell increases counter-clockwise for z<0
      phiCell *= -1;
      // This is not needed any more because we no currently do not calculate pad positions from IDs any longer
      //  now it is just important that we have the correct phiID range
      // LumiCal is (or not, if we fix it) rotated by pi around Z for negative side
      // phiCell += int(m_gmc.m_backwardRotationPhi/(2.0*M_PI)*_cellPhiMax+0.5);
    }
    // limit to range 0 to m_cellPhiMax-1
    if (phiCell >= m_cellPhiMax)
      phiCell -= m_cellPhiMax;
    if (phiCell < 0)
      phiCell += m_cellPhiMax;

    rCell = m_mydecoder->get(calHitIn.getCellID(), "r");     // from 0
    layer = m_mydecoder->get(calHitIn.getCellID(), "layer"); // from 1

    // Calculate internal cellID
    int cellId = cellIdZPR(layer, phiCell, rCell, arm);

    // skip this hit if the following conditions are met
    if (layer >= static_cast<int>(m_maxLayerToAnalyse) || layer < 0)
      continue;

    const auto locPos = rotateToLumiCal(calHitIn.getPosition());

    // streamlog_message(DEBUG2, std::stringstream p;
    //                   p << std::scientific << std::setprecision(3) << "\t Arm, CellId, Pos(x,y,z), hit energy [MeV]:
    //                   "
    //                     << std::setw(5) << arm << std::setw(13) << cellId << "\t (" << std::setw(13) << locPos[0]
    //                     << ", " << std::setw(13) << locPos[1] << ", " << std::setw(13) << locPos[2] << "), "
    //                     << 1000. * engyHit << "Layer/Phi/R/Arm" << std::setw(5) << layer << std::setw(5) << phiCell
    //                     << std::setw(5) << rCell << std::setw(5) << arm << std::endl;
    //                   , p.str(););

    // create a new LumiCalHit
    auto calHitNew = std::make_shared<LumiCalHit>();

    // write the parameters to the new LumiCalHit
    calHitNew->setCellID0(cellId);
    calHitNew->setEnergy(engyHit);
    calHitNew->setPosition(locPos);
    calHitNew->addHit(i);

    // add the LumiCalHit to a vector according to the detector
    // arm, and sum the total collected energy at either arm
    calHits[arm][layer].push_back(std::move(calHitNew));
    m_numHitsInArm[arm]++;
    m_totEngyArm[arm] += engyHit;
  } // for all simHits

  // streamlog_out(DEBUG4) << "Energy deposit: " << m_totEngyArm[-1] << "\t" << m_totEngyArm[1] << "\n"
  //                       << "Number of hits: " << m_numHitsInArm[-1] << "\t" << m_numHitsInArm[1] << std::endl;

  if (((m_numHitsInArm[-1] < m_clusterMinNumHits) || (m_totEngyArm[-1] < m_minClusterEngyGeV)) &&
      ((m_numHitsInArm[1] < m_clusterMinNumHits) || (m_totEngyArm[1] < m_minClusterEngyGeV))) {
    return {0, std::move(calohits)};
  } else {
    return {1, std::move(calohits)};
  }
}

edm4hep::CalorimeterHitCollection
LumiCalClustererClass::createCaloHitCollection(const edm4hep::SimCalorimeterHitCollection& input) const {
  m_alg->debug() << "Creating the CalorimeterHit collection with dummy digitization" << endmsg;

  const double calibrationFactor = getCalibrationFactor();

  auto caloHitCollection = edm4hep::CalorimeterHitCollection();
  for (const auto& simCaloHit : input) {
    auto calHit = caloHitCollection.create();

    calHit.setCellID(simCaloHit.getCellID());
    calHit.setEnergy(simCaloHit.getEnergy() * calibrationFactor);
    calHit.setPosition(simCaloHit.getPosition());
  }

  return caloHitCollection;
}
