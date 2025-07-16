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
#include "GaudiLumiCalClusterer.h"

#include <DD4hep/Detector.h>
#include <DD4hep/Readout.h>
#include <DD4hep/detail/ObjectsInterna.h>

#include <edm4hep/CalorimeterHitCollection.h>
#include <edm4hep/ClusterCollection.h>
#include <edm4hep/ReconstructedParticleCollection.h>

#include <string>
#include <utility>

GaudiLumiCalClusterer::GaudiLumiCalClusterer(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(name, svcLoc,
                       {
                           KeyValues("LumiCal_Collection", {"LumiCalCollection"}),
                       },
                       {
                           KeyValues("LumiCal_Hits", {"LumiCalHits"}),
                           KeyValues("LumiCal_Clusters", {"LumiCalClusters"}),
                           KeyValues("LumiCal_RecoParticles", {"LumiCalRecoParticles"}),
                       }),
      LumiCalClusterer(this) {}

StatusCode GaudiLumiCalClusterer::initialize() {

  m_geoSvc = serviceLocator()->service(m_geoSvcName);
  if (!m_geoSvc) {
    error() << "Unable to retrieve GeoSvc" << endmsg;
    return StatusCode::FAILURE;
  }
  const auto readout = m_geoSvc->getDetector()->readout("LumiCalCollection");
  const auto fieldDescription = readout->id.fieldDescription();

  LumiCalClusterer.createDecoder(fieldDescription);

  // Set parameters for the GlobalMethodsClass
  std::map<std::string, std::variant<int, float, std::string>> parameters;
  parameters["ZLayerPhiOffset"] = static_cast<float>(m_zLayerPhiOffset);
  parameters["EnergyCalibConst"] = static_cast<float>(m_EnergyCalibConst);
  parameters["LogWeigthConstant"] = static_cast<float>(m_logWeigthConstant);
  parameters["MinHitEnergy"] = static_cast<float>(m_minHitEnergy);
  parameters["MiddleEnergyHitBoundFrac"] = static_cast<float>(m_MiddleEnergyHitBoundFrac);
  parameters["ElementsPercentInShowerPeakLayer"] = static_cast<float>(m_ElementsPercentInShowerPeakLayer);
  parameters["ClusterMinNumHits"] = static_cast<int>(m_ClusterMinNumHits);
  parameters["MoliereRadius"] = static_cast<float>(m_rMoliere);
  parameters["MinClusterEngy"] = static_cast<float>(m_minClusterEngy);
  parameters["WeightingMethod"] = m_WeightingMethod;
  parameters["NumOfNearNeighbor"] = static_cast<int>(m_NumOfNearNeighbor);

  gmc.SetConstants(parameters);

  m_BeamCrossingAngle = gmc.GlobalParamD[GlobalMethodsClass::BeamCrossingAngle] / 2.;

  // printParameters();
  // /* --------------------------------------------------------------------------
  //    Print out Processor Parameters
  //    -------------------------------------------------------------------------- */
  // info() << "Global parameters for Processor:" << name() << "\t" << type() << endmsg;
  // gmc.PrintAllParameters();
  // info() << endmsg;

  LumiCalClusterer.init(gmc);
  LumiCalClusterer.setCutOnFiducialVolume(m_cutOnFiducialVolume);

  return StatusCode::SUCCESS;
}

std::tuple<edm4hep::CalorimeterHitCollection, edm4hep::ClusterCollection, edm4hep::ReconstructedParticleCollection>
GaudiLumiCalClusterer::operator()(const edm4hep::SimCalorimeterHitCollection& input) const {
  /* --------------------------------------------------------------------------
     create clusters using: LumiCalClustererClass
     -------------------------------------------------------------------------- */
  // TODO: Remove const_cast
  auto [status, calhits] = LumiCalClusterer.processEvent(const_cast<edm4hep::SimCalorimeterHitCollection&>(input));

  auto LCalClusterCol = edm4hep::ClusterCollection();

  auto LCalRPCol = edm4hep::ReconstructedParticleCollection();

  if (status == RETVAL::OK) {
    debug() << " Transfering reco results to LCalClusterCollection....." << endmsg;

    for (const int armNow : {-1, 1}) {
      const auto& pairIDCellsVector = LumiCalClusterer.m_superClusterIdToCellId.at(armNow);
      debug() << " Arm  " << std::setw(4) << armNow << "\t Number of clusters: " << pairIDCellsVector.size() << endmsg;

      for (const auto& pairIDCells : pairIDCellsVector) {
        const int clusterId = pairIDCells.first;
        LCCluster& thisClusterInfo =
            const_cast<LCCluster&>(LumiCalClusterer.m_superClusterIdClusterInfo.at(armNow).at(clusterId));
        thisClusterInfo.recalculatePositionFromHits(gmc);
        const auto& objectTuple = gmc.getLCIOObjects(thisClusterInfo, m_minClusterEngy, m_cutOnFiducialVolume, calhits);
        if (!std::get<0>(objectTuple).has_value())
          continue;

        LCalClusterCol.push_back(std::get<0>(objectTuple).value());
        LCalRPCol.push_back(std::get<1>(objectTuple).value());
      }
    }
  }

  return std::make_tuple(std::move(calhits), std::move(LCalClusterCol), std::move(LCalRPCol));
}
