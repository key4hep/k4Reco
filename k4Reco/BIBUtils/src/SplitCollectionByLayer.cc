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
#include "SplitCollectionByLayer.h"

#include <DDSegmentation/BitFieldCoder.h>

#include <cstddef>
#include <string>
#include <vector>

SplitCollectionByLayer::SplitCollectionByLayer(const std::string& name, ISvcLocator* svcLoc)
    : Transformer(name, svcLoc, KeyValues("InputCollection", {"VBTrackerHits"}),
                  KeyValues("OutputCollections", {"VBTrackerHitsInner", "VBTrackerHitsOuter"})) {}

StatusCode SplitCollectionByLayer::initialize() {
  m_geoSvc = serviceLocator()->service(m_geoSvcName);
  if (!m_geoSvc) {
    error() << "Unable to retrieve the GeoSvc" << endmsg;
    return StatusCode::FAILURE;
  }

  // One [start, end] layer interval must be given per output collection. The
  // number of output collections itself is enforced by the framework at the
  // first event, when the size of the returned vector is matched against the
  // OutputCollections list.
  if (m_startLayers.empty() || m_startLayers.size() != m_endLayers.size()) {
    error() << "StartLayers (" << m_startLayers.size() << ") and EndLayers (" << m_endLayers.size()
            << ") must be non-empty and have one entry per output collection" << endmsg;
    return StatusCode::FAILURE;
  }

  for (std::size_t i = 0; i < m_startLayers.size(); ++i) {
    if (m_startLayers[i] > m_endLayers[i]) {
      error() << "Output collection " << i << " has StartLayer " << m_startLayers[i] << " > EndLayer " << m_endLayers[i]
              << endmsg;
      return StatusCode::FAILURE;
    }
    info() << "Output collection " << i << " will collect layers [" << m_startLayers[i] << ", " << m_endLayers[i] << "]"
           << endmsg;
  }

  return StatusCode::SUCCESS;
}

std::vector<edm4hep::TrackerHitPlaneCollection>
SplitCollectionByLayer::operator()(const edm4hep::TrackerHitPlaneCollection& trackerHits) const {
  const std::string encoderString = m_geoSvc->constantAsString(m_encodingStringVariable.value());
  const dd4hep::DDSegmentation::BitFieldCoder bitFieldCoder(encoderString);

  const std::size_t nOutputs = m_startLayers.size();
  std::vector<edm4hep::TrackerHitPlaneCollection> outColls(nOutputs);
  for (auto& coll : outColls) {
    coll.setSubsetCollection();
  }

  for (const auto& hit : trackerHits) {
    const int layer = bitFieldCoder.get(hit.getCellID(), "layer");
    for (std::size_t i = 0; i < nOutputs; ++i) {
      if (layer >= m_startLayers[i] && layer <= m_endLayers[i]) {
        outColls[i].push_back(hit);
      }
    }
  }

  if (msgLevel(MSG::DEBUG)) {
    for (std::size_t i = 0; i < nOutputs; ++i) {
      debug() << outColls[i].size() << " hits routed to output collection " << i << endmsg;
    }
  }

  return outColls;
}

DECLARE_COMPONENT(SplitCollectionByLayer)
