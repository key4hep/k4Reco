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
#ifndef K4RECO_LUMICALHIT_H
#define K4RECO_LUMICALHIT_H 1

#include "ProjectionInfo.h"

#include <edm4hep/CalorimeterHit.h>

/** Class to hold position and energy information of CalorimeterHits in the LumiCal

 * To use of projected hits, which derive from multiple CalorimeterHits this
 * wrapper is necessary to assign the original CalorimeterHits to the created
 * cluster
 */

class LumiCalHit {
public:
  LumiCalHit() {}
  LumiCalHit(int cellIDProjection, ProjectionInfo const& projection)
      : m_cellID0(cellIDProjection), m_cellID1(projection.getCellIdHitZ()), m_energy(projection.getEnergy()),
        m_position{projection.getPosition()[0], projection.getPosition()[1], projection.getPosition()[2]} {
    for (size_t i = 0; i < projection.getCaloHits().size(); ++i) {
      m_caloHits.insert(i);
    }
  }

  LumiCalHit(const LumiCalHit& other) = delete;
  LumiCalHit& operator=(const LumiCalHit& other) = delete;
  LumiCalHit(LumiCalHit&& other) = default;
  LumiCalHit& operator=(LumiCalHit&& other) = default;
  ~LumiCalHit() = default;

private:
  int m_cellID0 = 0;
  int m_cellID1 = 0;
  double m_energy = 0.0;
  double m_position[3] = {0.0, 0.0, 0.0};
  /// original CalorimeterHits indexes (in the global coordinate system)
  std::set<size_t> m_caloHits{};

public:
  void addHit(const size_t index) { m_caloHits.insert(index); }

  const auto& getHits() const { return m_caloHits; }
  auto beginHits() const { return m_caloHits.begin(); }
  auto endHits() const { return m_caloHits.end(); }

  void setEnergy(double e) { m_energy = e; }
  void setCellID0(int i) { m_cellID0 = i; }
  void setCellID1(int i) { m_cellID1 = i; }
  void setPosition(double const* pos) {
    m_position[0] = pos[0];
    m_position[1] = pos[1];
    m_position[2] = pos[2];

  } // fixme memcopy??

  double getEnergy() const { return m_energy; }
  int getCellID0() const { return m_cellID0; }
  int getCellID1() const { return m_cellID1; }
  const double* getPosition() const { return m_position; }
};

#endif // K4RECO_LUMICALHIT_H
