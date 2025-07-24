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
#include "ProjectionInfo.h"
#include "LumiCalHit.h"

ProjectionInfo::ProjectionInfo() : m_energy(0.0), m_cellIdHitZ(0), newObject(true) {
  m_position[0] = 0.0;
  m_position[1] = 0.0;
  m_position[2] = 0.0;
}

ProjectionInfo::ProjectionInfo(CalHit const& calHit, int cellIdZ)
    : m_energy(calHit->getEnergy()), m_cellIdHitZ(cellIdZ), newObject(false) {
  m_position[0] = calHit->getPosition()[0];
  m_position[1] = calHit->getPosition()[1];
  m_position[2] = calHit->getPosition()[2];
}

void ProjectionInfo::addHit(CalHit const& calHit) { m_energy += calHit->getEnergy(); }
