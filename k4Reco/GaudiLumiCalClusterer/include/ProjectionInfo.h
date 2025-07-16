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
#ifndef K4RECO_PROJECTIONINFO_H
#define K4RECO_PROJECTIONINFO_H 1

#include <edm4hep/CalorimeterHit.h>

#include <memory>
#include <set>

class LumiCalHit;

using CalHit = std::shared_ptr<LumiCalHit>;

class ProjectionInfo {

public:
  ProjectionInfo();
  ProjectionInfo(CalHit const& calHit, int cellIdHitZ);
  ProjectionInfo(const ProjectionInfo& other) = delete;
  ProjectionInfo& operator=(const ProjectionInfo& other) = delete;
  ProjectionInfo(ProjectionInfo&& other) = default;
  ProjectionInfo& operator=(ProjectionInfo&& other) = default;
  ~ProjectionInfo() = default;

  void addHit(const CalHit& calHit);

  const double* getPosition() const { return position; }
  double getEnergy() const { return energy; }
  int getCellIdHitZ() const { return cellIdHitZ; }
  std::set<const edm4hep::CalorimeterHit*> const& getCaloHits() const { return hits; }

private:
  double energy;
  double position[3];
  int cellIdHitZ;
  std::set<const edm4hep::CalorimeterHit*> hits{};

public:
  bool newObject;
};

#endif // K4RECO_PROJECTIONINFO_H
