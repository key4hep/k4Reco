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
#ifndef K4RECO_VIRTUALCLUSTER_h
#define K4RECO_VIRTUALCLUSTER_h

#include <array>
#include <iostream>

class VirtualCluster {

public:
  VirtualCluster();
  VirtualCluster(double x, double y, double z);

  void clear();

  inline double getX() const { return m_position[0]; }
  inline double getY() const { return m_position[1]; }
  inline double getZ() const { return m_position[2]; }

  inline void setX(double x) { m_position[0] = x; }
  inline void setY(double y) { m_position[1] = y; }
  inline void setZ(double z) { m_position[2] = z; }

  inline const std::array<double, 3> getPosition() const { return m_position; }

  friend std::ostream& operator<<(std::ostream& o, const VirtualCluster& rhs);

private:
  std::array<double, 3> m_position; // Not sure the first parameter is Z or E, or Distance??
};

#endif // K4RECO_VIRTUALCLUSTER_h
