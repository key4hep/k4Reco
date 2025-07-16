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
#include "VirtualCluster.h"

#include <iomanip>
#include <iostream>

VirtualCluster::VirtualCluster() : _position{0.0, 0.0, 0.0} {}

VirtualCluster::VirtualCluster(double x, double y, double z) : _position{x, y, z} {}

void VirtualCluster::clear() {
  _position[0] = 0.0;
  _position[1] = 0.0;
  _position[2] = 0.0;
}

std::ostream& operator<<(std::ostream& o, const VirtualCluster& rhs) {
  o << std::setw(13) << rhs._position[0] << std::setw(13) << rhs._position[1] << std::setw(13) << rhs._position[2];
  return o;
}
