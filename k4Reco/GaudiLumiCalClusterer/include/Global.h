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
#ifndef K4RECO_GLOBAL_H
#define K4RECO_GLOBAL_H

#include <map>
#include <memory>
#include <vector>

class LumiCalHit;
class LCCluster;
class VirtualCluster;

using CalHit = std::shared_ptr<LumiCalHit>;
using VecCalHit = std::vector<CalHit>;
using VDouble = std::vector<double>;
using VInt = std::vector<int>;

using MapIntCalHit = std::map<int, CalHit>;
using MapIntLCCluster = std::map<int, LCCluster>;

using MapIntVCalHit = std::map<int, VecCalHit>;
using MapIntVDouble = std::map<int, VDouble>;
using MapIntVInt = std::map<int, VInt>;

using MapIntMapIntLCCluster = std::map<int, MapIntLCCluster>;
using MapIntMapIntVCalHit = std::map<int, MapIntVCalHit>;
using MapIntMapIntVDouble = std::map<int, MapIntVDouble>;
using MapIntMapIntVInt = std::map<int, MapIntVInt>;
using MapIntVirtualCluster = std::map<int, VirtualCluster>;
using MapIntDouble = std::map<int, double>;
using MapIntInt = std::map<int, int>;

using MapIntMapIntCalHit = std::map<int, MapIntCalHit>;

using VMapIntCalHit = std::vector<MapIntCalHit>;
using VMapIntInt = std::vector<MapIntInt>;
using VMapIntLCCluster = std::vector<MapIntLCCluster>;
using VMapIntVInt = std::vector<MapIntVInt>;
using VMapIntVirtualCluster = std::vector<MapIntVirtualCluster>;

using VVDouble = std::vector<VDouble>;

#endif // K4RECO_GLOBAL_H
