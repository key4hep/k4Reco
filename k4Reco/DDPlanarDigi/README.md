<!--
Copyright (c) 2020-2024 Key4hep-Project.

This file is part of Key4hep.
See https://key4hep.github.io/key4hep-doc/ for further info.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
-->
This is a port of the DDPlanarDigi processor. The original header and source
files can be found in
https://github.com/iLCSoft/MarlinTrkProcessors/blob/master/source/Digitisers/include/DDPlanarDigiProcessor.h
and
https://github.com/iLCSoft/MarlinTrkProcessors/blob/master/source/Digitisers/src/DDPlanarDigiProcessor.cc

It's very similar to the original one. The (possibly incomplete) list of differences would be:
- Different random number generator and seed (makes it impossible to compare on
  an event-by-event level). GSL was used in the processor, TRandom3 is being
  used in the algorithm that also uses the UniqueIDGenSvc to seed each event.
- The detector instance from DD4hep is obtained with GeoSvc in the algorithm.
- The algorithm has an extra option (`CellIDBits`) to only use a certain number
  of bits for the cell IDs that are obtained from the hits (see
  https://github.com/key4hep/k4Reco/pull/25).
