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
In this folder there are reimplementations of the other tracking processors that are used in the CLD reconstruction that are not ConformalTracking

* ClonesAndSplitTracksFinder

Originally ported from https://github.com/iLCSoft/MarlinTrkProcessors/blob/master/source/Refitting/include/ClonesAndSplitTracksFinder.h and https://github.com/iLCSoft/MarlinTrkProcessors/blob/master/source/Refitting/src/ClonesAndSplitTracksFinder.cc

Only the part that runs when `mergeSplitTracks` is false (default value in the CLD reconstruction) has been ported. The rest can be ported if there is demand.

Additional validation may be required, since in the simulations that have been done there is not overlap between hits.

* RefitFinal

Originally ported from
https://github.com/iLCSoft/MarlinTrkProcessors/blob/master/source/Refitting/include/RefitFinal.h
and
https://github.com/iLCSoft/MarlinTrkProcessors/blob/master/source/Refitting/src/RefitFinal.cc

The resulting `edm4hep::TrackMCParticleLinkCollection` need to be validated.

* TruthTrackFinder

Originally ported from
https://github.com/iLCSoft/MarlinTrkProcessors/blob/master/source/Refitting/include/TruthTrackFinder.h
and
https://github.com/iLCSoft/MarlinTrkProcessors/blob/master/source/Refitting/src/TruthTrackFinder.cc

Only the part that runs when `UseTruthInPrefit` is false (default value in the CLD reconstruction) has been ported. The rest can be ported if there on demand.

The member functions `getSubdetector` and `getLayer` have been implemented in `removeHitsSameLayer`.
