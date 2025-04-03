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
Reimplementation of the ConformalTracking processor from iLCSoft

- Conformal Tracking A few places have some value hardcoded from ILD, like
  `lcio::LCTrackerCellID::subdet()`. How to modify these is not clear since
  there are some functions that take them. If you want to know every place where
  those are, look for `ILD` in `ConformalTracking.cpp`.

  Some functions have been
  changed to return something instead of passing the value that is going to be
  returned by reference

- KDTrack
  Remove m_kalmanTrack that was never being set (and the corresponding code in ConformalTracking.cpp)
  Remove `m_nPoints` and use `m_clusters.size()` instead of keeping an additional counter

- HelixTrack
  Remove HelixTrack(const double* position, const double* p, double charge, double Bz);

- KDCluster
  Remove a few members and their corresponding setters and getters that were not used
