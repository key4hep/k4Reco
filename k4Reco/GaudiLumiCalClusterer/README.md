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
The `SKipNEvents` property is removed, since Gaudi takes care of that.

Only the parts that use DD4hep have been implemented

GlobalMethodsClass: Remove the members `_forwardCalo` and `_backwardCalo`

GaudiLumiCalClusterer: Originally `MarlinLumiCalClusterer`.
`TryMarlinLumiCalClusterer` has been moved inside `operator()`. The
`OutputManager` has not been implemented.

LumiCalClusterer: Remove templates where they are not needed. Removed
```
void setLumiCollectionName(std::string const& lumiNameNow)
void setLumiOutCollectionName(std::string const& name)
void getThetaPhiZCluster(MapIntCalHit const& calHitsCellId, VInt const& clusterIdToCellId, double totEngy, double * output);
template <class T> double posWeight(T const& calHit, double totEngy, GlobalMethodsClass::WeightingMethod_t method) const;
template <class T> double posWeight(T const& calHit, double totEngy, GlobalMethodsClass::WeightingMethod_t method, double logWeightConstNow) const;
double distance2DPolar(double * pos1, double * pos2 );
double getDistanceAroundCMWithEnergyPercent(LCCluster const& clusterCM, VInt const& clusterIdToCellId, MapIntCalHit const& calHitsCellId, double engyPercentage);
double getMoliereRadius(MapIntCalHit const& calHitsCellId, VInt const& clusterIdToCellId, LCCluster const& clusterCM );
```

GlobalMethodsClass: Removed
```
void thetaPhiCell(const int cellId, std::map<GlobalMethodsClass::Coordinate_t, double>& thetaPhiCell) const;
double toSignal(const double valNow) const;
static std::string getParameterName(Parameter_t par);
```

VirtualCluster: The class has been removed since it's a wrapper around a
3-vector, so `edm4hep::Vector3d` is used instead.

