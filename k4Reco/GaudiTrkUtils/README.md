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
# Reimplementation of Utilities from MarlinTrk and MarlinTrkUtils

The following files have been reimplemented from `MarlinTrk`:

- `MarlinDDKalTest`: Now called `GaudiDDKalTest`. Anything related to
  `MarlinTrk` has been removed. The rest is very similar, changing LCIO Tracker
  hits to EDM4hep tracker hits. A pointer to the current algorithm has been
  added to be able to use logging.

  The following have not been implemented:
  ```cpp
  MarlinTrk::IMarlinTrack* createTrack()
  std::string MarlinTrk::name()
  ```

- `MarlinDDKalTestTrack`: now called `GaudiDDKalTestTrack`. Similar changes as
  explained above.

  The following have not been implemented, since they have not been needed so far (some of these may be overloads):
  ```cpp
  void setMass(double mass)
  double getMass()
  int initialise(bool fitDirection)
  int smooth()
  int testChi2Increment(EVENT::TrackerHit* hit, double& chi2increment)
  int getTrackState(IMPL::TrackStateImpl& ts, double& chi2, int& ndf)
  int propagate(const edm4hep::Vector3d& point, edm4hep::TrackState& ts, double& chi2, int& ndf)
  int propagateToLayer(int layerID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID, int mode=modeClosest)
  int propagateToDetElement(int detElementID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode=modeClosest)
  int propagateToDetElement(int detEementID, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode=modeClosest)
  int propagateToDetElement(int detEementID, const TKalTrackSite& site, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode=modeClosest)
  int extrapolate(const Vector3D& point, IMPL::TrackStateImpl& ts, double& chi2, int& ndf)
  int extrapolate(const Vector3D& point, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf)
  int extrapolate(const Vector3D& point, const TKalTrackSite& site, IMPL::TrackStateImpl& ts, double& chi2, int& ndf)
  int extrapolateToLayer(int layerID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID, int mode=modeClosest)
  int extrapolateToLayer(int layerID, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID, int mode=modeClosest)
  int extrapolateToLayer(int layerID, const TKalTrackSite& site, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID, int mode=modeClosest)
  int extrapolateToDetElement(int detElementID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode=modeClosest)
  int extrapolateToDetElement(int detEementID, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode=modeClosest)
  int extrapolateToDetElement(int detEementID, const TKalTrackSite& site, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode=modeClosest)
  std::string toString()
  int intersectionWithLayer(int layerID, Vector3D& point, int& detElementID, int mode=modeClosest)
  int intersectionWithLayer(int layerID, EVENT::TrackerHit* hit, Vector3D& point, int& detElementID, int mode=modeClosest)
  int intersectionWithDetElement(int detElementID, Vector3D& point, int mode=modeClosest)
  int intersectionWithDetElement(int detElementID, EVENT::TrackerHit* hit, Vector3D& point, int mode=modeClosest)
  int intersectionWithDetElement(int detElementID, const TKalTrackSite& site, Vector3D& point, const DDVMeasLayer*& ml, int mode=modeClosest)
  ```

- `MarlinTrkUtils`: now called `GaudiTrkUtils`.
  - Add a function to set the streamlog level. This is mandatory for any
  algorithm that uses streamlog, since otherwise the maximum debug level will be
  in place without the messages printing anywhere but doing all the computations
  needed, leading to very bad performance.
  - Use references instead of pointers when possible
  - Add some members in the Gaudi implementation: the Gaudi algorithm (to be
    able to use logging from Gaudi), `GaudiDDKalTest`, a pointer to the geo
    service and where the encoding string is stored
