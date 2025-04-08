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
#ifndef K4RECO_PARAMETERS_H
#define K4RECO_PARAMETERS_H

#include <map>
#include <string>
#include <vector>

struct Parameters {
  using ParMap = std::map<std::string, double>;
  using StringVec = std::vector<std::string>;

public:
  Parameters(std::vector<int> const& collections, double maxCellAngle, double maxCellAngleRZ, double chi2cut,
             int minClustersOnTrack, double maxDistance, double maxSlopeZ, double highPTcut, bool highPTfit,
             bool onlyZSchi2cut, bool radialSearch, bool vertexToTracker, bool kalmanFitForward, int step, bool combine,
             bool build, bool extend, bool sortTracks)
      : m_collections(collections), m_maxCellAngle(maxCellAngle), m_maxCellAngleRZ(maxCellAngleRZ), m_chi2cut(chi2cut),
        m_minClustersOnTrack(minClustersOnTrack), m_maxDistance(maxDistance), m_maxSlopeZ(maxSlopeZ),
        m_highPTcut(highPTcut), m_highPTfit(highPTfit), m_onlyZSchi2cut(onlyZSchi2cut), m_radialSearch(radialSearch),
        m_vertexToTracker(vertexToTracker), m_kalmanFitForward(kalmanFitForward), m_step(step), m_combine(combine),
        m_build(build), m_extend(extend), m_sortTracks(sortTracks) {}

  Parameters(Parameters const&) = default;
  Parameters& operator=(Parameters const&) = delete;
  Parameters(Parameters&&) = default;
  Parameters& operator=(Parameters&&) = delete;
  ~Parameters() = default;

  std::vector<int> m_collections; /// which collections to combine
  double m_maxCellAngle;
  double m_maxCellAngleRZ;
  double m_chi2cut;
  int m_minClustersOnTrack;
  double m_maxDistance;
  double m_maxSlopeZ; /// Cut on the slope in the longitudinal plane for seeding
  double m_highPTcut;
  bool m_highPTfit;
  bool m_onlyZSchi2cut;
  bool m_radialSearch;
  bool m_vertexToTracker;
  bool m_kalmanFitForward;
  int m_step;
  bool m_combine;
  bool m_build;
  bool m_extend;
  bool m_sortTracks;
  double m_tightenStep = 1;

  const StringVec m_existingFunctions = {
      "CombineCollections",
      "ExtendTracks",
      "BuildNewTracks",
      "SortTracks",
  };
  const StringVec m_existingFlags = {
      "HighPTFit", "OnlyZSchi2cut", "RadialSearch", "VertexToTracker", "KalmanFitForward", "KalmanFitBackward",
  };
  const StringVec m_existingParameters = {
      "MaxCellAngle", "MaxCellAngleRZ", "Chi2Cut", "MinClustersOnTrack", "MaxDistance", "SlopeZRange", "HighPTCut",
  };

  void tighten() {
    double factor = (10.0 - m_tightenStep) / (10.0 - (m_tightenStep - 1.0));
    m_maxCellAngle *= factor;
    m_maxCellAngleRZ *= factor;
    m_chi2cut *= factor;
    m_tightenStep += 1.0;
  }

private:
  void check(StringVec const& values, StringVec const& options, std::string const& type);
  void check(ParMap const& values, StringVec const& options, std::string const& type);
};

#endif // K4RECO_PARAMETERS_H
