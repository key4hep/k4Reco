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

#ifndef K4RECO_GAUDITRKUTILS_H
#define K4RECO_GAUDITRKUTILS_H

#include "GaudiDDKalTest.h"
#include "GaudiDDKalTestTrack.h"

#include <streamlog/logstream.h>
#include <streamlog/streamlog.h>

#include <edm4hep/CovMatrix6f.h>
#include <edm4hep/Track.h>
#include <edm4hep/TrackState.h>
#include <edm4hep/TrackerHitPlane.h>

#include <DD4hep/BitFieldCoder.h>

#include <k4Interface/IGeoSvc.h>

#include <Gaudi/Algorithm.h>
#include <GaudiKernel/SmartIF.h>

#include <string>
#include <vector>

void setStreamlogOutputLevel(const Gaudi::Algorithm* thisAlg, streamlog::logscope* streamlogOutputLevel);

// Simple class to wrap functions from MarlinTrkUtils
// and hold a few objects that are needed for the functions
class GaudiTrkUtils {
public:
  GaudiTrkUtils() = delete;
  GaudiTrkUtils(const Gaudi::Algorithm* thisAlg, const GaudiDDKalTest& ddKalTest, SmartIF<IGeoSvc> geoSvc,
                const std::string encodingStringVariable)
      : m_thisAlg(thisAlg), m_ddkaltest(ddKalTest), m_geoSvc(geoSvc), m_encodingStringVariable(encodingStringVariable) {
  }

  GaudiTrkUtils(const GaudiTrkUtils&) = delete;
  GaudiTrkUtils& operator=(const GaudiTrkUtils&) = delete;
  GaudiTrkUtils(GaudiTrkUtils&&) = delete;
  GaudiTrkUtils& operator=(GaudiTrkUtils&&) = delete;

  int createFinalisedLCIOTrack(GaudiDDKalTestTrack& marlinTrk, const std::vector<const edm4hep::TrackerHit*>& hit_list,
                               edm4hep::MutableTrack& track, bool fit_direction,
                               const edm4hep::CovMatrix6f& initial_cov_for_prefit, float bfield_z,
                               double maxChi2Increment);

  int createFinalisedLCIOTrack(GaudiDDKalTestTrack& marlinTrk, const std::vector<const edm4hep::TrackerHit*>& hit_list,
                               edm4hep::MutableTrack& track, bool fit_direction, edm4hep::TrackState& pre_fit,
                               double maxChi2Increment);

  int createPrefit(const std::vector<const edm4hep::TrackerHit*>& hit_list, edm4hep::TrackState& pre_fit,
                   float bfield_z);

  int createFit(const std::vector<const edm4hep::TrackerHit*>& hit_list, GaudiDDKalTestTrack& marlinTrk,
                edm4hep::TrackState& pre_fit, bool fit_direction, double maxChi2Increment);

  int finaliseLCIOTrack(GaudiDDKalTestTrack& marlintrk, edm4hep::MutableTrack& track,
                        const std::vector<const edm4hep::TrackerHit*>& hit_list, bool fit_direction);

  void addHitNumbersToTrack(std::vector<int32_t>& subdetectorHitNumbers,
                            const std::vector<const edm4hep::TrackerHit*>& hit_list, bool hits_in_fit,
                            const dd4hep::DDSegmentation::BitFieldCoder& cellID_encoder) const;

  int createTrackStateAtCaloFace(GaudiDDKalTestTrack& marlintrk, edm4hep::TrackState& trkStateCalo,
                                 const edm4hep::TrackerHit* trkhit, bool tanL_is_positive);

private:
  const Gaudi::Algorithm* m_thisAlg;
  const GaudiDDKalTest& m_ddkaltest;
  SmartIF<IGeoSvc> m_geoSvc;
  const std::string m_encodingStringVariable;
};

#endif // K4RECO_GAUDITRKUTILS_H
