#
# Copyright (c) 2020-2024 Key4hep-Project.
#
# This file is part of Key4hep.
# See https://key4hep.github.io/key4hep-doc/ for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

# This steering file runs conformal tracking starting from a file with digitised collections
# It will not work with a file that has been obtained directly with ddsim

import os

from Gaudi.Configuration import INFO
from k4FWCore import ApplicationMgr
from k4FWCore import IOSvc
from Configurables import EventDataSvc
from Configurables import ClonesAndSplitTracksFinder
from Configurables import UniqueIDGenSvc
from Configurables import GeoSvc

id_service = UniqueIDGenSvc("UniqueIDGenSvc")

eds = EventDataSvc("EventDataSvc")

geoservice = GeoSvc("GeoSvc")
geoservice.detectors = [
    os.environ["K4GEO"] + "/FCCee/CLD/compact/CLD_o2_v07/CLD_o2_v07.xml"
]
geoservice.OutputLevel = INFO
geoservice.EnableGeant4Geo = False

iosvc = IOSvc()
iosvc.Input = "output_REC.edm4hep.root"
iosvc.Output = "output_clones_and_split_tracks_finder.root"

ClonesAndSplitTracksFinder = ClonesAndSplitTracksFinder(
    "ClonesAndSplitTracksFinder",
    InputTrackCollectionName=["SiTracksCT"],
    OutputTrackCollectionName=["GaudiSiTracks"],
    EnergyLossOn=True,
    MultipleScatteringOn=True,
    SmoothOn=False,
    extrapolateForward=True,
    maxSignificancePhi=3,
    maxSignificancePt=2,
    maxSignificanceTheta=3,
    mergeSplitTracks=False,
    minTrackPt=1,
    OutputLevel=INFO,
)

ApplicationMgr(
    TopAlg=[ClonesAndSplitTracksFinder],
    EvtSel="NONE",
    EvtMax=1,
    ExtSvc=[eds, geoservice],
    OutputLevel=INFO,
)
