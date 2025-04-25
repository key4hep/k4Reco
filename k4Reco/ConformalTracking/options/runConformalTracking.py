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
from Configurables import ConformalTracking
from Configurables import UniqueIDGenSvc
from Configurables import GeoSvc
from Configurables import RootHistSvc
from Configurables import Gaudi__Histograming__Sink__Root as RootHistoSink

from conformal_tracking_utils import configure_conformal_tracking_steps

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
iosvc.Output = "output_conformal_tracking.root"


tracking = ConformalTracking()
tracking.TrackerHitCollectionNames = [
    "VXDTrackerHits",
    "VXDEndcapTrackerHits",
    "ITrackerHits",
    "OTrackerHits",
    "ITrackerEndcapHits",
    "OTrackerEndcapHits",
]
tracking.RelationsNames = [
    "VXDTrackerHitRelations",
    "VXDEndcapTrackerHitRelations",
    "InnerTrackerBarrelHitsRelations",
    "InnerTrackerEndcapHitsRelations",
    "OuterTrackerBarrelHitsRelations",
    "OuterTrackerEndcapHitsRelations",
]
tracking.MCParticleCollectionName = ["MCParticles"]
tracking.SiTrackCollectionName = ["NewSiTracks"]

tracking.MainTrackerHitCollectionNames = [
    "ITrackerHits",
    "OTrackerHits",
    "ITrackerEndcapHits",
    "OTrackerEndcapHits",
]
tracking.VertexBarrelHitCollectionNames = ["VXDTrackerHits"]
tracking.VertexEndcapHitCollectionNames = ["VXDEndcapTrackerHits"]

# tracking.DebugHits = "DebugHits"
tracking.DebugPlots = True
tracking.DebugTiming = False
# tracking.MCParticleCollectionName = ["MCParticle"]
tracking.MaxHitInvertedFit = 0
tracking.MinClustersOnTrackAfterFit = 3
tracking.RetryTooManyTracks = False
# # tracking.SiTrackCollectionName = ["SiTracksCT"],
tracking.SortTreeResults = True
tracking.ThetaRange = 0.05
tracking.TooManyTracks = 100000
tracking.trackPurity = 0.7

CT_MAX_DIST = 0.05

# The keys (VXDBarrel, VXDEndcap...) are simply names and are not passed to ConformalTracking
parameters = {
        "VXDBarrel": {
            "collections": ["VXDTrackerHits"],
            "params": {
                "MaxCellAngle": 0.01,
                "MaxCellAngleRZ": 0.01,
                "Chi2Cut": 100,
                "MinClustersOnTrack": 4,
                "MaxDistance": CT_MAX_DIST,
                "SlopeZRange": 10.0,
                "HighPTCut": 10.0,
            },
            "flags": ["HighPTFit", "VertexToTracker"],
            "functions": ["CombineCollections", "BuildNewTracks"],
        },
        "VXDEncap": {
            "collections": ["VXDEndcapTrackerHits"],
            "params": {
                "MaxCellAngle": 0.01,
                "MaxCellAngleRZ": 0.01,
                "Chi2Cut": 100,
                "MinClustersOnTrack": 4,
                "MaxDistance": CT_MAX_DIST,
                "SlopeZRange": 10.0,
                "HighPTCut": 10.0,
            },
            "flags": ["HighPTFit", "VertexToTracker"],
            "functions": ["CombineCollections", "ExtendTracks"],
        },
        "LowerCellAngle1": {
            "collections": ["VXDTrackerHits", "VXDEndcapTrackerHits"],
            "params": {
                "MaxCellAngle": 0.05,
                "MaxCellAngleRZ": 0.05,
                "Chi2Cut": 100,
                "MinClustersOnTrack": 4,
                "MaxDistance": CT_MAX_DIST,
                "SlopeZRange": 10.0,
                "HighPTCut": 10.0,
            },
            "flags": ["HighPTFit", "VertexToTracker", "RadialSearch"],
            "functions": ["CombineCollections", "BuildNewTracks"],
        },
        "LowerCellAngle2": {
            "collections": [],
            "params": {
                "MaxCellAngle": 0.1,
                "MaxCellAngleRZ": 0.1,
                "Chi2Cut": 2000,
                "MinClustersOnTrack": 4,
                "MaxDistance": CT_MAX_DIST,
                "SlopeZRange": 10.0,
                "HighPTCut": 10.0,
            },
            "flags": ["HighPTFit", "VertexToTracker", "RadialSearch"],
            "functions": ["BuildNewTracks", "SortTracks"],
        },
        "Tracker": {
            "collections": ["ITrackerHits", "OTrackerHits", "ITrackerEndcapHits", "OTrackerEndcapHits"],
            "params": {
                "MaxCellAngle": 0.1,
                "MaxCellAngleRZ": 0.1,
                "Chi2Cut": 2000,
                "MinClustersOnTrack": 4,
                "MaxDistance": CT_MAX_DIST,
                "SlopeZRange": 10.0,
                "HighPTCut": 1.0,
            },
            "flags": ["HighPTFit", "VertexToTracker", "RadialSearch"],
            "functions": ["CombineCollections", "ExtendTracks"],
        },
        "Displaced": {
            "collections": ["VXDTrackerHits", "VXDEndcapTrackerHits", "ITrackerHits", "OTrackerHits", "ITrackerEndcapHits", "OTrackerEndcapHits"],
            "params": {
                "MaxCellAngle": 0.1,
                "MaxCellAngleRZ": 0.1,
                "Chi2Cut": 1000,
                "MinClustersOnTrack": 5,
                "MaxDistance": 0.015,
                "SlopeZRange": 10.0,
                "HighPTCut": 10.0,
            },
            "flags": ["OnlyZSchi2cut", "RadialSearch"],
            "functions": ["CombineCollections", "BuildNewTracks"],
        },
    }

configure_conformal_tracking_steps(tracking, parameters)

hps = RootHistSvc("HistogramPersistencySvc")
root_hist_svc = RootHistoSink("RootHistoSink")
root_hist_svc.FileName = "conformal_tracking_hist.root"

ApplicationMgr(
    TopAlg=[tracking],
    EvtSel="NONE",
    EvtMax=3,
    ExtSvc=[eds, geoservice, root_hist_svc],
    OutputLevel=INFO,
)
