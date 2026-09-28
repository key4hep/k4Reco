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
from Gaudi.Configuration import INFO
from k4FWCore import ApplicationMgr, IOSvc
from Configurables import EventDataSvc
from Configurables import SplitCollectionByLayer
from Configurables import GeoSvc
import os

geoservice = GeoSvc("GeoSvc")
geoservice.detectors = [
    os.environ["MUONCOLLIDER_GEO"]
    if "MUONCOLLIDER_GEO" in os.environ
    else os.environ["K4GEO"] + "/MuColl/MuColl_v1/MuColl_v1.xml"
]
geoservice.OutputLevel = INFO
geoservice.EnableGeant4Geo = False

# The number of output collections is arbitrary and follows the length of
# OutputCollections; StartLayers and EndLayers must have the same length. Each
# output collection i collects the hits whose layer is in [StartLayers[i],
# EndLayers[i]].
splitter = SplitCollectionByLayer("VBTrackerHitsByLayer")
splitter.InputCollection = ["VBTrackerHits"]
splitter.OutputCollections = ["VBTrackerHitsInner", "VBTrackerHitsOuter"]
splitter.StartLayers = [0, 4]
splitter.EndLayers = [3, 7]

iosvc = IOSvc()
iosvc.Input = "digi_output.edm4hep.root"
iosvc.Output = "split_by_layer_output.edm4hep.root"

ApplicationMgr(
    TopAlg=[splitter],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[EventDataSvc("EventDataSvc"), geoservice],
    OutputLevel=INFO,
)
