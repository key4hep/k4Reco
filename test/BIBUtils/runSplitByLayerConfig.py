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

# A valid SplitCollectionByLayer configuration, which the tests make invalid from
# the command line (e.g. --SplitCollectionByLayerConfig.StartLayers 0 1) to check
# that initialize() rejects it with the expected message.
import os

from Gaudi.Configuration import INFO
from k4FWCore import ApplicationMgr, IOSvc
from Configurables import EventDataSvc, GeoSvc
from Configurables import SplitCollectionByLayer

geoservice = GeoSvc("GeoSvc")
geoservice.detectors = [os.environ["K4GEO"] + "/MuColl/MAIA/compact/MAIA_v0/MAIA_v0.xml"]
geoservice.OutputLevel = INFO
geoservice.EnableGeant4Geo = False

splitter = SplitCollectionByLayer("SplitCollectionByLayerConfig")
splitter.InputCollection = "LayerTrackerHits"
splitter.OutputCollections = ["LayerHitsConfig"]
splitter.StartLayers = [0]
splitter.EndLayers = [3]

iosvc = IOSvc()
iosvc.Input = "bibutils_geo_input.edm4hep.root"
iosvc.Output = "bibutils_split_by_layer_config.edm4hep.root"

ApplicationMgr(
    TopAlg=[splitter],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[EventDataSvc("EventDataSvc")],
    OutputLevel=INFO,
)
