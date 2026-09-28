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
from Configurables import CaloHitSelector
from Configurables import GeoSvc
import os

# GeoSvc is needed to decode the calorimeter layer from the cellID.
geoservice = GeoSvc("GeoSvc")
geoservice.detectors = [
    os.environ["MUONCOLLIDER_GEO"]
    if "MUONCOLLIDER_GEO" in os.environ
    else os.environ["K4GEO"] + "/MuColl/MuColl_v1/MuColl_v1.xml"
]
geoservice.OutputLevel = INFO
geoservice.EnableGeant4Geo = False

selector = CaloHitSelector("MyEcalBarrelSelector")
selector.CaloHitCollectionName = ["EcalBarrelCollectionConed"]
selector.CaloRelationCollectionName = ["EcalBarrelRelationsSimConed"]
selector.GoodHitCollection = ["EcalBarrelCollectionSel"]
selector.GoodRelationCollection = ["EcalBarrelRelationsSimSel"]
# Path to the ROOT file holding the per-(theta, layer) threshold maps
# (th_2dmode_sym and stddev_sym). Alternatively set FlatThreshold > 0.
selector.ThresholdsFilePath = ""
selector.Nsigma = 0
selector.FlatThreshold = 0.0
selector.TimeWindowMin = -0.3
selector.TimeWindowMax = 0.3
selector.DoBIBsubtraction = False

iosvc = IOSvc()
iosvc.Input = "coned_output.edm4hep.root"
iosvc.Output = "selected_output.edm4hep.root"

ApplicationMgr(
    TopAlg=[selector],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[EventDataSvc("EventDataSvc")],
    OutputLevel=INFO,
)
