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
from Configurables import FilterConeHits
from Configurables import GeoSvc
from Configurables import RootHistSvc
from Configurables import Gaudi__Histograming__Sink__Root as RootHistoSink
import os

geoservice = GeoSvc("GeoSvc")
geoservice.detectors = [
    os.environ["MUONCOLLIDER_GEO"]
    if "MUONCOLLIDER_GEO" in os.environ
    else os.environ["K4GEO"] + "/MuColl/MuColl_v1/MuColl_v1.xml"
]
geoservice.OutputLevel = INFO
geoservice.EnableGeant4Geo = False

# One instance per tracker subdetector. The simulated hits are taken from the
# input reco-to-sim relations, so only the reco hit and relation collections
# are needed as inputs.
coner = FilterConeHits("VXDBarrelConer")
coner.MCParticleCollection = ["MCParticle"]
coner.TrackerHitInputCollections = ["VBTrackerHits"]
coner.TrackerHitInputRelations = ["VBTrackerHitsRelations"]
coner.TrackerHitOutputCollections = ["VBTrackerHitsConed"]
coner.TrackerSimHitOutputCollections = ["VertexBarrelCollectionConed"]
coner.TrackerHitOutputRelations = ["VBTrackerHitsRelationsConed"]
coner.Dist3DCut = 30.0
coner.DeltaRCut = -1.0
coner.FillHistograms = False

iosvc = IOSvc()
iosvc.Input = "digi_output.edm4hep.root"
iosvc.Output = "coned_output.edm4hep.root"

hps = RootHistSvc("HistogramPersistencySvc")
root_hist_svc = RootHistoSink("RootHistoSink")
root_hist_svc.FileName = "filterconehits_hist.root"

ApplicationMgr(
    TopAlg=[coner],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[EventDataSvc("EventDataSvc"), root_hist_svc],
    OutputLevel=INFO,
)
