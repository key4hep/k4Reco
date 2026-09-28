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
from Configurables import SplitCollectionByPolarAngle
from Configurables import RootHistSvc
from Configurables import Gaudi__Histograming__Sink__Root as RootHistoSink

# One instance per tracker collection. The simulated hits are taken from the
# input reco-to-sim relations, so only the reco hit and relation collections
# are needed as inputs.
splitter = SplitCollectionByPolarAngle("VXDBarrelSplitter")
splitter.TrackerHitInputCollections = ["VBTrackerHits"]
splitter.TrackerHitInputRelations = ["VBTrackerHitsRelations"]
splitter.TrackerHitOutputCollections = ["VBTrackerHitsSplit"]
splitter.TrackerSimHitOutputCollections = ["VertexBarrelCollectionSplit"]
splitter.TrackerHitOutputRelations = ["VBTrackerHitsRelationsSplit"]
splitter.PolarAngleLowerLimit = 50.0
splitter.PolarAngleUpperLimit = 130.0
splitter.FillHistograms = False

iosvc = IOSvc()
iosvc.Input = "digi_output.edm4hep.root"
iosvc.Output = "split_output.edm4hep.root"

hps = RootHistSvc("HistogramPersistencySvc")
root_hist_svc = RootHistoSink("RootHistoSink")
root_hist_svc.FileName = "splitcollectionbypolarangle_hist.root"

ApplicationMgr(
    TopAlg=[splitter],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[EventDataSvc("EventDataSvc"), root_hist_svc],
    OutputLevel=INFO,
)
