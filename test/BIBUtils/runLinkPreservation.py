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

# Runs the BIBUtils filters that need no geometry on the synthetic input
# written by makeLinkTestInput.py.
from Gaudi.Configuration import INFO
from k4FWCore import ApplicationMgr, IOSvc
from Configurables import EventDataSvc
from Configurables import CaloHitConeFilter, SplitCollectionByPolarAngle

coneFilter = CaloHitConeFilter("CaloHitConeFilter")
coneFilter.MCParticleCollectionName = "MCParticle"
coneFilter.CaloHitCollectionName = "CaloHits"
coneFilter.CaloRelationCollectionName = "CaloHitLinks"
coneFilter.GoodHitCollection = "CaloHitsConed"
coneFilter.GoodRelationCollection = "CaloHitLinksConed"
coneFilter.ConeWidth = 0.2

splitter = SplitCollectionByPolarAngle("SplitCollectionByPolarAngle")
splitter.TrackerHitInputCollections = "TrackerHits"
splitter.TrackerHitInputRelations = "TrackerHitLinks"
splitter.TrackerHitOutputCollections = "TrackerHitsSplit"
splitter.TrackerSimHitOutputCollections = "SimTrackerHitsSplit"
splitter.TrackerHitOutputRelations = "TrackerHitLinksSplit"
splitter.PolarAngleLowerLimit = 50.0
splitter.PolarAngleUpperLimit = 130.0

iosvc = IOSvc()
iosvc.Input = "bibutils_links_input.edm4hep.root"
iosvc.Output = "bibutils_links_output.edm4hep.root"

ApplicationMgr(
    TopAlg=[coneFilter, splitter],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[EventDataSvc("EventDataSvc")],
    OutputLevel=INFO,
)
