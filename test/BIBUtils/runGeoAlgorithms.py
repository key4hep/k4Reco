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

# Runs the BIBUtils algorithms that need the detector geometry on the input
# written by makeGeoTestInput.py. The expected outputs are listed in
# checkGeoAlgorithms.py.
import os

from Gaudi.Configuration import INFO
from k4FWCore import ApplicationMgr, IOSvc
from Configurables import EventDataSvc, GeoSvc
from Configurables import CaloHitSelector, SplitCollectionByLayer, TrackerHitHelixFilter

geoservice = GeoSvc("GeoSvc")
geoservice.detectors = [os.environ["K4GEO"] + "/MuColl/MAIA/compact/MAIA_v0/MAIA_v0.xml"]
geoservice.OutputLevel = INFO
geoservice.EnableGeant4Geo = False


def caloHitSelector(label, **properties):
    selector = CaloHitSelector(f"CaloHitSelector{label}")
    selector.CaloHitCollectionName = "CaloHits"
    selector.CaloRelationCollectionName = "CaloHitLinks"
    selector.GoodHitCollection = f"CaloHits{label}"
    selector.GoodRelationCollection = f"CaloHitLinks{label}"
    selector.ThresholdsFilePath = "bibutils_thresholds.root"
    selector.TimeWindowMin = -0.5
    selector.TimeWindowMax = 1.0
    for name, value in properties.items():
        setattr(selector, name, value)
    return selector


def trackerHitHelixFilter(label, **properties):
    helixFilter = TrackerHitHelixFilter(f"TrackerHitHelixFilter{label}")
    helixFilter.MCParticleCollection = "MCParticle"
    helixFilter.TrackerHitInputCollections = "TrackerHits"
    helixFilter.TrackerHitInputRelations = "TrackerHitLinks"
    helixFilter.TrackerHitOutputCollections = f"TrackerHits{label}"
    helixFilter.TrackerSimHitOutputCollections = f"SimTrackerHits{label}"
    helixFilter.TrackerHitOutputRelations = f"TrackerHitLinks{label}"
    helixFilter.ConeAroundStatus = [1]
    helixFilter.TrackerOuterRadius = 1500.0
    for name, value in properties.items():
        setattr(helixFilter, name, value)
    return helixFilter


algorithms = [
    caloHitSelector("Map", Nsigma=2),
    caloHitSelector("Flat", Nsigma=2, FlatThreshold=0.05),
    caloHitSelector("BIBSub", Nsigma=0, DoBIBsubtraction=True),
    trackerHitHelixFilter("Dist3D", Dist3DCut=30.0, DeltaRCut=-1.0),
    trackerHitHelixFilter("DeltaR", Dist3DCut=-1.0, DeltaRCut=0.05),
    SplitCollectionByLayer(
        "SplitCollectionByLayer",
        InputCollection="LayerTrackerHits",
        OutputCollections=["LayerHits0to3", "LayerHits2to5", "LayerHits7"],
        StartLayers=[0, 2, 7],
        EndLayers=[3, 5, 7],
    ),
]

iosvc = IOSvc()
iosvc.Input = "bibutils_geo_input.edm4hep.root"
iosvc.Output = "bibutils_geo_output.edm4hep.root"

ApplicationMgr(
    TopAlg=algorithms,
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[EventDataSvc("EventDataSvc")],
    OutputLevel=INFO,
)
