#
# Copyright (c) 2020-2023 Key4hep-Project.
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
from Configurables import ApplicationMgr, k4DataSvc, PodioInput, PodioOutput
from Configurables import DDPlanarDigi
from Configurables import GeoSvc
from Configurables import UniqueIDGenSvc
import os

id_service = UniqueIDGenSvc("UniqueIDGenSvc")

geoservice = GeoSvc("GeoSvc")
geoservice.detectors = [os.environ["K4GEO"]+"/CLIC/compact/CLIC_o3_v15/CLIC_o3_v15.xml"]
geoservice.OutputLevel = INFO
geoservice.EnableGeant4Geo = False

processor = DDPlanarDigi()
processor.SubDetectorName = "Vertex"
processor.IsStrip = False
processor.ResolutionU = [0.003, 0.003, 0.003, 0.003, 0.003, 0.003]
processor.ResolutionV = [0.003, 0.003, 0.003, 0.003, 0.003, 0.003]
processor.SimTrackerHitCollectionName = "VertexBarrelCollection"
processor.SimTrkHitRelCollection = "VXDTrackerHitRelations"
processor.TrackerHitCollectionName = "VXDTrackerHits"

data_svc = k4DataSvc("EventDataSvc")
data_svc.input = "input.root"

inp = PodioInput()
inp.collections = [
    "VertexBarrelCollection",
    "EventHeader",
]

out = PodioOutput("out")
out.filename = "planar_digi_histograms.root"

ApplicationMgr(TopAlg=[inp, processor, out],
               EvtSel="NONE",
               EvtMax=-1,
               ExtSvc=[data_svc],
               OutputLevel=INFO,
               )
