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

import os
print(os.getcwd())

from Gaudi.Configuration import INFO, WARNING
from k4FWCore import ApplicationMgr, IOSvc
from Configurables import EventDataSvc
from Configurables import DDPandoraPFANewAlgorithm
from Configurables import GeoSvc
from Configurables import UniqueIDGenSvc
from Configurables import MarlinProcessorWrapper
from FastJet.pandoraSettings import pandora

iosvc = IOSvc()
iosvc.Input = "output_REC.edm4hep.root"
iosvc.Output = "output_pandora.root"

id_service = UniqueIDGenSvc("UniqueIDGenSvc")

geoservice = GeoSvc("GeoSvc")
geoservice.detectors = [
    os.environ["K4GEO"] + "/FCCee/CLD/compact/CLD_o2_v07/CLD_o2_v07.xml"
]
geoservice.OutputLevel = INFO
geoservice.EnableGeant4Geo = False

marlinFastJet = MarlinProcessorWrapper("MyFastJetProcessor")
marlinFastJet.OutputLevel = WARNING
marlinFastJet.ProcessorType = "FastJetProcessor"
marlinFastJet.Parameters = {
    "algorithm": ["antikt_algorithm", "0.4"],
    "clusteringMode": ["Inclusive", "5.0"],
    "jetOut": ["JetOut"],
    "recParticleIn": ["GaudiPandoraPFOs"],
    "recParticleOut": ["PFOsFromJets"],
    "recombinationScheme": ["E_scheme"],
    "storeParticlesInJets": ["true"]
}

ApplicationMgr(
    TopAlg = [pandora, marlinFastJet],
    EvtSel = "NONE",
    EvtMax = 1,
    ExtSvc = [EventDataSvc("EventDataSvc")],
    OutputLevel=INFO,
)
