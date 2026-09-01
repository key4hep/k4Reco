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

# Run FastJetAlg on the PandoraPFOs of an existing REC file. The Marlin
# counterpart of this configuration lives in runMarlinFastJet.py and the two
# outputs are compared bit by bit in compareJets.py

from Gaudi.Configuration import INFO
from k4FWCore import ApplicationMgr, IOSvc
from Configurables import EventDataSvc
from Configurables import FastJetAlg

from k4FWCore.parseArgs import parser

parser.add_argument("--algorithm", default="antikt_algorithm", help="Jet algorithm")
parser.add_argument(
    "--algorithm-params", nargs="*", default=["0.4"], help="Parameters of the jet algorithm"
)
parser.add_argument("--clustering-mode", default="Inclusive", help="Clustering mode")
parser.add_argument(
    "--clustering-params", nargs="+", default=["5.0"], help="Parameters of the clustering mode"
)

args = parser.parse_known_args()[0]

iosvc = IOSvc()
iosvc.Input = "output_ttbar_REC.edm4hep.root"
iosvc.Output = "output_fastjet.root"

fastJet = FastJetAlg(
    "FastJet",
    algorithm=args.algorithm,
    algorithmParameters=[float(p) for p in args.algorithm_params],
    clusteringMode=args.clustering_mode,
    clusteringParams=[float(p) for p in args.clustering_params],
    jetOut="JetOut",
    recParticleIn="PandoraPFOs",
    recParticleOut="UsedPFOs",
    recombinationScheme="E_scheme",
    OutputLevel=INFO,
)

ApplicationMgr(
    TopAlg=[fastJet],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[EventDataSvc("EventDataSvc")],
    OutputLevel=INFO,
)
