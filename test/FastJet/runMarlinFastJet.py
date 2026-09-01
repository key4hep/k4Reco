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

# Marlin reference for runFastJet.py: the same FastJet configuration, run
# through the FastJetProcessor on the same PandoraPFOs input. Note that the
# CLD reconstruction does not produce jets by default (with Overlay disabled
# it only renames PandoraPFOs to PFOsFromJets), so the reference has to be
# produced here rather than taken from the CLDReconstruction output

from Gaudi.Configuration import INFO, WARNING
from k4FWCore import ApplicationMgr, IOSvc
from Configurables import EventDataSvc
from Configurables import MarlinProcessorWrapper

from k4MarlinWrapper.io_helpers import IOHandlerHelper
from k4FWCore.parseArgs import parser

parser.add_argument(
    "--inputfile", help="Input (REC) file", default="output_ttbar_REC.edm4hep.root"
)
parser.add_argument(
    "--outputfile", help="Output file after running FastJet", default="output_marlin_fastjet.root"
)
parser.add_argument("--algorithm", default="antikt_algorithm", help="Jet algorithm")
parser.add_argument(
    "--algorithm-params", nargs="*", default=["0.4"], help="Parameters of the jet algorithm"
)
parser.add_argument("--clustering-mode", default="Inclusive", help="Clustering mode")
parser.add_argument(
    "--clustering-params", nargs="+", default=["5.0"], help="Parameters of the clustering mode"
)

args = parser.parse_known_args()[0]

algList = []
iosvc = IOSvc()
io_handler = IOHandlerHelper(algList, iosvc)
io_handler.add_reader([args.inputfile])

marlinFastJet = MarlinProcessorWrapper("MyFastJetProcessor")
marlinFastJet.OutputLevel = WARNING
marlinFastJet.ProcessorType = "FastJetProcessor"
marlinFastJet.Parameters = {
    "algorithm": [args.algorithm] + args.algorithm_params,
    "clusteringMode": [args.clustering_mode] + args.clustering_params,
    "jetOut": ["JetOut"],
    "recParticleIn": ["PandoraPFOs"],
    "recParticleOut": ["UsedPFOs"],
    "recombinationScheme": ["E_scheme"],
    "storeParticlesInJets": ["true"],
}

algList.append(marlinFastJet)

io_handler.add_edm4hep_writer(args.outputfile)
io_handler.finalize_converters()

ApplicationMgr(
    TopAlg=algList,
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[EventDataSvc("EventDataSvc")],
    OutputLevel=INFO,
)
