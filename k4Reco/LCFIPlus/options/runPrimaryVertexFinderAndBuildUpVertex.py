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
from Configurables import PrimaryVertexFinder, BuildUpVertex, JetClustering, JetVertexRefiner
from Configurables import GeoSvc
from Configurables import UniqueIDGenSvc
import os

eds = EventDataSvc("EventDataSvc")
id_service = UniqueIDGenSvc("UniqueIDGenSvc")

# id_service = UniqueIDGenSvc("UniqueIDGenSvc")

# geoservice = GeoSvc("GeoSvc")
# geoservice.detectors = [os.environ["K4GEO"]+"/FCCee/CLD/compact/CLD_o2_v07/CLD_o2_v07.xml"]
# geoservice.OutputLevel = INFO
# geoservice.EnableGeant4Geo = False

iosvc = IOSvc()
iosvc.Input = "output_REC.edm4hep.root"
iosvc.Output = "output_primary_and_build_vertex.root"

PrimaryVertexFinder = PrimaryVertexFinder("PrimaryVertexFinder")

PrimaryVertexFinder.PFOCollection = ["PFOsFromJets"]
PrimaryVertexFinder.VertexNames = ["GAUDIPrimaryVertices"]

PrimaryVertexFinder.BeamspotConstraint = 1
PrimaryVertexFinder.BeamspotSmearing = False
PrimaryVertexFinder.MagneticField = 2.0
# The Marlin reference is CLDReconstruction.py's default --cms 240 setup.
PrimaryVertexFinder.BeamSizeX = 9.8e-3
PrimaryVertexFinder.BeamSizeY = 25.4e-6
PrimaryVertexFinder.BeamSizeZ = 0.65
PrimaryVertexFinder.Chi2Threshold = 25.0
PrimaryVertexFinder.TrackMaxD0 = 20.0
# PrimaryVertexFinder.TrackMaxInnermostHitRadius = 61.0
PrimaryVertexFinder.TrackMaxZ0 = 20.0
PrimaryVertexFinder.TrackMinFtdHits = 999999
PrimaryVertexFinder.TrackMinTpcHits = 999999
PrimaryVertexFinder.TrackMinTpcHitsMinPt = 999999
# PrimaryVertexFinder.TrackMinVtxFtdHits = 1
PrimaryVertexFinder.TrackMinVxdHits = 999999

BuildUpVertex = BuildUpVertex("BuildUpVertex")

BuildUpVertex.PFOCollection = ["PFOsFromJets"]
BuildUpVertex.PrimaryVertexCollectionName = ["GAUDIPrimaryVertices"]
BuildUpVertex.BuildUpVertexCollectionName = ["GAUDIBuildUpVertices"]
BuildUpVertex.V0VertexCollectionName = ["GAUDIBuildUpVertices_V0"]

BuildUpVertex.AVFTemperature = 5.0
BuildUpVertex.AssocIPTracks = 1
BuildUpVertex.AssocIPTracksChi2RatioSecToPri = 2.0
BuildUpVertex.AssocIPTracksMinDist = 0.0
BuildUpVertex.MassThreshold = 10.0
BuildUpVertex.MaxChi2ForDistOrder = 1.0
BuildUpVertex.MinDistFromIP = 0.3
BuildUpVertex.PrimaryChi2Threshold = 25.0
BuildUpVertex.SecondaryChi2Threshold = 9.0
BuildUpVertex.TrackMaxD0 = 10.0
BuildUpVertex.TrackMaxD0Err = 0.1
BuildUpVertex.TrackMaxZ0 = 20.0
BuildUpVertex.TrackMaxZ0Err = 0.1
BuildUpVertex.TrackMinFtdHits = 1
BuildUpVertex.TrackMinPt = 0.1
BuildUpVertex.TrackMinTpcHits = 1
BuildUpVertex.TrackMinTpcHitsMinPt = 999999
BuildUpVertex.TrackMinVxdFtdHits = 1
BuildUpVertex.TrackMinVxdHits = 1
BuildUpVertex.UseAVF = False
BuildUpVertex.UseV0Selection = 1
BuildUpVertex.BeamspotConstraint = True
BuildUpVertex.BeamspotSmearing = False

JetClustering = JetClustering("JetClustering")
JetClustering.PFOCollection = ["PFOsFromJets"]
JetClustering.PrimaryVertexCollectionName = ["GAUDIPrimaryVertices"]
JetClustering.InputVertexCollectionName = ["GAUDIBuildUpVertices"]
JetClustering.OutputJetCollectionName = ["GAUDIVertexJets"]
JetClustering.NJetsRequested = 2
JetClustering.YCut = 0.0
JetClustering.JetAlgorithm = "ValenciaVertex"
JetClustering.UseBeamJets = True
JetClustering.RParameter = 1.0
JetClustering.AlphaParameter = 1.0
JetClustering.BetaParameter = 1.0
JetClustering.GammaParameter = 1.0
JetClustering.YAddedForJetVertexVertex = 100.0
JetClustering.YAddedForJetVertexLepton = 0.0
JetClustering.YAddedForJetLeptonLepton = 100.0
JetClustering.MagneticField = 2.0

JetVertexRefiner = JetVertexRefiner("JetVertexRefiner")
JetVertexRefiner.PFOCollection = ["PFOsFromJets"]
JetVertexRefiner.InputJetCollectionName = ["GAUDIVertexJets"]
JetVertexRefiner.PrimaryVertexCollectionName = ["GAUDIPrimaryVertices"]
JetVertexRefiner.InputVertexCollectionName = ["GAUDIBuildUpVertices"]
JetVertexRefiner.V0VertexCollectionName = ["GAUDIBuildUpVertices_V0"]
JetVertexRefiner.OutputJetCollectionName = ["GAUDIRefinedVertexJets"]
JetVertexRefiner.OutputVertexCollectionName = ["GAUDIRefinedVertices"]


# VertexFinder.Algorithms = ["PrimaryVertexFinder", "BuildUpVertex"]
# VertexFinder.BeamSizeX = [str(BEAM_SPOT_SIZES[reco_args.cms][0])]
# VertexFinder.BeamSizeY = [str(BEAM_SPOT_SIZES[reco_args.cms][1])]
# VertexFinder.BeamSizeZ = [str(BEAM_SPOT_SIZES[reco_args.cms][2])]
# VertexFinder.MCPCollection = ["MCParticle"]
# VertexFinder.MCPFORelation = ["RecoMCTruthLink"]
# VertexFinder.MagneticField = ["2.0"]
# VertexFinder.PFOCollection = ["PFOsFromJets"]
# VertexFinder.PrintEventNumber = ["1"]
# VertexFinder.ReadSubdetectorEnergies = ["0"]
# VertexFinder.TrackHitOrdering = ["2"]
# VertexFinder.UpdateVertexRPDaughters = ["0"]
# VertexFinder.UseMCP = ["0"]

# VertexFinder.AlgorithmParameters = {
#                            # "BuildUpVertex.AVFTemperature": ["5.0"],
#                            # "BuildUpVertex.AssocIPTracks": ["1"],
#                            # "BuildUpVertex.AssocIPTracksChi2RatioSecToPri": ["2.0"],
#                            # "BuildUpVertex.AssocIPTracksMinDist": ["0."],
#                            # "BuildUpVertex.MassThreshold": ["10."],
#                            # "BuildUpVertex.MaxChi2ForDistOrder": ["1.0"],
#                            # "BuildUpVertex.MinDistFromIP": ["0.3"],
#                            # "BuildUpVertex.PrimaryChi2Threshold": ["25."],
#                            # "BuildUpVertex.SecondaryChi2Threshold": ["9."],
#                            # "BuildUpVertex.TrackMaxD0": ["10."],
#                            # "BuildUpVertex.TrackMaxD0Err": ["0.1"],
#                            # "BuildUpVertex.TrackMaxZ0": ["20."],
#                            # "BuildUpVertex.TrackMaxZ0Err": ["0.1"],
#                            # "BuildUpVertex.TrackMinFtdHits": ["1"],
#                            # "BuildUpVertex.TrackMinPt": ["0.1"],
#                            # "BuildUpVertex.TrackMinTpcHits": ["1"],
#                            # "BuildUpVertex.TrackMinTpcHitsMinPt": ["999999"],
#                            # "BuildUpVertex.TrackMinVxdFtdHits": ["1"],
#                            # "BuildUpVertex.TrackMinVxdHits": ["1"],
#                            # "BuildUpVertex.UseAVF": ["false"],
#                            # "BuildUpVertex.UseV0Selection": ["1"],
#                            # "BuildUpVertex.V0VertexCollectionName": ["BuildUpVertices_V0"],
#                            # "BuildUpVertexCollectionName": ["BuildUpVertices"],
# }

ApplicationMgr(TopAlg=[PrimaryVertexFinder,
                       BuildUpVertex,
                       JetClustering,
                       JetVertexRefiner
                       ],
               EvtSel="NONE",
               EvtMax=3,
               ExtSvc=[eds, id_service],
               OutputLevel=INFO,
               )
