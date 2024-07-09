from Gaudi.Configuration import INFO

from k4FWCore import ApplicationMgr
from k4FWCore import IOSvc
from Configurables import EventDataSvc
from Configurables import OverlayTiming
from Configurables import UniqueIDGenSvc

id_service = UniqueIDGenSvc("UniqueIDGenSvc")

eds = EventDataSvc("EventDataSvc")

iosvc = IOSvc()
# iosvc.input = "input.root"
iosvc.input = "/home/juanmi/Key4hep/Algorithm-validation/Overlay/signal.root"
iosvc.output = "output_overlay.root"

# inp.collections = [
#     "EventHeader",
#     "MCParticle",
#     "VertexBarrelCollection",
#     "VertexEndcapCollection",
#     "HCalRingCollection",
# ]

overlay = OverlayTiming()
overlay.MCParticles = ["MCParticle"]
overlay.SimTrackerHits = ["VertexBarrelCollection", "VertexEndcapCollection"]
overlay.SimTrackerHitNames = ["NewVertexBarrelCollection", "NewVertexEndcapCollection"]
overlay.SimCalorimeterHits = ["HCalRingCollection"]
overlay.SimCalorimeterHitNames = ["NewHCalRingCollection"]
overlay.CaloHitContributions = ["CaloHitContributionsCollection"]
overlay.OutputSimTrackerHits = ["NewVertexBarrelCollection", "NewVertexEndcapCollection"]
overlay.OutputSimCalorimeterHits = ["NewHCalRingCollection"]
overlay.OutputCaloHitContributions = ["NewCaloHitCollection"]
# overlay.StartBackgroundEventIndex = 0
overlay.BackgroundFileNames = [
    ["/home/juanmi/Key4hep/Algorithm-validation/Overlay/background1.root"],
    ["/home/juanmi/Key4hep/Algorithm-validation/Overlay/background2.root"],
]
overlay.TimeWindows = {"MCParticle": [0, 23.5], "VertexBarrelCollection": [0, 23.5], "VertexEndcapCollection": [0, 23.5], "HCalRingCollection": [0, 23.5]}

ApplicationMgr(TopAlg=[overlay],
               EvtSel="NONE",
               EvtMax=10,
               ExtSvc=[eds],
               OutputLevel=INFO,
               )
