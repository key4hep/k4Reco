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
overlay.SimCalorimeterHits = ["HCalRingCollection"]
overlay.CaloHitContributions = ["CaloHitContributionsCollection"]
overlay.OutputSimTrackerHits = ["NewVertexBarrelCollection", "NewVertexEndcapCollection"]
overlay.OutputSimCalorimeterHits = ["NewHCalRingCollection"]
overlay.OutputCaloHitContributions = ["NewCaloHitCollection"]
# overlay.StartBackgroundEventIndex = 0
overlay.BackgroundFileNames = [
    ["/home/juanmi/Key4hep/Algorithm-validation/Overlay/background1.root"],
    ["/home/juanmi/Key4hep/Algorithm-validation/Overlay/background2.root"],
]
# overlay.MCParticles = "MCParticles"
# overlay.filterTimeMin = -0.25
# overlay.filterTimeMax = 23.25
# Supported formats:
#   <collection name>: []  << all objects
#   <collection name>: [t_max]  << all objects with time < t_max
#   <collection name>: [t_min, t_max]  << all objects with time between t_min and t_max

# overlay.inputCollections = {
#     "MCParticles": [],
#     "VertexBarrelCollection": [],
#     "HCalRingCollection": [],

# }

ApplicationMgr(TopAlg=[overlay],
               EvtSel="NONE",
               EvtMax=1,
               ExtSvc=[eds],
               OutputLevel=INFO,
               )
