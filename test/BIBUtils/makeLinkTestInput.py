#!/usr/bin/env python
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

# Writes a one-event file with reco hits carrying several weighted links to sim
# hits, to check that the BIBUtils filters keep all the links of the hits they
# accept. Hits are identified by their cellID in checkLinkPreservation.py.
import argparse

import ROOT
import edm4hep
from podio import Frame, root_io

ROOT.gInterpreter.Declare(
    "#include <edm4hep/CaloHitSimCaloHitLinkCollection.h>\n"
    "#include <edm4hep/TrackerHitSimTrackerHitLinkCollection.h>\n"
)

parser = argparse.ArgumentParser(description="Write the input for the BIBUtils link tests")
parser.add_argument("--output", default="bibutils_links_input.edm4hep.root", help="Output file")
args = parser.parse_args()

mcParticles = edm4hep.MCParticleCollection()
particle = mcParticles.create()
particle.setGeneratorStatus(1)
particle.setPDG(13)
particle.setCharge(-1.0)
particle.setMomentum(edm4hep.Vector3d(10.0, 0.0, 0.0))

# Calorimeter: (cellID, position, [(sim cellID, weight), ...]).
# The cone filter keeps hits within 0.2 rad of the particle (along +x).
caloHitDefs = [
    (1, (1000.0, 0.0, 0.0), [(101, 0.7), (102, 0.3)]),  # in the cone, two links
    (2, (1000.0, 50.0, 0.0), [(102, 1.0)]),  # in the cone, shares a sim hit with 1
    (3, (0.0, 1000.0, 0.0), [(103, 1.0)]),  # outside the cone
    (4, (1000.0, -50.0, 0.0), []),  # in the cone, no link
]
simCaloHits = edm4hep.SimCalorimeterHitCollection()
simCaloByID = {}
for simID in (101, 102, 103):
    simHit = simCaloHits.create()
    simHit.setCellID(simID)
    simCaloByID[simID] = simHit

caloHits = edm4hep.CalorimeterHitCollection()
caloLinks = edm4hep.CaloHitSimCaloHitLinkCollection()
for cellID, pos, links in caloHitDefs:
    hit = caloHits.create()
    hit.setCellID(cellID)
    hit.setPosition(edm4hep.Vector3f(*pos))
    for simID, weight in links:
        link = caloLinks.create()
        link.setFrom(hit)
        link.setTo(simCaloByID[simID])
        link.setWeight(weight)

# Tracker: the polar-angle split keeps hits with 50 deg < theta < 130 deg, and
# drops hits without links since it also writes out their sim hits.
trackerHitDefs = [
    (11, (100.0, 0.0, 0.0), [(111, 0.6), (112, 0.4)]),  # theta = 90 deg, two links
    (12, (0.0, 100.0, 0.0), [(112, 1.0)]),  # theta = 90 deg, shares a sim hit with 11
    (13, (10.0, 0.0, 100.0), [(113, 1.0)]),  # theta ~ 6 deg, outside the window
    (14, (0.0, -100.0, 0.0), []),  # theta = 90 deg, no link
]
simTrackerHits = edm4hep.SimTrackerHitCollection()
simTrackerByID = {}
for simID in (111, 112, 113):
    simHit = simTrackerHits.create()
    simHit.setCellID(simID)
    simTrackerByID[simID] = simHit

trackerHits = edm4hep.TrackerHitPlaneCollection()
trackerLinks = edm4hep.TrackerHitSimTrackerHitLinkCollection()
for cellID, pos, links in trackerHitDefs:
    hit = trackerHits.create()
    hit.setCellID(cellID)
    hit.setPosition(edm4hep.Vector3d(*pos))
    for simID, weight in links:
        link = trackerLinks.create()
        link.setFrom(hit)
        link.setTo(simTrackerByID[simID])
        link.setWeight(weight)

frame = Frame()
frame.put(mcParticles, "MCParticle")
frame.put(simCaloHits, "SimCaloHits")
frame.put(caloHits, "CaloHits")
frame.put(caloLinks, "CaloHitLinks")
frame.put(simTrackerHits, "SimTrackerHits")
frame.put(trackerHits, "TrackerHits")
frame.put(trackerLinks, "TrackerHitLinks")

# The file is closed by podio when the interpreter exits.
writer = root_io.Writer(args.output)
writer.write_frame(frame, "events")
