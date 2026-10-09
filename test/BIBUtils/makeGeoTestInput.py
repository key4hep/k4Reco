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

# Writes the input for the BIBUtils algorithms that need the detector geometry
# (run by runGeoAlgorithms.py, checked by checkGeoAlgorithms.py): a one-event file
# whose cellIDs and helices use the encodings and field of the geometry, and a
# small CaloHitSelector threshold-map file. Hits are identified by their type field.
import argparse
import math
import os
import re
import xml.etree.ElementTree as ET

import ROOT
import edm4hep
from podio import Frame, root_io

ROOT.gInterpreter.Declare(
    "#include <edm4hep/CaloHitSimCaloHitLinkCollection.h>\n"
    "#include <edm4hep/TrackerHitSimTrackerHitLinkCollection.h>\n"
)

parser = argparse.ArgumentParser(description="Write the input for the BIBUtils geometry tests")
parser.add_argument(
    "--compact",
    default=os.path.join(os.environ.get("K4GEO", ""), "MuColl/MAIA/compact/MAIA_v0/MAIA_v0.xml"),
    help="Compact file of the geometry",
)
parser.add_argument("--output", default="bibutils_geo_input.edm4hep.root", help="Output event file")
parser.add_argument("--thresholds", default="bibutils_thresholds.root", help="Output threshold-map file")
args = parser.parse_args()

# The constants are read directly from the compact file: building the full
# geometry with DD4hep only to read them would take tens of seconds.
compactRoot = ET.parse(args.compact).getroot()
compactConstants = {constant.get("name"): constant.get("value") for constant in compactRoot.iter("constant")}

ROOT.gSystem.Load("libDDCore")


def layerEncoder(encodingName):
    """Returns a function giving a cellID with only the layer field set."""
    coder = ROOT.dd4hep.DDSegmentation.BitFieldCoder(compactConstants[encodingName])
    offset = coder[coder.index("layer")].offset()

    def encode(layer):
        cellID = layer << offset
        assert coder.get(cellID, "layer") == layer
        return cellID

    return encode


# ---------------------------------------------------------------------------
# CaloHitSelector
#
# Threshold maps: 2 bins in the polar angle folded onto [0, pi/2] (iTheta = 1 below
# 45 deg, 2 above) times 4 layers (0-3), with
#   mode   = 0.02 * iTheta + 0.004 * layer   [GeV]
#   stddev = 0.001 * (iTheta + layer)        [GeV]
# so that the thresholds of any two cells differ by at least 6 MeV.
#
# Selector instances (see runGeoAlgorithms.py), all with time window (-0.5, 1.0) ns:
#   Map:    threshold = mode + 2 * stddev
#   Flat:   FlatThreshold = 0.05 GeV, which overrides the map
#   BIBSub: Nsigma = 0 with BIB subtraction, i.e. kept if E - mode > mode
# ---------------------------------------------------------------------------


def caloMode(iTheta, layer):
    return 0.02 * iTheta + 0.004 * layer


def caloStddev(iTheta, layer):
    return 0.001 * (iTheta + layer)


thresholdFile = ROOT.TFile(args.thresholds, "RECREATE")
modeMap = ROOT.TH2D("th_2dmode_sym", "", 2, 0.0, math.pi / 2, 4, 0.0, 4.0)
stddevMap = ROOT.TH2D("stddev_sym", "", 2, 0.0, math.pi / 2, 4, 0.0, 4.0)
for iTheta in (1, 2):
    for layer in range(4):
        modeMap.SetBinContent(iTheta, layer + 1, caloMode(iTheta, layer))
        stddevMap.SetBinContent(iTheta, layer + 1, caloStddev(iTheta, layer))
thresholdFile.Write()
thresholdFile.Close()

# (type, layer, theta [deg], energy [GeV], time [ns], [(sim cellID, weight), ...])
# The comments give the threshold of each instance and the expected outcome.
caloHitDefs = [
    # Map: own cell (iTheta 1, layer 0) threshold 0.022.
    (1, 0, 30.0, 0.0225, 0.0, [(1001, 0.7), (1002, 0.3)]),  # Map kept; Flat, BIBSub (0.04) dropped
    (2, 0, 30.0, 0.0215, 0.0, [(1003, 1.0)]),  # all dropped
    # Layer lookup: (iTheta 1, layer 2) threshold 0.034; layers 0, 1 and 3 give 0.022, 0.028, 0.040.
    (3, 2, 30.0, 0.0345, 0.0, [(1004, 1.0)]),  # Map kept, would be dropped if read as layer 3
    (4, 2, 30.0, 0.0335, 0.0, [(1005, 1.0)]),  # Map dropped, would be kept if read as layer 0 or 1
    # Theta bin: (iTheta 2, layer 1) threshold 0.050; iTheta 1 gives 0.028.
    (5, 1, 70.0, 0.0505, 0.0, [(1006, 1.0)]),  # Map, Flat kept
    (6, 1, 70.0, 0.0495, 0.0, [(1007, 1.0)]),  # all dropped; Map would keep it with iTheta 1
    # Folding: theta = 110 deg must use the 70 deg bin; unfolded it would fall outside the map (threshold 0).
    (7, 1, 110.0, 0.0495, 0.0, [(1008, 1.0)]),  # all dropped
    (8, 1, 110.0, 0.0505, 0.0, [(1009, 1.0)]),  # Map, Flat kept
    (9, 3, 150.0, 0.0405, 0.0, [(1010, 1.0)]),  # Map kept (folded to 30 deg: 0.040; 0.062 at 70 deg)
    # Time window (-0.5, 1.0) ns, both edges excluded; the energy passes every instance.
    (10, 0, 30.0, 0.2, -0.5, [(1011, 1.0)]),  # all dropped
    (11, 0, 30.0, 0.2, -0.45, [(1012, 1.0)]),  # all kept
    (12, 0, 30.0, 0.2, 0.95, []),  # all kept, no link
    (13, 0, 30.0, 0.2, 1.0, [(1013, 1.0)]),  # all dropped
    # BIB subtraction: (iTheta 1, layer 0) has mode 0.02, so BIBSub needs E > 0.04.
    (14, 0, 30.0, 0.03, 0.0, [(1014, 1.0)]),  # Map kept, BIBSub dropped (kept without the subtraction)
    (15, 0, 30.0, 0.041, 0.0, [(1015, 1.0)]),  # Map, BIBSub kept
    # Flat override: (iTheta 2, layer 3) map threshold 0.062, flat 0.05.
    (16, 3, 70.0, 0.055, 0.0, [(1016, 1.0)]),  # Flat kept, Map dropped
]

encodeCaloLayer = layerEncoder("GlobalCalorimeterReadoutID")
caloRadius = 2000.0  # [mm], only the polar angle matters

simCaloHits = edm4hep.SimCalorimeterHitCollection()
simCaloByID = {}
for _, _, _, _, _, links in caloHitDefs:
    for simID, _ in links:
        if simID not in simCaloByID:
            simHit = simCaloHits.create()
            simHit.setCellID(simID)
            simCaloByID[simID] = simHit

caloHits = edm4hep.CalorimeterHitCollection()
caloLinks = edm4hep.CaloHitSimCaloHitLinkCollection()
for hitType, layer, thetaDeg, energy, time, links in caloHitDefs:
    theta = math.radians(thetaDeg)
    hit = caloHits.create()
    hit.setType(hitType)
    hit.setCellID(encodeCaloLayer(layer))
    hit.setPosition(edm4hep.Vector3f(caloRadius * math.sin(theta), 0.0, caloRadius * math.cos(theta)))
    hit.setEnergy(energy)
    hit.setTime(time)
    for simID, weight in links:
        link = caloLinks.create()
        link.setFrom(hit)
        link.setTo(simCaloByID[simID])
        link.setWeight(weight)

# ---------------------------------------------------------------------------
# TrackerHitHelixFilter
#
# Hits are placed on trajectories obtained by integrating the Lorentz force, so
# that they do not depend on the helix parametrisation used by the algorithm.
# Filter instances (see runGeoAlgorithms.py), both with ConeAroundStatus = [1] and
# TrackerOuterRadius = 1500 mm:
#   Dist3D: hits closer than 30 mm to the helix are kept
#   DeltaR: hits seen at an angle below 0.05 rad from the helix are kept
# ---------------------------------------------------------------------------

# Field at the origin, which is what the algorithm takes from the geometry.
solenoid = next(field for field in compactRoot.iter("field") if field.get("type") == "solenoid")
fieldMatch = re.fullmatch(r"\s*([-+0-9.eE]+)\s*\*\s*tesla\s*", solenoid.get("inner_field"))
assert fieldMatch, f"Unexpected format of the solenoid field: {solenoid.get('inner_field')}"
bField = float(fieldMatch.group(1))  # [T], along +z

# Speed of light in GeV / (T mm): dp/ds = charge * kC * (p_hat x B) with p in GeV, s in mm.
kC = 0.299792458e-3


def trajectory(vertex, momentum, charge, step=0.5, maxLength=20000.0):
    """Yields (position, momentum) along the trajectory, every `step` mm of path
    (backwards if step < 0), integrating dp/ds = charge * kC * (p_hat x B) with RK4."""

    def derivatives(mom):
        p = math.sqrt(sum(c * c for c in mom))
        return [c / p for c in mom], [charge * kC * bField * mom[1] / p, -charge * kC * bField * mom[0] / p, 0.0]

    pos, mom = list(vertex), list(momentum)
    for _ in range(int(maxLength / abs(step))):
        kPos, kMom = [None] * 4, [None] * 4
        kPos[0], kMom[0] = derivatives(mom)
        for k in (1, 2, 3):
            f = step if k == 3 else 0.5 * step
            kPos[k], kMom[k] = derivatives([mom[i] + f * kMom[k - 1][i] for i in range(3)])
        for i in range(3):
            pos[i] += step / 6.0 * (kPos[0][i] + 2 * kPos[1][i] + 2 * kPos[2][i] + kPos[3][i])
            mom[i] += step / 6.0 * (kMom[0][i] + 2 * kMom[1][i] + 2 * kMom[2][i] + kMom[3][i])
        yield list(pos), list(mom)
    raise RuntimeError("Trajectory did not reach the requested point")


def rho(pos):
    return math.hypot(pos[0], pos[1])


def pointAtRadius(vertex, momentum, charge, radius, step=0.5):
    """First point of the trajectory at a transverse distance >= radius from the z axis."""
    return next((pos, mom) for pos, mom in trajectory(vertex, momentum, charge, step) if rho(pos) >= radius)


def pointOnReturnArc(vertex, momentum, charge, radius):
    """First point at a transverse distance <= radius after the trajectory has turned back towards the axis."""
    previous = 0.0
    turnedBack = False
    for pos, mom in trajectory(vertex, momentum, charge):
        turnedBack = turnedBack or rho(pos) < previous
        previous = rho(pos)
        if turnedBack and rho(pos) <= radius:
            return pos, mom


def displaced(pointAndMomentum, distance):
    """Point moved by `distance` in the transverse plane, perpendicular to the trajectory."""
    pos, mom = pointAndMomentum
    pT = math.hypot(mom[0], mom[1])
    return [pos[0] - distance * mom[1] / pT, pos[1] + distance * mom[0] / pT, pos[2]]


def fromPolar(pT, phiDeg, pz):
    return [pT * math.cos(math.radians(phiDeg)), pT * math.sin(math.radians(phiDeg)), pz]


origin = [0.0, 0.0, 0.0]
# P1: the particle the filter must build the cone around. Its radius in a 5 T field
# is about 1.33 m, so it leaves the 1.5 m tracker cylinder and comes back.
muonMom = fromPolar(2.0, 30.0, 4.0)
# P2: a neutral particle, P3: a particle with a generator status not in ConeAroundStatus.
photonMom = fromPolar(2.0, 200.0, 0.0)
status2Mom = fromPolar(2.0, 110.0, 0.0)

mcParticles = edm4hep.MCParticleCollection()
for pdg, charge, status, mom in ((13, -1.0, 1, muonMom), (22, 0.0, 1, photonMom), (-13, 1.0, 2, status2Mom)):
    particle = mcParticles.create()
    particle.setPDG(pdg)
    particle.setCharge(charge)
    particle.setGeneratorStatus(status)
    particle.setMomentum(edm4hep.Vector3d(*mom))

muonAt = {radius: pointAtRadius(origin, muonMom, -1.0, radius) for radius in (100.0, 500.0, 700.0, 1000.0)}
returnPoint = pointOnReturnArc(origin, muonMom, -1.0, 1000.0)[0]
# The return-arc hit must be dropped by the 1.5 m clipping alone, not by the hemisphere cut.
assert sum(returnPoint[i] * muonMom[i] for i in range(3)) > 0.0

# (type, position, [(sim cellID, weight), ...]); the comments give the expected
# outcome for the Dist3D and DeltaR instances.
trackerHitDefs = [
    (101, muonAt[100.0][0], [(2001, 0.6), (2002, 0.4)]),  # on the trajectory: kept by both
    (102, muonAt[500.0][0], [(2002, 1.0)]),  # kept by both, shares a sim hit with 101
    (103, muonAt[1000.0][0], [(2003, 1.0)]),  # kept by both
    (104, displaced(muonAt[1000.0], 50.0), [(2004, 1.0)]),  # 50 mm, angle ~0.02: DeltaR kept, Dist3D dropped
    (105, displaced(muonAt[100.0], 20.0), [(2005, 1.0)]),  # 20 mm, angle ~0.09: Dist3D kept, DeltaR dropped
    # On the trajectory of the opposite charge: hundreds of mm away at r = 1000 mm.
    (106, pointAtRadius(origin, muonMom, 1.0, 1000.0)[0], [(2006, 1.0)]),  # dropped by both
    # On the straight line along the momentum, ~375 mm from the helix at r = 1000 mm.
    (107, [c * 1000.0 / math.hypot(muonMom[0], muonMom[1]) for c in muonMom], [(2007, 1.0)]),  # dropped
    (108, pointAtRadius(origin, muonMom, -1.0, 200.0, step=-0.5)[0], [(2008, 1.0)]),  # behind the vertex: dropped
    (109, returnPoint, [(2009, 1.0)]),  # back inside after leaving the 1.5 m cylinder: dropped
    # Along the photon: dropped. This only checks that the photon adds no cone; even
    # without the neutral-particle skip, a charge-0 helix would not pass near this hit.
    (110, [c * 500.0 / 2.0 for c in photonMom], [(2010, 1.0)]),
    (111, pointAtRadius(origin, status2Mom, 1.0, 500.0)[0], [(2011, 1.0)]),  # on the status-2 particle: dropped
    (112, muonAt[700.0][0], []),  # on the trajectory, but without links: dropped
]

simTrackerHits = edm4hep.SimTrackerHitCollection()
simTrackerByID = {}
for _, _, links in trackerHitDefs:
    for simID, _ in links:
        if simID not in simTrackerByID:
            simHit = simTrackerHits.create()
            simHit.setCellID(simID)
            simTrackerByID[simID] = simHit

trackerHits = edm4hep.TrackerHitPlaneCollection()
trackerLinks = edm4hep.TrackerHitSimTrackerHitLinkCollection()
for hitType, pos, links in trackerHitDefs:
    hit = trackerHits.create()
    hit.setType(hitType)
    hit.setPosition(edm4hep.Vector3d(*pos))
    for simID, weight in links:
        link = trackerLinks.create()
        link.setFrom(hit)
        link.setTo(simTrackerByID[simID])
        link.setWeight(weight)

frame = Frame()
frame.put(simCaloHits, "SimCaloHits")
frame.put(caloHits, "CaloHits")
frame.put(caloLinks, "CaloHitLinks")
frame.put(mcParticles, "MCParticle")
frame.put(simTrackerHits, "SimTrackerHits")
frame.put(trackerHits, "TrackerHits")
frame.put(trackerLinks, "TrackerHitLinks")

# The file is closed by podio when the interpreter exits.
writer = root_io.Writer(args.output)
writer.write_frame(frame, "events")
