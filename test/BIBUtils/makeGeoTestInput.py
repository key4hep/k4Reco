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
# whose cellIDs are encoded with the encodings of the geometry, and a small
# CaloHitSelector threshold-map file. Hits are identified by their type field.
import argparse
import math
import os
import xml.etree.ElementTree as ET

import ROOT
import edm4hep
from podio import Frame, root_io

ROOT.gInterpreter.Declare("#include <edm4hep/CaloHitSimCaloHitLinkCollection.h>\n")

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
compactConstants = {
    constant.get("name"): constant.get("value") for constant in ET.parse(args.compact).getroot().iter("constant")
}

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

frame = Frame()
frame.put(simCaloHits, "SimCaloHits")
frame.put(caloHits, "CaloHits")
frame.put(caloLinks, "CaloHitLinks")

# The file is closed by podio when the interpreter exits.
writer = root_io.Writer(args.output)
writer.write_frame(frame, "events")
