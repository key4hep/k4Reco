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

# Checks the output of runLinkPreservation.py: the filters must keep every link
# of the hits they accept, with its original weight, and write each linked sim
# hit only once. Objects are identified by cellID (see makeLinkTestInput.py).
import argparse
import sys

from podio.reading import get_reader

parser = argparse.ArgumentParser(description="Check the links written by the BIBUtils filters")
parser.add_argument("--input", default="bibutils_links_output.edm4hep.root", help="Output of the filters")
args = parser.parse_args()

frame = get_reader(args.input).get("events")[0]


def cellIDs(collection):
    return sorted(obj.getCellID() for obj in collection)


def links(collection):
    return sorted((link.getFrom().getCellID(), link.getTo().getCellID(), round(link.getWeight(), 5)) for link in collection)


expected = {
    "CaloHitsConed hits": ([1, 2, 4], cellIDs(frame.get("CaloHitsConed"))),
    "CaloHitLinksConed links": (
        [(1, 101, 0.7), (1, 102, 0.3), (2, 102, 1.0)],
        links(frame.get("CaloHitLinksConed")),
    ),
    "TrackerHitsSplit hits": ([11, 12], cellIDs(frame.get("TrackerHitsSplit"))),
    "SimTrackerHitsSplit sim hits": ([111, 112], cellIDs(frame.get("SimTrackerHitsSplit"))),
    "TrackerHitLinksSplit links": (
        [(11, 111, 0.6), (11, 112, 0.4), (12, 112, 1.0)],
        links(frame.get("TrackerHitLinksSplit")),
    ),
}

failed = False
for label, (want, got) in expected.items():
    if want != got:
        print(f"FAIL {label}: expected {want}, got {got}")
        failed = True

if failed:
    sys.exit(1)
print("All BIBUtils link checks passed")
