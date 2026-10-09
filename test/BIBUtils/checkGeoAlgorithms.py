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

# Checks the output of runGeoAlgorithms.py. The hits each algorithm instance must
# keep are listed below, identified by their type field; the reasons are given
# next to the hit definitions in makeGeoTestInput.py. Every link of a kept hit
# must be copied with its weight.
import argparse
import sys

from podio.reading import get_reader

parser = argparse.ArgumentParser(description="Check the output of the BIBUtils geometry tests")
parser.add_argument("--input", default="bibutils_geo_input.edm4hep.root", help="Input of the algorithms")
parser.add_argument("--output", default="bibutils_geo_output.edm4hep.root", help="Output of the algorithms")
args = parser.parse_args()

inputFrame = get_reader(args.input).get("events")[0]
outputFrame = get_reader(args.output).get("events")[0]


def hitTypes(collection):
    return sorted(hit.getType() for hit in collection)


def links(collection):
    return sorted(
        (link.getFrom().getType(), link.getTo().getCellID(), round(link.getWeight(), 5)) for link in collection
    )


def expectedLinks(inputLinks, keptTypes):
    return [link for link in links(inputLinks) if link[0] in keptTypes]


expectedCaloHits = {
    "Map": [1, 3, 5, 8, 9, 11, 12, 14, 15],
    "Flat": [5, 8, 11, 12, 16],
    "BIBSub": [11, 12, 15],
}

checks = {}
for label, kept in expectedCaloHits.items():
    checks[f"CaloHits{label} hits"] = (kept, hitTypes(outputFrame.get(f"CaloHits{label}")))
    checks[f"CaloHitLinks{label} links"] = (
        expectedLinks(inputFrame.get("CaloHitLinks"), kept),
        links(outputFrame.get(f"CaloHitLinks{label}")),
    )

failed = False
for label, (want, got) in checks.items():
    if want != got:
        print(f"FAIL {label}: expected {want}, got {got}")
        failed = True

if failed:
    sys.exit(1)
print("All BIBUtils geometry checks passed")
