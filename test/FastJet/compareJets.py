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

# A simple script to compare the jets from the Gaudi and Marlin output
import argparse
from podio.reading import get_reader

parser = argparse.ArgumentParser(description="Compare jets from Gaudi and Marlin")
parser.add_argument(
    "--gaudi-file", default="output_fastjet_ttbar.root", help="Gaudi output file"
)
parser.add_argument(
    "--marlin-file", default="output_marlin_fastjet_ttbar.root", help="Marlin output file"
)
parser.add_argument("--gaudi-jets", default="JetOut", help="Gaudi jet collection")
parser.add_argument("--marlin-jets", default="JetOut", help="Marlin jet collection")
parser.add_argument(
    "--gaudi-constituents", default="UsedPFOs", help="Gaudi jet constituent collection"
)
parser.add_argument(
    "--marlin-constituents", default="UsedPFOs", help="Marlin jet constituent collection"
)

args = parser.parse_args()

# Everything a ReconstructedParticle carries. Only energy, mass, momentum and the
# particle relation are actually set on the jets, but comparing all of them also
# catches anything that ends up set on one side and not on the other
MEMBERS = [
    "PDG",
    "Energy",
    "Momentum",
    "ReferencePoint",
    "Charge",
    "Mass",
    "GoodnessOfPID",
    "CovMatrix",
]
MULTI_RELATIONS = ["Clusters", "Tracks", "Particles"]
SINGLE_RELATIONS = ["DecayVertex"]


def multi_relation(particle, relation):
    # The related objects can only be compared through their index within the
    # collection they live in, the collection IDs differ between the two files
    return [elem.id().index for elem in getattr(particle, f"get{relation}")()]


def single_relation(particle, relation):
    related = getattr(particle, f"get{relation}")()
    return related.id().index if related.isAvailable() else None


def compare(particles_gaudi, particles_marlin, what, event):
    print(f"Checking event {event} with {len(particles_gaudi)} {what}")
    assert len(particles_gaudi) == len(
        particles_marlin
    ), f"Number of {what} differ in event {event}: {len(particles_gaudi)} vs {len(particles_marlin)}"

    for j, (gaudi, marlin) in enumerate(zip(particles_gaudi, particles_marlin)):
        for member in MEMBERS:
            value_gaudi = getattr(gaudi, f"get{member}")()
            value_marlin = getattr(marlin, f"get{member}")()
            assert (
                value_gaudi == value_marlin
            ), f"{member} differ for {what} {j} in event {event}: {value_gaudi} vs {value_marlin}"

        for relation in MULTI_RELATIONS:
            value_gaudi = multi_relation(gaudi, relation)
            value_marlin = multi_relation(marlin, relation)
            assert (
                value_gaudi == value_marlin
            ), f"{relation} differ for {what} {j} in event {event}: {value_gaudi} vs {value_marlin}"

        for relation in SINGLE_RELATIONS:
            value_gaudi = single_relation(gaudi, relation)
            value_marlin = single_relation(marlin, relation)
            assert (
                value_gaudi == value_marlin
            ), f"{relation} differ for {what} {j} in event {event}: {value_gaudi} vs {value_marlin}"


reader_gaudi = get_reader(args.gaudi_file)
reader_marlin = get_reader(args.marlin_file)

events_gaudi = reader_gaudi.get("events")
events_marlin = reader_marlin.get("events")

assert len(events_gaudi) == len(
    events_marlin
), f"Number of events differ: {len(events_gaudi)} vs {len(events_marlin)}"

for i, frame_gaudi in enumerate(events_gaudi):
    frame_marlin = events_marlin[i]
    compare(
        frame_gaudi.get(args.gaudi_jets), frame_marlin.get(args.marlin_jets), "jets", i
    )
    compare(
        frame_gaudi.get(args.gaudi_constituents),
        frame_marlin.get(args.marlin_constituents),
        "jet constituents",
        i,
    )
