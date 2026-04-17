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
import argparse
from podio.reading import get_reader

def main(args) -> int:
    gaudiReader = get_reader(args.gaudi_file)
    marlinReader = get_reader(args.marlin_file)

    gaudiEvents = gaudiReader.get("events")
    marlinEvents = marlinReader.get("events")

    for i, gaudiFrame in enumerate(gaudiEvents):
        marlinFrame = marlinEvents[i]

        gaudiJets = gaudiFrame.get(args.gaudi_jets)
        marlinJets = marlinFrame.get(args.marlin_jets)

        eventStr = f"Event #{i:0>4}"

        assert( len(gaudiJets) == len(marlinJets)), f"{eventStr} Number of jets differ, Gaudi: {len(gaudiJets)}, Marlin: {len(marlinJets)}"

        for j, (gaudiJet, marlinJet) in enumerate(zip(gaudiJets, marlinJets)):
            jetStr = f"Jet #{j:0>4}"
            assert(
                gaudiJet.getEnergy() == marlinJet.getEnergy()
            ), f"{eventStr} {jetStr} Jet energy discrepancy, Gaudi {gaudiJet.getEnergy()}, Marlin: {marlinJet.getEnergy()}"

            assert(
                gaudiJet.getMass() == marlinJet.getMass()
            ), f"{eventStr} {jetStr} Jet mass discrepancy, Gaudi {gaudiJet.getMass()}, Marlin: {marlinJet.getMass()}"

            assert(
                gaudiJet.getMomentum().x == marlinJet.getMomentum().x
            ), f"{eventStr} {jetStr} Jet momentum discrepancy, x component, Gaudi {gaudiJet.getMomentum().x}, Marlin: {marlinJet.getMomentum().x}"

            assert(
                gaudiJet.getMomentum().y == marlinJet.getMomentum().y
            ), f"{eventStr} {jetStr} Jet momentum discrepancy, y component, Gaudi {gaudiJet.getMomentum().y}, Marlin: {marlinJet.getMomentum().y}"

            assert(
                gaudiJet.getMomentum().z == marlinJet.getMomentum().z
            ), f"{eventStr} {jetStr} Jet momentum discrepancy, z component, Gaudi {gaudiJet.getMomentum().z}, Marlin: {marlinJet.getMomentum().z}"


    print("Comparison succeeded!")
    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = "Compare Gaudi FastJet jets to Marlin"
    )

    parser.add_argument(
        "--gaudi-file",
        default = "output_pandoraPlusFastJet_ttbar.root",
        help="File containing the gaudi created fast jet outputs"
    )

    parser.add_argument(
        "--marlin-file",
        default = "output_marlinPandoraPlusFastJet_ttbar.root",
        help = "File containing the marlin created fast jet outputs"
    )

    parser.add_argument(
        "--gaudi-jets",
        default = "JetOut",
        help = "Gaudi jet collection"
    )

    parser.add_argument(
        "--marlin-jets",
        default = "JetOut",
        help = "Marlin jet collection"
    )

    args = parser.parse_args()

    returnCode = main(args)
    exit(returnCode)
