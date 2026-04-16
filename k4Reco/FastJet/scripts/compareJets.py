#!/usr/bin/env python
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
        marlinJets = gaudiFrame.get(args.marlin_jets)

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
            ), f"{eventStr} {jetStr} Jet momentum discrepancy, y component, Gaudi {gaudiJet.getMomentum().x}, Marlin: {marlinJet.getMomentum().x}"

            assert(
                gaudiJet.getMomentum().z == marlinJet.getMomentum().z
            ), f"{eventStr} {jetStr} Jet momentum discrepancy, z component, Gaudi {gaudiJet.getMomentum().x}, Marlin: {marlinJet.getMomentum().x}"
            

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
