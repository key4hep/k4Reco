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

# A simple script to compare the hits from the Gaudi and Marlin output
import argparse
from podio.root_io import Reader

parser = argparse.ArgumentParser(description="Compare hits from Gaudi and Marlin")
parser.add_argument(
    "--gaudi-file", default="output_gaudi_lumical_clustering.root", help="Gaudi output file"
)
parser.add_argument(
    "--marlin-file", default="output_lumical_REC.edm4hep.root", help="Marlin output file"
)

parser.add_argument(
    "--gaudi-hit-collections",
    default=[
        "GaudiLumiCalHits",
    ], nargs="+",
    help="Gaudi hits collection",
)

parser.add_argument(
    "--marlin-hit-collections",
    default=[
        "LumiCalHits",
    ], nargs="+", help="Marlin hits collection"
)

parser.add_argument(
    "--gaudi-cluster-collections",
    default=[
        "GaudiLumiCalClusters",
    ], nargs="+",
    help="Gaudi hits collection",
)

parser.add_argument(
    "--marlin-cluster-collections",
    default=[
        "LumiCalClusters",
    ], nargs="+", help="Marlin hits collection"
)

parser.add_argument(
    "--gaudi-recoparticle-collections",
    default=[
        "GaudiLumiCalRecoParticles",
    ], nargs="+",
    help="Gaudi hits collection",
)

parser.add_argument(
    "--marlin-recoparticle-collections",
    default=[
        "LumiCalRecoParticles",
    ], nargs="+", help="Marlin hits collection"
)

args = parser.parse_args()

reader_gaudi = Reader(args.gaudi_file)
reader_marlin = Reader(args.marlin_file)

events_gaudi = reader_gaudi.get("events")
events_marlin = reader_marlin.get("events")

for i, frame_gaudi in enumerate(events_gaudi):
    frame_marlin = events_marlin[i]
    for collection_marlin, collection_gaudi in zip(args.marlin_hit_collections, args.gaudi_hit_collections):
        print(f'Checking collection "{collection_marlin}" (Marlin) and "{collection_gaudi}" (Gaudi)')
        hits_gaudi = frame_gaudi.get(collection_gaudi)
        hits_marlin = frame_marlin.get(collection_marlin)
        print(f"Checking event {i} with {len(hits_gaudi)} hits")
        assert len(hits_gaudi) == len(
            hits_marlin
        ), f"Number of hits differ: {len(hits_gaudi)} vs {len(hits_marlin)}"
        for j, (hit_gaudi, hit_marlin) in enumerate(zip(hits_gaudi, hits_marlin)):
            print(f"Checking hit {j}")
            for attr in [
                "CellID",
                "Energy",
                "EnergyError",
                "Time",
                "Position",
                "Type",
            ]:
                assert (
                    getattr(hit_gaudi, f"get{attr}")() == getattr(hit_marlin, f"get{attr}")()
                ), f"{attr} differ for hit {j}: {getattr(hit_gaudi, f'get{attr}')()} vs {getattr(hit_marlin, f'get{attr}')()}"

    for collection_marlin, collection_gaudi in zip(args.marlin_cluster_collections, args.gaudi_cluster_collections):
        print(f'Checking collection "{collection_marlin}" (Marlin) and "{collection_gaudi}" (Gaudi)')
        clusters_gaudi = frame_gaudi.get(collection_gaudi)
        clusters_marlin = frame_marlin.get(collection_marlin)
        print(f"Checking event {i} with {len(clusters_gaudi)} clusters")
        assert len(clusters_gaudi) == len(
            clusters_marlin
        ), f"Number of clusters differ: {len(clusters_gaudi)} vs {len(clusters_marlin)}"
        for j, (hit_gaudi, hit_marlin) in enumerate(zip(clusters_gaudi, clusters_marlin)):
            print(f"Checking cluster {j}")
            for attr in [
                "Type",
                "Energy",
                "EnergyError",
                "Position",
                "PositionError",
                "ITheta",
                "Phi",
                "DirectionError",
            ]:
                assert (
                    getattr(hit_gaudi, f"get{attr}")() == getattr(hit_marlin, f"get{attr}")()
                ), f"{attr} differ for hit {j}: {getattr(hit_gaudi, f'get{attr}')()} vs {getattr(hit_marlin, f'get{attr}')()}"

            for vmember in [
                "ShapeParameters",
                "SubdetectorEnergies",
            ]:
                assert (
                    list(getattr(hit_gaudi, f"get{vmember}")()) == list(getattr(hit_marlin, f"get{vmember}")())
                ), f"{vmember} differ for hit {j}: {getattr(hit_gaudi, f'get{vmember}')()} vs {getattr(hit_marlin, f'get{vmember}')()}"

            for relation in [
                "Clusters",
                "Hits",
            ]:
                assert (
                    [elem.id().index for elem in getattr(hit_gaudi, f"get{relation}")()] ==
                    [elem.id().index for elem in getattr(hit_marlin, f"get{relation}")()]
                ), f"{relation} differ for hit {j}: {getattr(hit_gaudi, f'get{relation}')()} vs {getattr(hit_marlin, f'get{relation}')()}"

    for collection_marlin, collection_gaudi in zip(args.marlin_recoparticle_collections, args.gaudi_recoparticle_collections):
        print(f'Checking collection "{collection_marlin}" (Marlin) and "{collection_gaudi}" (Gaudi)')
        clusters_gaudi = frame_gaudi.get(collection_gaudi)
        clusters_marlin = frame_marlin.get(collection_marlin)
        print(f"Checking event {i} with {len(clusters_gaudi)} reconstructed particles")
        assert len(clusters_gaudi) == len(
            clusters_marlin
        ), f"Number of clusters differ: {len(clusters_gaudi)} vs {len(clusters_marlin)}"
        for j, (hit_gaudi, hit_marlin) in enumerate(zip(clusters_gaudi, clusters_marlin)):
            print(f"Checking reconstructed particle {j}")

            for attr in [
                "PDG",
                "Energy",
                "Momentum",
                "ReferencePoint",
                "Charge",
                "Mass",
                "GoodnessOfPID",
                "CovMatrix",
            ]:
                assert (
                    getattr(hit_gaudi, f"get{attr}")() == getattr(hit_marlin, f"get{attr}")()
                ), f"{attr} differ for hit {j}: {getattr(hit_gaudi, f'get{attr}')()} vs {getattr(hit_marlin, f'get{attr}')()}"


            for relation in [
                "DecayVertex",
            ]:
                assert (
                    getattr(hit_gaudi, f"get{relation}")().id().index ==
                    getattr(hit_marlin, f"get{relation}")().id().index
                ), f"{relation} differ for hit {j}: {getattr(hit_gaudi, f'get{relation}')()} vs {getattr(hit_marlin, f'get{relation}')()}"

            for relation in [
                "Clusters",
                "Tracks",
                "Particles",
            ]:
                assert (
                    [elem.id().index for elem in getattr(hit_gaudi, f"get{relation}")()] ==
                    [elem.id().index for elem in getattr(hit_marlin, f"get{relation}")()]
                ), f"{relation} differ for hit {j}: {getattr(hit_gaudi, f'get{relation}')()} vs {getattr(hit_marlin, f'get{relation}')()}"
