#!/usr/bin/env python3
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
"""Compare EDM4hep-native LCFIPlus output against the Marlin reference."""
import argparse
import math
from podio.root_io import Reader, RNTupleReader

p = argparse.ArgumentParser()
p.add_argument("--gaudi-file", default="output_lcfiplus_native.root")
p.add_argument("--marlin-file", default="output_lcfiplus_marlin_REC.edm4hep.root")
p.add_argument("--verbose", action="store_true", help="print per-stage summaries before comparing")
args = p.parse_args()

def summary(frame, vertices, jets):
    vs = frame.get(vertices)
    js = frame.get(jets)
    return (len(vs), len(js), sum(v.getChi2() for v in vs),
            sum(j.getEnergy() for j in js), sum(len(j.getParticles()) for j in js))

def stage_summary(frame, vertices, jets):
    vs = frame.get(vertices) if vertices else ()
    js = frame.get(jets) if jets else ()
    return (len(vs), len(js), sum(v.getChi2() for v in vs),
            sum(j.getEnergy() for j in js), sum(len(j.getParticles()) for j in js))

def vertex_memberships(frame, vertices):
    return [(round(vertex.getChi2(), 6), [particle.getObjectID().index for particle in vertex.getParticles()])
            for vertex in frame.get(vertices)]

def vertex_position(vertex):
    position = vertex.getPosition()
    return (position.x, position.y, position.z)

def vertex_covariance(vertex):
    covariance = vertex.getCovMatrix()
    return tuple(covariance[index] for index in range(6))

def compare_vertices(event, native_frame, marlin_frame, native_collection, marlin_collection):
    """Require an exact object-by-object match for one EDM4hep vertex collection."""
    native_vertices = native_frame.get(native_collection)
    marlin_vertices = marlin_frame.get(marlin_collection)
    assert len(native_vertices) == len(marlin_vertices), (
        f"event {event} {native_collection}/{marlin_collection}: "
        f"vertex count {len(native_vertices)} != {len(marlin_vertices)}"
    )

    for index, (native_vertex, marlin_vertex) in enumerate(zip(native_vertices, marlin_vertices)):
        # edm4hep.yaml: Vertex has scalar members type, chi2, ndf, position,
        # covMatrix and algorithmType; parameters; and particles. cppyy does
        # not provide value equality for Vector3f/CovMatrix3f, so compare their
        # YAML-defined scalar components explicitly.
        for attribute in ("Type", "Chi2", "Ndf", "AlgorithmType"):
            native_value = getattr(native_vertex, f"get{attribute}")()
            marlin_value = getattr(marlin_vertex, f"get{attribute}")()
            assert native_value == marlin_value, (
                f"event {event} {native_collection}/{marlin_collection} vertex {index}: "
                f"{attribute} {native_value} != {marlin_value}"
            )

        for attribute, native_value, marlin_value in (
            ("Position", vertex_position(native_vertex), vertex_position(marlin_vertex)),
            ("CovMatrix", vertex_covariance(native_vertex), vertex_covariance(marlin_vertex)),
        ):
            assert native_value == marlin_value, (
                f"event {event} {native_collection}/{marlin_collection} vertex {index}: "
                f"{attribute} {native_value} != {marlin_value}"
            )

        native_parameters = list(native_vertex.getParameters())
        marlin_parameters = list(marlin_vertex.getParameters())
        assert native_parameters == marlin_parameters, (
            f"event {event} {native_collection}/{marlin_collection} vertex {index}: "
            f"parameters {native_parameters} != {marlin_parameters}"
        )

        native_particles = [particle.getObjectID().index for particle in native_vertex.getParticles()]
        marlin_particles = [particle.getObjectID().index for particle in marlin_vertex.getParticles()]
        assert native_particles == marlin_particles, (
            f"event {event} {native_collection}/{marlin_collection} vertex {index}: "
            f"particle indices {native_particles} != {marlin_particles}"
        )

def read_events(filename):
    """Accept either podio's RNTuple output or its legacy ROOT-tree output."""
    try:
        return RNTupleReader([filename]).get("events")
    except Exception:
        return Reader([filename]).get("events")

native = read_events(args.gaudi_file)
marlin = read_events(args.marlin_file)
assert len(native) == len(marlin), "event count differs"
errors = []
for event, (n, m) in enumerate(zip(native, marlin)):
    for native_vertices, marlin_vertices in (
        ("GAUDIPrimaryVertices", "PrimaryVertices"),
        ("GAUDIBuildUpVertices", "BuildUpVertices"),
        ("GAUDIBuildUpVertices_V0", "BuildUpVertices_V0"),
        ("GAUDIRefinedVertices", "RefinedVertices"),
    ):
        compare_vertices(event, n, m, native_vertices, marlin_vertices)

    if args.verbose:
        stages = (("primary", "GAUDIPrimaryVertices", "PrimaryVertices", None, None),
                  ("build-up", "GAUDIBuildUpVertices", "BuildUpVertices", None, None),
                  ("V0", "GAUDIBuildUpVertices_V0", "BuildUpVertices_V0", None, None),
                  ("jet", None, None, "GAUDIVertexJets", "VertexJets"),
                  ("refined", "GAUDIRefinedVertices", "RefinedVertices", "GAUDIRefinedVertexJets", "RefinedVertexJets"))
        for label, native_vertices, marlin_vertices, native_jets, marlin_jets in stages:
            native_stage = stage_summary(n, native_vertices, native_jets)
            marlin_stage = stage_summary(m, marlin_vertices, marlin_jets)
            print(f"event {event} {label}: native={native_stage} marlin={marlin_stage}")
        print(f"event {event} primary memberships: native={vertex_memberships(n, 'GAUDIPrimaryVertices')} marlin={vertex_memberships(m, 'PrimaryVertices')}")
        print(f"event {event} build-up memberships: native={vertex_memberships(n, 'GAUDIBuildUpVertices')} marlin={vertex_memberships(m, 'BuildUpVertices')}")
    a = summary(n, "GAUDIRefinedVertices", "GAUDIRefinedVertexJets")
    b = summary(m, "RefinedVertices", "RefinedVertexJets")
    if a[:2] != b[:2]:
        errors.append(f"event {event}: vertex/jet multiplicities {a[:2]} != {b[:2]}")
    for native_value, marlin_value, label in zip(a[2:], b[2:], ("vertex chi2", "jet energy", "jet constituents")):
        if not math.isclose(native_value, marlin_value, rel_tol=2e-4, abs_tol=2e-5):
            errors.append(f"event {event}: {label} {native_value} != {marlin_value}")
assert not errors, "\n".join(errors)
