/*
 * Copyright (c) 2020-2024 Key4hep-Project.
 *
 * This file is part of Key4hep.
 * See https://key4hep.github.io/key4hep-doc/ for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#pragma once
#include <edm4hep/ReconstructedParticleCollection.h>
#include <edm4hep/VertexCollection.h>
#include <k4FWCore/Transformer.h>
struct JetVertexRefiner final
    : k4FWCore::MultiTransformer<std::tuple<edm4hep::ReconstructedParticleCollection, edm4hep::VertexCollection>(
          const edm4hep::ReconstructedParticleCollection&, const edm4hep::ReconstructedParticleCollection&,
          const edm4hep::VertexCollection&, const edm4hep::VertexCollection&, const edm4hep::VertexCollection&)> {
  JetVertexRefiner(const std::string&, ISvcLocator*);
  std::tuple<edm4hep::ReconstructedParticleCollection, edm4hep::VertexCollection>
  operator()(const edm4hep::ReconstructedParticleCollection&, const edm4hep::ReconstructedParticleCollection&,
             const edm4hep::VertexCollection&, const edm4hep::VertexCollection&,
             const edm4hep::VertexCollection&) const override;
};
