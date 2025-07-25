#[[
Copyright (c) 2020-2024 Key4hep-Project.

This file is part of Key4hep.
See https://key4hep.github.io/key4hep-doc/ for further info.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
]]
cmake_minimum_required(VERSION 3.23)

project(k4Reco LANGUAGES CXX)

# please keep this layout for version setting, used by the automatic tagging script
set(${PROJECT_NAME}_VERSION_MAJOR 0)
set(${PROJECT_NAME}_VERSION_MINOR 2)
set(${PROJECT_NAME}_VERSION_PATCH 1)

set(${PROJECT_NAME}_VERSION "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}.${${PROJECT_NAME}_VERSION_PATCH}")

option(BUILD_TRACKING "Build conformal tracking" ON)

find_package(podio 1.2.99 REQUIRED)
find_package(EDM4HEP REQUIRED)
find_package(DD4hep REQUIRED)
find_package(Gaudi REQUIRED)
find_package(k4FWCore REQUIRED)
find_package(ROOT REQUIRED COMPONENTS RIO Tree)
find_package(k4SimGeant4 REQUIRED)

include(cmake/Key4hepConfig.cmake)

include(GNUInstallDirs)

add_subdirectory(k4Reco)
include(CTest)
if(BUILD_TESTING)
  add_subdirectory(test)
endif()

add_subdirectory(python)

include(cmake/CreateProjectConfig.cmake)
