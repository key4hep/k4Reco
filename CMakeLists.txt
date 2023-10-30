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
CMAKE_MINIMUM_REQUIRED(VERSION 3.12)

project(k4Reco)

# please keep this layout for version setting, used by the automatic tagging script
set(PACKAGE_VERSION_MAJOR 0)
set(PACKAGE_VERSION_MINOR 0)
set(PACKAGE_VERSION_PATCH 0)

set(PACKAGE_VERSION "${PACKAGE_VERSION_MAJOR}.${PACKAGE_VERSION_MINOR}.${PACKAGE_VERSION_PATCH}")

find_package(ROOT REQUIRED COMPONENTS RIO Tree)
find_package(EDM4HEP REQUIRED)
find_package(k4FWCore REQUIRED)
find_package(Gaudi REQUIRED)
find_package(DD4hep REQUIRED)
find_package(k4SimGeant4 REQUIRED)

include(GNUInstallDirs)

# Set up C++ Standard
# ``-DCMAKE_CXX_STANDARD=<standard>`` when invoking CMake
set(CMAKE_CXX_STANDARD 17 CACHE STRING "")

if(NOT CMAKE_CXX_STANDARD MATCHES "14|17|20")
  message(FATAL_ERROR "Unsupported C++ standard: ${CMAKE_CXX_STANDARD}")
endif()

include(CTest)

add_subdirectory(k4Reco)
add_subdirectory(test)

# include(cmake/CreateProjectConfig.cmake)
