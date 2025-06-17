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

set(sources
  DDPlanarDigi/components/DDPlanarDigi.cpp
  Overlay/components/OverlayTiming.cpp
)

set(libs
        podio::podioIO
        EDM4HEP::edm4hep
        k4FWCore::k4FWCore
        k4FWCore::k4Interface
        Gaudi::GaudiKernel
        DD4hep::DDRec
        DD4hep::DDCore
        ROOT::Core
        ROOT::MathCore
        ROOT::Tree
)

if(BUILD_TRACKING)

  find_package(LCIO REQUIRED)
  find_package(ILCUTIL COMPONENTS streamlog REQUIRED)
  find_package(KalTest REQUIRED)
  find_package(DDKalTest REQUIRED)

  gaudi_add_library(GaudiTrkUtils
    SOURCES
    GaudiTrkUtils/src/GaudiDDKalTest.cpp
    GaudiTrkUtils/src/GaudiDDKalTestTrack.cpp
    LINK PRIVATE
    LCIO::lcio EDM4HEP::edm4hep DD4hep::DDRec Gaudi::GaudiKernel
    ${KalTest_LIBRARIES}
    ${DDKalTest_LIBRARIES}
  )

  # Needed for the warnings in KalTest
  target_compile_options(GaudiTrkUtils PRIVATE "-Wno-error=deprecated-copy")

  target_sources(GaudiTrkUtils
    PUBLIC FILE_SET headers
    TYPE HEADERS FILES
    GaudiTrkUtils/include/GaudiDDKalTest.h
    GaudiTrkUtils/include/GaudiDDKalTestTrack.h
    BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/GaudiTrkUtils/include
  )

  target_include_directories(GaudiTrkUtils PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/GaudiTrkUtils/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE
    ${KalTest_INCLUDE_DIRS}
    ${DDKalTest_INCLUDE_DIRS}
    ${streamlog_INCLUDE_DIRS}
  )

  # Expected by KalTest
  target_compile_definitions(GaudiTrkUtils PRIVATE __NOT0__)

  list(APPEND sources
       ConformalTracking/components/ConformalTracking.cpp
       ConformalTracking/src/KDTrack.cpp
       ConformalTracking/src/kdtree2.cpp
       ConformalTracking/src/KDTree.cpp
       ConformalTracking/src/HelixTrack.cpp
       ConformalTracking/src/GaudiTrkUtils.cpp
       Tracking/src/ClonesAndSplitTracksFinder.cpp
       Tracking/src/RefitFinal.cpp
       Tracking/src/TruthTrackFinder.cpp
  )

  list(APPEND libs
       ${KalTest_LIBRARIES}
       ${DDKalTest_LIBRARIES}
       GaudiTrkUtils
  )
endif()


gaudi_add_module(k4RecoPlugins
                 SOURCES ${sources}
                 LINK ${libs}
)

# Common include directories
set(k4RecoPlugins_public_includes
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_include_directories(k4RecoPlugins PUBLIC ${k4RecoPlugins_public_includes})

if(BUILD_TRACKING)
  target_include_directories(k4RecoPlugins
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/ConformalTracking/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/Tracking/include>
    PRIVATE
    ${k4RecoPlugins_public_includes}
    ${KalTest_INCLUDE_DIRS}
    ${DDKalTest_INCLUDE_DIRS}
    ${streamlog_INCLUDE_DIRS}
  )
endif()

install(TARGETS GaudiTrkUtils k4RecoPlugins
  EXPORT k4RecoTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  FILE_SET headers DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/k4Reco"
  INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  COMPONENT dev)

