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

if(BUILD_CONFORMAL_TRACKING)

  find_package(LCIO REQUIRED)
  find_package(KalTest REQUIRED)
  find_package(DDKalTest REQUIRED)

  list(APPEND sources
       ConformalTracking/components/ConformalTracking.cpp
       ConformalTracking/src/GaudiDDKalTest.cpp
       ConformalTracking/src/GaudiDDKalTestTrack.cpp
       ConformalTracking/src/KDTrack.cpp
       ConformalTracking/src/kdtree2.cpp
       ConformalTracking/src/KDTree.cpp
       ConformalTracking/src/HelixTrack.cpp
       ConformalTracking/src/GaudiTrkUtils.cpp
  )

  list(APPEND libs
       LCIO::lcio
       ${KalTest_LIBRARIES}
       ${DDKalTest_LIBRARIES}
  )


  # ADD_DEFINITIONS( -D __NOT0__ )
  # include_directories("${KalTest_INCLUDE_DIRS};${DDKalTest_INCLUDE_DIRS}")
  # link_libraries("${KalTest_LIBRARIES};${DDKalTest_LIBRARIES}")

endif()


gaudi_add_module(k4RecoPlugins
                 SOURCES ${sources}
                 LINK
                 ${libs}
)
target_include_directories(k4RecoPlugins PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/ConformalTracking/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

if(BUILD_CONFORMAL_TRACKING)
  # Expected by KalTest
  target_compile_definitions(k4RecoPlugins PRIVATE __NOT0__)

  target_include_directories(k4RecoPlugins PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/ConformalTracking/include>
    ${KalTest_INCLUDE_DIRS}
    ${DDKalTest_INCLUDE_DIRS}
    ${streamlog_INCLUDE_DIRS}
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

endif()

install(TARGETS k4RecoPlugins
  EXPORT ${CMAKE_PROJECT_NAME}Targets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/@{CMAKE_PROJECT_NAME}"
  COMPONENT dev)
