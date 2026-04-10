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


find_package(k4geo REQUIRED)

function(set_test_env testname)
  set_property(TEST ${testname} APPEND PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${PROJECT_BINARY_DIR}:${PROJECT_BINARY_DIR}/${PROJECT_NAME}:${PROJECT_BINARY_DIR}/${PROJECT_NAME}/genConfDir/${PROJECT_NAME}:${PROJECT_BINARY_DIR}/test/k4FWCoreTest:${PROJECT_BINARY_DIR}/test/k4FWCoreTest/genConfDir/k4FWCoreTest:$<$<TARGET_EXISTS:ROOT::Core>:$<TARGET_FILE_DIR:ROOT::Core>>:$<$<TARGET_EXISTS:EDM4HEP::edm4hep>:$<TARGET_FILE_DIR:EDM4HEP::edm4hep>>:$<$<TARGET_EXISTS:podio::podio>:$<TARGET_FILE_DIR:podio::podio>>:$ENV{LD_LIBRARY_PATH}")
  set_property(TEST ${testname} APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${PROJECT_SOURCE_DIR}/python:${PROJECT_BINARY_DIR}/${PROJECT_NAME}/genConfDir:/${PROJECT_NAME}/genConfDir:$ENV{PYTHONPATH}")
endfunction()

add_test(NAME "clone_CLDConfig"
  COMMAND bash -c "if [ ! -d CLDConfig ]; then git clone https://github.com/key4hep/CLDConfig.git --depth 1; fi && cp -r CLDConfig/CLDConfig/* . && sed -i -e '/Overlay.Overlay/d' -e '/HighLevelReco.RecoMCTruthLink/d' CLDReconstruction.py"
)
set_test_env("clone_CLDConfig")
set_tests_properties("clone_CLDConfig" PROPERTIES FIXTURES_SETUP CLDConfig_fixture)

add_test(NAME "run_ddsim"
  COMMAND bash -c "ddsim --compactFile $K4GEO/FCCee/CLD/compact/CLD_o2_v07/CLD_o2_v07.xml -G --numberOfEvents 3 --outputFile sim.edm4hep.root --gun.energy '10*GeV' --gun.particle 'mu-' --gun.distribution uniform --gun.multiplicity 10"
  )
set_test_env("run_ddsim")
set_tests_properties("run_ddsim" PROPERTIES FIXTURES_SETUP SimOutput_fixture)

# If this is modified, make sure that this simulations produces clusters and reconstructed particles in the LumiCal
# as this doesn't always happen with the current setting (but does for the current seed)
add_test(NAME "run_ddsim_LumiCal"
  COMMAND bash -c "ddsim --compactFile $K4GEO/FCCee/CLD/compact/CLD_o2_v07/CLD_o2_v07.xml -G --numberOfEvents 3 --outputFile sim_lumical.edm4hep.root --gun.energy '50*GeV' --gun.particle 'gamma' --gun.distribution uniform --gun.multiplicity 10 --gun.thetaMin='50*mrad' --gun.thetaMax='100*mrad'"
  )
set_test_env("run_ddsim_LumiCal")
set_tests_properties("run_ddsim_LumiCal" PROPERTIES FIXTURES_SETUP SimLumiCalOutput_fixture)

add_test(NAME "run_marlin_wrapper_truth_tracking"
  COMMAND k4run CLDReconstruction.py --inputFiles sim.edm4hep.root --trackingOnly --truthTracking --outputBasename output_truth_tracking
)
set_test_env("run_marlin_wrapper_truth_tracking")
set_tests_properties("run_marlin_wrapper_truth_tracking" PROPERTIES
  FIXTURES_REQUIRED "CLDConfig_fixture;SimOutput_fixture"
  FIXTURES_SETUP TruthTrackingOutput_fixture)

add_test(NAME "run_wrapper"
  COMMAND k4run CLDReconstruction.py --inputFiles sim.edm4hep.root --trackingOnly
)
set_test_env("run_wrapper")
set_tests_properties("run_wrapper" PROPERTIES
  FIXTURES_REQUIRED "CLDConfig_fixture;SimOutput_fixture"
  FIXTURES_SETUP WrapperOutput_fixture)

add_test(NAME "run_wrapper_LumiCal"
  COMMAND k4run CLDReconstruction.py --inputFiles sim_lumical.edm4hep.root --outputBasename output_lumical
)
set_test_env("run_wrapper_LumiCal")
set_tests_properties("run_wrapper_LumiCal" PROPERTIES
  FIXTURES_REQUIRED "CLDConfig_fixture;SimLumiCalOutput_fixture"
  FIXTURES_SETUP WrapperLumiCalOutput_fixture)

add_test(NAME "run_GaudiLumiCalClusterer"
  COMMAND k4run ${PROJECT_SOURCE_DIR}/k4Reco/GaudiLumiCalClusterer/options/runGaudiLumiCalClusterer.py --IOSvc.Input output_lumical_REC.edm4hep.root
)
set_test_env("run_GaudiLumiCalClusterer")
set_tests_properties("run_GaudiLumiCalClusterer" PROPERTIES
  FIXTURES_REQUIRED "WrapperLumiCalOutput_fixture"
  FIXTURES_SETUP GaudiLumiCalOutput_fixture)

add_test(NAME "compare_LumiCalClusterer"
  COMMAND "${PROJECT_SOURCE_DIR}/k4Reco/GaudiLumiCalClusterer/scripts/compare-lumical.py"
)
set_test_env("compare_LumiCalClusterer")
set_tests_properties("compare_LumiCalClusterer" PROPERTIES FIXTURES_REQUIRED "WrapperLumiCalOutput_fixture;GaudiLumiCalOutput_fixture")

if(BUILD_TRACKING)

  add_test(NAME "run_ConformalTracking"
    COMMAND k4run ${PROJECT_SOURCE_DIR}/k4Reco/ConformalTracking/options/runConformalTracking.py
  )
  set_test_env("run_ConformalTracking")
  set_tests_properties("run_ConformalTracking" PROPERTIES
    FIXTURES_REQUIRED "WrapperOutput_fixture"
    FIXTURES_SETUP ConformalTrackingOutput_fixture)

  add_test(NAME "compare_ConformalTracking"
    COMMAND "${PROJECT_SOURCE_DIR}/k4Reco/ConformalTracking/scripts/compare-tracks.py"
  )
  set_test_env("compare_ConformalTracking")
  set_tests_properties("compare_ConformalTracking" PROPERTIES FIXTURES_REQUIRED "ConformalTrackingOutput_fixture")

  add_test(NAME "run_ClonesAndSplitTracksFinder"
    COMMAND k4run ${PROJECT_SOURCE_DIR}/k4Reco/Tracking/options/runClonesAndSplitTracksFinder.py
  )
  set_test_env("run_ClonesAndSplitTracksFinder")
  set_tests_properties("run_ClonesAndSplitTracksFinder" PROPERTIES
    FIXTURES_REQUIRED "WrapperOutput_fixture"
    FIXTURES_SETUP ClonesAndSplitTracksFinderOutput_fixture)

  add_test(NAME "compare_ClonesAndSplitTracksFinder"
    COMMAND "${PROJECT_SOURCE_DIR}/k4Reco/ConformalTracking/scripts/compare-tracks.py" --gaudi-file output_clones_and_split_tracks_finder.root --marlin-tracks SiTracks --gaudi-tracks GaudiSiTracks
  )
  set_test_env("compare_ClonesAndSplitTracksFinder")
  set_tests_properties("compare_ClonesAndSplitTracksFinder" PROPERTIES FIXTURES_REQUIRED "ClonesAndSplitTracksFinderOutput_fixture")

  add_test(NAME "run_RefitFinal"
    COMMAND k4run ${PROJECT_SOURCE_DIR}/k4Reco/Tracking/options/runRefitFinal.py
  )
  set_test_env("run_RefitFinal")
  set_tests_properties("run_RefitFinal" PROPERTIES
    FIXTURES_REQUIRED "WrapperOutput_fixture"
    FIXTURES_SETUP RefitFinalOutput_fixture)

  add_test(NAME "compare_RefitFinal"
    COMMAND "${PROJECT_SOURCE_DIR}/k4Reco/ConformalTracking/scripts/compare-tracks.py" --gaudi-file output_refit_final.root --marlin-tracks SiTracks_Refitted --gaudi-tracks GaudiSiTracks_Refitted
  )
  set_test_env("compare_RefitFinal")
  set_tests_properties("compare_RefitFinal" PROPERTIES FIXTURES_REQUIRED "RefitFinalOutput_fixture")

  add_test(NAME "run_TruthTrackFinder"
    COMMAND k4run ${PROJECT_SOURCE_DIR}/k4Reco/Tracking/options/runTruthTrackFinder.py
  )
  set_test_env("run_TruthTrackFinder")
  set_tests_properties("run_TruthTrackFinder" PROPERTIES
    FIXTURES_REQUIRED "TruthTrackingOutput_fixture"
    FIXTURES_SETUP TruthTrackFinderOutput_fixture)

  add_test(NAME "compare_TruthTrackFinder"
    COMMAND "${PROJECT_SOURCE_DIR}/k4Reco/ConformalTracking/scripts/compare-tracks.py" --marlin-file output_truth_tracking_REC.edm4hep.root --gaudi-file output_truth_track_finder.root --marlin-tracks SiTracks --gaudi-tracks GaudiSiTracks
  )
  set_test_env("compare_TruthTrackFinder")
  set_tests_properties("compare_TruthTrackFinder" PROPERTIES FIXTURES_REQUIRED "TruthTrackFinderOutput_fixture")

endif()
