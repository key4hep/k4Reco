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
  COMMAND bash -c "git clone https://github.com/key4hep/CLDConfig.git --depth 1 && cp -r CLDConfig/CLDConfig/* . && sed -i -e '/Overlay.Overlay/d' -e '/HighLevelReco.RecoMCTruthLink/d' CLDReconstruction.py"
)
set_test_env("clone_CLDConfig")

add_test(NAME "run_ddsim"
  COMMAND bash -c "ddsim --compactFile $K4GEO/FCCee/CLD/compact/CLD_o2_v07/CLD_o2_v07.xml -G --numberOfEvents 3 --outputFile sim.edm4hep.root --gun.energy '10*GeV' --gun.particle 'mu-' --gun.distribution uniform --gun.multiplicity 10"
  )
set_test_env("run_ddsim")

add_test(NAME "run_marlin_wrapper"
  COMMAND k4run CLDReconstruction.py --inputFiles sim.edm4hep.root --trackingOnly
)
set_test_env("run_marlin_wrapper")
set_tests_properties("run_marlin_wrapper" PROPERTIES DEPENDS "clone_CLDConfig;run_ddsim")

add_test(NAME "run_marlin_wrapper with Truth Tracking"
  COMMAND k4run CLDReconstruction.py --inputFiles sim.edm4hep.root --trackingOnly --truthTracking --outputBasename output_truth_tracking
)
set_test_env("run_marlin_wrapper with Truth Tracking")
set_tests_properties("run_marlin_wrapper with Truth Tracking" PROPERTIES DEPENDS "clone_CLDConfig;run_ddsim")

if(BUILD_TRACKING)

  add_test(NAME "run_ConformalTracking"
    COMMAND k4run ${PROJECT_SOURCE_DIR}/k4Reco/ConformalTracking/options/runConformalTracking.py
  )
  set_test_env("run_ConformalTracking")
  set_tests_properties("run_ConformalTracking" PROPERTIES DEPENDS "run_marlin_wrapper")

  add_test(NAME "compare_ConformalTracking"
    COMMAND "${PROJECT_SOURCE_DIR}/k4Reco/ConformalTracking/scripts/compare-tracks.py"
  )
  set_test_env("compare_ConformalTracking")
  set_tests_properties("compare_ConformalTracking" PROPERTIES DEPENDS "run_ConformalTracking")


  add_test(NAME "run_ClonesAndSplitTracksFinder"
    COMMAND k4run ${PROJECT_SOURCE_DIR}/k4Reco/Tracking/options/runClonesAndSplitTracksFinder.py
  )
  set_test_env("run_ClonesAndSplitTracksFinder")
  set_tests_properties("run_ClonesAndSplitTracksFinder" PROPERTIES DEPENDS "run_marlin_wrapper")

  add_test(NAME "compare_ClonesAndSplitTracksFinder"
    COMMAND "${PROJECT_SOURCE_DIR}/k4Reco/ConformalTracking/scripts/compare-tracks.py" --gaudi-file output_clones_and_split_tracks_finder.root --marlin-tracks SiTracks --gaudi-tracks GaudiSiTracks
  )
  set_test_env("compare_ClonesAndSplitTracksFinder")
  set_tests_properties("compare_ClonesAndSplitTracksFinder" PROPERTIES DEPENDS "run_ConformalTracking")


  add_test(NAME "run_RefitFinal"
    COMMAND k4run ${PROJECT_SOURCE_DIR}/k4Reco/Tracking/options/runRefitFinal.py
  )
  set_test_env("run_RefitFinal")
  set_tests_properties("run_RefitFinal" PROPERTIES DEPENDS "run_marlin_wrapper")

  add_test(NAME "compare_RefitFinal"
    COMMAND "${PROJECT_SOURCE_DIR}/k4Reco/ConformalTracking/scripts/compare-tracks.py" --gaudi-file output_refit_final.root --marlin-tracks SiTracks_Refitted --gaudi-tracks GaudiSiTracks_Refitted
  )
  set_test_env("compare_RefitFinal")
  set_tests_properties("compare_RefitFinal" PROPERTIES DEPENDS "run_ConformalTracking")

  add_test(NAME "run_TruthTrackFinder"
    COMMAND k4run ${PROJECT_SOURCE_DIR}/k4Reco/Tracking/options/runTruthTrackFinder.py
  )
  set_test_env("run_TruthTrackFinder")
  set_tests_properties("run_TruthTrackFinder" PROPERTIES DEPENDS "run_marlin_wrapper with Truth Tracking")

  add_test(NAME "compare_TruthTrackFinder"
    COMMAND "${PROJECT_SOURCE_DIR}/k4Reco/ConformalTracking/scripts/compare-tracks.py" --marlin-file output_truth_tracking_REC.edm4hep.root --gaudi-file output_truth_track_finder.root --marlin-tracks SiTracks --gaudi-tracks GaudiSiTracks
  )
  set_test_env("compare_TruthTrackFinder")
  set_tests_properties("compare_TruthTrackFinder" PROPERTIES DEPENDS "run_ConformalTracking")

endif()
