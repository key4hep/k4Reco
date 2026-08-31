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

add_test(NAME "clone_CLDConfig"
  COMMAND bash -c "if [ ! -d CLDConfig ]; then git clone https://github.com/key4hep/CLDConfig.git --depth 1; fi && cp -r CLDConfig/CLDConfig/* . && sed -i -e '/Overlay.Overlay/d' -e '/HighLevelReco.RecoMCTruthLink/d' CLDReconstruction.py"
)
set_tests_properties("clone_CLDConfig" PROPERTIES FIXTURES_SETUP CLDConfig)

add_test(NAME "run_ddsim"
  COMMAND bash -c "ddsim --compactFile $K4GEO/FCCee/CLD/compact/CLD_o2_v07/CLD_o2_v07.xml -G --numberOfEvents 3 --outputFile sim.edm4hep.root --gun.energy '10*GeV' --gun.particle 'mu-' --gun.distribution uniform --gun.multiplicity 10"
  )
set_tests_properties("run_ddsim" PROPERTIES FIXTURES_SETUP SimOutput)

# If this is modified, make sure that this simulations produces clusters and reconstructed particles in the LumiCal
# as this doesn't always happen with the current setting (but does for the current seed)
add_test(NAME "run_ddsim_LumiCal"
  COMMAND bash -c "ddsim --compactFile $K4GEO/FCCee/CLD/compact/CLD_o2_v07/CLD_o2_v07.xml -G --numberOfEvents 3 --outputFile sim_lumical.edm4hep.root --gun.energy '50*GeV' --gun.particle 'gamma' --gun.distribution uniform --gun.multiplicity 10 --gun.thetaMin='50*mrad' --gun.thetaMax='100*mrad'"
  )
set_tests_properties("run_ddsim_LumiCal" PROPERTIES FIXTURES_SETUP SimLumiCalOutput)

add_test(NAME "run_marlin_wrapper_truth_tracking"
  COMMAND k4run CLDReconstruction.py --inputFiles sim.edm4hep.root --trackingOnly --truthTracking --outputBasename output_truth_tracking
)
set_tests_properties("run_marlin_wrapper_truth_tracking" PROPERTIES
  FIXTURES_REQUIRED "CLDConfig;SimOutput"
  FIXTURES_SETUP TruthTrackingOutput)

add_test(NAME "run_wrapper"
  COMMAND k4run CLDReconstruction.py --inputFiles sim.edm4hep.root --trackingOnly
)
set_tests_properties("run_wrapper" PROPERTIES
  FIXTURES_REQUIRED "CLDConfig;SimOutput"
  FIXTURES_SETUP WrapperOutput)

add_test(NAME "run_wrapper_LumiCal"
  COMMAND k4run CLDReconstruction.py --inputFiles sim_lumical.edm4hep.root --outputBasename output_lumical
)
set_tests_properties("run_wrapper_LumiCal" PROPERTIES
  FIXTURES_REQUIRED "CLDConfig;SimLumiCalOutput"
  FIXTURES_SETUP WrapperLumiCalOutput)

add_test(NAME "run_GaudiLumiCalClusterer"
  COMMAND k4run ${PROJECT_SOURCE_DIR}/k4Reco/GaudiLumiCalClusterer/options/runGaudiLumiCalClusterer.py --IOSvc.Input output_lumical_REC.edm4hep.root
)
set_tests_properties("run_GaudiLumiCalClusterer" PROPERTIES
  FIXTURES_REQUIRED "WrapperLumiCalOutput"
  FIXTURES_SETUP GaudiLumiCalOutput)

add_test(NAME "compare_LumiCalClusterer"
  COMMAND "${PROJECT_SOURCE_DIR}/k4Reco/GaudiLumiCalClusterer/scripts/compare-lumical.py"
)
set_tests_properties("compare_LumiCalClusterer" PROPERTIES FIXTURES_REQUIRED "WrapperLumiCalOutput;GaudiLumiCalOutput")

if(BUILD_TRACKING)

  add_test(NAME "run_ddsim_ClonesAndSplitTracksFinder"
    COMMAND bash -c "ddsim --compactFile $K4GEO/FCCee/CLD/compact/CLD_o2_v07/CLD_o2_v07.xml -G --numberOfEvents 3 --outputFile sim_clones_and_split_tracks_finder.edm4hep.root --gun.energy '10*GeV' --gun.particle 'mu-' --gun.distribution uniform --gun.multiplicity 30 --gun.thetaMin='89*deg' --gun.thetaMax='91*deg' --gun.phiMin='89*deg' --gun.phiMax='91*deg'"
  )
  set_tests_properties("run_ddsim_ClonesAndSplitTracksFinder" PROPERTIES FIXTURES_SETUP ClonesAndSplitTracksFinderSimOutput)

  add_test(NAME "run_marlin_wrapper_ClonesAndSplitTracksFinder"
    COMMAND k4run CLDReconstruction.py --inputFiles sim_clones_and_split_tracks_finder.edm4hep.root --trackingOnly --outputBasename output_clones_and_split_tracks_finder
  )
  set_tests_properties("run_marlin_wrapper_ClonesAndSplitTracksFinder" PROPERTIES
    FIXTURES_REQUIRED "CLDConfig;ClonesAndSplitTracksFinderSimOutput"
    FIXTURES_SETUP ClonesAndSplitTracksFinderWrapperOutput)

  add_test(NAME "run_ConformalTracking"
    COMMAND k4run ${PROJECT_SOURCE_DIR}/k4Reco/ConformalTracking/options/runConformalTracking.py
  )
  set_tests_properties("run_ConformalTracking" PROPERTIES
    FIXTURES_REQUIRED "WrapperOutput"
    FIXTURES_SETUP ConformalTrackingOutput)

  add_test(NAME "compare_ConformalTracking"
    COMMAND "${PROJECT_SOURCE_DIR}/k4Reco/ConformalTracking/scripts/compare-tracks.py"
  )
  set_tests_properties("compare_ConformalTracking" PROPERTIES FIXTURES_REQUIRED "ConformalTrackingOutput")

  add_test(NAME "run_ClonesAndSplitTracksFinder"
    COMMAND k4run ${PROJECT_SOURCE_DIR}/k4Reco/Tracking/options/runClonesAndSplitTracksFinder.py --IOSvc.Input output_clones_and_split_tracks_finder_REC.edm4hep.root
  )
  set_tests_properties("run_ClonesAndSplitTracksFinder" PROPERTIES
    FIXTURES_REQUIRED "ClonesAndSplitTracksFinderWrapperOutput"
    FIXTURES_SETUP ClonesAndSplitTracksFinderOutput)

  add_test(NAME "compare_ClonesAndSplitTracksFinder"
    COMMAND "${PROJECT_SOURCE_DIR}/k4Reco/ConformalTracking/scripts/compare-tracks.py" --marlin-file output_clones_and_split_tracks_finder_REC.edm4hep.root --gaudi-file output_clones_and_split_tracks_finder.root --marlin-tracks SiTracks --gaudi-tracks GaudiSiTracks
  )
  set_tests_properties("compare_ClonesAndSplitTracksFinder" PROPERTIES FIXTURES_REQUIRED "ClonesAndSplitTracksFinderOutput")

  add_test(NAME "run_RefitFinal"
    COMMAND k4run ${PROJECT_SOURCE_DIR}/k4Reco/Tracking/options/runRefitFinal.py
  )
  set_tests_properties("run_RefitFinal" PROPERTIES
    FIXTURES_REQUIRED "WrapperOutput"
    FIXTURES_SETUP RefitFinalOutput)

  add_test(NAME "compare_RefitFinal"
    COMMAND "${PROJECT_SOURCE_DIR}/k4Reco/ConformalTracking/scripts/compare-tracks.py" --gaudi-file output_refit_final.root --marlin-tracks SiTracks_Refitted --gaudi-tracks GaudiSiTracks_Refitted
  )
  set_tests_properties("compare_RefitFinal" PROPERTIES FIXTURES_REQUIRED "RefitFinalOutput")

  add_test(NAME "run_TruthTrackFinder"
    COMMAND k4run ${PROJECT_SOURCE_DIR}/k4Reco/Tracking/options/runTruthTrackFinder.py
  )
  set_tests_properties("run_TruthTrackFinder" PROPERTIES
    FIXTURES_REQUIRED "TruthTrackingOutput"
    FIXTURES_SETUP TruthTrackFinderOutput)

  add_test(NAME "compare_TruthTrackFinder"
    COMMAND "${PROJECT_SOURCE_DIR}/k4Reco/ConformalTracking/scripts/compare-tracks.py" --marlin-file output_truth_tracking_REC.edm4hep.root --gaudi-file output_truth_track_finder.root --marlin-tracks SiTracks --gaudi-tracks GaudiSiTracks
  )
  set_tests_properties("compare_TruthTrackFinder" PROPERTIES FIXTURES_REQUIRED "TruthTrackFinderOutput")

endif()

add_test( NAME "run_ddsim_ttbar"
  COMMAND bash -c "ddsim --compactFile $K4GEO/FCCee/CLD/compact/CLD_o2_v07/CLD_o2_v07.xml --numberOfEvents 3 --outputFile sim_ttbar.edm4hep.root --inputFile CLDConfig/test/yyxyev_000.stdhep"
)
set_tests_properties("run_ddsim_ttbar" PROPERTIES
  FIXTURES_REQUIRED "CLDConfig"
  FIXTURES_SETUP SimOutputTTbar
)

add_test(NAME "run_wrapper_ttbar"
  COMMAND k4run CLDReconstruction.py --inputFiles sim_ttbar.edm4hep.root --outputBasename output_ttbar
)
set_tests_properties("run_wrapper_ttbar" PROPERTIES
  FIXTURES_REQUIRED "CLDConfig;SimOutputTTbar"
  FIXTURES_SETUP WrapperTTbar
)

add_test(NAME "run_FastJet_ttbar"
  COMMAND k4run ${CMAKE_CURRENT_SOURCE_DIR}/FastJet/runFastJet.py --IOSvc.Input output_ttbar_REC.edm4hep.root --IOSvc.Output output_fastjet_ttbar.root
)
set_tests_properties("run_FastJet_ttbar" PROPERTIES
  FIXTURES_REQUIRED "WrapperTTbar"
  FIXTURES_SETUP GaudiFastJet
)

add_test(NAME "run_MarlinFastJet_ttbar"
  COMMAND k4run ${CMAKE_CURRENT_SOURCE_DIR}/FastJet/runMarlinFastJet.py --inputfile output_ttbar_REC.edm4hep.root --outputfile output_marlin_fastjet_ttbar.root
)
set_tests_properties("run_MarlinFastJet_ttbar" PROPERTIES
  FIXTURES_REQUIRED "WrapperTTbar"
  FIXTURES_SETUP MarlinFastJet
)

add_test(NAME "compare_FastJet"
  COMMAND "${PROJECT_SOURCE_DIR}/k4Reco/FastJet/scripts/compareJets.py" --marlin-file output_marlin_fastjet_ttbar.root --gaudi-file output_fastjet_ttbar.root --marlin-jets JetOut --gaudi-jets JetOut
)
set_tests_properties("compare_FastJet" PROPERTIES
  FIXTURES_REQUIRED "GaudiFastJet;MarlinFastJet"
  FIXTURES_SETUP GaudiMarlinFastJetComparison
)

get_filename_component(CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)

set(test_environment "\
LD_LIBRARY_PATH=\
${PROJECT_BINARY_DIR}:\
${PROJECT_BINARY_DIR}/${PROJECT_NAME}:\
${CMAKE_CURRENT_BINARY_DIR}:\
${CMAKE_CURRENT_BINARY_DIR}/genConfDir/${CURRENT_DIR_NAME}:\
$ENV{LD_LIBRARY_PATH};\
\
GAUDI_PLUGIN_PATH=
${PROJECT_BINARY_DIR}:\
${PROJECT_BINARY_DIR}/${PROJECT_NAME}:\
${CMAKE_CURRENT_BINARY_DIR}:\
${CMAKE_CURRENT_BINARY_DIR}/genConfDir/${CURRENT_DIR_NAME}:\
$ENV{GAUDI_PLUGIN_PATH}:\
$ENV{LD_LIBRARY_PATH};\
\
PYTHONPATH=\
${PROJECT_SOURCE_DIR}/python:\
${CMAKE_CURRENT_SOURCE_DIR}:\
${PROJECT_BINARY_DIR}/${PROJECT_NAME}/genConfDir:\
${CMAKE_CURRENT_BINARY_DIR}/genConfDir:\
$ENV{PYTHONPATH};\
"
)
get_property(test_names DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY TESTS)
set_tests_properties(${test_names} PROPERTIES ENVIRONMENT "${test_environment}")
