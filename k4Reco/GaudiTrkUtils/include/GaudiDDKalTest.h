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
#ifndef K4RECO_GAUDIDDKALTEST_H
#define K4RECO_GAUDIDDKALTEST_H

#include <DDSegmentation/BitFieldCoder.h>

#include <TVector3.h>

// Remove when these warnings are fixed in KalTest
// #pragma GCC diagnostic push
// #pragma GCC diagnostic ignored "-Wdeprecated-copy"
// #pragma GCC diagnostic ignored "-Woverloaded-virtual"
// #include <kaltest/TKalDetCradle.h>
// #pragma GCC diagnostic pop

#include <edm4hep/TrackerHit.h>

#include <cstdint>
#include <map>
#include <vector>

class DDKalDetector;
class DDVMeasLayer;
class THelicalTrack;
class TVKalDetector;
class TKalDetCradle;
namespace edm4hep {
class TrackerHitPlane;
}
namespace Gaudi {
class Algorithm;
}

class DDCylinderMeasLayer;

/** Interface to KaltTest Kalman fitter - instantiates and holds the detector geometry.
 */
class GaudiDDKalTest {
public:
  friend class GaudiDDKalTestTrack;

  //   // define some configuration constants
  //   static const bool FitBackward   = kIterBackward;
  //   static const bool FitForward    = kIterForward;
  //   static const bool OrderOutgoing = true;
  //   static const bool OrderIncoming = false;

  GaudiDDKalTest() = delete;
  GaudiDDKalTest(const Gaudi::Algorithm* algorithm);
  GaudiDDKalTest(const GaudiDDKalTest&) = delete;
  GaudiDDKalTest const& operator=(const GaudiDDKalTest&) = delete;
  GaudiDDKalTest(GaudiDDKalTest&&) = delete;
  GaudiDDKalTest const& operator=(GaudiDDKalTest&&) = delete;
  ~GaudiDDKalTest();

  //   /** Sets the specified option ( one of the constants defined in IMarlinTrkSystem::CFG )
  //     *  to the given value. Override here to re-configure E-loss and QMS
  //     *  after the initialization.
  //     */
  //   virtual void setOption(unsigned CFGOption, bool val);

  /** initialise track fitter system */
  void init();

  // Copy the encoder
  void setEncoder(const dd4hep::DDSegmentation::BitFieldCoder& encoder) { m_encoder = encoder; }

private:
  /** take multiple scattering into account during the fit */
  void includeMultipleScattering(bool on);

  /** take energy loss into account during the fit */
  void includeEnergyLoss(bool on);

  /** Store active measurement module IDs for a given TVKalDetector needed for navigation  */
  void storeActiveMeasurementModuleIDs(const TVKalDetector* detector);

  /** Store active measurement module IDs needed for navigation  */
  std::vector<const DDVMeasLayer*> getSensitiveMeasurementModules(const std::uint64_t detElementID) const;

  /** Store active measurement module IDs needed for navigation  */
  std::vector<const DDVMeasLayer*> getSensitiveMeasurementModulesForLayer(std::uint64_t layerID) const;

  // find the measurment layer for a given det element ID and point in space
  const DDVMeasLayer* findMeasLayer(const std::uint64_t detElementID, const TVector3& point) const;

  // get the last layer crossed by the helix when extrapolating from the present position to the pca to point
  const DDVMeasLayer* getLastMeasLayer(const THelicalTrack& helix, const TVector3& point) const;

  // find the measurement layer for a given hit
  const DDVMeasLayer* findMeasLayer(const edm4hep::TrackerHit& trkhit) const;

  const DDCylinderMeasLayer* getIPLayer() const { return m_ipLayer; }

private:
  bool is_initialised = false;
  const DDCylinderMeasLayer* m_ipLayer = nullptr;
  std::unique_ptr<TKalDetCradle> m_det; // the detector cradle, pointer to be able to forward declare
  std::multimap<int, const DDVMeasLayer*> m_active_measurement_modules{};
  std::multimap<int, const DDVMeasLayer*> m_active_measurement_modules_by_layer{};
  std::vector<DDKalDetector*> m_detectors{};

  // Originally not present in MarlinDDKalTest, this is needed to be able to use
  // logging from Gaudi
  const Gaudi::Algorithm* m_thisAlg;

  dd4hep::DDSegmentation::BitFieldCoder m_encoder;
};

#endif
