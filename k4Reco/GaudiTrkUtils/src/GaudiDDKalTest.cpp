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
#include "GaudiDDKalTest.h"

#include <DD4hep/DetElement.h>
#include <DD4hep/Detector.h>
#include <DDRec/Material.h>
#include <DDRec/SurfaceManager.h>

#include <GaudiKernel/MsgStream.h>
#include <kaltest/TKalDetCradle.h>
#include <kaltest/TVKalDetector.h>
#include <kaltest/TVSurface.h>

#include <DDKalTest/DDCylinderMeasLayer.h>
#include <DDKalTest/DDKalDetector.h>
#include <DDKalTest/DDVMeasLayer.h>

#include <edm4hep/TrackerHitPlane.h>

#include <Gaudi/Algorithm.h>

#include <TVector3.h>

#include <iostream>
#include <map>
#include <stdexcept>

GaudiDDKalTest::GaudiDDKalTest(const Gaudi::Algorithm* algorithm) : m_thisAlg(algorithm) {
  m_det.reset(new TKalDetCradle());
  m_det->SetOwner(true); // takes care of deleting subdetector in the end ...

  // this->registerOptions();
}

GaudiDDKalTest::~GaudiDDKalTest() {
  for (auto* ddKalDet : m_detectors) {
    delete ddKalDet;
  }
}

// void GaudiDDKalTest::setOption(unsigned CFGOption, bool val) {
//   IMarlinTrkSystem::setOption(CFGOption, val);

//   switch (CFGOption) {
//     case IMarlinTrkSystem::CFG::useQMS:
//       this->includeMultipleScattering(val);
//       break;
//     case IMarlinTrkSystem::CFG::usedEdx:
//       this->includeEnergyLoss(val);
//       break;
//       // IMarlinTrkSystem::CFG::useSmoothing handled directly in MarlinDDKalTestTrack
//   }
// }

void GaudiDDKalTest::init() {
  // TODO: Don't hardcode the options
  this->includeMultipleScattering(true);
  this->includeEnergyLoss(true);
  // this->sm

  // std::cout << " -------------------------------------------------------------------------------- " << endmsg;
  // std::cout << "  GaudiDDKalTest::init() called with the following options :                     " << endmsg;
  // std::cout << this->getOptions();
  // std::cout << " -------------------------------------------------------------------------------- " << endmsg;

  if (is_initialised) {
    m_thisAlg->warning() << "  GaudiDDKalTest::init()  - already initialized - only options are set .. " << endmsg;
    return;
  }

  m_thisAlg->info() << "GaudiDDKalTest::init()  - initializing  " << endmsg;

  dd4hep::Detector& lcdd = dd4hep::Detector::getInstance();

  double minS = 1.e99;
  DDCylinderMeasLayer* ipLayer = nullptr;

  // for the tracking we get all tracking detectors and all passive detectors (beam pipe,...)

  std::vector<dd4hep::DetElement> detectors = lcdd.detectors("tracker");
  const std::vector<dd4hep::DetElement>& passiveDets = lcdd.detectors("passive");
  const std::vector<dd4hep::DetElement>& calos = lcdd.detectors("calorimeter");

  detectors.reserve(detectors.size() + passiveDets.size() + calos.size());

  std::copy(passiveDets.begin(), passiveDets.end(), std::back_inserter(detectors));

  for (const auto& det : calos) {
    std::string name = det.name();
    std::transform(name.begin(), name.end(), name.begin(), ::tolower);
    if (name.find("ecal") != std::string::npos) {
      detectors.push_back(det);
    }
  }

  for (const auto& det : detectors) {
    m_thisAlg->debug() << "  GaudiDDKalTest::init() - creating DDKalDetector for : " << det.name() << endmsg;

    m_detectors.push_back(new DDKalDetector(det));
    DDKalDetector* kalDet = m_detectors.back();

    this->storeActiveMeasurementModuleIDs(kalDet);

    m_det->Install(*kalDet);

    Int_t nLayers = kalDet->GetEntriesFast();

    // --- keep the cylinder meas layer with smallest sorting policy (radius) as ipLayer
    // fixme: this should be implemented in a more explicit way ...
    for (int i = 0; i < nLayers; ++i) {
      const TVSurface* tvs = static_cast<const TVSurface*>(kalDet->At(i));

      double s = tvs->GetSortingPolicy();
      if (s < minS && dynamic_cast<DDCylinderMeasLayer*>(kalDet->At(i))) {
        minS = s;
        ipLayer = dynamic_cast<DDCylinderMeasLayer*>(kalDet->At(i));
      }
    }
  }

  if (ipLayer) {
    m_ipLayer = ipLayer;

    m_thisAlg->debug() << " GaudiDDKalTest: install IP layer at radius : " << minS << endmsg;
  }
  //-------------------------------------------------------------------------------

  m_det->Close(); // close the cradle
  // done in Close()    m_det->Sort() ;           // sort meas. layers from inside to outside

  m_thisAlg->debug() << "  GaudiDDKalTest - number of layers = " << m_det->GetEntriesFast() << endmsg;

  // if (streamlog_level(DEBUG)) {
  //   lcio::BitField64 bf(UTIL::LCTrackerCellID::encoding_string());

  //   for (unsigned i = 0, N = m_det->GetEntriesFast(); i < N; ++i) {
  //     DDVMeasLayer* ml = dynamic_cast<DDVMeasLayer*>(m_det->At(i));

  //     bf.setValue(ml->getLayerID());

  //     TVSurface* s = dynamic_cast<TVSurface*>(m_det->At(i));

  //     streamlog_out(DEBUG) << " *** meas. layer : " << bf.valueString() << "  sorting: " << s->GetSortingPolicy()
  //                          << endmsg;
  //   }
  // }

  is_initialised = true;
}

void GaudiDDKalTest::includeMultipleScattering(bool msOn) {
  m_thisAlg->debug() << "  **** GaudiDDKalTest::includeMultipleScattering( " << msOn << " ) called " << endmsg;

  if (msOn) {
    m_det->SwitchOnMS();
  } else {
    m_det->SwitchOffMS();
  }
}

void GaudiDDKalTest::includeEnergyLoss(bool energyLossOn) {
  m_thisAlg->debug() << "  **** GaudiDDKalTest::includeEnergyLoss( " << energyLossOn << " ) called " << endmsg;

  if (energyLossOn) {
    m_det->SwitchOnDEDX();
  } else {
    m_det->SwitchOffDEDX();
  }
}

std::vector<const DDVMeasLayer*> GaudiDDKalTest::getSensitiveMeasurementModules(const std::uint64_t moduleID) const {
  std::vector<const DDVMeasLayer*> measmodules;

  auto ii = this->m_active_measurement_modules.equal_range(moduleID); // set the first and last entry in ii;

  for (auto& it = ii.first; it != ii.second; ++it) {
    //      std::cout<<"Key = "<<it->first<<"    Value = "<<it->second << endmsg ;
    measmodules.push_back(it->second);
  }
  return measmodules;
}

std::vector<const DDVMeasLayer*> GaudiDDKalTest::getSensitiveMeasurementModulesForLayer(std::uint64_t layerID) const {
  std::vector<const DDVMeasLayer*> measmodules;

  m_thisAlg->debug() << "GaudiDDKalTest::getSensitiveMeasurementModulesForLayer: layerID = " << layerID << endmsg;

  //  for(it = m_active_measurement_modules_by_layer.begin(); it != m_active_measurement_modules_by_layer.end(); ++it) {
  //    streamlog_out( DEBUG0 ) << "Key = "<< ttdecodeILD(it->first) <<"    Value = "<<it->second << endmsg ;
  //  }

  // set the module and sensor bit ranges to zero as these are not used in the map
  // TODO: Pass the encoding string from Gaudi
  // TODO: Check if correct
  m_encoder.set(layerID, "module", 0);
  m_encoder.set(layerID, "sensor", 0);
  layerID = m_encoder.lowWord(layerID);

  auto ii = this->m_active_measurement_modules_by_layer.equal_range(layerID); // set the first and last entry in ii;

  for (auto& it = ii.first; it != ii.second; ++it) {
    //    streamlog_out( DEBUG0 ) <<"Key = "<< it->first <<"    Value = "<<it->second << endmsg ;
    measmodules.push_back(it->second);
  }
  return measmodules;
}

void GaudiDDKalTest::storeActiveMeasurementModuleIDs(const TVKalDetector* detector) {
  Int_t nLayers = detector->GetEntriesFast();

  for (int i = 0; i < nLayers; ++i) {
    const auto* ml = dynamic_cast<const DDVMeasLayer*>(detector->At(i));

    if (!ml) {
      std::stringstream errorMsg;
      errorMsg << "GaudiDDKalTest::storeActiveMeasurementLayerIDs dynamic_cast to DDVMeasLayer* failed ";
      throw std::runtime_error(errorMsg.str());
    }

    if (ml->IsActive()) {
      // then get all the sensitive element id's assosiated with this DDVMeasLayer and store them in the map
      auto it = ml->getCellIDs().begin();

      while (it != ml->getCellIDs().end()) {
        int sensitive_element_id = *it;
        this->m_active_measurement_modules.insert(std::pair<int, const DDVMeasLayer*>(sensitive_element_id, ml));
        ++it;
      }

      int subdet_layer_id = ml->getLayerID();

      this->m_active_measurement_modules_by_layer.insert(std::pair<int, const DDVMeasLayer*>(subdet_layer_id, ml));

      // std::cout << "GaudiDDKalTest::storeActiveMeasurementLayerIDs added active layer with "
      //                       << " LayerID = " << subdet_layer_id << " and DetElementIDs  ";

      // for (it = ml->getCellIDs().begin(); it != ml->getCellIDs().end(); ++it) {
      //   std::cout << " : " << *it;
      // }

      // std::cout << endmsg;
    }
  }
}

const DDVMeasLayer* GaudiDDKalTest::getLastMeasLayer(const THelicalTrack& hel, TVector3 const& point) const {
  THelicalTrack helix = hel;

  double deflection_to_point = 0;
  helix.MoveTo(point, deflection_to_point, nullptr, nullptr);

  bool isfwd = ((helix.GetKappa() > 0 && deflection_to_point < 0) || (helix.GetKappa() <= 0 && deflection_to_point > 0))
                   ? true
                   : false;

  const int mode = isfwd ? -1 : +1;

  //  streamlog_out( DEBUG4 ) << "  GaudiDDKalTest - getLastMeasLayer deflection to point = " << deflection_to_point <<
  //  " kappa = " << helix.GetKappa()  << "  mode = " << mode << endmsg ; streamlog_out( DEBUG4 ) << " Point to move
  //  to:" << endmsg; point.Print();

  const TVSurface* ml_retval = nullptr;
  double min_deflection = DBL_MAX;

  for (int i = 0; i < m_det->GetEntriesFast(); ++i) {
    const auto* sfp = static_cast<const TVSurface*>(m_det->At(i)); // surface at destination

    double defection_angle = 0;
    TVector3 crossing_point;
    const int does_cross = sfp->CalcXingPointWith(helix, crossing_point, defection_angle, mode);

    if (does_cross) {
      const double deflection = std::abs(deflection_to_point - defection_angle);

      if (deflection < min_deflection) {
        //      streamlog_out( DEBUG4 ) << "  GaudiDDKalTest - crossing found for suface = " << ml.GetMLName()
        //                              << endmsg
        //                              << "  min_deflection = " << min_deflection
        //                              << "  deflection = " << deflection
        //                              << "  deflection angle = " << defection_angle
        //                              << endmsg
        //                              << " x = " << crossing_point.X()
        //                              << " y = " << crossing_point.Y()
        //                              << " z = " << crossing_point.Z()
        //                              << " r = " << crossing_point.Perp()
        //                              << endmsg ;

        min_deflection = deflection;
        ml_retval = sfp;
      }
    }
  }

  return dynamic_cast<const DDVMeasLayer*>(ml_retval);
}

const DDVMeasLayer* GaudiDDKalTest::findMeasLayer(const edm4hep::TrackerHitPlane& trkhit) const {
  const TVector3 hit_pos(trkhit.getPosition()[0], trkhit.getPosition()[1], trkhit.getPosition()[2]);

  return this->findMeasLayer(trkhit.getCellID(), hit_pos);
}

const DDVMeasLayer* GaudiDDKalTest::findMeasLayer(const std::uint64_t detElementID, const TVector3& point) const {
  const DDVMeasLayer* ml = nullptr;

  // search for the list of measurement layers associated with this CellID
  const auto meas_modules = this->getSensitiveMeasurementModules(detElementID);

  if (meas_modules.size() == 0) { // no measurement layers found

    throw std::runtime_error("GaudiDDKalTest::findMeasLayer module id " + std::to_string(detElementID) +
                             " not found in active measurement modules");

  } else if (meas_modules.size() == 1) { // one to one mapping

    ml = meas_modules[0];

  } else { // layer has been split

    bool surf_found = false;

    // loop over the measurement layers associated with this CellID and find the correct one using the position of the
    // hit
    for (const auto& module : meas_modules) {
      const auto* surf = dynamic_cast<const TVSurface*>(module);

      if (!surf) {
        throw std::runtime_error("GaudiDDKalTest::findMeasLayer dynamic_cast failed for surface type: moduleID = " +
                                 std::to_string(detElementID));
      }

      const bool hit_on_surface = surf->IsOnSurface(point);

      if (!surf_found && hit_on_surface) {
        ml = module;
        surf_found = true;

      } else if (surf_found && hit_on_surface) { // only one surface should be found, if not throw

        std::stringstream errorMsg;
        errorMsg << "GaudiDDKalTest::findMeasLayer point found to be on two surfaces: moduleID = " << detElementID;
        throw std::runtime_error(errorMsg.str());
      }
    }
    if (!surf_found) { // print out debug info
      m_thisAlg->debug() << "GaudiDDKalTest::findMeasLayer point not found to be on any surface matching moduleID = "
                         << detElementID << ": x = " << point.x() << " y = " << point.y() << " z = " << point.z()
                         << endmsg;
    } else {
      m_thisAlg->debug() << "GaudiDDKalTest::findMeasLayer point found to be on surface matching moduleID = "
                         << detElementID << ": x = " << point.x() << " y = " << point.y() << " z = " << point.z()
                         << endmsg;
    }
  }

  return ml;
}
