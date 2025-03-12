#ifndef K4RECO_GAUDIDDKALTEST_H
#define K4RECO_GAUDIDDKALTEST_H

#include <TVector3.h>

#include <kaltest/TKalDetCradle.h>

#include <cstdint>
#include <map>
#include <vector>

class DDKalDetector;
class DDVMeasLayer;
class THelicalTrack;
class TVKalDetector;
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
  GaudiDDKalTest(const GaudiDDKalTest&)                  = delete;
  GaudiDDKalTest const& operator=(const GaudiDDKalTest&) = delete;
  GaudiDDKalTest(GaudiDDKalTest&&)                       = delete;
  GaudiDDKalTest const& operator=(GaudiDDKalTest&&)      = delete;
  ~GaudiDDKalTest();

  //   /** Sets the specified option ( one of the constants defined in IMarlinTrkSystem::CFG )
  //     *  to the given value. Override here to re-configure E-loss and QMS
  //     *  after the initialization.
  //     */
  //   virtual void setOption(unsigned CFGOption, bool val);

  /** initialise track fitter system */
  void init();

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
  const DDVMeasLayer* findMeasLayer(const edm4hep::TrackerHitPlane& trkhit) const;

  const DDCylinderMeasLayer* getIPLayer() const { return m_ipLayer; }

private:
  bool                                    is_initialised = false;
  const DDCylinderMeasLayer*              m_ipLayer      = nullptr;
  TKalDetCradle                           m_det{};  // the detector cradle
  std::multimap<int, const DDVMeasLayer*> m_active_measurement_modules{};
  std::multimap<int, const DDVMeasLayer*> m_active_measurement_modules_by_layer{};
  std::vector<DDKalDetector*>             m_detectors{};

  const Gaudi::Algorithm* m_thisAlg;
};

#endif
