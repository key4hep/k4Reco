Reimplementation of the ConformalTracking processor from iLCSoft, as well as
utilities from MarlinTrk and MarlinTrkUtils.

The following files have been reimplemented from `MarlinTrk`:
- MarlinDDKalTest: Now called GaudiDDKalTest. Anything related to `MarlinTrk`
  has been removed. The rest is very similar, changing LCIO Tracker hits to
  EDM4hep tracker hits.

  Remove MarlinTrk::IMarlinTrack* createTrack()
         std::string MarlinTrk::name()

- MarlinDDKalTestTrack: now called GaudiDDKalTestTrack. 

  Remove void setMass(double mass)
         double getMass() 
         int initialise( bool fitDirection );
         int smooth() ;
         int testChi2Increment( EVENT::TrackerHit* hit, double& chi2increment ) ;
         int getTrackState( IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) ;
         int propagate( const edm4hep::Vector3d& point, edm4hep::TrackState& ts, double& chi2, int& ndf ) ;
         int propagateToLayer( int layerID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID, int mode=modeClosest )  ;
         int propagateToDetElement( int detElementID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode=modeClosest ) ;
         int propagateToDetElement( int detEementID, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode=modeClosest ) ;
         int propagateToDetElement( int detEementID, const TKalTrackSite& site, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode=modeClosest ) ;
         int extrapolate( const Vector3D& point, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) ;
         int extrapolate( const Vector3D& point, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) ;
         int extrapolate( const Vector3D& point, const TKalTrackSite& site, IMPL::TrackStateImpl& ts, double& chi2, int& ndf ) ;
         int extrapolateToLayer( int layerID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID, int mode=modeClosest )  ;
         int extrapolateToLayer( int layerID, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID, int mode=modeClosest )  ;
         int extrapolateToLayer( int layerID, const TKalTrackSite& site, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int& detElementID, int mode=modeClosest )  ;
         int extrapolateToDetElement( int detElementID, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode=modeClosest ) ;
         int extrapolateToDetElement( int detEementID, EVENT::TrackerHit* hit, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode=modeClosest ) ;
         int extrapolateToDetElement( int detEementID, const TKalTrackSite& site, IMPL::TrackStateImpl& ts, double& chi2, int& ndf, int mode=modeClosest ) ;
         std::string toString() ;
         int intersectionWithLayer( int layerID, Vector3D& point, int& detElementID, int mode=modeClosest )  ;
         int intersectionWithLayer( int layerID, EVENT::TrackerHit* hit, Vector3D& point, int& detElementID, int mode=modeClosest )  ;
         int intersectionWithDetElement( int detElementID, Vector3D& point, int mode=modeClosest )  ;
         int intersectionWithDetElement( int detElementID, EVENT::TrackerHit* hit, Vector3D& point, int mode=modeClosest )  ;
         int intersectionWithDetElement( int detElementID, const TKalTrackSite& site, Vector3D& point, const DDVMeasLayer*& ml, int mode=modeClosest ) ;


- KDTrack
  Remove m_kalmanTrack that was never being set (and the corresponding code in ConformalTracking.cpp)
  Remove `m_nPoints` and use `m_clusters.size()` instead of keeping an additional counter

- HelixTrack
  Remove HelixTrack(const double* position, const double* p, double charge, double Bz);
