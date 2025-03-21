In this folder there are reimplementations of the other tracking processors that are used in the CLD reconstruction that are not ConformalTracking

* ClonesAndSplitTracksFinder

Originally ported from https://github.com/key4hep/k4Reco/blob/conformal-tracking/k4Reco/Tracking/include/ClonesAndSplitTracksFinder.h and https://github.com/key4hep/k4Reco/blob/conformal-tracking/k4Reco/Tracking/src/ClonesAndSplitTracksFinder.cpp. 

Only the part that runs when `mergeSplitTracks` is false (default value in the CLD reconstruction) has been ported. The rest can be ported if there on demand.

Additional validation may be required, since in the simulations that have been done there is not overlap between hits.

* RefitFinal

Originally ported from
https://github.com/iLCSoft/MarlinTrkProcessors/blob/master/source/Refitting/include/RefitFinal.h
and
https://github.com/iLCSoft/MarlinTrkProcessors/blob/master/source/Refitting/src/RefitFinal.cc

The resulting `edm4hep::TrackMCParticleLinkCollection` need to be validated.

* TruthTrackFinder

Originally ported from
https://github.com/iLCSoft/MarlinTrkProcessors/blob/master/source/Refitting/include/TruthTrackFinder.h
and
https://github.com/iLCSoft/MarlinTrkProcessors/blob/master/source/Refitting/src/TruthTrackFinder.cc 

Only the part that runs when `UseTruthInPrefit` is false (default value in the CLD reconstruction) has been ported. The rest can be ported if there on demand.

The member functions `getSubdetector` and `getLayer` have been implemented in `removeHitsSameLayer`.
