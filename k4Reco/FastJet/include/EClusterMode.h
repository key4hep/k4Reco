#ifndef EClusterMode_h
#define EClusterMode_h 1

#include <iostream>


// The enum, name and value of the enum for the Cluster Mode
enum EClusterMode {
	NONE = 0,
	FJ_exclusive_yCut = 1,		// exclusive clustering mode implemented in FastJet
	FJ_exclusive_nJets = 2,		// exclusive clustering mode implemented in FastJet
	FJ_inclusive = 4,		// inclusive "-"
	OWN_inclusiveIteration = 8	// use FJ inclusive Clustering, but iterate until we have the desired number of jets
};

namespace {

std::ostream& operator<<(std::ostream& out, EClusterMode& m) {
  switch (m) {
  case OWN_inclusiveIteration:
    out << "InclusiveIterativeNJets"; break;
  case FJ_inclusive:
    out << "Inclusive"; break;
  case FJ_exclusive_nJets:
    out << "ExclusiveNJets"; break;
  case FJ_exclusive_yCut:
    out << "ExclusiveYCut"; break;
  default:
    out << "unknown"; break;
  }
  return out;
}

}
#endif // EClusterMode_h