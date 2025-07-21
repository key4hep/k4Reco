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

#include "Global.h"
#include "LCCluster.h"
#include "LumiCalClusterer.h"
#include "LumiCalHit.h"

#include <TH1.h>

#include <cmath>
#include <map>
#include <string>
#include <utility>
#include <vector>

void LumiCalClustererClass::energyCorrections(MapIntVInt& superClusterIdToCellId,
                                              MapIntVDouble& superClusterIdToCellEngy, MapIntLCCluster& superClusterCM,
                                              const MapIntCalHit& calHitsCellIdGlobal) {
  std::map<int, std::vector<int>>::iterator superClusterIdToCellIdIterator;

  std::vector<int> cellIdV;
  std::vector<double> cellEngyV;

  std::map<int, double> engyCorrectionCellId, engyCorrectionEngy;
  std::map<int, double>::iterator engyCorrectionCellIdIterator;

  double distanceNow, distanceAB, distanceAC, distanceBC, engyNow;
  double correctionFactor;
  int cellIdHit;
  int superClusterId;

  int maxEngySuperClusterId(0);
  double maxEngyCluster, engyClusterNow;

  std::string hisName;
  int numBins1;
  double hisRange1[2];

  hisRange1[0] = -m_moliereRadius * 5;
  hisRange1[1] = m_moliereRadius * 5;
  numBins1 = int(hisRange1[1] - hisRange1[0]); // 1 mm bin width

  // internal temporary histograms that are not to be written out
  hisName = "leftSideLargeHis";
  TH1F leftSideLargeHisH(hisName.c_str(), hisName.c_str(), numBins1, hisRange1[0], hisRange1[1]);

  hisName = "rightSideSmallHis";
  TH1F rightSideSmallHisH(hisName.c_str(), hisName.c_str(), numBins1, hisRange1[0], hisRange1[1]);

  hisName = "correctionRatio";
  TH1F correctionRatioH(hisName.c_str(), hisName.c_str(), numBins1, hisRange1[0], hisRange1[1]);

  /* --------------------------------------------------------------------------
     Find the reconstructed clusters with the most energy
     -------------------------------------------------------------------------- */
  maxEngyCluster = 0.;
  superClusterIdToCellIdIterator = superClusterIdToCellId.begin();
  for (size_t superClusterNow = 0; superClusterNow < superClusterIdToCellId.size();
       superClusterNow++, superClusterIdToCellIdIterator++) {
    superClusterId = (int)(*superClusterIdToCellIdIterator).first; // Id of cluster

    engyClusterNow = superClusterCM[superClusterId].getE();
    if (maxEngyCluster < engyClusterNow) {
      maxEngyCluster = engyClusterNow;
      maxEngySuperClusterId = superClusterId;
    }
  }

  /* --------------------------------------------------------------------------
     fill correction histograms with large cluster projected hits at
     negative distanceAC, and small cluster hits at positive distanceAB
     (assumeing that there is no contribution from the small cluster
     at distanceAC<0 and mixing of the two clusters at distanceAC>0)
     -------------------------------------------------------------------------- */
  superClusterIdToCellIdIterator = superClusterIdToCellId.begin();
  double pos1[2] = {0.0, 0.0}, pos2[2] = {0.0, 0.0};
  for (size_t superClusterNow = 0; superClusterNow < superClusterIdToCellId.size();
       superClusterNow++, superClusterIdToCellIdIterator++) {
    superClusterId = (int)(*superClusterIdToCellIdIterator).first; // Id of cluster

    if (superClusterId == maxEngySuperClusterId) {
      pos1[0] = superClusterCM[superClusterId].getX();
      pos1[1] = superClusterCM[superClusterId].getY();
    }
    if (superClusterId != maxEngySuperClusterId) {
      pos2[0] = superClusterCM[superClusterId].getX();
      pos2[1] = superClusterCM[superClusterId].getY();
    }
  }

  distanceAB = std::hypot(pos1[0] - pos2[0], pos1[1] - pos2[1]);

  superClusterIdToCellIdIterator = superClusterIdToCellId.begin();
  for (size_t superClusterNow = 0; superClusterNow < superClusterIdToCellId.size();
       superClusterNow++, superClusterIdToCellIdIterator++) {
    superClusterId = (*superClusterIdToCellIdIterator).first; // Id of cluster

    for (size_t hitNow = 0; hitNow < superClusterIdToCellId[superClusterId].size(); hitNow++) {
      cellIdHit = superClusterIdToCellId[superClusterId][hitNow];

      const auto& thisHit = calHitsCellIdGlobal.at(cellIdHit);
      double pos3[2] = {thisHit->getPosition()[0], thisHit->getPosition()[1]};

      engyNow = superClusterIdToCellEngy[superClusterId][hitNow];

      distanceNow = (pos3[0] - pos1[0]) * (pos2[0] - pos1[0]) + (pos3[1] - pos1[1]) * (pos2[1] - pos1[1]);
      distanceNow /= distanceAB * distanceAB;

      // distanceNow point of the tangent from point pos3[] to the line connecting the two CMs
      pos3[0] = pos1[0] + distanceNow * (pos2[0] - pos1[0]);
      pos3[1] = pos1[1] + distanceNow * (pos2[1] - pos1[1]);

      distanceAC = std::hypot(pos1[0] - pos3[0], pos1[1] - pos3[1]);
      distanceBC = std::hypot(pos2[0] - pos3[0], pos2[1] - pos3[1]);

      // distanceNow == 0 if point A is in between points C and B
      distanceNow = fabs(distanceBC - distanceAC - distanceAB) / distanceBC;

      if (distanceNow < 1e-7 && superClusterId == maxEngySuperClusterId)
        leftSideLargeHisH.Fill(distanceAC, engyNow);

      if (distanceNow > 1e-7 && superClusterId != maxEngySuperClusterId)
        rightSideSmallHisH.Fill(distanceAC, engyNow);
    }
  }

  /* --------------------------------------------------------------------------
     fill the correctionRatioH histogram with correction ratios
     -------------------------------------------------------------------------- */
  int nBinsX = rightSideSmallHisH.GetNbinsX();
  for (int binNowX = 0; binNowX < nBinsX; binNowX++) {

    distanceNow = leftSideLargeHisH.GetBinCenter(binNowX);
    double engyLargeNow = leftSideLargeHisH.GetBinContent(binNowX);
    double engySmallNow = rightSideSmallHisH.GetBinContent(binNowX);
    double deltaEngy = engySmallNow - engyLargeNow;
    double engyRatio = deltaEngy / engySmallNow;
    /* (BP) seems engyNow is for nothing ?
    if(deltaEngy > 0)	engyNow = deltaEngy;
    else	        engyNow = engySmallNow;

    engyNow = GlobalMethodsClass::SignalGevConversion(GlobalMethodsClass::Signal_to_GeV, engyNow);
    */
    if (engyRatio > 0)
      correctionRatioH.Fill(distanceNow, engyRatio);
  }

  /* --------------------------------------------------------------------------
     decrease the energy of hits from the small cluster and store the changes
     in order to increase the energy of the large cluster later on
     -------------------------------------------------------------------------- */
  superClusterIdToCellIdIterator = superClusterIdToCellId.begin();
  for (size_t superClusterNow = 0; superClusterNow < superClusterIdToCellId.size();
       superClusterNow++, superClusterIdToCellIdIterator++) {
    superClusterId = (int)(*superClusterIdToCellIdIterator).first; // Id of cluster

    for (size_t hitNow = 0; hitNow < superClusterIdToCellId[superClusterId].size(); hitNow++) {
      cellIdHit = superClusterIdToCellId[superClusterId][hitNow];

      const auto& thisHit = calHitsCellIdGlobal.at(cellIdHit);
      double pos3[2] = {thisHit->getPosition()[0], thisHit->getPosition()[1]};

      engyNow = superClusterIdToCellEngy[superClusterId][hitNow];

      distanceNow = (pos3[0] - pos1[0]) * (pos2[0] - pos1[0]) + (pos3[1] - pos1[1]) * (pos2[1] - pos1[1]);
      distanceNow /= distanceAB * distanceAB;

      // distanceNow point of the tangent from point pos3[] to the line connecting the two CMs
      pos3[0] = pos1[0] + distanceNow * (pos2[0] - pos1[0]);
      pos3[1] = pos1[1] + distanceNow * (pos2[1] - pos1[1]);

      distanceAC = std::hypot(pos1[0] - pos3[0], pos1[1] - pos3[1]);
      distanceBC = std::hypot(pos2[0] - pos3[0], pos2[1] - pos3[1]);

      // distanceNow == 0 if point A is in between points C and B
      distanceNow = fabs(distanceBC - distanceAC - distanceAB) / distanceBC;

      if (distanceNow > 1e-7) {
        int binNow = correctionRatioH.Fill(distanceAC, 0);
        correctionFactor = correctionRatioH.GetBinContent(binNow);

        if (correctionFactor > 0 && superClusterId != maxEngySuperClusterId) {
          engyNow = engyNow * correctionFactor;

          // store the cell id that is changed for modigying the lareg cluster later
          engyCorrectionCellId[cellIdHit] = correctionFactor;
          engyCorrectionEngy[cellIdHit] = engyNow * (1 / correctionFactor - 1);

          // modify the energy of the small cluster
          superClusterIdToCellEngy[superClusterId][hitNow] = engyNow;
        }
      }
    }
  }

  /* --------------------------------------------------------------------------
     increase the energy of the large cluster according to what was decreased
     from the small cluster
     -------------------------------------------------------------------------- */
  superClusterIdToCellIdIterator = superClusterIdToCellId.begin();
  for (size_t superClusterNow = 0; superClusterNow < superClusterIdToCellId.size();
       superClusterNow++, superClusterIdToCellIdIterator++) {
    superClusterId = (int)(*superClusterIdToCellIdIterator).first; // Id of cluster

    if (superClusterId != maxEngySuperClusterId)
      continue;

    // go over all exisiting cells in the cluster and increase their energy
    for (size_t hitNow = 0; hitNow < superClusterIdToCellId[superClusterId].size(); hitNow++) {
      cellIdHit = superClusterIdToCellId[superClusterId][hitNow];

      if (engyCorrectionCellId[cellIdHit] > 0) {
        correctionFactor = 2 - engyCorrectionCellId[cellIdHit];
        superClusterIdToCellEngy[superClusterId][hitNow] *= correctionFactor;
        engyCorrectionCellId[cellIdHit] = 0.;
      }
    }

    // go over all the new cells and add them to the cluster
    engyCorrectionCellIdIterator = engyCorrectionCellId.begin();
    for (size_t hitNow = 0; hitNow < engyCorrectionCellId.size(); hitNow++, engyCorrectionCellIdIterator++) {
      cellIdHit = (int)(*engyCorrectionCellIdIterator).first;

      if (engyCorrectionCellId[cellIdHit] > 0) {
        engyNow = engyCorrectionEngy[cellIdHit];

        superClusterIdToCellId[superClusterId].push_back(cellIdHit);
        superClusterIdToCellEngy[superClusterId].push_back(engyNow);
      }
    }
  }

  /* --------------------------------------------------------------------------
     verbosity
     -------------------------------------------------------------------------- */
#if _MCPARTICLE_CLUSTER_DEBUG == 1
  cout << endl << coutBlue << "Original Super Clusters:  " << coutDefault << endl;

  superClusterIdToCellIdIterator = superClusterIdToCellId.begin();
  numSuperClusters = superClusterIdToCellId.size();
  for (int superClusterNow = 0; superClusterNow < numSuperClusters;
       superClusterNow++, superClusterIdToCellIdIterator++) {
    superClusterId = (int)(*superClusterIdToCellIdIterator).first; // Id of cluster

    engyNow = superClusterCM[superClusterId][0];

    cout << "\t Id " << superClusterId << "  \t  energy(signal,GeV) = ( " << engyNow << " , "
         << engySignalGeV(engyNow, GlobalMethodsClass::Signal_to_GeV) << " )  \t pos(x,y) =  ( "
         << superClusterCM[superClusterId][1] << " , " << superClusterCM[superClusterId][2] << " )" << endl;
  }
#endif

  /* --------------------------------------------------------------------------
     compute the total energy and center of mass of the superClusters
     according to the corrected energy vectors
     -------------------------------------------------------------------------- */
  superClusterCM.clear();
  superClusterIdToCellIdIterator = superClusterIdToCellId.begin();
  for (size_t superClusterNow = 0; superClusterNow < superClusterIdToCellId.size();
       superClusterNow++, superClusterIdToCellIdIterator++) {
    superClusterId = (int)(*superClusterIdToCellIdIterator).first; // Id of cluster

    cellIdV = superClusterIdToCellId[superClusterId];
    cellEngyV = superClusterIdToCellEngy[superClusterId];

    // initialize the energy/position std::vector for new clusters only
    superClusterCM[superClusterId] = LCCluster();

    // calculate/update the energy/position of the CM
    calculateEngyPosCM_EngyV(cellIdV, cellEngyV, calHitsCellIdGlobal, superClusterCM, superClusterId, m_methodCM);
  }
  cellIdV.clear();
  cellEngyV.clear();

  /* --------------------------------------------------------------------------
     verbosity
     -------------------------------------------------------------------------- */
#if _MCPARTICLE_CLUSTER_DEBUG == 1
  cout << endl << coutBlue << "Fixed Super Clusters:  (distanceAB =  " << distanceAB << ")" << coutDefault << endl;

  superClusterIdToCellIdIterator = superClusterIdToCellId.begin();
  numSuperClusters = superClusterIdToCellId.size();
  for (int superClusterNow = 0; superClusterNow < numSuperClusters;
       superClusterNow++, superClusterIdToCellIdIterator++) {
    superClusterId = (int)(*superClusterIdToCellIdIterator).first; // Id of cluster

    engyNow = superClusterCM[superClusterId][0];

    cout << "\t Id " << superClusterId << "  \t  energy(signal,GeV) = ( " << engyNow << " , "
         << engySignalGeV(engyNow, GlobalMethodsClass::Signal_to_GeV) << " )  \t pos(x,y) =  ( "
         << superClusterCM[superClusterId][1] << " , " << superClusterCM[superClusterId][2] << " )" << endl;
  }
#endif

  return;
}
