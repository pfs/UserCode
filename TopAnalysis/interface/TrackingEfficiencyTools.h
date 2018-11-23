#ifndef _correction_tools_h_
#define _correction_tools_h_

#include <vector>
#include <string>

#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"
#include "TGraphAsymmErrors.h"

//apply tracking efficiency
std::map<TString, std::map<TString, std::vector<double> > > getTrackingEfficiencyMap(TString era);
void applyTrackingEfficiencySF(MiniEvent_t &ev, double sf, double minEta = -2.4, double maxEta = 2.4);
void applyEtaDepTrackingEfficiencySF(MiniEvent_t &ev, std::vector<double> sfs, std::vector<double> etaBins);

#endif
