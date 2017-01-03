#ifndef _TopCharmedMesonAnalysis_h_
#define _TopCharmedMesonAnalysis_h_

#include "TString.h"
#include "TH1F.h"
#include "TopLJets2015/TopAnalysis/interface/ObjectTools.h"

void RunTopCharmedMesonAnalysis(TString filename,
				TString outname,
				Int_t channelSelection, 
				Int_t chargeSelection, 
				FlavourSplitting flavourSplitting,
				TH1F *normH, 
				Bool_t runSysts,
				TString era,
				Bool_t debug=false);

#endif
