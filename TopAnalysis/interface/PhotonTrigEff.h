#ifndef _photontrigeff_h_
#define _photontrigeff_h_

#include "TLorentzVector.h"
#include "TopLJets2015/TopAnalysis/interface/ObjectTools.h"
#include "TopLJets2015/TopAnalysis/interface/SelectionTools.h"

void RunPhotonTrigEff(TString filename,
                     TString outname,
                      TH1F *normH, 
                      TH1F *genPU,
                      TString era,
                      bool debug=false); 

#endif
