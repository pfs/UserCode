#ifndef _VBFVectorBoson_h_
#define _VBFVectorBoson_h_

#include "TLorentzVector.h"
#include "TopLJets2015/TopAnalysis/interface/ObjectTools.h"
#include "TopLJets2015/TopAnalysis/interface/SelectionTools.h"

//Vector boson will be either Z or photon at the moment

void RunVBFVectorBoson(TString filename,
                       TString outname,
                       Int_t channelSelection, 
                       Int_t chargeSelection, 
                       TH1F *normH, 
                       TH1F *genPU,
                       TString era,
                       Bool_t debug=false);

#endif
