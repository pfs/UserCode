#ifndef _topue_h_
#define _topue_h_

#include "TH1.h"
#include "TString.h"
#include "TFile.h"

#include <map>

struct TopUE_t
{
  Int_t run,event,lumi,cat;
  Int_t nvtx;
  Int_t nw;
  Float_t weight[400];
  Int_t passSel,gen_passSel;
  Int_t nj[9],  nb[9];
  Int_t gen_nj, gen_nb;
  Float_t ptpos[5],  phipos[5],  ptll[5],  phill[5],  mll[5],  sumpt[5],  dphill[5];
  Float_t gen_ptpos, gen_phipos, gen_ptll, gen_phill, gen_mll, gen_sumpt, gen_dphill;
  Float_t ptttbar[9],     phittbar[9];
  Float_t parton_ptttbar, parton_phittbar;
  Float_t gen_ptttbar,    gen_phittbar;  
  Int_t n,gen_n,isInBFlags[5000];
  Float_t pt[5000],     eta[5000],     phi[5000];
  Float_t gen_pt[5000], gen_eta[5000], gen_phi[5000]; 
};

void createTopUETree(TTree *t,TopUE_t &tue);
void resetTopUE(TopUE_t &tue);
void RunTopUE(TString filename,
	      TString outname,
	      Int_t channelSelection, 
	      Int_t chargeSelection, 
	      FlavourSplitting flavourSplitting,
	      TH1F *normH, 
	      Bool_t runSysts,
	      TString era);
#endif
