#ifndef _MttbarAnalyzer_h_
#define _MttbarAnalyzer_h_

#include "TLorentzVector.h"
#include "TopLJets2015/TopAnalysis/interface/ObjectTools.h"
#include "TopLJets2015/TopAnalysis/interface/SelectionTools.h"

struct MttbarSummary_t
{
  //new
  float mttbar_chi,mttbar_deltaR,mttbar_random;

  float weight;
  float nvtx;
  float rho;
  float l_pt,l_eta,l_phi,l_m;
  float met_pt, met_phi;
  int nj;
  float nj_f,nb,nl;
  float j_pt[100],j_eta[100],j_phi[100],j_m[100],j_csv[100];
  float j_nch[100],j_PtD[100],j_PtDs[100],j_width[100],j_tau21[100],j_tau32[100],j_tau43[100],j_zg[100];
  float ue_nch,ue_chsumpt,ue_chsumpz;
  float mttbar,gen_mttbar;
};


void RunMttbarAnalyzer(TString filename,
                     TString outname,
                     Int_t channelSelection, 
                     Int_t chargeSelection, 
                     TH1F *normH, 
                     TString era,
                     Bool_t debug=false);

void createMttbarSummaryTree(TTree *t,MttbarSummary_t &summary);


#endif
