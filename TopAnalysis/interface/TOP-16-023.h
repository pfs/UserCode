#ifndef _top16023_h_
#define _top16023_h_

#include "TopLJets2015/TopAnalysis/interface/TOP-16-006.h"

struct LJEvent_t
{
  Float_t w;
  Int_t nj,ngj,ngp,nb;
  Bool_t j_btag[20];
  Float_t gp_pt[20],gp_eta[20],gp_phi[20],gp_m[20];
  Float_t j_pt[20],j_eta[20],j_phi[20],j_m[20];
  Float_t gj_pt[20],gj_eta[20],gj_phi[20],gj_m[20];
  Int_t l_id;
  Float_t l_pt, l_eta, l_phi, l_m;
  Float_t gl_pt, gl_eta, gl_phi, gl_m;
};

bool sortDijetKeys(std::pair< std::pair<int,int>, float > a,
		   std::pair< std::pair<int,int>, float > b);
std::pair<int,int> getDijetsSystemCandidate(std::vector<TLorentzVector> &lightJets);
void RunTop16023(TString filename,
		     TString outname,
		     Int_t channelSelection, 
		     Int_t chargeSelection, 
		     FlavourSplitting flavourSplitting,
		     TH1F *normH, 
		     Bool_t runSysts,
		     TString era);

#endif
