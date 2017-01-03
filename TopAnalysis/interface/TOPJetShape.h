#ifndef _topjetshape_h_
#define _topjetshape_h_

#include "TLorentzVector.h"
#include "TopLJets2015/TopAnalysis/interface/ObjectTools.h"

void RunTopJetShape(TString filename,
				TString outname,
				Int_t channelSelection, 
				Int_t chargeSelection, 
				FlavourSplitting flavourSplitting,
				TH1F *normH, 
				Bool_t runSysts,
				TString era,
				Bool_t debug=false);

struct TopJetShapeEvent_t
{
  Int_t nw, nl, ngl, nj, ngj;
  Int_t gen_sel, reco_sel;
  
  Float_t weight[10];
  
  Float_t l_pt[5], l_eta[5], l_phi[5], l_m[5];
  Int_t l_id[5];
  
  Float_t gl_pt[5], gl_eta[5], gl_phi[5], gl_m[5];
  Int_t gl_id[5];
  
  Float_t j_pt[50], j_eta[50], j_phi[50], j_m[50];
  Int_t j_flavor[50], j_overlap[50], j_gj[50];
  
  Float_t gj_pt[50], gj_eta[50], gj_phi[50], gj_m[50];
  Int_t gj_flavor[50], gj_overlap[50], gj_j[50];
  
  Float_t met_pt,met_phi;
  
  Float_t j_mult_charged[50];
  Float_t j_mult_puppi[50];
  Float_t j_mult_all[50];
  Float_t j_width_charged[50];
  Float_t j_width_puppi[50];
  Float_t j_width_all[50];
  Float_t j_ptd_charged[50];
  Float_t j_ptd_puppi[50];
  Float_t j_ptd_all[50];
  Float_t j_ecc_charged[50];
  Float_t j_ecc_puppi[50];
  Float_t j_ecc_all[50];
  Float_t j_tau21_charged[50];
  Float_t j_tau21_puppi[50];
  Float_t j_tau21_all[50];
  Float_t j_tau32_charged[50];
  Float_t j_tau32_puppi[50];
  Float_t j_tau32_all[50];
  Float_t j_tau43_charged[50];
  Float_t j_tau43_puppi[50];
  Float_t j_tau43_all[50];
  Float_t j_zg_charged[50];
  Float_t j_zg_puppi[50];
  Float_t j_zg_all[50];
  Float_t j_zgxdr_charged[50];
  Float_t j_zgxdr_puppi[50];
  Float_t j_zgxdr_all[50];
  Float_t j_zgdr_charged[50];
  Float_t j_zgdr_puppi[50];
  Float_t j_zgdr_all[50];
  Float_t j_ga_width_charged[50];
  Float_t j_ga_width_puppi[50];
  Float_t j_ga_width_all[50];
  Float_t j_ga_lha_charged[50];
  Float_t j_ga_lha_puppi[50];
  Float_t j_ga_lha_all[50];
  Float_t j_ga_thrust_charged[50];
  Float_t j_ga_thrust_puppi[50];
  Float_t j_ga_thrust_all[50];
  Float_t j_c1_02_charged[50];
  Float_t j_c1_02_puppi[50];
  Float_t j_c1_02_all[50];
  Float_t j_c1_05_charged[50];
  Float_t j_c1_05_puppi[50];
  Float_t j_c1_05_all[50];
  Float_t j_c1_10_charged[50];
  Float_t j_c1_10_puppi[50];
  Float_t j_c1_10_all[50];
  Float_t j_c1_20_charged[50];
  Float_t j_c1_20_puppi[50];
  Float_t j_c1_20_all[50];
  Float_t j_c2_02_charged[50];
  Float_t j_c2_02_puppi[50];
  Float_t j_c2_02_all[50];
  Float_t j_c2_05_charged[50];
  Float_t j_c2_05_puppi[50];
  Float_t j_c2_05_all[50];
  Float_t j_c2_10_charged[50];
  Float_t j_c2_10_puppi[50];
  Float_t j_c2_10_all[50];
  Float_t j_c2_20_charged[50];
  Float_t j_c2_20_puppi[50];
  Float_t j_c2_20_all[50];
  Float_t j_c3_02_charged[50];
  Float_t j_c3_02_puppi[50];
  Float_t j_c3_02_all[50];
  Float_t j_c3_05_charged[50];
  Float_t j_c3_05_puppi[50];
  Float_t j_c3_05_all[50];
  Float_t j_c3_10_charged[50];
  Float_t j_c3_10_puppi[50];
  Float_t j_c3_10_all[50];
  Float_t j_c3_20_charged[50];
  Float_t j_c3_20_puppi[50];
  Float_t j_c3_20_all[50];

  Float_t gj_mult_charged[50];
  Float_t gj_mult_puppi[50];
  Float_t gj_mult_all[50];
  Float_t gj_width_charged[50];
  Float_t gj_width_puppi[50];
  Float_t gj_width_all[50];
  Float_t gj_ptd_charged[50];
  Float_t gj_ptd_puppi[50];
  Float_t gj_ptd_all[50];
  Float_t gj_ecc_charged[50];
  Float_t gj_ecc_puppi[50];
  Float_t gj_ecc_all[50];
  Float_t gj_tau21_charged[50];
  Float_t gj_tau21_puppi[50];
  Float_t gj_tau21_all[50];
  Float_t gj_tau32_charged[50];
  Float_t gj_tau32_puppi[50];
  Float_t gj_tau32_all[50];
  Float_t gj_tau43_charged[50];
  Float_t gj_tau43_puppi[50];
  Float_t gj_tau43_all[50];
  Float_t gj_zg_charged[50];
  Float_t gj_zg_puppi[50];
  Float_t gj_zg_all[50];
  Float_t gj_zgxdr_charged[50];
  Float_t gj_zgxdr_puppi[50];
  Float_t gj_zgxdr_all[50];
  Float_t gj_zgdr_charged[50];
  Float_t gj_zgdr_puppi[50];
  Float_t gj_zgdr_all[50];
  Float_t gj_ga_width_charged[50];
  Float_t gj_ga_width_puppi[50];
  Float_t gj_ga_width_all[50];
  Float_t gj_ga_lha_charged[50];
  Float_t gj_ga_lha_puppi[50];
  Float_t gj_ga_lha_all[50];
  Float_t gj_ga_thrust_charged[50];
  Float_t gj_ga_thrust_puppi[50];
  Float_t gj_ga_thrust_all[50];
  Float_t gj_c1_02_charged[50];
  Float_t gj_c1_02_puppi[50];
  Float_t gj_c1_02_all[50];
  Float_t gj_c1_05_charged[50];
  Float_t gj_c1_05_puppi[50];
  Float_t gj_c1_05_all[50];
  Float_t gj_c1_10_charged[50];
  Float_t gj_c1_10_puppi[50];
  Float_t gj_c1_10_all[50];
  Float_t gj_c1_20_charged[50];
  Float_t gj_c1_20_puppi[50];
  Float_t gj_c1_20_all[50];
  Float_t gj_c2_02_charged[50];
  Float_t gj_c2_02_puppi[50];
  Float_t gj_c2_02_all[50];
  Float_t gj_c2_05_charged[50];
  Float_t gj_c2_05_puppi[50];
  Float_t gj_c2_05_all[50];
  Float_t gj_c2_10_charged[50];
  Float_t gj_c2_10_puppi[50];
  Float_t gj_c2_10_all[50];
  Float_t gj_c2_20_charged[50];
  Float_t gj_c2_20_puppi[50];
  Float_t gj_c2_20_all[50];
  Float_t gj_c3_02_charged[50];
  Float_t gj_c3_02_puppi[50];
  Float_t gj_c3_02_all[50];
  Float_t gj_c3_05_charged[50];
  Float_t gj_c3_05_puppi[50];
  Float_t gj_c3_05_all[50];
  Float_t gj_c3_10_charged[50];
  Float_t gj_c3_10_puppi[50];
  Float_t gj_c3_10_all[50];
  Float_t gj_c3_20_charged[50];
  Float_t gj_c3_20_puppi[50];
  Float_t gj_c3_20_all[50];

};

double calcGA(double beta, double kappa, Jet jet, bool includeNeutrals = false, bool usePuppi = false, double ptcut = 1.0);
void createTopJetShapeEventTree(TTree *t,TopJetShapeEvent_t &tjsev);
void resetTopJetShapeEvent(TopJetShapeEvent_t &tjsev);

double getMult(Jet jet, bool includeNeutrals = false, bool usePuppi = false, double ptcut = 1.0);
double getPtD(Jet jet, bool includeNeutrals = false, bool usePuppi = false, double ptcut = 1.0);
double getWidth(Jet jet, bool includeNeutrals = false, bool usePuppi = false, double ptcut = 1.0);
double getEcc(Jet jet, bool includeNeutrals = false, bool usePuppi = false, double ptcut = 1.0);
double getTau(int N, int M, Jet jet, bool includeNeutrals = false, bool usePuppi = false, double ptcut = 1.0);
double getC(int N, double beta, Jet jet, bool includeNeutrals = false, bool usePuppi = false, double ptcut = 1.0);
std::vector<double> getZg(Jet jet, bool includeNeutrals = false, bool usePuppi = false, double ptcut = 1.0);

double deltaR(TLorentzVector v1, TLorentzVector v2) { return v1.DeltaR(v2); }
double mapAngleMPiToPi(double phi) { return TVector2::Phi_mpi_pi(phi); }
#endif
