#ifndef _minievent_h_
#define _minievent_h_

#include "TTree.h"

struct MiniEvent_t
{
  MiniEvent_t()
  {
    g_nw=0; g_npsw=0; ng=0; ngtop=0; 
    ngamma=0; nl=0; nj=0; 
  }

  Bool_t isData;
  UInt_t run,lumi;
  ULong64_t event;

  //gen level event
  Int_t g_id1, g_id2;
  Float_t g_x1, g_x2, g_qscale;
  Int_t g_pu,g_putrue;
  Int_t g_nw,g_npsw,g_nup;
  Float_t g_w[500],g_psw[500];
  Int_t ng,ngtop;
  Int_t g_id[500],g_bid[500],g_tagCtrs[500];
  Bool_t g_isSemiLepBhad[500];
  Float_t g_pt[500],g_eta[500],g_phi[500],g_m[500],g_xb[500],g_xbp[500]; 
  Int_t gtop_id[25];
  Float_t gtop_pt[25],gtop_eta[25],gtop_phi[25],gtop_m[25]; 
  Int_t g_nchPV;
  Float_t g_sumPVChPt,g_sumPVChPz,g_sumPVChHt;

  //reco level event
  Int_t nvtx;
  Int_t triggerBits,addTriggerBits;
  Int_t zeroBiasPS;
  Float_t rho;

  //leptons
  Int_t nl;
  Bool_t l_isPromptFinalState[50], l_isDirectPromptTauDecayProductFinalState[50];
  Int_t l_id[50],l_charge[50],l_pid[50],l_g[200];
  Float_t l_pt[50],l_eta[50],l_phi[50], l_mass[50], l_highpt[50],
    l_scaleUnc1[50], l_scaleUnc2[50], l_scaleUnc3[50], l_scaleUnc4[50], l_scaleUnc5[50], l_scaleUnc6[50], l_scaleUnc7[50],   
    l_miniIso[50], l_chargedHadronIso[50], l_relIso[50], l_ip3d[50], l_ip3dsig[50],l_mva[50],l_mvaCats[50];

  Int_t ngamma;
  Bool_t gamma_isPromptFinalState[50];
  Int_t gamma_pid[50],gamma_idFlags[50],gamma_g[50];
  Float_t gamma_pt[50],gamma_eta[50],gamma_phi[50], 
    gamma_scaleUnc1[50],gamma_scaleUnc2[50],gamma_scaleUnc3[50],gamma_scaleUnc4[50],gamma_scaleUnc5[50],gamma_scaleUnc6[50],gamma_scaleUnc7[50],
    gamma_mva[50], gamma_mvaCats[50],
    gamma_chargedHadronIso[50],gamma_neutralHadronIso[50],gamma_photonIso[50],gamma_hoe[50],gamma_r9[50],gamma_sieie[50];

  Int_t nj;
  Float_t j_pt[200],j_eta[200],j_phi[200],j_mass[200],j_area[200],j_rawsf[200];
  Float_t j_jerUp[200],j_jerDn[200],j_jecUp[30][200],j_jecDn[30][200];
  Float_t j_csv[200],j_deepcsv[200],j_pumva[200],j_emf[200],j_qg[200];
  Float_t j_c2_00[200],j_c2_02[200],j_c2_05[200],j_c2_10[200],j_c2_20[200];
  Float_t j_zg[200],j_mult[200],j_gaptd[200],j_gawidth[200],j_gathrust[200],j_tau32[200],j_tau21[200];
  Float_t j_vtxmass[200],j_vtx3DVal[200],j_vtx3DSig[200],j_vtxpx[200],j_vtxpy[200],j_vtxpz[200];
  Bool_t j_btag[200];
  Int_t j_vtxNtracks[200],j_flav[200],j_id[200],j_pid[200],j_hadflav[200],j_g[200];

  //met
  Float_t met_pt,met_phi,met_sig;
  Float_t met_ptShifted[14],met_phiShifted[14];
  Int_t met_filterBits;
  
  //event energy fluxes (PF-based)
  Int_t nchPV;
  Float_t sumPVChPt,sumPVChPz,sumPVChHt;
  Int_t nPFCands[8],nPFChCands[8];
  Float_t sumPFHt[8],sumPFEn[8],sumPFPz[8],sumPFChHt[8],sumPFChEn[8],sumPFChPz[8];

  //CTPPS protons
  Short_t nfwdtrk,fwdtrk_pot[50],fwdtrk_method[50];
  Float_t fwdtrk_ex[50],fwdtrk_ey[50],fwdtrk_ez[50],fwdtrk_y[50],fwdtrk_chisqnorm[50],fwdtrk_xi[50],fwdtrk_t[50];

  //these are crazy variables for the cross check
  Int_t nrawmu;
  Short_t rawmu_pt[200],rawmu_eta[200],rawmu_phi[200];
  Int_t rawmu_pid[200];
};

void createMiniEventTree(TTree *t,MiniEvent_t &ev,Int_t njecUncs=0);
void attachToMiniEventTree(TTree *t, MiniEvent_t &ev,bool full=false);

#endif
