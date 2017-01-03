#ifndef _minievent_h_
#define _minievent_h_

#include "TTree.h"

struct MiniEvent_t
{
  MiniEvent_t()
  {
    ttbar_nw=0;
    ng=0; ngtop=0; ngpf=0;
    nl=0; nj=0; nmet=0; npf=0;
  }

  Bool_t isData;
  Int_t run,event,lumi;


  //gen level event
  Int_t pu,putrue;
  Int_t ttbar_nw, ttbar_allmepartons, ttbar_matchmepartons;
  Float_t ttbar_w[500];
  Int_t ng,ngjets,ngbjets,ngleptons,ngtop,ngpf;
  Int_t g_id[500];
  Float_t g_pt[500],g_eta[500],g_phi[500],g_m[500]; 
  Int_t gtop_id[15];
  Float_t gtop_pt[15],gtop_eta[15],gtop_phi[15],gtop_m[15]; 
  Int_t gpf_id[5000],gpf_c[5000],gpf_g[5000];
  Float_t gpf_pt[5000],gpf_eta[5000],gpf_phi[5000],gpf_m[5000];

  //reco level event
  Int_t nvtx;
  Int_t muTrigger,elTrigger;
  Float_t rho;
  Int_t nl,nleptons;
  Bool_t isPromptFinalState[50], isDirectPromptTauDecayProductFinalState[50];
  Int_t l_id[50],l_charge[50],l_pid[50],l_g[200];
  Float_t l_pt[50],l_eta[50],l_phi[50], l_mass[50], l_scaleUnc[50], l_miniIso[50], l_chargedHadronIso[50], l_relIso[50], l_ip3d[50], l_ip3dsig[50];

  Int_t nj;
  Float_t j_pt[200],j_eta[200],j_phi[200],j_mass[200],j_area[200],j_rawsf[200];
  Float_t j_csv[200],j_cvsl[200],j_cvsb[200],j_vtxmass[200],j_vtx3DVal[200],j_vtx3DSig[200],j_puid[200],j_vtxpx[200],j_vtxpy[200],j_vtxpz[200];
  Bool_t j_btag[200];
  Int_t j_vtxNtracks[200],j_flav[200],j_pid[200],j_hadflav[200],j_g[200];

  //met
  Int_t nmet;
  Float_t met_pt[10],met_phi[10];

  //PF candidates
  Int_t npf,pf_j[5000];
  Int_t pf_id[5000],pf_c[5000];
  Float_t pf_pt[5000],pf_eta[5000],pf_phi[5000],pf_m[5000],pf_puppiWgt[5000];
};

void createMiniEventTree(TTree *t,MiniEvent_t &ev);
void attachToMiniEventTree(TTree *t, MiniEvent_t &ev,bool full=false);

#endif
