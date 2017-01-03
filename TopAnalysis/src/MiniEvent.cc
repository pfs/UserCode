#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"

//
void createMiniEventTree(TTree *t,MiniEvent_t &ev)
{
  //event header
  t->Branch("isData",     &ev.isData,     "isData/O");
  t->Branch("run",       &ev.run,       "run/I");
  t->Branch("event",     &ev.event,     "event/I");
  t->Branch("lumi",      &ev.lumi,      "lumi/I");

  //generator level event
  t->Branch("pu",      &ev.pu,      "pu/I");
  t->Branch("putrue",      &ev.putrue,      "putrue/I");
  t->Branch("ttbar_nw",        &ev.ttbar_nw,        "ttbar_nw/I");
  t->Branch("ttbar_allmepartons",        &ev.ttbar_allmepartons,        "ttbar_allmepartons/I");
  t->Branch("ttbar_matchmepartons",        &ev.ttbar_matchmepartons,        "ttbar_matchmepartons/I");
  t->Branch("ttbar_w",        ev.ttbar_w,        "ttbar_w[ttbar_nw]/F");

  //gen event (jets and dressed leptons)
  t->Branch("ngjets",       &ev.ngjets,       "ngjets/I");
  t->Branch("ngbjets",       &ev.ngbjets,       "ngbjets/I");
  t->Branch("ngleptons",       &ev.ngleptons,       "ngleptons/I");
  t->Branch("ng",       &ev.ng,       "ng/I");
  t->Branch("g_id",      ev.g_id,     "g_id[ng]/I");
  t->Branch("g_pt",      ev.g_pt,     "g_pt[ng]/F");
  t->Branch("g_eta",     ev.g_eta,    "g_eta[ng]/F");
  t->Branch("g_phi",     ev.g_phi,    "g_phi[ng]/F");
  t->Branch("g_m",       ev.g_m,      "g_m[ng]/F");

  //top (lastCopy and pseudo-top)
  t->Branch("ngtop",     &ev.ngtop,      "ngtop/I");
  t->Branch("gtop_id",    ev.gtop_id,    "gtop_id[ngtop]/I");
  t->Branch("gtop_pt",    ev.gtop_pt,    "gtop_pt[ngtop]/F");
  t->Branch("gtop_eta",   ev.gtop_eta,   "gtop_eta[ngtop]/F");
  t->Branch("gtop_phi",   ev.gtop_phi,   "gtop_phi[ngtop]/F");
  t->Branch("gtop_m",     ev.gtop_m,     "gtop_m[ngtop]/F");

  //final state
  t->Branch("ngpf",       &ev.ngpf,       "ngpf/I");
  t->Branch("gpf_id",      ev.gpf_id,     "gpf_id[ngpf]/I");
  t->Branch("gpf_c",       ev.gpf_c,      "gpf_c[ngpf]/I");
  t->Branch("gpf_g",       ev.gpf_g,      "gpf_g[ngpf]/I");
  t->Branch("gpf_pt",      ev.gpf_pt,     "gpf_pt[ngpf]/F");
  t->Branch("gpf_eta",     ev.gpf_eta,    "gpf_eta[ngpf]/F");
  t->Branch("gpf_phi",     ev.gpf_phi,    "gpf_phi[ngpf]/F");
  t->Branch("gpf_m",       ev.gpf_m,      "gpf_m[ngpf]/F");

  //reco level event
  t->Branch("nvtx",      &ev.nvtx,      "nvtx/I");
  t->Branch("rho",      &ev.rho,      "rho/F");
  t->Branch("muTrigger",        &ev.muTrigger,        "muTrigger/I");
  t->Branch("elTrigger",        &ev.elTrigger,        "elTrigger/I");

  t->Branch("nleptons", &ev.nleptons, "nleptons/I");
  t->Branch("nl", &ev.nl, "nl/I");
  t->Branch("isPromptFinalState",                         ev.isPromptFinalState,        "isPromptFinalState[nl]/O");
  t->Branch("isDirectPromptTauDecayProductFinalState",    ev.isDirectPromptTauDecayProductFinalState,        "isDirectPromptTauDecayProductFinalState[nl]/O");
  t->Branch("l_id",       ev.l_id,      "l_id[nl]/I");
  t->Branch("l_pid",      ev.l_pid,     "l_pid[nl]/I");
  t->Branch("l_g",        ev.l_g,       "l_g[nl]/I");
  t->Branch("l_charge",   ev.l_charge,  "l_charge[nl]/I");
  t->Branch("l_pt",       ev.l_pt,      "l_pt[nl]/F");
  t->Branch("l_eta",      ev.l_eta,     "l_eta[nl]/F");
  t->Branch("l_phi",      ev.l_phi,     "l_phi[nl]/F");
  t->Branch("l_mass",     ev.l_mass,    "l_mass[nl]/F");
  t->Branch("l_scaleUnc",       ev.l_scaleUnc,      "l_scaleUnc[nl]/F");
  t->Branch("l_chargedHadronIso", ev.l_chargedHadronIso, "l_chargedHadronIso[nl]/F");
  t->Branch("l_miniIso",          ev.l_miniIso,          "l_miniIso[nl]/F");
  t->Branch("l_relIso",           ev.l_relIso,           "l_relIso[nl]/F");
  t->Branch("l_ip3d",             ev.l_ip3d,             "l_ip3d[nl]/F");
  t->Branch("l_ip3dsig",          ev.l_ip3dsig,          "l_ip3dsig[nl]/F");

  //jet info
  t->Branch("nj",        &ev.nj,        "nj/I");
  t->Branch("j_g",        ev.j_g,   "j_g[nj]/I");
  t->Branch("j_area",     ev.j_area,      "j_area[nj]/F");
  t->Branch("j_rawsf",    ev.j_rawsf,      "j_rawsf[nj]/F");
  t->Branch("j_pt",       ev.j_pt,      "j_pt[nj]/F");
  t->Branch("j_eta",      ev.j_eta,     "j_eta[nj]/F");
  t->Branch("j_phi",      ev.j_phi,     "j_phi[nj]/F");
  t->Branch("j_mass",     ev.j_mass,     "j_mass[nj]/F");
  t->Branch("j_csv",      ev.j_csv,     "j_csv[nj]/F");
  t->Branch("j_btag",     ev.j_btag,    "j_btag[nj]/O");
  t->Branch("j_csvl",     ev.j_cvsl,     "j_cvsl[nj]/F");
  t->Branch("j_cvsb",     ev.j_cvsb,     "j_cvsb[nj]/F");
  t->Branch("j_vtxpx",    ev.j_vtxpx, "j_vtxpx[nj]/F");
  t->Branch("j_vtxpy",    ev.j_vtxpy, "j_vtxpy[nj]/F");
  t->Branch("j_vtxpz",    ev.j_vtxpz, "j_vtxpz[nj]/F");
  t->Branch("j_vtxmass",  ev.j_vtxmass, "j_vtxmass[nj]/F");
  t->Branch("j_vtxNtracks",  ev.j_vtxNtracks, "j_vtxNtracks[nj]/I");
  t->Branch("j_vtx3DVal",    ev.j_vtx3DVal, "j_vtx3DVal[nj]/F");
  t->Branch("j_vtx3DSig",    ev.j_vtx3DSig, "j_vtx3DSig[nj]/F");
  t->Branch("j_puid",        ev.j_puid,    "j_puid[nj]/F");  
  t->Branch("j_flav",        ev.j_flav,    "j_flav[nj]/I");
  t->Branch("j_hadflav",     ev.j_hadflav,    "j_hadflav[nj]/I");
  t->Branch("j_pid",         ev.j_pid,     "j_pid[nj]/I");

  //pf candidates (only charged if outside jets)
  t->Branch("npf",        &ev.npf,         "npf/I");
  t->Branch("pf_j",        ev.pf_j,        "pf_j[npf]/I");
  t->Branch("pf_id",       ev.pf_id,       "pf_id[npf]/I");
  t->Branch("pf_c",        ev.pf_c,        "pf_c[npf]/I");
  t->Branch("pf_pt",       ev.pf_pt,       "pf_pt[npf]/F");
  t->Branch("pf_eta",      ev.pf_eta,      "pf_eta[npf]/F");
  t->Branch("pf_phi",      ev.pf_phi,      "pf_phi[npf]/F");
  t->Branch("pf_m",        ev.pf_m,        "pf_m[npf]/F");
  t->Branch("pf_puppiWgt", ev.pf_puppiWgt, "pf_puppiWgt[npf]/F");

  //MET
  t->Branch("nmet",      &ev.nmet,     "nmet/I");
  t->Branch("met_pt",     ev.met_pt,   "met_pt[nmet]/F");
  t->Branch("met_phi",    ev.met_phi,  "met_phi[nmet]/F");
}

//
void attachToMiniEventTree(TTree *t,MiniEvent_t &ev,bool full)
{
  //event header
  t->SetBranchAddress("isData",     &ev.isData);
  t->SetBranchAddress("run",       &ev.run);
  t->SetBranchAddress("event",     &ev.event);
  t->SetBranchAddress("lumi",      &ev.lumi);

  //generator level event
  t->SetBranchAddress("pu",      &ev.pu);
  t->SetBranchAddress("putrue",      &ev.putrue);
  t->SetBranchAddress("ttbar_nw",        &ev.ttbar_nw);
  t->SetBranchAddress("ttbar_allmepartons",        &ev.ttbar_allmepartons);
  t->SetBranchAddress("ttbar_matchmepartons",        &ev.ttbar_matchmepartons);
  t->SetBranchAddress("ttbar_w",        ev.ttbar_w);

  //gen event (jets and dressed leptons)
  t->SetBranchAddress("ngjets",       &ev.ngjets);
  t->SetBranchAddress("ngbjets",       &ev.ngbjets);
  t->SetBranchAddress("ngleptons",       &ev.ngleptons);
  t->SetBranchAddress("ng",       &ev.ng);
  t->SetBranchAddress("g_id",      ev.g_id);
  t->SetBranchAddress("g_pt",      ev.g_pt);
  t->SetBranchAddress("g_eta",     ev.g_eta);
  t->SetBranchAddress("g_phi",     ev.g_phi);
  t->SetBranchAddress("g_m",       ev.g_m);

  //top (lastCopy and pseudo-top)
  t->SetBranchAddress("ngtop",     &ev.ngtop);
  t->SetBranchAddress("gtop_id",    ev.gtop_id);
  t->SetBranchAddress("gtop_pt",    ev.gtop_pt);
  t->SetBranchAddress("gtop_eta",   ev.gtop_eta);
  t->SetBranchAddress("gtop_phi",   ev.gtop_phi);
  t->SetBranchAddress("gtop_m",     ev.gtop_m);

  //final state
  if(full)
    {
      t->SetBranchAddress("ngpf",       &ev.ngpf);
      t->SetBranchAddress("gpf_id",      ev.gpf_id);
      t->SetBranchAddress("gpf_c",       ev.gpf_c);
      t->SetBranchAddress("gpf_g",       ev.gpf_g); 
      t->SetBranchAddress("gpf_pt",      ev.gpf_pt);
      t->SetBranchAddress("gpf_eta",     ev.gpf_eta);
      t->SetBranchAddress("gpf_phi",     ev.gpf_phi);
      t->SetBranchAddress("gpf_m",       ev.gpf_m);
    }

  //reco level event
  t->SetBranchAddress("nvtx",      &ev.nvtx);
  t->SetBranchAddress("muTrigger",        &ev.muTrigger);
  t->SetBranchAddress("elTrigger",        &ev.elTrigger);

  t->SetBranchAddress("nleptons", &ev.nleptons);
  t->SetBranchAddress("nl", &ev.nl);
  t->SetBranchAddress("isPromptFinalState",                         ev.isPromptFinalState);
  t->SetBranchAddress("isDirectPromptTauDecayProductFinalState",    ev.isDirectPromptTauDecayProductFinalState);
  t->SetBranchAddress("l_id",       ev.l_id);
  t->SetBranchAddress("l_pid",      ev.l_pid);
  t->SetBranchAddress("l_g",        ev.l_g);
  t->SetBranchAddress("l_charge",   ev.l_charge);
  t->SetBranchAddress("l_pt",       ev.l_pt);
  t->SetBranchAddress("l_eta",      ev.l_eta);
  t->SetBranchAddress("l_phi",      ev.l_phi);
  t->SetBranchAddress("l_mass",     ev.l_mass);
  t->SetBranchAddress("l_scaleUnc",       ev.l_scaleUnc);
  t->SetBranchAddress("l_chargedHadronIso", ev.l_chargedHadronIso);
  t->SetBranchAddress("l_miniIso",          ev.l_miniIso);
  t->SetBranchAddress("l_relIso",           ev.l_relIso);
  t->SetBranchAddress("l_ip3d",             ev.l_ip3d);
  t->SetBranchAddress("l_ip3dsig",          ev.l_ip3dsig);

  //jet info
  t->SetBranchAddress("nj",        &ev.nj);
  t->SetBranchAddress("j_g",        ev.j_g);
  t->SetBranchAddress("j_area",     ev.j_area);
  t->SetBranchAddress("j_rawsf",    ev.j_rawsf);
  t->SetBranchAddress("j_pt",       ev.j_pt);
  t->SetBranchAddress("j_eta",      ev.j_eta);
  t->SetBranchAddress("j_phi",      ev.j_phi);
  t->SetBranchAddress("j_mass",     ev.j_mass);
  t->SetBranchAddress("j_csv",      ev.j_csv);
  t->SetBranchAddress("j_btag",     ev.j_btag);
  t->SetBranchAddress("j_csvl",     ev.j_cvsl);
  t->SetBranchAddress("j_cvsb",     ev.j_cvsb);
  t->SetBranchAddress("j_vtxpx",    ev.j_vtxpx);
  t->SetBranchAddress("j_vtxpy",    ev.j_vtxpy);
  t->SetBranchAddress("j_vtxpz",    ev.j_vtxpz);
  t->SetBranchAddress("j_vtxmass",  ev.j_vtxmass);
  t->SetBranchAddress("j_vtxNtracks",  ev.j_vtxNtracks);
  t->SetBranchAddress("j_vtx3DVal",    ev.j_vtx3DVal);
  t->SetBranchAddress("j_vtx3DSig",    ev.j_vtx3DSig);
  t->SetBranchAddress("j_puid",        ev.j_puid);
  t->SetBranchAddress("j_flav",        ev.j_flav);
  t->SetBranchAddress("j_hadflav",     ev.j_hadflav);
  t->SetBranchAddress("j_pid",         ev.j_pid);

  //pf candidates (only charged if outside jets)
  if(full)
    {
      t->SetBranchAddress("npf",        &ev.npf);
      t->SetBranchAddress("pf_j",        ev.pf_j);
      t->SetBranchAddress("pf_id",       ev.pf_id);
      t->SetBranchAddress("pf_c",        ev.pf_c);
      t->SetBranchAddress("pf_pt",       ev.pf_pt);
      t->SetBranchAddress("pf_eta",      ev.pf_eta);
      t->SetBranchAddress("pf_phi",      ev.pf_phi);
      t->SetBranchAddress("pf_m",        ev.pf_m);
      t->SetBranchAddress("pf_puppiWgt", ev.pf_puppiWgt);
    }

  //MET
  t->SetBranchAddress("nmet",      &ev.nmet);
  t->SetBranchAddress("met_pt",    ev.met_pt);
  t->SetBranchAddress("met_phi",   ev.met_phi);
}
