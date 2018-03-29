#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"

//
void createMiniEventTree(TTree *t,MiniEvent_t &ev)
{
  //event header
  t->Branch("isData",    &ev.isData,   "isData/O");
  t->Branch("run",       &ev.run,      "run/i");
  t->Branch("event",     &ev.event,    "event/i");
  t->Branch("lumi",      &ev.lumi,     "lumi/i");

  //generator level information
  t->Branch("g_pu",      &ev.g_pu,     "g_pu/I");
  t->Branch("g_putrue",  &ev.g_putrue, "g_putrue/I");
  t->Branch("g_id1",     &ev.g_id1,    "g_id1/I");
  t->Branch("g_id2",     &ev.g_id2,    "g_id2/I");
  t->Branch("g_x1",      &ev.g_x1,     "g_x1/F");
  t->Branch("g_x2",      &ev.g_x2,     "g_x2/F");
  t->Branch("g_qscale",  &ev.g_qscale, "g_qscale/F");
  t->Branch("g_nw",      &ev.g_nw,     "g_nw/I");
  t->Branch("g_w",        ev.g_w,      "g_w[g_nw]/F");

  //gen event (jets and dressed leptons)
  t->Branch("ng",       &ev.ng,       "ng/I");
  t->Branch("g_id",      ev.g_id,     "g_id[ng]/I");
  t->Branch("g_bid",     ev.g_bid,    "g_bid[ng]/I");
  t->Branch("g_tagCtrs",     ev.g_tagCtrs,    "g_tagCtrs[ng]/I");
  t->Branch("g_isSemiLepBhad",     ev.g_isSemiLepBhad,    "g_isSemiLepBhad[ng]/O");
  t->Branch("g_xb",      ev.g_xb,     "g_xb[ng]/F");  
  t->Branch("g_xbp",      ev.g_xbp,   "g_xbp[ng]/F");
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

  //final state particles
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
  t->Branch("triggerBits",   &ev.triggerBits,        "triggerBits/I");

  //leptons
  t->Branch("nl", &ev.nl, "nl/I");
  t->Branch("l_isPromptFinalState",                         ev.l_isPromptFinalState,                       "l_isPromptFinalState[nl]/O");
  t->Branch("l_isDirectPromptTauDecayProductFinalState",    ev.l_isDirectPromptTauDecayProductFinalState,  "l_isDirectPromptTauDecayProductFinalState[nl]/O");
  t->Branch("l_id",       ev.l_id,      "l_id[nl]/I");
  t->Branch("l_pid",      ev.l_pid,     "l_pid[nl]/I");
  t->Branch("l_g",        ev.l_g,       "l_g[nl]/I");
  t->Branch("l_charge",   ev.l_charge,  "l_charge[nl]/I");
  t->Branch("l_mva",       ev.l_mva,      "l_mva[nl]/F");
  t->Branch("l_pt",       ev.l_pt,      "l_pt[nl]/F");
  t->Branch("l_eta",      ev.l_eta,     "l_eta[nl]/F");
  t->Branch("l_phi",      ev.l_phi,     "l_phi[nl]/F");
  t->Branch("l_mass",     ev.l_mass,    "l_mass[nl]/F");
  t->Branch("l_scaleUnc",         ev.l_scaleUnc,         "l_scaleUnc[nl]/F");
  t->Branch("l_smearUnc",         ev.l_smearUnc,         "l_smearUnc[nl]/F");
  t->Branch("l_chargedHadronIso", ev.l_chargedHadronIso, "l_chargedHadronIso[nl]/F");
  t->Branch("l_miniIso",          ev.l_miniIso,          "l_miniIso[nl]/F");
  t->Branch("l_relIso",           ev.l_relIso,           "l_relIso[nl]/F");
  t->Branch("l_ip3d",             ev.l_ip3d,             "l_ip3d[nl]/F");
  t->Branch("l_ip3dsig",          ev.l_ip3dsig,          "l_ip3dsig[nl]/F");

  //photons
  t->Branch("ngamma",                   &ev.ngamma,                   "ngamma/I");
  t->Branch("gamma_isPromptFinalState",  ev.gamma_isPromptFinalState, "gamma_isPromptFinalState[ngamma]/O");
  t->Branch("gamma_passElectronVeto",    ev.gamma_passElectronVeto,   "gamma_passElectronVeto[ngamma]/O");
  t->Branch("gamma_hasPixelSeed",        ev.gamma_hasPixelSeed,       "gamma_hasPixelSeed[ngamma]/O");
  t->Branch("gamma_pid",                 ev.gamma_pid,                "gamma_pid[ngamma]/I");
  t->Branch("gamma_g",                   ev.gamma_g,                  "gamma_g[ngamma]/I");
  t->Branch("gamma_mva",                 ev.gamma_mva,                "gamma_mva[ngamma]/F");
  t->Branch("gamma_pt",                  ev.gamma_pt,                 "gamma_pt[ngamma]/F");
  t->Branch("gamma_eta",                 ev.gamma_eta,                "gamma_eta[ngamma]/F");
  t->Branch("gamma_phi",                 ev.gamma_phi,                "gamma_phi[ngamma]/F");  
  t->Branch("gamma_scaleUnc",            ev.gamma_scaleUnc,           "gamma_scaleUnc[ngamma]/F");
  t->Branch("gamma_smearUnc",            ev.gamma_smearUnc,           "gamma_smearUnc[ngamma]/F");
  t->Branch("gamma_chargedHadronIso",    ev.gamma_chargedHadronIso,   "gamma_chargedHadronIso[ngamma]/F");
  t->Branch("gamma_neutralHadronIso",    ev.gamma_neutralHadronIso,   "gamma_neutralHadronIso[ngamma]/F");
  t->Branch("gamma_photonIso",           ev.gamma_photonIso,          "gamma_photonIso[ngamma]/F");
  t->Branch("gamma_hoe",                 ev.gamma_hoe,                "gamma_hoe[ngamma]/F");
  t->Branch("gamma_r9",                  ev.gamma_r9,                 "gamma_r9[ngamma]/F");
  t->Branch("gamma_sieie",               ev.gamma_sieie,              "gamma_sieie[ngamma]/F");

  //jet info
  t->Branch("nj",        &ev.nj,        "nj/I");
  t->Branch("j_g",        ev.j_g,       "j_g[nj]/I");
  t->Branch("j_area",     ev.j_area,    "j_area[nj]/F");
  t->Branch("j_rawsf",    ev.j_rawsf,   "j_rawsf[nj]/F");
  t->Branch("j_pt",       ev.j_pt,      "j_pt[nj]/F");
  t->Branch("j_eta",      ev.j_eta,     "j_eta[nj]/F");
  t->Branch("j_phi",      ev.j_phi,     "j_phi[nj]/F");
  t->Branch("j_mass",     ev.j_mass,    "j_mass[nj]/F");
  t->Branch("j_pumva",    ev.j_pumva,   "j_pumva[nj]/F");
  t->Branch("j_id",       ev.j_id,      "j_id[nj]/I");
  t->Branch("j_csv",      ev.j_csv,     "j_csv[nj]/F");
  t->Branch("j_btag",     ev.j_btag,    "j_btag[nj]/O");
  t->Branch("j_c1_00",    ev.j_c1_00,   "j_c1_00[nj]/F");
  t->Branch("j_c1_02",    ev.j_c1_02,   "j_c1_02[nj]/F");
  t->Branch("j_c1_05",    ev.j_c1_05,   "j_c1_05[nj]/F");
  t->Branch("j_c1_10",    ev.j_c1_10,   "j_c1_10[nj]/F");
  t->Branch("j_c1_20",    ev.j_c1_20,   "j_c1_20[nj]/F");
  t->Branch("j_c2_00",    ev.j_c2_00,   "j_c2_00[nj]/F");
  t->Branch("j_c2_02",    ev.j_c2_02,   "j_c2_02[nj]/F");
  t->Branch("j_c2_05",    ev.j_c2_05,   "j_c2_05[nj]/F");
  t->Branch("j_c2_10",    ev.j_c2_10,   "j_c2_10[nj]/F");
  t->Branch("j_c2_20",    ev.j_c2_20,   "j_c2_20[nj]/F");
  t->Branch("j_c3_00",    ev.j_c3_00,   "j_c3_00[nj]/F");
  t->Branch("j_c3_02",    ev.j_c3_02,   "j_c3_02[nj]/F");
  t->Branch("j_c3_05",    ev.j_c3_05,   "j_c3_05[nj]/F");
  t->Branch("j_c3_10",    ev.j_c3_10,   "j_c3_10[nj]/F");
  t->Branch("j_c3_20",    ev.j_c3_20,   "j_c3_20[nj]/F");
  t->Branch("j_zg",       ev.j_zg,      "j_zg[nj]/F");
  t->Branch("j_mult",     ev.j_mult,    "j_mult[nj]/F");
  t->Branch("j_gaptd",    ev.j_gaptd,   "j_gaptd[nj]/F");
  t->Branch("j_gawidth",  ev.j_gawidth, "j_gawidth[nj]/F");
  t->Branch("j_gathrust", ev.j_gathrust,"j_gathrust[nj]/F");
  t->Branch("j_tau32",    ev.j_tau32,   "j_tau32[nj]/F");
  t->Branch("j_tau21",    ev.j_tau21,   "j_tau21[nj]/F");
  t->Branch("j_deepcsv",  ev.j_deepcsv, "j_deepcsv[nj]/F");
  t->Branch("j_vtxpx",    ev.j_vtxpx,   "j_vtxpx[nj]/F");
  t->Branch("j_vtxpy",    ev.j_vtxpy,   "j_vtxpy[nj]/F");
  t->Branch("j_vtxpz",    ev.j_vtxpz,   "j_vtxpz[nj]/F");
  t->Branch("j_vtxmass",  ev.j_vtxmass, "j_vtxmass[nj]/F");
  t->Branch("j_vtxNtracks",  ev.j_vtxNtracks, "j_vtxNtracks[nj]/I");
  t->Branch("j_vtx3DVal",    ev.j_vtx3DVal,   "j_vtx3DVal[nj]/F");
  t->Branch("j_vtx3DSig",    ev.j_vtx3DSig,   "j_vtx3DSig[nj]/F");
  t->Branch("j_flav",        ev.j_flav,       "j_flav[nj]/I");
  t->Branch("j_hadflav",     ev.j_hadflav,    "j_hadflav[nj]/I");
  t->Branch("j_pid",         ev.j_pid,        "j_pid[nj]/I");

  //pf candidates (only charged if outside jets)
  t->Branch("npf",        &ev.npf,         "npf/I");
  t->Branch("pf_j",        ev.pf_j,        "pf_j[npf]/I");
  t->Branch("pf_id",       ev.pf_id,       "pf_id[npf]/I");
  t->Branch("pf_c",        ev.pf_c,        "pf_c[npf]/I");
  t->Branch("pf_pt",       ev.pf_pt,       "pf_pt[npf]/F");
  t->Branch("pf_eta",      ev.pf_eta,      "pf_eta[npf]/F");
  t->Branch("pf_phi",      ev.pf_phi,      "pf_phi[npf]/F");
  t->Branch("pf_m",        ev.pf_m,        "pf_m[npf]/F");
  t->Branch("pf_dxy",      ev.pf_dxy,      "pf_dxy[npf]/F");
  t->Branch("pf_dz",       ev.pf_dz,       "pf_dz[npf]/F");
  //t->Branch("pf_dxyUnc",   ev.pf_dxyUnc,   "pf_dxyUnc[npf]/F");
  //t->Branch("pf_dzUnc",    ev.pf_dzUnc,    "pf_dzUnc[npf]/F");
  //t->Branch("pf_vtxRef",   ev.pf_vtxRef,   "pf_vtxRef[npf]/I");
  //t->Branch("pf_pvAssoc",  ev.pf_pvAssoc,  "pf_pvAssoc[npf]/I");
  t->Branch("pf_puppiWgt", ev.pf_puppiWgt, "pf_puppiWgt[npf]/F");

  //MET
  t->Branch("nmet",       &ev.nmet,       "nmet/I");
  t->Branch("met_pt",      ev.met_pt,     "met_pt[nmet]/F");
  t->Branch("met_phi",     ev.met_phi,    "met_phi[nmet]/F");
  t->Branch("met_sig",     ev.met_sig,    "met_sig[nmet]/F");
  t->Branch("met_filterBits", &ev.met_filterBits, "met_filterBits/I");

  //CTPPS local tracks
  t->Branch("nfwdtrk",      &ev.nfwdtrk,       "nfwdtrk/I");
  t->Branch("fwdtrk_arm",    ev.fwdtrk_arm,    "fwdtrk_arm[nfwdtrk]/I");
  t->Branch("fwdtrk_station",ev.fwdtrk_station,"fwdtrk_station[nfwdtrk]/I");
  t->Branch("fwdtrk_pot",    ev.fwdtrk_pot,    "fwdtrk_pot[nfwdtrk]/I");
  t->Branch("fwdtrk_x",      ev.fwdtrk_x,      "fwdtrk_x[nfwdtrk]/F");
  t->Branch("fwdtrk_x_unc",  ev.fwdtrk_x_unc,  "fwdtrk_x_unc[nfwdtrk]/F");
  t->Branch("fwdtrk_y",      ev.fwdtrk_y,      "fwdtrk_y[nfwdtrk]/F");
  t->Branch("fwdtrk_y_unc",  ev.fwdtrk_y_unc,  "fwdtrk_y_unc[nfwdtrk]/F");
}

//
void attachToMiniEventTree(TTree *t,MiniEvent_t &ev,bool full)
{
  //event header
  t->SetBranchAddress("isData",    &ev.isData);
  t->SetBranchAddress("run",       &ev.run);
  t->SetBranchAddress("event",     &ev.event);
  t->SetBranchAddress("lumi",      &ev.lumi);

  //generator level event
  t->SetBranchAddress("g_pu",      &ev.g_pu);
  t->SetBranchAddress("g_putrue",  &ev.g_putrue);
  t->SetBranchAddress("g_id1",     &ev.g_id1);
  t->SetBranchAddress("g_id2",     &ev.g_id2);
  t->SetBranchAddress("g_x1",      &ev.g_x1);
  t->SetBranchAddress("g_x2",      &ev.g_x2);
  t->SetBranchAddress("g_qscale",  &ev.g_qscale);
  t->SetBranchAddress("g_nw",      &ev.g_nw);
  t->SetBranchAddress("g_w",       ev.g_w);

  //gen event (jets and dressed leptons)
  t->SetBranchAddress("ng",       &ev.ng);
  t->SetBranchAddress("g_id",      ev.g_id);
  t->SetBranchAddress("g_tagCtrs", ev.g_tagCtrs);
  t->SetBranchAddress("g_bid",     ev.g_bid);
  t->SetBranchAddress("g_isSemiLepBhad", ev.g_isSemiLepBhad);
  t->SetBranchAddress("g_xb",      ev.g_xb);
  t->SetBranchAddress("g_xbp",     ev.g_xbp);
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
  t->SetBranchAddress("rho",      &ev.rho);
  t->SetBranchAddress("triggerBits",        &ev.triggerBits);
  
  //lepton info
  t->SetBranchAddress("nl", &ev.nl);
  t->SetBranchAddress("l_isPromptFinalState",                         ev.l_isPromptFinalState);
  t->SetBranchAddress("l_isDirectPromptTauDecayProductFinalState",    ev.l_isDirectPromptTauDecayProductFinalState);
  t->SetBranchAddress("l_mva",      ev.l_mva);
  t->SetBranchAddress("l_id",       ev.l_id);
  t->SetBranchAddress("l_pid",      ev.l_pid);
  t->SetBranchAddress("l_g",        ev.l_g);
  t->SetBranchAddress("l_charge",   ev.l_charge);
  t->SetBranchAddress("l_pt",       ev.l_pt);
  t->SetBranchAddress("l_eta",      ev.l_eta);
  t->SetBranchAddress("l_phi",      ev.l_phi);
  t->SetBranchAddress("l_mass",     ev.l_mass);
  t->SetBranchAddress("l_scaleUnc",         ev.l_scaleUnc);
  t->SetBranchAddress("l_smearUnc",         ev.l_smearUnc);
  t->SetBranchAddress("l_chargedHadronIso", ev.l_chargedHadronIso);
  t->SetBranchAddress("l_miniIso",          ev.l_miniIso);
  t->SetBranchAddress("l_relIso",           ev.l_relIso);
  t->SetBranchAddress("l_ip3d",             ev.l_ip3d);
  t->SetBranchAddress("l_ip3dsig",          ev.l_ip3dsig);

  //photon info
  t->SetBranchAddress("ngamma", &ev.ngamma);
  t->SetBranchAddress("gamma_isPromptFinalState",  ev.gamma_isPromptFinalState);
  t->SetBranchAddress("gamma_passElectronVeto",    ev.gamma_passElectronVeto);
  t->SetBranchAddress("gamma_hasPixelSeed",        ev.gamma_hasPixelSeed);
  t->SetBranchAddress("gamma_pid",                 ev.gamma_pid);
  t->SetBranchAddress("gamma_g",                   ev.gamma_g);
  t->SetBranchAddress("gamma_mva",                 ev.gamma_mva);
  t->SetBranchAddress("gamma_pt",                  ev.gamma_pt);
  t->SetBranchAddress("gamma_eta",                 ev.gamma_eta);
  t->SetBranchAddress("gamma_phi",                 ev.gamma_phi);
  t->SetBranchAddress("gamma_scaleUnc",            ev.gamma_scaleUnc);
  t->SetBranchAddress("gamma_smearUnc",            ev.gamma_smearUnc);
  t->SetBranchAddress("gamma_chargedHadronIso",    ev.gamma_chargedHadronIso);
  t->SetBranchAddress("gamma_neutralHadronIso",    ev.gamma_neutralHadronIso);
  t->SetBranchAddress("gamma_photonIso",           ev.gamma_photonIso);
  t->SetBranchAddress("gamma_hoe",                 ev.gamma_hoe);
  t->SetBranchAddress("gamma_r9",                  ev.gamma_r9);
  t->SetBranchAddress("gamma_sieie",               ev.gamma_sieie);

  //jet info
  t->SetBranchAddress("nj",        &ev.nj);
  t->SetBranchAddress("j_g",        ev.j_g);
  t->SetBranchAddress("j_area",     ev.j_area);
  t->SetBranchAddress("j_rawsf",    ev.j_rawsf);
  t->SetBranchAddress("j_pt",       ev.j_pt);
  t->SetBranchAddress("j_eta",      ev.j_eta);
  t->SetBranchAddress("j_phi",      ev.j_phi);
  t->SetBranchAddress("j_mass",     ev.j_mass);
  t->SetBranchAddress("j_pumva",    ev.j_pumva);
  t->SetBranchAddress("j_csv",      ev.j_csv);
  t->SetBranchAddress("j_btag",     ev.j_btag);
  t->SetBranchAddress("j_deepcsv",     ev.j_deepcsv);  
  t->SetBranchAddress("j_vtxpx",    ev.j_vtxpx);
  t->SetBranchAddress("j_vtxpy",    ev.j_vtxpy);
  t->SetBranchAddress("j_vtxpz",    ev.j_vtxpz);
  t->SetBranchAddress("j_vtxmass",  ev.j_vtxmass);
  t->SetBranchAddress("j_vtxNtracks",  ev.j_vtxNtracks);
  t->SetBranchAddress("j_vtx3DVal",    ev.j_vtx3DVal);
  t->SetBranchAddress("j_vtx3DSig",    ev.j_vtx3DSig);
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
      t->SetBranchAddress("pf_dxy",      ev.pf_dxy);
      t->SetBranchAddress("pf_dz",       ev.pf_dz);
      //t->SetBranchAddress("pf_dxyUnc",   ev.pf_dxyUnc);
      //t->SetBranchAddress("pf_dzUnc",    ev.pf_dzUnc);
      //t->SetBranchAddress("pf_pvAssoc",  ev.pf_pvAssoc);
      t->SetBranchAddress("pf_puppiWgt", ev.pf_puppiWgt);
    }

  //MET
  t->SetBranchAddress("nmet",      &ev.nmet);
  t->SetBranchAddress("met_pt",    ev.met_pt);
  t->SetBranchAddress("met_phi",   ev.met_phi);
  t->SetBranchAddress("met_sig",   ev.met_sig);
  t->SetBranchAddress("met_filterBits", &ev.met_filterBits);

  //CTPPS local tracks
  t->SetBranchAddress("nfwdtrk",       &ev.nfwdtrk);
  t->SetBranchAddress("fwdtrk_arm",    ev.fwdtrk_arm);
  t->SetBranchAddress("fwdtrk_station",ev.fwdtrk_station);
  t->SetBranchAddress("fwdtrk_pot",    ev.fwdtrk_pot);
  t->SetBranchAddress("fwdtrk_x",      ev.fwdtrk_x);
  t->SetBranchAddress("fwdtrk_x_unc",  ev.fwdtrk_x_unc);
  t->SetBranchAddress("fwdtrk_y",      ev.fwdtrk_y);
  t->SetBranchAddress("fwdtrk_y_unc",  ev.fwdtrk_y_unc);
}
