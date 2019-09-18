#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"

//
void createMiniEventTree(TTree *t,MiniEvent_t &ev,Int_t njecUncs)
{
  //event header
  t->Branch("isData",    &ev.isData,   "isData/O");
  t->Branch("run",       &ev.run,      "run/i");
  t->Branch("event",     &ev.event,    "event/l");
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
  t->Branch("g_npsw",    &ev.g_npsw,   "g_npsw/I");
  t->Branch("g_nup",     &ev.g_nup,    "g_nup/I");
  t->Branch("g_w",        ev.g_w,      "g_w[g_nw]/F");
  t->Branch("g_psw",      ev.g_psw,    "g_psw[g_npsw]/F");

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

  t->Branch("g_nchPV",       &ev.g_nchPV,      "g_nchPV/I");
  t->Branch("g_sumPVChPt",   &ev.g_sumPVChPt,  "g_sumPVChPt/F");
  t->Branch("g_sumPVChPz",   &ev.g_sumPVChPz,  "g_sumPVChPz/F");
  t->Branch("g_sumPVChHt",   &ev.g_sumPVChHt,  "g_sumPVChHt/F");

  //top (lastCopy and pseudo-top)
  t->Branch("ngtop",     &ev.ngtop,      "ngtop/I");
  t->Branch("gtop_id",    ev.gtop_id,    "gtop_id[ngtop]/I");
  t->Branch("gtop_pt",    ev.gtop_pt,    "gtop_pt[ngtop]/F");
  t->Branch("gtop_eta",   ev.gtop_eta,   "gtop_eta[ngtop]/F");
  t->Branch("gtop_phi",   ev.gtop_phi,   "gtop_phi[ngtop]/F");
  t->Branch("gtop_m",     ev.gtop_m,     "gtop_m[ngtop]/F");

  //reco level event
  t->Branch("nvtx",          &ev.nvtx,        "nvtx/I");
  t->Branch("rho",           &ev.rho,         "rho/F");
  t->Branch("triggerBits",   &ev.triggerBits, "triggerBits/I");
  t->Branch("addTriggerBits",   &ev.addTriggerBits, "addTriggerBits/I");

  t->Branch("zeroBiasPS",   &ev.zeroBiasPS, "zeroBiasPS/I");

  //leptons
  t->Branch("nl", &ev.nl, "nl/I");
  t->Branch("l_isPromptFinalState",                         ev.l_isPromptFinalState,                       "l_isPromptFinalState[nl]/O");
  t->Branch("l_isDirectPromptTauDecayProductFinalState",    ev.l_isDirectPromptTauDecayProductFinalState,  "l_isDirectPromptTauDecayProductFinalState[nl]/O");
  t->Branch("l_id",       ev.l_id,      "l_id[nl]/I");
  t->Branch("l_pid",      ev.l_pid,     "l_pid[nl]/I");
  t->Branch("l_g",        ev.l_g,       "l_g[nl]/I");
  t->Branch("l_charge",   ev.l_charge,  "l_charge[nl]/I");
  t->Branch("l_mva",      ev.l_mva,     "l_mva[nl]/F");
  t->Branch("l_mvaCats",  ev.l_mvaCats, "l_mvaCats[nl]/F");
  t->Branch("l_highpt",   ev.l_highpt,  "l_highpt[nl]/F");
  t->Branch("l_pt",       ev.l_pt,      "l_pt[nl]/F");
  t->Branch("l_eta",      ev.l_eta,     "l_eta[nl]/F");
  t->Branch("l_phi",      ev.l_phi,     "l_phi[nl]/F");
  t->Branch("l_mass",     ev.l_mass,    "l_mass[nl]/F");
  t->Branch("l_scaleUnc1",        ev.l_scaleUnc1,        "l_scaleUnc1[nl]/F");
  t->Branch("l_scaleUnc2",        ev.l_scaleUnc2,        "l_scaleUnc2[nl]/F");
  t->Branch("l_scaleUnc3",        ev.l_scaleUnc3,        "l_scaleUnc3[nl]/F");
  t->Branch("l_scaleUnc4",        ev.l_scaleUnc4,        "l_scaleUnc4[nl]/F");
  t->Branch("l_scaleUnc5",        ev.l_scaleUnc5,        "l_scaleUnc5[nl]/F");
  t->Branch("l_scaleUnc6",        ev.l_scaleUnc6,        "l_scaleUnc6[nl]/F");
  t->Branch("l_scaleUnc7",        ev.l_scaleUnc7,        "l_scaleUnc7[nl]/F");
  t->Branch("l_chargedHadronIso", ev.l_chargedHadronIso, "l_chargedHadronIso[nl]/F");
  t->Branch("l_miniIso",          ev.l_miniIso,          "l_miniIso[nl]/F");
  t->Branch("l_relIso",           ev.l_relIso,           "l_relIso[nl]/F");
  t->Branch("l_ip3d",             ev.l_ip3d,             "l_ip3d[nl]/F");
  t->Branch("l_ip3dsig",          ev.l_ip3dsig,          "l_ip3dsig[nl]/F");

  //photons
  t->Branch("ngamma",                   &ev.ngamma,                   "ngamma/I");
  t->Branch("gamma_isPromptFinalState",  ev.gamma_isPromptFinalState, "gamma_isPromptFinalState[ngamma]/O");
  t->Branch("gamma_pid",                 ev.gamma_pid,                "gamma_pid[ngamma]/I");
  t->Branch("gamma_idFlags",             ev.gamma_idFlags,                "gamma_idFlags[ngamma]/I");
  t->Branch("gamma_g",                   ev.gamma_g,                  "gamma_g[ngamma]/I");
  t->Branch("gamma_mva",                 ev.gamma_mva,                "gamma_mva[ngamma]/F");
  t->Branch("gamma_mvaCats",             ev.gamma_mvaCats,            "gamma_mvaCats[ngamma]/F");
  t->Branch("gamma_pt",                  ev.gamma_pt,                 "gamma_pt[ngamma]/F");
  t->Branch("gamma_eta",                 ev.gamma_eta,                "gamma_eta[ngamma]/F");
  t->Branch("gamma_phi",                 ev.gamma_phi,                "gamma_phi[ngamma]/F");  
  t->Branch("gamma_scaleUnc1",           ev.gamma_scaleUnc1,          "gamma_scaleUnc1[ngamma]/F");
  t->Branch("gamma_scaleUnc2",           ev.gamma_scaleUnc2,          "gamma_scaleUnc2[ngamma]/F");
  t->Branch("gamma_scaleUnc3",           ev.gamma_scaleUnc3,          "gamma_scaleUnc3[ngamma]/F");
  t->Branch("gamma_scaleUnc4",           ev.gamma_scaleUnc4,          "gamma_scaleUnc4[ngamma]/F");
  t->Branch("gamma_scaleUnc5",           ev.gamma_scaleUnc5,          "gamma_scaleUnc5[ngamma]/F");
  t->Branch("gamma_scaleUnc6",           ev.gamma_scaleUnc6,          "gamma_scaleUnc6[ngamma]/F");
  t->Branch("gamma_scaleUnc7",           ev.gamma_scaleUnc7,          "gamma_scaleUnc7[ngamma]/F");
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
  t->Branch("j_jerUp",    ev.j_jerUp,   "j_jerUp[nj]/F");
  t->Branch("j_jerDn",    ev.j_jerDn,   "j_jerDn[nj]/F");
  for(int i=0; i<njecUncs; i++) {
    t->Branch(Form("j_jecUp%d",i),    ev.j_jecUp[i],   Form("j_jecUp%d[nj]/F",i));
    t->Branch(Form("j_jecDn%d",i),    ev.j_jecDn[i],   Form("j_jecDn%d[nj]/F",i));
  }
  t->Branch("j_rawsf",    ev.j_rawsf,   "j_rawsf[nj]/F");
  t->Branch("j_pt",       ev.j_pt,      "j_pt[nj]/F");
  t->Branch("j_eta",      ev.j_eta,     "j_eta[nj]/F");
  t->Branch("j_phi",      ev.j_phi,     "j_phi[nj]/F");
  t->Branch("j_mass",     ev.j_mass,    "j_mass[nj]/F");
  t->Branch("j_pumva",    ev.j_pumva,   "j_pumva[nj]/F");
  t->Branch("j_id",       ev.j_id,      "j_id[nj]/I");
  t->Branch("j_csv",      ev.j_csv,     "j_csv[nj]/F");
  t->Branch("j_btag",     ev.j_btag,    "j_btag[nj]/O");
  t->Branch("j_emf",      ev.j_emf,     "j_emf[nj]/F");
  t->Branch("j_qg",       ev.j_qg,      "j_qg[nj]/F");
  t->Branch("j_c2_00",    ev.j_c2_00,   "j_c2_00[nj]/F");
  t->Branch("j_c2_02",    ev.j_c2_02,   "j_c2_02[nj]/F");
  t->Branch("j_c2_05",    ev.j_c2_05,   "j_c2_05[nj]/F");
  t->Branch("j_c2_10",    ev.j_c2_10,   "j_c2_10[nj]/F");
  t->Branch("j_c2_20",    ev.j_c2_20,   "j_c2_20[nj]/F");
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

  //pf sums
  t->Branch("nchPV",        &ev.nchPV,         "nchPV/I");
  t->Branch("sumPVChPt",    &ev.sumPVChPt,     "sumPVChPt/F");
  t->Branch("sumPVChPz",    &ev.sumPVChPz,     "sumPVChPz/F");
  t->Branch("sumPVChHt",    &ev.sumPVChHt,     "sumPVChHt/F");
  t->Branch("nPFCands",     ev.nPFCands,       "nPFCands[8]/I");
  t->Branch("sumPFEn",      ev.sumPFEn,       "sumPFEn[8]/F");
  t->Branch("sumPFPz",      ev.sumPFPz,       "sumPFPz[8]/F");
  t->Branch("sumPFHt",      ev.sumPFHt,       "sumPFHt[8]/F");
  t->Branch("nPFChCands",   ev.nPFChCands,    "nPFChCands[8]/I");
  t->Branch("sumPFChEn",    ev.sumPFChEn,     "sumPFChEn[8]/F");
  t->Branch("sumPFChPz",    ev.sumPFChPz,     "sumPFChPz[8]/F");
  t->Branch("sumPFChHt",    ev.sumPFChHt,     "sumPFChHt[8]/F");

  //MET
  t->Branch("met_pt",      &ev.met_pt,     "met_pt/F");
  t->Branch("met_phi",     &ev.met_phi,    "met_phi/F");
  t->Branch("met_sig",     &ev.met_sig,    "met_sig/F");
  t->Branch("met_ptShifted",   ev.met_ptShifted,    "met_ptShifted[14]/F");
  t->Branch("met_phiShifted",   ev.met_phiShifted,    "met_phiShifted[14]/F");
  t->Branch("met_filterBits", &ev.met_filterBits, "met_filterBits/I");

  //CTPPS local tracks
  t->Branch("nfwdtrk",         &ev.nfwdtrk,          "nfwdtrk/S");
  t->Branch("fwdtrk_pot",       ev.fwdtrk_pot,       "fwdtrk_pot[nfwdtrk]/S");
  t->Branch("fwdtrk_method",    ev.fwdtrk_method,    "fwdtrk_method[nfwdtrk]/S");
  t->Branch("fwdtrk_ex",        ev.fwdtrk_ex,        "fwdtrk_ex[nfwdtrk]/F");
  t->Branch("fwdtrk_ey",        ev.fwdtrk_ey,        "fwdtrk_ey[nfwdtrk]/F");
  t->Branch("fwdtrk_ez",        ev.fwdtrk_ez,        "fwdtrk_ez[nfwdtrk]/F");
  t->Branch("fwdtrk_y",         ev.fwdtrk_y,         "fwdtrk_y[nfwdtrk]/F");
  t->Branch("fwdtrk_chisqnorm", ev.fwdtrk_chisqnorm, "fwdtrk_chisqnorm[nfwdtrk]/F");
  t->Branch("fwdtrk_xi",        ev.fwdtrk_xi,        "fwdtrk_xi[nfwdtrk]/F");
  t->Branch("fwdtrk_t",         ev.fwdtrk_t,         "fwdtrk_t[nfwdtrk]/F");

  t->Branch("nrawmu", &ev.nrawmu, "nrawmu/I");
  t->Branch("rawmu_pt", ev.rawmu_pt, "rawmu_pt[nrawmu]/S");
  t->Branch("rawmu_eta", ev.rawmu_eta, "rawmu_eta[nrawmu]/S");
  t->Branch("rawmu_phi", ev.rawmu_phi, "rawmu_phi[nrawmu]/S");
  t->Branch("rawmu_pid", ev.rawmu_pid, "rawmu_pid[nrawmu]/I");
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
  t->SetBranchAddress("g_npsw",      &ev.g_npsw);
  t->SetBranchAddress("g_psw",       ev.g_psw);
  t->SetBranchAddress("g_nup",     &ev.g_nup);

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

  t->SetBranchAddress("g_nchPV",       &ev.g_nchPV);
  t->SetBranchAddress("g_sumPVChPt",   &ev.g_sumPVChPt);
  t->SetBranchAddress("g_sumPVChPz",   &ev.g_sumPVChPz);
  t->SetBranchAddress("g_sumPVChHt",   &ev.g_sumPVChHt);


  //top (lastCopy and pseudo-top)
  t->SetBranchAddress("ngtop",     &ev.ngtop);
  t->SetBranchAddress("gtop_id",    ev.gtop_id);
  t->SetBranchAddress("gtop_pt",    ev.gtop_pt);
  t->SetBranchAddress("gtop_eta",   ev.gtop_eta);
  t->SetBranchAddress("gtop_phi",   ev.gtop_phi);
  t->SetBranchAddress("gtop_m",     ev.gtop_m);

  //reco level event
  t->SetBranchAddress("nvtx",      &ev.nvtx);
  t->SetBranchAddress("rho",      &ev.rho);
  t->SetBranchAddress("triggerBits",        &ev.triggerBits);
  t->SetBranchAddress("addTriggerBits",        &ev.addTriggerBits);
  if(t->FindBranch("zeroBiasPS"))
    t->SetBranchAddress("zeroBiasPS",   &ev.zeroBiasPS);

  //lepton info
  t->SetBranchAddress("nl", &ev.nl);
  t->SetBranchAddress("l_isPromptFinalState",                         ev.l_isPromptFinalState);
  t->SetBranchAddress("l_isDirectPromptTauDecayProductFinalState",    ev.l_isDirectPromptTauDecayProductFinalState);
  t->SetBranchAddress("l_mva",      ev.l_mva);
  t->SetBranchAddress("l_mvaCats",  ev.l_mvaCats);
  t->SetBranchAddress("l_id",       ev.l_id);
  t->SetBranchAddress("l_pid",      ev.l_pid);
  t->SetBranchAddress("l_g",        ev.l_g);
  t->SetBranchAddress("l_charge",   ev.l_charge);
  t->SetBranchAddress("l_highpt",   ev.l_highpt);
  t->SetBranchAddress("l_pt",       ev.l_pt);
  t->SetBranchAddress("l_eta",      ev.l_eta);
  t->SetBranchAddress("l_phi",      ev.l_phi);
  t->SetBranchAddress("l_mass",     ev.l_mass);
  t->SetBranchAddress("l_scaleUnc1",        ev.l_scaleUnc1);
  t->SetBranchAddress("l_scaleUnc2",        ev.l_scaleUnc2);
  t->SetBranchAddress("l_scaleUnc3",        ev.l_scaleUnc3);
  t->SetBranchAddress("l_scaleUnc4",        ev.l_scaleUnc4);
  t->SetBranchAddress("l_scaleUnc5",        ev.l_scaleUnc5);
  t->SetBranchAddress("l_scaleUnc6",        ev.l_scaleUnc6);
  t->SetBranchAddress("l_scaleUnc7",        ev.l_scaleUnc7);
  t->SetBranchAddress("l_chargedHadronIso", ev.l_chargedHadronIso);
  t->SetBranchAddress("l_miniIso",          ev.l_miniIso);
  t->SetBranchAddress("l_relIso",           ev.l_relIso);
  t->SetBranchAddress("l_ip3d",             ev.l_ip3d);
  t->SetBranchAddress("l_ip3dsig",          ev.l_ip3dsig);

  //photon info
  t->SetBranchAddress("ngamma", &ev.ngamma);
  t->SetBranchAddress("gamma_isPromptFinalState",  ev.gamma_isPromptFinalState);
  t->SetBranchAddress("gamma_pid",                 ev.gamma_pid);
  t->SetBranchAddress("gamma_idFlags",                 ev.gamma_idFlags);
  t->SetBranchAddress("gamma_g",                   ev.gamma_g);
  t->SetBranchAddress("gamma_mva",                 ev.gamma_mva);
  t->SetBranchAddress("gamma_mvaCats",             ev.gamma_mvaCats);
  t->SetBranchAddress("gamma_pt",                  ev.gamma_pt);
  t->SetBranchAddress("gamma_eta",                 ev.gamma_eta);
  t->SetBranchAddress("gamma_phi",                 ev.gamma_phi);
  t->SetBranchAddress("gamma_scaleUnc1",           ev.gamma_scaleUnc1);
  t->SetBranchAddress("gamma_scaleUnc2",           ev.gamma_scaleUnc2);
  t->SetBranchAddress("gamma_scaleUnc3",           ev.gamma_scaleUnc3);
  t->SetBranchAddress("gamma_scaleUnc4",           ev.gamma_scaleUnc4);
  t->SetBranchAddress("gamma_scaleUnc5",           ev.gamma_scaleUnc5);
  t->SetBranchAddress("gamma_scaleUnc6",           ev.gamma_scaleUnc6);
  t->SetBranchAddress("gamma_scaleUnc7",           ev.gamma_scaleUnc7);
  t->SetBranchAddress("gamma_chargedHadronIso",    ev.gamma_chargedHadronIso);
  t->SetBranchAddress("gamma_neutralHadronIso",    ev.gamma_neutralHadronIso);
  t->SetBranchAddress("gamma_photonIso",           ev.gamma_photonIso);
  t->SetBranchAddress("gamma_hoe",                 ev.gamma_hoe);
  t->SetBranchAddress("gamma_r9",                  ev.gamma_r9);
  t->SetBranchAddress("gamma_sieie",               ev.gamma_sieie);

  //jet info
  t->SetBranchAddress("nj",        &ev.nj);
  t->SetBranchAddress("j_g",        ev.j_g);
  t->SetBranchAddress("j_jerUp",    ev.j_jerUp);
  t->SetBranchAddress("j_jerDn",    ev.j_jerDn);
  for(int i=0; i<29; i++) {
    t->SetBranchAddress(Form("j_jecUp%d",i), ev.j_jecUp[i]);
    t->SetBranchAddress(Form("j_jecDn%d",i), ev.j_jecDn[i]);
  }

  t->SetBranchAddress("j_area",     ev.j_area);
  t->SetBranchAddress("j_rawsf",    ev.j_rawsf);
  t->SetBranchAddress("j_pt",       ev.j_pt);
  t->SetBranchAddress("j_eta",      ev.j_eta);
  t->SetBranchAddress("j_phi",      ev.j_phi);
  t->SetBranchAddress("j_mass",     ev.j_mass);
  t->SetBranchAddress("j_pumva",    ev.j_pumva);
  t->SetBranchAddress("j_id",       ev.j_id);
  t->SetBranchAddress("j_csv",      ev.j_csv);
  t->SetBranchAddress("j_btag",     ev.j_btag);
  t->SetBranchAddress("j_emf",      ev.j_emf);
  t->SetBranchAddress("j_qg",       ev.j_qg);
  t->SetBranchAddress("j_c2_00",    ev.j_c2_00);
  t->SetBranchAddress("j_c2_02",    ev.j_c2_02);
  t->SetBranchAddress("j_c2_05",    ev.j_c2_05);
  t->SetBranchAddress("j_c2_10",    ev.j_c2_10);
  t->SetBranchAddress("j_c2_20",    ev.j_c2_20);
  t->SetBranchAddress("j_zg",       ev.j_zg);
  t->SetBranchAddress("j_mult",     ev.j_mult);
  t->SetBranchAddress("j_gaptd",    ev.j_gaptd);
  t->SetBranchAddress("j_gawidth",  ev.j_gawidth);
  t->SetBranchAddress("j_gathrust", ev.j_gathrust);
  t->SetBranchAddress("j_tau32",    ev.j_tau32);
  t->SetBranchAddress("j_tau21",    ev.j_tau21);
  t->SetBranchAddress("j_deepcsv",  ev.j_deepcsv);  
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

  //PF Sums
  t->SetBranchAddress("nchPV",        &ev.nchPV);
  t->SetBranchAddress("sumPVChPt",    &ev.sumPVChPt);
  t->SetBranchAddress("sumPVChPz",    &ev.sumPVChPz);
  t->SetBranchAddress("sumPVChHt",    &ev.sumPVChHt);
  t->SetBranchAddress("nPFCands",     ev.nPFCands);
  t->SetBranchAddress("sumPFEn",      ev.sumPFEn);
  t->SetBranchAddress("sumPFPz",      ev.sumPFPz);
  t->SetBranchAddress("sumPFHt",      ev.sumPFHt);
  t->SetBranchAddress("nPFChCands",   ev.nPFChCands);
  t->SetBranchAddress("sumPFChEn",    ev.sumPFChEn);
  t->SetBranchAddress("sumPFChPz",    ev.sumPFChPz);
  t->SetBranchAddress("sumPFChHt",    ev.sumPFChHt);

  //MET
  t->SetBranchAddress("met_pt",    &ev.met_pt);
  t->SetBranchAddress("met_phi",   &ev.met_phi);
  t->SetBranchAddress("met_sig",   &ev.met_sig);
  t->SetBranchAddress("met_ptShifted",   ev.met_ptShifted);
  t->SetBranchAddress("met_phiShifted",   ev.met_phiShifted);
  t->SetBranchAddress("met_filterBits", &ev.met_filterBits);

  //CTPPS local tracks
  t->SetBranchAddress("nfwdtrk",         &ev.nfwdtrk);
  t->SetBranchAddress("fwdtrk_pot",       ev.fwdtrk_pot);
  t->SetBranchAddress("fwdtrk_method",    ev.fwdtrk_method);
  t->SetBranchAddress("fwdtrk_ex",        ev.fwdtrk_ex);
  t->SetBranchAddress("fwdtrk_ey",        ev.fwdtrk_ey);
  t->SetBranchAddress("fwdtrk_ez",        ev.fwdtrk_ez);
  t->SetBranchAddress("fwdtrk_y",         ev.fwdtrk_y);
  t->SetBranchAddress("fwdtrk_chisqnorm", ev.fwdtrk_chisqnorm);
  t->SetBranchAddress("fwdtrk_xi",        ev.fwdtrk_xi);
  t->SetBranchAddress("fwdtrk_t",         ev.fwdtrk_t);

  //
  t->SetBranchAddress("nrawmu",   &ev.nrawmu);
  t->SetBranchAddress("rawmu_pt",  ev.rawmu_pt);
  t->SetBranchAddress("rawmu_eta", ev.rawmu_eta);
  t->SetBranchAddress("rawmu_phi", ev.rawmu_phi);
  t->SetBranchAddress("rawmu_pid", ev.rawmu_pid);
}
