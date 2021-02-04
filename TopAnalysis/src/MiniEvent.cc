#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"

void createMiniEventTree(TTree *tree, MiniEvent_t & ev, Int_t njecUncs)
{
  //event header
  tree -> Branch("isData",                                    & ev.isData,                                      "isData/O"                                          );
  tree -> Branch("run",                                       & ev.run,                                         "run/i"                                             );
  tree -> Branch("event",                                     & ev.event,                                       "event/l"                                           );
  tree -> Branch("event",                                     & ev.event,                                       "event/L"                                           );
  tree -> Branch("lumi",                                      & ev.lumi,                                        "lumi/i"                                            );
  tree -> Branch("beamXangle",                                & ev.beamXangle,                                  "beamXangle/F"                                      );
  tree -> Branch("instLumi",                                  & ev.instLumi,                                    "instLumi/F"                                        );

  //generator level information
  tree -> Branch("g_pu",                                      & ev.g_pu,                                        "g_pu/I"                                            );
  tree -> Branch("g_putrue",                                  & ev.g_putrue,                                    "g_putrue/I"                                        );
  tree -> Branch("g_id1",                                     & ev.g_id1,                                       "g_id1/I"                                           );
  tree -> Branch("g_id2",                                     & ev.g_id2,                                       "g_id2/I"                                           );
  tree -> Branch("g_x1",                                      & ev.g_x1,                                        "g_x1/F"                                            );
  tree -> Branch("g_x2",                                      & ev.g_x2,                                        "g_x2/F"                                            );
  tree -> Branch("g_qscale",                                  & ev.g_qscale,                                    "g_qscale/F"                                        );
  tree -> Branch("g_nw",                                      & ev.g_nw,                                        "g_nw/I"                                            );
  tree -> Branch("g_npsw",                                    & ev.g_npsw,                                      "g_npsw/I"                                          );
  tree -> Branch("g_nup",                                     & ev.g_nup,                                       "g_nup/I"                                           );
  tree -> Branch("g_w",                                         ev.g_w,                                         "g_w[g_nw]/F"                                       );
  tree -> Branch("g_psw",                                       ev.g_psw,                                       "g_psw[g_npsw]/F"                                   );

  //gen event (jets and dressed leptons)
  tree -> Branch("ng",                                        & ev.ng,                                          "ng/I"                                              );
  tree -> Branch("g_id",                                        ev.g_id,                                        "g_id[ng]/I"                                        );
  tree -> Branch("g_bid",                                       ev.g_bid,                                       "g_bid[ng]/I"                                       );
  tree -> Branch("g_tagCtrs",                                   ev.g_tagCtrs,                                   "g_tagCtrs[ng]/I"                                   );
  tree -> Branch("g_isSemiLepBhad",                             ev.g_isSemiLepBhad,                             "g_isSemiLepBhad[ng]/O"                             );
  tree -> Branch("g_xb",                                        ev.g_xb,                                        "g_xb[ng]/F"                                        );
  tree -> Branch("g_xbp",                                       ev.g_xbp,                                       "g_xbp[ng]/F"                                       );
  tree -> Branch("g_pt",                                        ev.g_pt,                                        "g_pt[ng]/F"                                        );
  tree -> Branch("g_eta",                                       ev.g_eta,                                       "g_eta[ng]/F"                                       );
  tree -> Branch("g_phi",                                       ev.g_phi,                                       "g_phi[ng]/F"                                       );
  tree -> Branch("g_m",                                         ev.g_m,                                         "g_m[ng]/F"                                         );

  tree -> Branch("g_nchPV",                                   & ev.g_nchPV,                                     "g_nchPV/I"                                         );
  tree -> Branch("g_sumPVChPt",                               & ev.g_sumPVChPt,                                 "g_sumPVChPt/F"                                     );
  tree -> Branch("g_sumPVChPz",                               & ev.g_sumPVChPz,                                 "g_sumPVChPz/F"                                     );
  tree -> Branch("g_sumPVChHt",                               & ev.g_sumPVChHt,                                 "g_sumPVChHt/F"                                     );
 
  //top (lastCopy and pseudo-top)
  tree -> Branch("ngtop",                                     & ev.ngtop,                                       "ngtop/I"                                           );
  tree -> Branch("gtop_id",                                     ev.gtop_id,                                     "gtop_id[ngtop]/I"                                  );
  tree -> Branch("gtop_pt",                                     ev.gtop_pt,                                     "gtop_pt[ngtop]/F"                                  );
  tree -> Branch("gtop_eta",                                    ev.gtop_eta,                                    "gtop_eta[ngtop]/F"                                 );
  tree -> Branch("gtop_phi",                                    ev.gtop_phi,                                    "gtop_phi[ngtop]/F"                                 );
  tree -> Branch("gtop_m",                                      ev.gtop_m,                                      "gtop_m[ngtop]/F"                                   );

  //final state particles                                                                                                                                                                                   
  tree -> Branch("ngpf",                                      & ev.ngpf,                                        "ngpf/I"                                            );
  tree -> Branch("gpf_id",                                      ev.gpf_id,                                      "gpf_id[ngpf]/I"                                    );
  tree -> Branch("gpf_c",                                       ev.gpf_c,                                       "gpf_c[ngpf]/I"                                     );
  tree -> Branch("gpf_g",                                       ev.gpf_g,                                       "gpf_g[ngpf]/I"                                     );
  tree -> Branch("gpf_pt",                                      ev.gpf_pt,                                      "gpf_pt[ngpf]/F"                                    );
  tree -> Branch("gpf_eta",                                     ev.gpf_eta,                                     "gpf_eta[ngpf]/F"                                   );
  tree -> Branch("gpf_phi",                                     ev.gpf_phi,                                     "gpf_phi[ngpf]/F"                                   );
  tree -> Branch("gpf_m",                                       ev.gpf_m,                                       "gpf_m[ngpf]/F"                                     );


  //reco level event
  tree -> Branch("nvtx",                                       & ev.nvtx,                                       "nvtx/I"                                            );
  tree -> Branch("rho",                                        & ev.rho,                                        "rho/F"                                             );
  tree -> Branch("triggerBits",                                & ev.triggerBits,                                "triggerBits/I"                                     );
  tree -> Branch("addTriggerBits",                             & ev.addTriggerBits,                             "addTriggerBits/I"                                  );

  tree -> Branch("zeroBiasPS",                                 & ev.zeroBiasPS,                                 "zeroBiasPS/I"                                      );

  //leptons
  tree -> Branch("nl",                                         & ev.nl,                                         "nl/I"                                              );
  tree -> Branch("l_isPromptFinalState",                         ev.l_isPromptFinalState,                       "l_isPromptFinalState[nl]/O"                        );
  tree -> Branch("l_isDirectPromptTauDecayProductFinalState",    ev.l_isDirectPromptTauDecayProductFinalState,  "l_isDirectPromptTauDecayProductFinalState[nl]/O"   );
  tree -> Branch("l_id",                                         ev.l_id,                                       "l_id[nl]/I"                                        );
  tree -> Branch("l_pid",                                        ev.l_pid,                                      "l_pid[nl]/I"                                       );
  tree -> Branch("l_g",                                          ev.l_g,                                        "l_g[nl]/I"                                         );
  tree -> Branch("l_charge",                                     ev.l_charge,                                   "l_charge[nl]/I"                                    );
  tree -> Branch("l_mva",                                        ev.l_mva,                                      "l_mva[nl]/F"                                       );
  tree -> Branch("l_mvaCats",                                    ev.l_mvaCats,                                  "l_mvaCats[nl]/F"                                   );
  tree -> Branch("l_highpt",                                     ev.l_highpt,                                   "l_highpt[nl]/F"                                    );
  tree -> Branch("l_pt",                                         ev.l_pt,                                       "l_pt[nl]/F"                                        );
  tree -> Branch("l_eta",                                        ev.l_eta,                                      "l_eta[nl]/F"                                       );
  tree -> Branch("l_phi",                                        ev.l_phi,                                      "l_phi[nl]/F"                                       );
  tree -> Branch("l_mass",                                       ev.l_mass,                                     "l_mass[nl]/F"                                      );
  tree -> Branch("l_scaleUnc1",                                  ev.l_scaleUnc1,                                "l_scaleUnc1[nl]/F"                                 );
  tree -> Branch("l_scaleUnc2",                                  ev.l_scaleUnc2,                                "l_scaleUnc2[nl]/F"                                 );
  tree -> Branch("l_scaleUnc3",                                  ev.l_scaleUnc3,                                "l_scaleUnc3[nl]/F"                                 );
  tree -> Branch("l_scaleUnc4",                                  ev.l_scaleUnc4,                                "l_scaleUnc4[nl]/F"                                 );
  tree -> Branch("l_scaleUnc5",                                  ev.l_scaleUnc5,                                "l_scaleUnc5[nl]/F"                                 );
  tree -> Branch("l_scaleUnc6",                                  ev.l_scaleUnc6,                                "l_scaleUnc6[nl]/F"                                 );
  tree -> Branch("l_scaleUnc7",                                  ev.l_scaleUnc7,                                "l_scaleUnc7[nl]/F"                                 );
  tree -> Branch("l_chargedHadronIso",                           ev.l_chargedHadronIso,                         "l_chargedHadronIso[nl]/F"                          );
  tree -> Branch("l_miniIso",                                    ev.l_miniIso,                                  "l_miniIso[nl]/F"                                   );
  tree -> Branch("l_relIso",                                     ev.l_relIso,                                   "l_relIso[nl]/F"                                    );
  tree -> Branch("l_ip3d",                                       ev.l_ip3d,                                     "l_ip3d[nl]/F"                                      );
  tree -> Branch("l_ip3dsig",                                    ev.l_ip3dsig,                                  "l_ip3dsig[nl]/F"                                   );

  //photons
  tree -> Branch("ngamma",                                     & ev.ngamma,                                     "ngamma/I"                                          );
  tree -> Branch("gamma_isPromptFinalState",                     ev.gamma_isPromptFinalState,                   "gamma_isPromptFinalState[ngamma]/O"                );
  tree -> Branch("gamma_pid",                                    ev.gamma_pid,                                  "gamma_pid[ngamma]/I"                               );
  tree -> Branch("gamma_idFlags",                                ev.gamma_idFlags,                              "gamma_idFlags[ngamma]/I"                           );
  tree -> Branch("gamma_g",                                      ev.gamma_g,                                    "gamma_g[ngamma]/I"                                 );
  tree -> Branch("gamma_mva",                                    ev.gamma_mva,                                  "gamma_mva[ngamma]/F"                               );
  tree -> Branch("gamma_mvaCats",                                ev.gamma_mvaCats,                              "gamma_mvaCats[ngamma]/F"                           );
  tree -> Branch("gamma_pt",                                     ev.gamma_pt,                                   "gamma_pt[ngamma]/F"                                );
  tree -> Branch("gamma_eta",                                    ev.gamma_eta,                                  "gamma_eta[ngamma]/F"                               );
  tree -> Branch("gamma_phi",                                    ev.gamma_phi,                                  "gamma_phi[ngamma]/F"                               );
  tree -> Branch("gamma_scaleUnc1",                              ev.gamma_scaleUnc1,                            "gamma_scaleUnc1[ngamma]/F"                         );
  tree -> Branch("gamma_scaleUnc2",                              ev.gamma_scaleUnc2,                            "gamma_scaleUnc2[ngamma]/F"                         );
  tree -> Branch("gamma_scaleUnc3",                              ev.gamma_scaleUnc3,                            "gamma_scaleUnc3[ngamma]/F"                         );
  tree -> Branch("gamma_scaleUnc4",                              ev.gamma_scaleUnc4,                            "gamma_scaleUnc4[ngamma]/F"                         );
  tree -> Branch("gamma_scaleUnc5",                              ev.gamma_scaleUnc5,                            "gamma_scaleUnc5[ngamma]/F"                         );
  tree -> Branch("gamma_scaleUnc6",                              ev.gamma_scaleUnc6,                            "gamma_scaleUnc6[ngamma]/F"                         );
  tree -> Branch("gamma_scaleUnc7",                              ev.gamma_scaleUnc7,                            "gamma_scaleUnc7[ngamma]/F"                         );
  tree -> Branch("gamma_chargedHadronIso",                       ev.gamma_chargedHadronIso,                     "gamma_chargedHadronIso[ngamma]/F"                  );
  tree -> Branch("gamma_neutralHadronIso",                       ev.gamma_neutralHadronIso,                     "gamma_neutralHadronIso[ngamma]/F"                  );
  tree -> Branch("gamma_photonIso",                              ev.gamma_photonIso,                            "gamma_photonIso[ngamma]/F"                         );
  tree -> Branch("gamma_hoe",                                    ev.gamma_hoe,                                  "gamma_hoe[ngamma]/F"                               );
  tree -> Branch("gamma_r9",                                     ev.gamma_r9,                                   "gamma_r9[ngamma]/F"                                );
  tree -> Branch("gamma_sieie",                                  ev.gamma_sieie,                                "gamma_sieie[ngamma]/F"                             );

  //jet info
  tree -> Branch("nj",                                         & ev.nj,                                         "nj/I"                                              );
  tree -> Branch("j_g",                                          ev.j_g,                                        "j_g[nj]/I"                                         );
  tree -> Branch("j_area",                                       ev.j_area,                                     "j_area[nj]/F"                                      );
  tree -> Branch("j_jerUp",                                      ev.j_jerUp,                                    "j_jerUp[nj]/F"                                     );
  tree -> Branch("j_jerDn",                                      ev.j_jerDn,                                    "j_jerDn[nj]/F"                                     );
  for(unsigned int ind = 0; ind < (unsigned int ) njecUncs; ind ++) 
    {
      tree -> Branch(Form("j_jecUp%d", ind),                         ev.j_jecUp[ind],                               Form("j_jecUp%d[nj]/F", ind)                        );
      tree -> Branch(Form("j_jecDn%d", ind),                         ev.j_jecDn[ind],                               Form("j_jecDn%d[nj]/F", ind)                        );
    }
  tree -> Branch("j_rawsf",                                      ev.j_rawsf,                                    "j_rawsf[nj]/F"                                     );
  tree -> Branch("j_pt",                                         ev.j_pt,                                       "j_pt[nj]/F"                                        );
  tree -> Branch("j_eta",                                        ev.j_eta,                                      "j_eta[nj]/F"                                       );
  tree -> Branch("j_phi",                                        ev.j_phi,                                      "j_phi[nj]/F"                                       );
  tree -> Branch("j_mass",                                       ev.j_mass,                                     "j_mass[nj]/F"                                      );
  tree -> Branch("j_pumva",                                      ev.j_pumva,                                    "j_pumva[nj]/F"                                     );
  tree -> Branch("j_id",                                         ev.j_id,                                       "j_id[nj]/I"                                        );
  tree -> Branch("j_csv",                                        ev.j_csv,                                      "j_csv[nj]/F"                                       );
  tree -> Branch("j_btag",                                       ev.j_btag,                                     "j_btag[nj]/O"                                      );
  tree -> Branch("j_emf",                                        ev.j_emf,                                      "j_emf[nj]/F"                                       );
  tree -> Branch("j_qg",                                         ev.j_qg,                                       "j_qg[nj]/F"                                        );
  tree -> Branch("j_c2_00",                                      ev.j_c2_00,                                    "j_c2_00[nj]/F"                                     );
  tree -> Branch("j_c2_02",                                      ev.j_c2_02,                                    "j_c2_02[nj]/F"                                     );
  tree -> Branch("j_c2_05",                                      ev.j_c2_05,                                    "j_c2_05[nj]/F"                                     );
  tree -> Branch("j_c2_10",                                      ev.j_c2_10,                                    "j_c2_10[nj]/F"                                     );
  tree -> Branch("j_c2_20",                                      ev.j_c2_20,                                    "j_c2_20[nj]/F"                                     );
  tree -> Branch("j_zg",                                         ev.j_zg,                                       "j_zg[nj]/F"                                        );
  tree -> Branch("j_mult",                                       ev.j_mult,                                     "j_mult[nj]/F"                                      );
  tree -> Branch("j_gaptd",                                      ev.j_gaptd,                                    "j_gaptd[nj]/F"                                     );
  tree -> Branch("j_gawidth",                                    ev.j_gawidth,                                  "j_gawidth[nj]/F"                                   );
  tree -> Branch("j_gathrust",                                   ev.j_gathrust,                                 "j_gathrust[nj]/F"                                  );
  tree -> Branch("j_tau32",                                      ev.j_tau32,                                    "j_tau32[nj]/F"                                     );
  tree -> Branch("j_tau21",                                      ev.j_tau21,                                    "j_tau21[nj]/F"                                     );
  tree -> Branch("j_deepcsv",                                    ev.j_deepcsv,                                  "j_deepcsv[nj]/F"                                   );
  tree -> Branch("j_probc",                                      ev.j_probc,                                    "j_probc[nj]/F"                                     );
  tree -> Branch("j_probudsg",                                   ev.j_probudsg,                                 "j_probudsg[nj]/F"                                  );
  tree -> Branch("j_probb",                                      ev.j_probb,                                    "j_probb[nj]/F"                                     );
  tree -> Branch("j_probbb",                                     ev.j_probbb,                                   "j_probbb[nj]/F"                                    );
  tree -> Branch("j_CvsL",                                       ev.j_CvsL,                                     "j_CvsL[nj]/F"                                      );
  tree -> Branch("j_CvsB",                                       ev.j_CvsB,                                     "j_CvsB[nj]/F"                                      );
  tree -> Branch("j_vtxpx",                                      ev.j_vtxpx,                                    "j_vtxpx[nj]/F"                                     );
  tree -> Branch("j_vtxpy",                                      ev.j_vtxpy,                                    "j_vtxpy[nj]/F"                                     );
  tree -> Branch("j_vtxpz",                                      ev.j_vtxpz,                                    "j_vtxpz[nj]/F"                                     );
  tree -> Branch("j_vtxmass",                                    ev.j_vtxmass,                                  "j_vtxmass[nj]/F"                                   );
  tree -> Branch("j_vtxNtracks",                                 ev.j_vtxNtracks,                               "j_vtxNtracks[nj]/I"                                );
  tree -> Branch("j_vtx3DVal",                                   ev.j_vtx3DVal,                                 "j_vtx3DVal[nj]/F"                                  );
  tree -> Branch("j_vtx3DSig",                                   ev.j_vtx3DSig,                                 "j_vtx3DSig[nj]/F"                                  );
  tree -> Branch("j_flav",                                       ev.j_flav,                                     "j_flav[nj]/I"                                      );
  tree -> Branch("j_hadflav",                                    ev.j_hadflav,                                  "j_hadflav[nj]/I"                                   );
  tree -> Branch("j_pid",                                        ev.j_pid,                                      "j_pid[nj]/I"                                       );

  //pf candidates (only charged if outside jets)                                                                                                                                                            
  tree -> Branch("npf",                                        & ev.npf,                                        "npf/I"                                             );
  tree -> Branch("pf_j",                                         ev.pf_j,                                       "pf_j[npf]/I"                                       );
  tree -> Branch("pf_id",                                        ev.pf_id,                                      "pf_id[npf]/I"                                      );
  tree -> Branch("pf_c",                                         ev.pf_c,                                       "pf_c[npf]/I"                                       );
  tree -> Branch("pf_pt",                                        ev.pf_pt,                                      "pf_pt[npf]/F"                                      );
  tree -> Branch("pf_eta",                                       ev.pf_eta,                                     "pf_eta[npf]/F"                                     );
  tree -> Branch("pf_phi",                                       ev.pf_phi,                                     "pf_phi[npf]/F"                                     );
  tree -> Branch("pf_m",                                         ev.pf_m,                                       "pf_m[npf]/F"                                       );


  //pf sums
  tree -> Branch("nchPV",                                      & ev.nchPV,                                      "nchPV/I"                                           );
  tree -> Branch("sumPVChPt",                                  & ev.sumPVChPt,                                  "sumPVChPt/F"                                       );
  tree -> Branch("sumPVChPz",                                  & ev.sumPVChPz,                                  "sumPVChPz/F"                                       );
  tree -> Branch("sumPVChHt",                                  & ev.sumPVChHt,                                  "sumPVChHt/F"                                       );
  tree -> Branch("nPFCands",                                     ev.nPFCands,                                   "nPFCands[8]/I"                                     );
  tree -> Branch("sumPFEn",                                      ev.sumPFEn,                                    "sumPFEn[8]/F"                                      );
  tree -> Branch("sumPFPz",                                      ev.sumPFPz,                                    "sumPFPz[8]/F"                                      );
  tree -> Branch("sumPFHt",                                      ev.sumPFHt,                                    "sumPFHt[8]/F"                                      );
  tree -> Branch("nPFChCands",                                   ev.nPFChCands,                                 "nPFChCands[8]/I"                                   );
  tree -> Branch("sumPFChEn",                                    ev.sumPFChEn,                                  "sumPFChEn[8]/F"                                    );
  tree -> Branch("sumPFChPz",                                    ev.sumPFChPz,                                  "sumPFChPz[8]/F"                                    );
  tree -> Branch("sumPFChHt",                                    ev.sumPFChHt,                                  "sumPFChHt[8]/F"                                    );

  //MET
  tree -> Branch("nmet",                                       & ev.nmet,                                       "nmet/I"                                            );
  tree -> Branch("met_pt",                                     & ev.met_pt,                                     "met_pt/F"                                          );
  tree -> Branch("met_phi",                                    & ev.met_phi,                                    "met_phi/F"                                         );
  tree -> Branch("met_sig",                                    & ev.met_sig,                                    "met_sig/F"                                         );
  tree -> Branch("met_ptShifted",                                ev.met_ptShifted,                              "met_ptShifted[14]/F"                               );
  tree -> Branch("met_phiShifted",                               ev.met_phiShifted,                             "met_phiShifted[14]/F"                              );
  tree -> Branch("met_filterBits",                             & ev.met_filterBits,                             "met_filterBits/I"                                  );

  //CTPPS local tracks
  tree -> Branch("nppstrk",                                    & ev.nppstrk,                                    "nppstrk/S"                                         );
  tree -> Branch("ppstrk_pot",                                   ev.ppstrk_pot,                                 "ppstrk_pot[nppstrk]/S"                             );
  tree -> Branch("ppstrk_x",                                     ev.ppstrk_x,                                   "ppstrk_x[nppstrk]/F"                               );
  tree -> Branch("ppstrk_y",                                     ev.ppstrk_y,                                   "ppstrk_y[nppstrk]/F"                               );
  tree -> Branch("ppstrk_xUnc",                                  ev.ppstrk_xUnc,                                "ppstrk_xUnc[nppstrk]/F"                            );
  tree -> Branch("ppstrk_yUnc",                                  ev.ppstrk_yUnc,                                "ppstrk_yUnc[nppstrk]/F"                            );
  tree -> Branch("ppstrk_tx",                                    ev.ppstrk_tx,                                  "ppstrk_tx[nppstrk]/F"                              );
  tree -> Branch("ppstrk_ty",                                    ev.ppstrk_ty,                                  "ppstrk_ty[nppstrk]/F"                              );
  tree -> Branch("ppstrk_txUnc",                                 ev.ppstrk_txUnc,                               "ppstrk_txUnc[nppstrk]/F"                           );
  tree -> Branch("ppstrk_tyUnc",                                 ev.ppstrk_tyUnc,                               "ppstrk_tyUnc[nppstrk]/F"                           );
  tree -> Branch("ppstrk_recoInfo",                              ev.ppstrk_recoInfo,                            "ppstrk_recoInfo[nppstrk]/S"                        );
  tree -> Branch("ppstrk_chisqnorm",                             ev.ppstrk_chisqnorm,                           "ppstrk_chisqnorm[nppstrk]/F"                       );
  tree -> Branch("ppstrk_t",                                     ev.ppstrk_t,                                   "ppstrk_t[nppstrk]/F"                               );
  tree -> Branch("ppstrk_tUnc",                                  ev.ppstrk_tUnc,                                "ppstrk_tUnc[nppstrk]/F"                            );
  tree -> Branch("nfwdtrk",                                    & ev.nfwdtrk,                                    "nfwdtrk/S"                                         );
  tree -> Branch("fwdtrk_pot",                                   ev.fwdtrk_pot,                                 "fwdtrk_pot[nfwdtrk]/S"                             );
  tree -> Branch("fwdtrk_method",                                ev.fwdtrk_method,                              "fwdtrk_method[nfwdtrk]/S"                          );
  tree -> Branch("fwdtrk_recoInfo",                              ev.fwdtrk_recoInfo,                            "fwdtrk_recoInfo[nfwdtrk]/S"                        );
  tree -> Branch("fwdtrk_thetax",                                ev.fwdtrk_thetax,                              "fwdtrk_thetax[nfwdtrk]/F"                          );
  tree -> Branch("fwdtrk_thetay",                                ev.fwdtrk_thetay,                              "fwdtrk_thetay[nfwdtrk]/F"                          );
  tree -> Branch("fwdtrk_x",                                     ev.fwdtrk_x,                                   "fwdtrk_x[nfwdtrk]/F"                               );
  tree -> Branch("fwdtrk_y",                                     ev.fwdtrk_y,                                   "fwdtrk_y[nfwdtrk]/F"                               );
  tree -> Branch("fwdtrk_tx",                                    ev.fwdtrk_tx,                                  "fwdtrk_tx[nfwdtrk]/F"                              );
  tree -> Branch("fwdtrk_ty",                                    ev.fwdtrk_ty,                                  "fwdtrk_ty[nfwdtrk]/F"                              );
  tree -> Branch("fwdtrk_vx",                                    ev.fwdtrk_vx,                                  "fwdtrk_vx[nfwdtrk]/F"                              );
  tree -> Branch("fwdtrk_vy",                                    ev.fwdtrk_vy,                                  "fwdtrk_vy[nfwdtrk]/F"                              );
  tree -> Branch("fwdtrk_vz",                                    ev.fwdtrk_vz,                                  "fwdtrk_vz[nfwdtrk]/F"                              );
  tree -> Branch("fwdtrk_time",                                  ev.fwdtrk_time,                                "fwdtrk_time[nfwdtrk]/F"                            );
  tree -> Branch("fwdtrk_timeError",                             ev.fwdtrk_timeError,                           "fwdtrk_timeError[nfwdtrk]/F"                       );
  tree -> Branch("fwdtrk_chisqnorm",                             ev.fwdtrk_chisqnorm,                           "fwdtrk_chisqnorm[nfwdtrk]/F"                       );
  tree -> Branch("fwdtrk_xi",                                    ev.fwdtrk_xi,                                  "fwdtrk_xi[nfwdtrk]/F"                              );
  tree -> Branch("fwdtrk_xiError",                               ev.fwdtrk_xiError,                             "fwdtrk_xiError[nfwdtrk]/F"                         );
  tree -> Branch("fwdtrk_t",                                     ev.fwdtrk_t,                                   "fwdtrk_t[nfwdtrk]/F"                               );

  tree -> Branch("nrawmu",                                     & ev.nrawmu,                                     "nrawmu/I"                                          );
  tree -> Branch("rawmu_pt",                                     ev.rawmu_pt,                                   "rawmu_pt[nrawmu]/S"                                );
  tree -> Branch("rawmu_eta",                                    ev.rawmu_eta,                                  "rawmu_eta[nrawmu]/S"                               );
  tree -> Branch("rawmu_phi",                                    ev.rawmu_phi,                                  "rawmu_phi[nrawmu]/S"                               );
  tree -> Branch("rawmu_pid",                                    ev.rawmu_pid,                                  "rawmu_pid[nrawmu]/I"                               );
}

//
void attachToMiniEventTree(TTree * tree, MiniEvent_t & ev, bool full)
{
  //event header
  tree -> SetBranchAddress("isData",                                      & ev.isData                                     );
  tree -> SetBranchAddress("run",                                         & ev.run                                        );
  tree -> SetBranchAddress("event",                                       & ev.event                                      );
  tree -> SetBranchAddress("event",                                       & ev.event                                      );
  tree -> SetBranchAddress("lumi",                                        & ev.lumi                                       );
  // tree -> SetBranchAddress("beamXangle",                                  & ev.beamXangle                                 );
  // tree -> SetBranchAddress("instLumi",                                    & ev.instLumi                                   );

  //generator level event
  tree -> SetBranchAddress("g_pu",                                        & ev.g_pu                                       );
  // tree -> SetBranchAddress("g_putrue",                                    & ev.g_putrue                                   );
  tree -> SetBranchAddress("g_id1",                                       & ev.g_id1                                      );
  tree -> SetBranchAddress("g_id2",                                       & ev.g_id2                                      );
  tree -> SetBranchAddress("g_x1",                                        & ev.g_x1                                       );
  tree -> SetBranchAddress("g_x2",                                        & ev.g_x2                                       );
  tree -> SetBranchAddress("g_qscale",                                    & ev.g_qscale                                   );
  tree -> SetBranchAddress("g_nw",                                        & ev.g_nw                                       );
  tree -> SetBranchAddress("g_w",                                           ev.g_w                                        );
  // tree -> SetBranchAddress("g_npsw",                                      & ev.g_npsw                                     );
  // tree -> SetBranchAddress("g_psw",                                         ev.g_psw                                      );
  // tree -> SetBranchAddress("g_nup",                                       & ev.g_nup                                      );

  //gen event (jets and dressed leptons)
  tree -> SetBranchAddress("ng",                                          & ev.ng                                         );
  tree -> SetBranchAddress("g_id",                                          ev.g_id                                       );
  tree -> SetBranchAddress("g_tagCtrs",                                     ev.g_tagCtrs                                  );
  tree -> SetBranchAddress("g_bid",                                         ev.g_bid                                      );
  tree -> SetBranchAddress("g_isSemiLepBhad",                               ev.g_isSemiLepBhad                            );
  tree -> SetBranchAddress("g_xb",                                          ev.g_xb                                       );
  tree -> SetBranchAddress("g_xbp",                                         ev.g_xbp                                      );
  tree -> SetBranchAddress("g_pt",                                          ev.g_pt                                       );
  tree -> SetBranchAddress("g_eta",                                         ev.g_eta                                      );
  tree -> SetBranchAddress("g_phi",                                         ev.g_phi                                      );
  tree -> SetBranchAddress("g_m",                                           ev.g_m                                        );

  tree -> SetBranchAddress("ngpf",                                        & ev.ngpf                                       );
  tree -> SetBranchAddress("gpf_id",                                        ev.gpf_id                                     );
  tree -> SetBranchAddress("gpf_c",                                         ev.gpf_c                                      );
  tree -> SetBranchAddress("gpf_g",                                         ev.gpf_g                                      );
  tree -> SetBranchAddress("gpf_pt",                                        ev.gpf_pt                                     );
  tree -> SetBranchAddress("gpf_eta",                                       ev.gpf_eta                                    );
  tree -> SetBranchAddress("gpf_phi",                                       ev.gpf_phi                                    );
  tree -> SetBranchAddress("gpf_m",                                         ev.gpf_m                                      );



  // tree -> SetBranchAddress("g_nchPV",                                     & ev.g_nchPV                                    );
  // tree -> SetBranchAddress("g_sumPVChPt",                                 & ev.g_sumPVChPt                                );
  // tree -> SetBranchAddress("g_sumPVChPz",                                 & ev.g_sumPVChPz                                );
  // tree -> SetBranchAddress("g_sumPVChHt",                                 & ev.g_sumPVChHt                                );


  //top (lastCopy and pseudo-top)
  tree -> SetBranchAddress("ngtop",                                       & ev.ngtop                                      );
  tree -> SetBranchAddress("gtop_id",                                       ev.gtop_id                                    );
  tree -> SetBranchAddress("gtop_pt",                                       ev.gtop_pt                                    );
  tree -> SetBranchAddress("gtop_eta",                                      ev.gtop_eta                                   );
  tree -> SetBranchAddress("gtop_phi",                                      ev.gtop_phi                                   );
  tree -> SetBranchAddress("gtop_m",                                        ev.gtop_m                                     );

  //reco level event
  tree -> SetBranchAddress("nvtx",                                        & ev.nvtx                                       );
  tree -> SetBranchAddress("rho",                                         & ev.rho                                        );
  tree -> SetBranchAddress("triggerBits",                                 & ev.triggerBits                                );
  // tree -> SetBranchAddress("addTriggerBits",                              & ev.addTriggerBits                             );
  if(tree -> FindBranch("zeroBiasPS"))
    tree -> SetBranchAddress("zeroBiasPS",                                  & ev.zeroBiasPS                                 );

  //lepton info
  tree -> SetBranchAddress("nl",                                          & ev.nl                                         );
  tree -> SetBranchAddress("l_isPromptFinalState",                          ev.l_isPromptFinalState                       );
  tree -> SetBranchAddress("l_isDirectPromptTauDecayProductFinalState",     ev.l_isDirectPromptTauDecayProductFinalState  );
  // tree -> SetBranchAddress("l_mva",                                         ev.l_mva                                      );
  // tree -> SetBranchAddress("l_mvaCats",                                     ev.l_mvaCats                                  );
  tree -> SetBranchAddress("l_id",                                          ev.l_id                                       );
  tree -> SetBranchAddress("l_pid",                                         ev.l_pid                                      );
  tree -> SetBranchAddress("l_g",                                           ev.l_g                                        );
  tree -> SetBranchAddress("l_charge",                                      ev.l_charge                                   );
  // tree -> SetBranchAddress("l_highpt",                                      ev.l_highpt                                   );
  tree -> SetBranchAddress("l_pt",                                          ev.l_pt                                       );
  tree -> SetBranchAddress("l_eta",                                         ev.l_eta                                      );
  tree -> SetBranchAddress("l_phi",                                         ev.l_phi                                      );
  tree -> SetBranchAddress("l_mass",                                        ev.l_mass                                     );
  // tree -> SetBranchAddress("l_scaleUnc1",                                   ev.l_scaleUnc1                                );
  // tree -> SetBranchAddress("l_scaleUnc2",                                   ev.l_scaleUnc2                                );
  // tree -> SetBranchAddress("l_scaleUnc3",                                   ev.l_scaleUnc3                                );
  // tree -> SetBranchAddress("l_scaleUnc4",                                   ev.l_scaleUnc4                                );
  // tree -> SetBranchAddress("l_scaleUnc5",                                   ev.l_scaleUnc5                                );
  // tree -> SetBranchAddress("l_scaleUnc6",                                   ev.l_scaleUnc6                                );
  // tree -> SetBranchAddress("l_scaleUnc7",                                   ev.l_scaleUnc7                                );
  tree -> SetBranchAddress("l_chargedHadronIso",                            ev.l_chargedHadronIso                         );
  tree -> SetBranchAddress("l_miniIso",                                     ev.l_miniIso                                  );
  tree -> SetBranchAddress("l_relIso",                                      ev.l_relIso                                   );
  tree -> SetBranchAddress("l_ip3d",                                        ev.l_ip3d                                     );
  tree -> SetBranchAddress("l_ip3dsig",                                     ev.l_ip3dsig                                  );

  //photon info
  // tree -> SetBranchAddress("ngamma",                                      & ev.ngamma                                     );
  // tree -> SetBranchAddress("gamma_isPromptFinalState",                      ev.gamma_isPromptFinalState                   );
  // tree -> SetBranchAddress("gamma_pid",                                     ev.gamma_pid                                  );
  // tree -> SetBranchAddress("gamma_idFlags",                                 ev.gamma_idFlags                              );
  // tree -> SetBranchAddress("gamma_g",                                       ev.gamma_g                                    );
  // tree -> SetBranchAddress("gamma_mva",                                     ev.gamma_mva                                  );
  // tree -> SetBranchAddress("gamma_mvaCats",                                 ev.gamma_mvaCats                              );
  // tree -> SetBranchAddress("gamma_pt",                                      ev.gamma_pt                                   );
  // tree -> SetBranchAddress("gamma_eta",                                     ev.gamma_eta                                  );
  // tree -> SetBranchAddress("gamma_phi",                                     ev.gamma_phi                                  );
  // tree -> SetBranchAddress("gamma_scaleUnc1",                               ev.gamma_scaleUnc1                            );
  // tree -> SetBranchAddress("gamma_scaleUnc2",                               ev.gamma_scaleUnc2                            );
  // tree -> SetBranchAddress("gamma_scaleUnc3",                               ev.gamma_scaleUnc3                            );
  // tree -> SetBranchAddress("gamma_scaleUnc4",                               ev.gamma_scaleUnc4                            );
  // tree -> SetBranchAddress("gamma_scaleUnc5",                               ev.gamma_scaleUnc5                            );
  // tree -> SetBranchAddress("gamma_scaleUnc6",                               ev.gamma_scaleUnc6                            );
  // tree -> SetBranchAddress("gamma_scaleUnc7",                               ev.gamma_scaleUnc7                            );
  // tree -> SetBranchAddress("gamma_chargedHadronIso",                        ev.gamma_chargedHadronIso                     );
  // tree -> SetBranchAddress("gamma_neutralHadronIso",                        ev.gamma_neutralHadronIso                     );
  // tree -> SetBranchAddress("gamma_photonIso",                               ev.gamma_photonIso                            );
  // tree -> SetBranchAddress("gamma_hoe",                                     ev.gamma_hoe                                  );
  // tree -> SetBranchAddress("gamma_r9",                                      ev.gamma_r9                                   );
  // tree -> SetBranchAddress("gamma_sieie",                                   ev.gamma_sieie                                );

  //jet info
  tree -> SetBranchAddress("nj",                                          & ev.nj                                         );
  tree -> SetBranchAddress("j_g",                                           ev.j_g                                        );
  // tree -> SetBranchAddress("j_jerUp",                                       ev.j_jerUp                                    );
  // tree -> SetBranchAddress("j_jerDn",                                       ev.j_jerDn                                    );
  // for(unsigned char ind = 0; ind < 29; ind ++) 
  //   {
  //     tree -> SetBranchAddress(Form("j_jecUp%d", ind),                          ev.j_jecUp[ind]                               );
  //     tree -> SetBranchAddress(Form("j_jecDn%d", ind),                          ev.j_jecDn[ind]                               );
  //   }

  tree -> SetBranchAddress("j_area",                                        ev.j_area                                     );
  tree -> SetBranchAddress("j_rawsf",                                       ev.j_rawsf                                    );
  tree -> SetBranchAddress("j_pt",                                          ev.j_pt                                       );
  tree -> SetBranchAddress("j_eta",                                         ev.j_eta                                      );
  tree -> SetBranchAddress("j_phi",                                         ev.j_phi                                      );
  tree -> SetBranchAddress("j_mass",                                        ev.j_mass                                     );
  // tree -> SetBranchAddress("j_pumva",                                       ev.j_pumva                                    );
  // tree -> SetBranchAddress("j_id",                                          ev.j_id                                       );
  tree -> SetBranchAddress("j_csv",                                         ev.j_csv                                      );
  tree -> SetBranchAddress("j_btag",                                        ev.j_btag                                     );
  // tree -> SetBranchAddress("j_emf",                                         ev.j_emf                                      );
  // tree -> SetBranchAddress("j_qg",                                          ev.j_qg                                       );
  // tree -> SetBranchAddress("j_c2_00",                                       ev.j_c2_00                                    );
  // tree -> SetBranchAddress("j_c2_02",                                       ev.j_c2_02                                    );
  // tree -> SetBranchAddress("j_c2_05",                                       ev.j_c2_05                                    );
  // tree -> SetBranchAddress("j_c2_10",                                       ev.j_c2_10                                    );
  // tree -> SetBranchAddress("j_c2_20",                                       ev.j_c2_20                                    );
  // tree -> SetBranchAddress("j_zg",                                          ev.j_zg                                       );
  // tree -> SetBranchAddress("j_mult",                                        ev.j_mult                                     );
  // tree -> SetBranchAddress("j_gaptd",                                       ev.j_gaptd                                    );
  // tree -> SetBranchAddress("j_gawidth",                                     ev.j_gawidth                                  );
  // tree -> SetBranchAddress("j_gathrust",                                    ev.j_gathrust                                 );
  // tree -> SetBranchAddress("j_tau32",                                       ev.j_tau32                                    );
  // tree -> SetBranchAddress("j_tau21",                                       ev.j_tau21                                    );
  // tree -> SetBranchAddress("j_deepcsv",                                     ev.j_deepcsv                                  );
  // tree -> SetBranchAddress("j_probc",                                       ev.j_probc                                    );
  // tree -> SetBranchAddress("j_probudsg",                                    ev.j_probudsg                                 );
  // tree -> SetBranchAddress("j_probb",                                       ev.j_probb                                    );
  // tree -> SetBranchAddress("j_probbb",                                      ev.j_probbb                                   );
  // tree -> SetBranchAddress("j_CvsL",                                        ev.j_CvsL                                     );
  // tree -> SetBranchAddress("j_CvsB",                                        ev.j_CvsB                                     );
  tree -> SetBranchAddress("j_vtxpx",                                       ev.j_vtxpx                                    );
  tree -> SetBranchAddress("j_vtxpy",                                       ev.j_vtxpy                                    );
  tree -> SetBranchAddress("j_vtxpz",                                       ev.j_vtxpz                                    );
  tree -> SetBranchAddress("j_vtxmass",                                     ev.j_vtxmass                                  );
  tree -> SetBranchAddress("j_vtxNtracks",                                  ev.j_vtxNtracks                               );
  tree -> SetBranchAddress("j_vtx3DVal",                                    ev.j_vtx3DVal                                 );
  tree -> SetBranchAddress("j_vtx3DSig",                                    ev.j_vtx3DSig                                 );
  tree -> SetBranchAddress("j_flav",                                        ev.j_flav                                     );
  tree -> SetBranchAddress("j_hadflav",                                     ev.j_hadflav                                  );
  tree -> SetBranchAddress("j_pid",                                         ev.j_pid                                      );

  //PF Sums
  // tree -> SetBranchAddress("nchPV",                                       & ev.nchPV                                      );
  // tree -> SetBranchAddress("sumPVChPt",                                   & ev.sumPVChPt                                  );
  // tree -> SetBranchAddress("sumPVChPz",                                   & ev.sumPVChPz                                  );
  // tree -> SetBranchAddress("sumPVChHt",                                   & ev.sumPVChHt                                  );
  // tree -> SetBranchAddress("nPFCands",                                      ev.nPFCands                                   );
  // tree -> SetBranchAddress("sumPFEn",                                       ev.sumPFEn                                    );
  // tree -> SetBranchAddress("sumPFPz",                                       ev.sumPFPz                                    );
  // tree -> SetBranchAddress("sumPFHt",                                       ev.sumPFHt                                    );
  // tree -> SetBranchAddress("nPFChCands",                                    ev.nPFChCands                                 );
  // tree -> SetBranchAddress("sumPFChEn",                                     ev.sumPFChEn                                  );
  // tree -> SetBranchAddress("sumPFChPz",                                     ev.sumPFChPz                                  );
  // tree -> SetBranchAddress("sumPFChHt",                                     ev.sumPFChHt                                  );

  tree -> SetBranchAddress("npf",                                        & ev.npf                                          );
  tree -> SetBranchAddress("pf_j",                                         ev.pf_j                                         );
  tree -> SetBranchAddress("pf_id",                                        ev.pf_id                                        );
  tree -> SetBranchAddress("pf_c",                                         ev.pf_c                                         );
  tree -> SetBranchAddress("pf_pt",                                        ev.pf_pt                                        );
  tree -> SetBranchAddress("pf_eta",                                       ev.pf_eta                                       );
  tree -> SetBranchAddress("pf_phi",                                       ev.pf_phi                                       );
  tree -> SetBranchAddress("pf_m",                                         ev.pf_m                                         );
  tree -> SetBranchAddress("pf_dxy",                                       ev.pf_dxy                                       );
  tree -> SetBranchAddress("pf_dz",                                        ev.pf_dz                                        );
  //t->SetBranchAddress("pf_dxyUnc",   ev.pf_dxyUnc);                                                                                                                                                  
  //t->SetBranchAddress("pf_dzUnc",    ev.pf_dzUnc);                                                                                                                                                   
  //t->SetBranchAddress("pf_pvAssoc",  ev.pf_pvAssoc);                                                                                                                                                 
  tree -> SetBranchAddress("pf_puppiWgt",                                  ev.pf_puppiWgt                                  );


  //MET
  // tree -> SetBranchAddress("met_pt",                                      & ev.met_pt                                     );
  // tree -> SetBranchAddress("met_phi",                                     & ev.met_phi                                    );
  // tree -> SetBranchAddress("met_sig",                                     & ev.met_sig                                    );
  tree -> SetBranchAddress("met_pt",                                      & ev.met_pt                                     );
  tree -> SetBranchAddress("met_phi",                                     & ev.met_phi                                    );
  tree -> SetBranchAddress("met_sig",                                     & ev.met_sig                                    );
  // tree -> SetBranchAddress("met_ptShifted",                                 ev.met_ptShifted                              );
  // tree -> SetBranchAddress("met_phiShifted",                                ev.met_phiShifted                             );
  tree -> SetBranchAddress("met_filterBits",                              & ev.met_filterBits                             );

  //CTPPS local tracks
  // tree -> SetBranchAddress("nppstrk",                                     & ev.nppstrk                                    );
  // tree -> SetBranchAddress("ppstrk_pot",                                    ev.ppstrk_pot                                 );
  // tree -> SetBranchAddress("ppstrk_x",                                      ev.ppstrk_x                                   );
  // tree -> SetBranchAddress("ppstrk_y",                                      ev.ppstrk_y                                   );
  // tree -> SetBranchAddress("ppstrk_xUnc",                                   ev.ppstrk_xUnc                                );
  // tree -> SetBranchAddress("ppstrk_yUnc",                                   ev.ppstrk_yUnc                                );
  // tree -> SetBranchAddress("ppstrk_tx",                                     ev.ppstrk_tx                                  );
  // tree -> SetBranchAddress("ppstrk_ty",                                     ev.ppstrk_ty                                  );
  // tree -> SetBranchAddress("ppstrk_txUnc",                                  ev.ppstrk_txUnc                               );
  // tree -> SetBranchAddress("ppstrk_tyUnc",                                  ev.ppstrk_tyUnc                               );
  // tree -> SetBranchAddress("ppstrk_recoInfo",                               ev.ppstrk_recoInfo                            );
  // tree -> SetBranchAddress("ppstrk_chisqnorm",                              ev.ppstrk_chisqnorm                           );
  // tree -> SetBranchAddress("ppstrk_t",                                      ev.ppstrk_t                                   );
  // tree -> SetBranchAddress("ppstrk_tUnc",                                   ev.ppstrk_tUnc                                );
  // tree -> SetBranchAddress("nfwdtrk",                                     & ev.nfwdtrk                                    );
  // tree -> SetBranchAddress("fwdtrk_pot",                                    ev.fwdtrk_pot                                 );
  // tree -> SetBranchAddress("fwdtrk_method",                                 ev.fwdtrk_method                              );
  // tree -> SetBranchAddress("fwdtrk_recoInfo",                               ev.fwdtrk_recoInfo                            );
  // tree -> SetBranchAddress("fwdtrk_thetax",                                 ev.fwdtrk_thetax                              );
  // tree -> SetBranchAddress("fwdtrk_thetay",                                 ev.fwdtrk_thetay                              );
  // tree -> SetBranchAddress("fwdtrk_x",                                      ev.fwdtrk_x                                   );
  // tree -> SetBranchAddress("fwdtrk_y",                                      ev.fwdtrk_y                                   );
  // tree -> SetBranchAddress("fwdtrk_tx",                                     ev.fwdtrk_tx                                  );
  // tree -> SetBranchAddress("fwdtrk_ty",                                     ev.fwdtrk_ty                                  );
  // tree -> SetBranchAddress("fwdtrk_vx",                                     ev.fwdtrk_vx                                  );
  // tree -> SetBranchAddress("fwdtrk_vy",                                     ev.fwdtrk_vy                                  );
  // tree -> SetBranchAddress("fwdtrk_vz",                                     ev.fwdtrk_vz                                  );
  // tree -> SetBranchAddress("fwdtrk_time",                                   ev.fwdtrk_time                                );
  // tree -> SetBranchAddress("fwdtrk_timeError",                              ev.fwdtrk_timeError                           );
  // tree -> SetBranchAddress("fwdtrk_chisqnorm",                              ev.fwdtrk_chisqnorm                           );
  // tree -> SetBranchAddress("fwdtrk_xi",                                     ev.fwdtrk_xi                                  );
  // tree -> SetBranchAddress("fwdtrk_xiError",                                ev.fwdtrk_xiError                             );
  // tree -> SetBranchAddress("fwdtrk_t",                                      ev.fwdtrk_t                                   );

  //
  // tree -> SetBranchAddress("nrawmu",                                      & ev.nrawmu                                     );
  // tree -> SetBranchAddress("rawmu_pt",                                      ev.rawmu_pt                                   );
  // tree -> SetBranchAddress("rawmu_eta",                                     ev.rawmu_eta                                  );
  // tree -> SetBranchAddress("rawmu_phi",                                     ev.rawmu_phi                                  );
  // tree -> SetBranchAddress("rawmu_pid",                                     ev.rawmu_pid                                  );
}

void MiniEvent_t::ls()
{
  printf("npf %u\n", npf);
  for (int ind = 0; ind < npf; ind++)
    {
      // printf("%u: %f %f %f %f\n", ind, 
    }

}
