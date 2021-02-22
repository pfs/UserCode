#include "TopLJets2015/TopAnalysis/interface/SelectionTools.h"
#include "DataFormats/MuonReco/interface/Muon.h"

#include <iostream>
using namespace std;
TString SelectionTool::flagFinalState_old(MiniEvent_t & ev, std::vector<Particle> preselleptons) 
{

  //clear vectors
  leptons_.clear(); 
  vetoLeptons_.clear();
  jets_.clear();
  jet_indices_.clear();
  gen_jet_indices_.clear();
  //if no set of pre-selected leptons has been passed, use standard top selections
  if (preselleptons.size() == 0) preselleptons = getTopFlaggedLeptons(ev);

  //decide the channel based on the lepton multiplicity and set lepton collections
  std::vector<Particle> passllid( getLeptons(preselleptons, PASSLLID) ), passlid( getLeptons(preselleptons, PASSLID) );
  TString chTag("");
  if (passllid.size() >= 2)
    {
      const int ch( abs(passllid[0].id() * passllid[1].id()) );
      if      (ch == 11*13) chTag = "EM";
      else if (ch == 13*13) chTag = "MM";
      else if (ch == 11*11) chTag = "EE";
      leptons_ = passllid;
    }
  else if (passlid.size() == 1)
    {
      const int ch(abs(passlid[0].id()) );
      if      (ch == 13) chTag = "M";
      else if (ch == 11) chTag = "E";
      leptons_ = passlid;
      vetoLeptons_ = getLeptons(preselleptons, SelectionTool::PASSLVETO, 0.0, 99.0, & leptons_);
    }

  //select jets based on the leptons
  jets_ = getGoodJets_old(ev, 30.0, 2.4, leptons_);

  //build the met
  met_.SetPtEtaPhiM( ev.met_pt, 0.0, ev.met_phi, 0.0 );


  //check if triggers have fired
  bool hasETrigger(  hasTriggerBit("HLT_Ele32_eta2p1_WPTight_Gsf_v",                       ev.triggerBits) );
  bool hasMTrigger(  hasTriggerBit("HLT_IsoMu24_v",                                        ev.triggerBits) || 
                     hasTriggerBit("HLT_IsoTkMu24_v",                                      ev.triggerBits) );
  bool hasEMTrigger( hasTriggerBit("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",     ev.triggerBits) ||
     	             hasTriggerBit("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",    ev.triggerBits) ||
		     hasTriggerBit("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v", ev.triggerBits) ||
		     hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",    ev.triggerBits) ||
		     hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", ev.triggerBits) );
  bool hasMMTrigger( hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",                ev.triggerBits) ||
		     hasTriggerBit("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",              ev.triggerBits) );
  bool hasEETrigger( hasTriggerBit("HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v",              ev.triggerBits) ||
		     hasTriggerBit("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",          ev.triggerBits) );
  if(ev.isData)
    {
      hasEMTrigger = false;
      if(ev.run <= 280385)
	{
	  hasEMTrigger |= hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",     ev.triggerBits);
	  hasEMTrigger |= hasTriggerBit("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",      ev.triggerBits);
	}
      if(ev.run >= 278273 && ev.run <= 280385)
	{
	  hasEMTrigger |= hasTriggerBit("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",     ev.triggerBits);
	}
      if(ev.run >= 278273)
	{
	  hasEMTrigger |= hasTriggerBit("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",   ev.triggerBits);
	  hasEMTrigger |= hasTriggerBit("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",  ev.triggerBits);
	  hasEMTrigger |= hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",  ev.triggerBits);
	}
    }


  //check consistency with data
  if (chTag == "EM")
    {
      if (!hasEMTrigger && !hasETrigger && !hasMTrigger)                            chTag = "";
      if (isDoubleEGPD_      || isDoubleMuonPD_)                                    chTag = "";
      if (isSingleElectronPD_ && (hasEMTrigger || !hasETrigger))                    chTag = "";
      if (isSingleMuonPD_    && (hasEMTrigger  || hasETrigger || !hasMTrigger))     chTag = "";
      if (isMuonEGPD_        && !hasEMTrigger)                                      chTag = "";
    }
  if (chTag == "EE")
    {
      if (!hasEETrigger && !hasETrigger)                                            chTag = "";
      if (isMuonEGPD_ || isSingleMuonPD_ || isDoubleMuonPD_)                        chTag = "";
      if (isSingleElectronPD_ && (hasEETrigger || !hasETrigger) )                   chTag = "";
      if (isDoubleEGPD_      && !hasEETrigger)                                      chTag = "";
    }
  if (chTag == "MM")
    {
      if (!hasMMTrigger && !hasMTrigger)                                            chTag = "";
      if (isMuonEGPD_ || isSingleElectronPD_ || isDoubleEGPD_)                      chTag = "";
      if (isSingleMuonPD_ && (hasMMTrigger || !hasMTrigger) )                       chTag = "";
      if (isDoubleMuonPD_ && !hasMMTrigger)                                         chTag = "";
    }
  if (chTag == "M")
    {
      if (!hasMTrigger)                                                             chTag = "";
      if (isMuonEGPD_ || isDoubleMuonPD_ || isDoubleEGPD_ || isSingleElectronPD_)   chTag = "";
      if (isSingleMuonPD_ && !hasMTrigger)                                          chTag = "";
    }
  if (chTag == "E")
    {
      if (!hasETrigger)                                                             chTag = "";
      if (isMuonEGPD_ || isDoubleMuonPD_ || isDoubleEGPD_ || isSingleMuonPD_)       chTag = "";
      if (isSingleElectronPD_ && !hasETrigger)                                      chTag = "";
    }
      

  if (debug_) cout << "[flagFinalState] chTag=" << chTag << endl
  		   << "\t Pre-selection lepton mult." << preselleptons.size() << endl
		   << "\t 2l cands=" << passllid.size() << " 1l cands=" << passlid.size() << endl
		   << "\t Trig bits. e=" << hasETrigger << " m=" << hasMTrigger << " em=" << hasEMTrigger << " mm=" << hasMMTrigger << " ee=" << hasEETrigger << endl
		   << "\t Sel mult. l=" << leptons_.size() << " vl=" << vetoLeptons_.size() << " j=" << jets_.size() << endl;

  //all done
  return chTag;
}


std::vector<Jet> SelectionTool::getGoodJets_old(MiniEvent_t & ev, double minPt, double maxEta, std::vector<Particle> leptons) 
{
  std::vector<Jet> jets;
  
  for (int kind = 0; kind < ev.nj; kind ++) 
    {
      TLorentzVector jp4;
      jp4.SetPtEtaPhiM(ev.j_pt[kind], ev.j_eta[kind], ev.j_phi[kind], ev.j_mass[kind]);
      //flavor based on b tagging

      int flavor = 0;
      if (ev.j_btag[kind]) 
	{
	  flavor = 5;
	}
      //cross clean with leptons
      bool overlapsWithLepton(false);
      for (auto & lepton : leptons) 
	{
	  if (jp4.DeltaR(lepton.p4()) < 0.4) 
	    overlapsWithLepton = true;
	}
      if (overlapsWithLepton) 
	continue;
    
      //jet kinematic selection
      if (jp4.Pt() < minPt || abs(jp4.Eta()) > maxEta) 
	continue;

    
      Jet jet(jp4, flavor, kind);
      jet.setCSV(ev.j_csv[kind]);
      jet.setPartonFlavor(ev.j_flav[kind]);
      //fill jet constituents
      for (int prtcl = 0; prtcl < ev.npf; prtcl ++) 
	{
	  if (ev.pf_j[prtcl] == kind) 
	    {
	      TLorentzVector pp4;
	      pp4.SetPtEtaPhiM(ev.pf_pt[prtcl], ev.pf_eta[prtcl], ev.pf_phi[prtcl], ev.pf_m[prtcl]);
	      jet.addParticle(Particle(pp4, ev.pf_c[prtcl], ev.pf_id[prtcl], 0, prtcl, ev.pf_puppiWgt[prtcl]));
	      if (ev.pf_c[prtcl] != 0) 
		jet.addTrack(pp4, ev.pf_id[prtcl]);
	    }
	}

      if(debug_) cout << "Jet #" << jets.size() 
		      << " pt=" << jp4.Pt() << " eta=" << jp4.Eta() << " csv=" << ev.j_csv[kind] << endl;

    
    
      jets.push_back(jet);
      jet_indices_.push_back(kind);
    }
  
  //additional jet-jet information
  for (unsigned int iind = 0; iind < jets.size(); iind ++) 
    {
      for (unsigned int jind = iind + 1; jind < jets.size(); jind ++) 
	{
	  //flag jet-jet overlaps
	  if (jets[iind].p4().DeltaR(jets[jind].p4()) < 0.8) 
	    {
	      jets[iind].setOverlap(1);
	      jets[jind].setOverlap(1);
	    }
	  //flag non-b jets as part of W boson candidates: flavor 0->1
	  if (jets[iind].flavor() == 5 or jets[jind].flavor() == 5) continue;
	  TLorentzVector wCand = jets[iind].p4() + jets[jind].p4();
	  if (abs(wCand.M() - 80.4) < 15.0) 
	    {
	      jets[iind].setFlavor(1);
	      jets[jind].setFlavor(1);
	    }
	}
    }
  
  return jets;
}
