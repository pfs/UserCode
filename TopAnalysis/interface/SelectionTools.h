#ifndef _selection_tools_h_
#define _selection_tools_h_

#include <vector>

#include "TString.h"

#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TopLJets2015/TopAnalysis/interface/ObjectTools.h"

std::vector<Particle> getGoodLeptons(MiniEvent_t ev, double tightMinPt = 30., double tightMaxEta = 2.1, bool vetoLoose = false, double looseMinPt = 15., double looseMaxEta = 2.5) {
  std::vector<Particle> leptons;
  
  for (int il=0; il<ev.nl; il++) {
	  bool passTightKin(ev.l_pt[il]>tightMinPt && fabs(ev.l_eta[il])<tightMaxEta);
    bool passLooseKin(ev.l_pt[il]>looseMinPt && fabs(ev.l_eta[il])<looseMaxEta);
	  bool passTightId (ev.l_id[il]==13 ? (ev.l_pid[il]>>1)&0x1 : (ev.l_pid[il]>>2)&0x1);
    bool passTightIso(ev.l_id[il]==13 ?  ev.l_relIso[il]<0.15 : (ev.l_pid[il]>>1)&0x1);
    bool passLooseIso(ev.l_id[il]==13 ?  ev.l_relIso[il]<0.25 : (ev.l_pid[il]   )&0x1);
	  if(passTightKin && passTightId && passTightIso) {
      TLorentzVector lp4;
      lp4.SetPtEtaPhiM(ev.l_pt[il],ev.l_eta[il],ev.l_phi[il],ev.l_mass[il]);
      leptons.push_back(Particle(lp4, ev.l_charge[il], ev.l_id[il]));
    }
    else if(vetoLoose && passLooseKin && passLooseIso) {
      return {};
    }
	}
  
  return leptons;
}

std::vector<Particle> getGenLeptons(MiniEvent_t ev, double tightMinPt = 30., double tightMaxEta = 2.1, bool vetoLoose = false, double looseMinPt = 15., double looseMaxEta = 2.5) {
  std::vector<Particle> leptons;
  
  //loop over leptons from pseudotop producer
  for (int i = 0; i < ev.ng; i++) {
    if (abs(ev.g_id[i])>10) {
      bool passTightKin(ev.g_pt[i]>tightMinPt && fabs(ev.g_eta[i])<tightMaxEta);
      bool passLooseKin(ev.g_pt[i]>looseMinPt && fabs(ev.g_eta[i])<looseMaxEta);
      
      if(passTightKin) {
        TLorentzVector lp4;
        lp4.SetPtEtaPhiM(ev.g_pt[i],ev.g_eta[i],ev.g_phi[i],ev.g_m[i]);
        leptons.push_back(Particle(lp4, -ev.g_id[i]/abs(ev.g_id[i]), ev.g_id[i]));
      }
      else if(vetoLoose && passLooseKin) {
        return {};
      }
    }
  }
  
  return leptons;
}

std::vector<Jet> getGoodJets(MiniEvent_t ev, double minPt = 30., double maxEta = 2.4, std::vector<Particle> leptons = {}) {
  std::vector<Jet> jets;
  
  for (int k=0; k<ev.nj; k++) {
    TLorentzVector jp4;
    jp4.SetPtEtaPhiM(ev.j_pt[k],ev.j_eta[k],ev.j_phi[k],ev.j_mass[k]);

    //cross clean with leptons
    bool overlapsWithLepton(false);
    for (auto& lepton : leptons) {
      if(jp4.DeltaR(lepton.p4())<0.4) overlapsWithLepton=true;
    }
    if(overlapsWithLepton) continue;
    
    //jet kinematic selection
    if(jp4.Pt() < minPt || abs(jp4.Eta()) > maxEta) continue;

    //flavor based on b tagging
    int flavor = 0;
    if (ev.j_btag[k]) {
      flavor = 5;
    }
    
    Jet jet(jp4, flavor, k);
    
    //fill jet constituents
    for (int p = 0; p < ev.npf; p++) {
      if (ev.pf_j[p] == k) {
        TLorentzVector pp4;
        pp4.SetPtEtaPhiM(ev.pf_pt[p],ev.pf_eta[p],ev.pf_phi[p],ev.pf_m[p]);
        jet.addParticle(Particle(pp4, ev.pf_c[p], ev.pf_id[p], ev.pf_puppiWgt[p]));
        if (ev.pf_c[p] != 0) jet.addTrack(pp4, ev.pf_id[p]);
      }
    }
    
    jets.push_back(jet);
  }
  
  //additional jet-jet information
  for (unsigned int i = 0; i < jets.size(); i++) {
    for (unsigned int j = i+1; j < jets.size(); j++) {
      //flag jet-jet overlaps
      if (jets[i].p4().DeltaR(jets[j].p4()) < 0.8) {
        jets[i].setOverlap(1);
        jets[j].setOverlap(1);
      }
      //flag non-b jets as part of W boson candidates: flavor 0->1
      if (jets[i].flavor()==5 or jets[j].flavor()==5) continue;
      TLorentzVector wCand = jets[i].p4() + jets[j].p4();
      if (abs(wCand.M()-80.4) < 15.) {
        jets[i].setFlavor(1);
        jets[j].setFlavor(1);
      }
    }
  }
  
  return jets;
}

std::vector<Jet> getGenJets(MiniEvent_t ev, double minPt = 30., double maxEta = 2.4, std::vector<Particle> leptons = {}) {
  std::vector<Jet> jets;
  
  for (int i = 0; i < ev.ng; i++) {
    if (abs(ev.g_id[i])<10) {
      TLorentzVector jp4;
      jp4.SetPtEtaPhiM(ev.g_pt[i],ev.g_eta[i],ev.g_phi[i],ev.g_m[i]);

      //cross clean with leptons
      bool overlapsWithLepton(false);
      for (auto& lepton : leptons) {
        if(jp4.DeltaR(lepton.p4())<0.4) overlapsWithLepton=true;
      }
      if(overlapsWithLepton) continue;
      
      //jet kinematic selection
      if(jp4.Pt() < minPt || abs(jp4.Eta()) > maxEta) continue;

      //flavor
      int flavor = ev.g_id[i];
      
      Jet jet(jp4, flavor, i);
      
      //fill jet constituents
      for (int p = 0; p < ev.ngpf; p++) {
        if (ev.gpf_g[p] == i) {
          TLorentzVector pp4;
          pp4.SetPtEtaPhiM(ev.gpf_pt[p],ev.gpf_eta[p],ev.gpf_phi[p],ev.gpf_m[p]);
          jet.addParticle(Particle(pp4, ev.gpf_c[p], ev.gpf_id[p], 1.));
          if (ev.gpf_c[p] != 0) jet.addTrack(pp4, ev.gpf_id[p]);
        }
      }
      
      jets.push_back(jet);
    }
  }
  
  //additional jet-jet information
  for (unsigned int i = 0; i < jets.size(); i++) {
    for (unsigned int j = i+1; j < jets.size(); j++) {
      //flag jet-jet overlaps
      if (jets[i].p4().DeltaR(jets[j].p4()) < 0.8) {
        jets[i].setOverlap(1);
        jets[j].setOverlap(1);
      }
      //flag non-b jets as part of W boson candidates: flavor 0->1
      if (jets[i].flavor()==5 or jets[j].flavor()==5) continue;
      TLorentzVector wCand = jets[i].p4() + jets[j].p4();
      if (abs(wCand.M()-80.4) < 15.) {
        jets[i].setFlavor(1);
        jets[j].setFlavor(1);
      }
    }
  }
  
  return jets;
}

TString getChannel(MiniEvent_t ev, TString filename, std::vector<Particle> leptons) {
  bool requireEETriggers(false);
  if(ev.isData && filename.Contains("DoubleEG"))       requireEETriggers=true;
  bool requireMMTriggers(false);
  if(ev.isData && filename.Contains("DoubleMuon"))     requireMMTriggers=true;
  bool requireEMTriggers(false);
  if(ev.isData && filename.Contains("MuonEG"))         requireEMTriggers=true;
  bool requireETriggers(false);
  if(ev.isData && filename.Contains("SingleElectron")) requireETriggers=true;
  bool requireMTriggers(false);
  if(ev.isData && filename.Contains("SingleMuon"))     requireMTriggers=true;
  
  //check if triggers have fired
  bool hasEETrigger(requireEETriggers && ((ev.elTrigger>>1)&0x1)!=0 || ((ev.elTrigger>>4)&0x1)!=0);
  bool hasMMTrigger(requireMMTriggers && ((ev.muTrigger>>2)&0x3)!=0);
  bool hasEMTrigger(requireEMTriggers && ((ev.elTrigger>>2)&0x3)!=0);
  bool hasETrigger (requireETriggers  && ((ev.elTrigger>>0)&0x1)!=0);
  bool hasMTrigger (requireMTriggers  && ((ev.muTrigger>>0)&0x3)!=0);
  if(!ev.isData) { 
    hasEETrigger = true;
    hasMMTrigger = true;
    hasEMTrigger = true;
    hasETrigger  = true;
    hasMTrigger  = true;
  }

  //decide the channel
  TString chTag("");
  if(leptons.size()>=2) {
    if      (abs(leptons[0].id()*leptons[1].id())==11*13 && hasEMTrigger) chTag = "EM";
    else if (abs(leptons[0].id()*leptons[1].id())==13*13 && hasMMTrigger) chTag = "MM";
    else if (abs(leptons[0].id()*leptons[1].id())==11*11 && hasEETrigger) chTag = "EE";
  }
  if(leptons.size()==1) {
    if      (abs(leptons[0].id())==13 && hasMTrigger) chTag = "M";
    else if (abs(leptons[0].id())==11 && hasETrigger) chTag = "E";
  }
  
  return chTag;
}

TString getGenChannel(std::vector<Particle> leptons) {
  //decide the channel
  TString chTag("");
  if(leptons.size()>=2) {
    if      (abs(leptons[0].id()*leptons[1].id())==11*13) chTag = "EM";
    else if (abs(leptons[0].id()*leptons[1].id())==13*13) chTag = "MM";
    else if (abs(leptons[0].id()*leptons[1].id())==11*11) chTag = "EE";
  }
  if(leptons.size()==1) {
    if      (abs(leptons[0].id())==13) chTag = "M";
    else if (abs(leptons[0].id())==11) chTag = "E";
  }
  
  return chTag;
}

#endif
