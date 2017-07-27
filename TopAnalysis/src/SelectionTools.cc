#include "TopLJets2015/TopAnalysis/interface/SelectionTools.h"

#include <iostream>

using namespace std;

SelectionTool::SelectionTool(TString dataset,bool debug,TH1 *triggerList) :
  isSingleElectronPD_(dataset.Contains("SingleElectron")), 
  isSingleMuonPD_(dataset.Contains("SingleMuon")), 
  isDoubleEGPD_(dataset.Contains("DoubleEG")), 
  isDoubleMuonPD_(dataset.Contains("DoubleMuon")), 
  isMuonEGPD_(dataset.Contains("MuonEG"))
{
  if(triggerList!=0)
    for(int xbin=0; xbin<triggerList->GetNbinsX(); xbin++)
      triggerBits_[ triggerList->GetXaxis()->GetBinLabel(xbin+1) ] = xbin;
}

//
bool SelectionTool::passMETFilters(MiniEvent_t &ev){  
  if(ev.isData) return ev.met_filterBits==0xff;
  else          return ((ev.met_filterBits&0xf)==0xf) && ((ev.met_filterBits>>5)==0x7);
}

//
std::vector<Particle> SelectionTool::getTopFlaggedLeptons(MiniEvent_t &ev){
  std::vector<Particle> leptons;
  for (int il=0; il<ev.nl; il++) {

    float pt(ev.l_pt[il]);
    float eta(fabs(ev.l_eta[il]));
    int pid(ev.l_pid[il]);
    float relIso(ev.l_relIso[il]);

    int topLeptonQualityFlagsWord(0);
    if(abs(ev.l_id[il])==11)
      {
	if( pt>20 && eta<2.4 && ((pid>>7) &0x1))                                     topLeptonQualityFlagsWord |= (0x1 << PASSLLID);
	if( pt>30 && eta<2.1 && ((pid>>4) &0x1))                                     topLeptonQualityFlagsWord |= (0x1 << PASSLID);
	if( pt>15 && eta<2.4 && ((pid>>2) &0x1))                                     topLeptonQualityFlagsWord |= (0x1 << PASSLVETO);
	if( pt>26 && eta<2.1 && ((pid>>5) &0x1) && ((pid>>4) &0x1)==0 && relIso>0.4) topLeptonQualityFlagsWord |= (0x1 << PASSLIDNONISO);
      }
    else
      {
	if( pt>20 && eta<2.4 && ((pid>>4) &0x1) && relIso<0.15)  topLeptonQualityFlagsWord |= (0x1 << PASSLLID);
	if( pt>30 && eta<2.1 && ((pid>>4) &0x1) && relIso<0.15)  topLeptonQualityFlagsWord |= (0x1 << PASSLID);
	if( pt>15 && eta<2.4 && ((pid>>1) &0x1) && relIso<0.25)  topLeptonQualityFlagsWord |= (0x1 << PASSLVETO);
	if( pt>26 && eta<2.1 && ((pid>>4) &0x1) && relIso>0.25)  topLeptonQualityFlagsWord |= (0x1 << PASSLIDNONISO);
      }

    if(debug_) cout << "Lepton #" << il << " id=" << ev.l_id[il] 
		    << " pt=" << pt << " eta=" << eta << " relIso=" << relIso 
		    << " charge=" << ev.l_charge[il]
		    << " flag=0x" << std::hex << topLeptonQualityFlagsWord << std::dec << endl;

    if(topLeptonQualityFlagsWord==0) continue;

    TLorentzVector lp4;
    lp4.SetPtEtaPhiM(ev.l_pt[il],ev.l_eta[il],ev.l_phi[il],ev.l_mass[il]);
    leptons.push_back(Particle(lp4, ev.l_charge[il], ev.l_id[il], topLeptonQualityFlagsWord, il, 1.0) );
  }
  
  return leptons;
}

//
std::vector<Particle> SelectionTool::getLeptons(std::vector<Particle> &leptons,int qualBit,double minPt, double maxEta,std::vector<Particle> *vetoParticles){
  std::vector<Particle> selLeptons;

  for(size_t i =0; i<leptons.size(); i++)
    {
      //check quality flag
      if( !leptons[i].hasQualityFlag(qualBit) ) continue;

      //check kinematics
      if(leptons[i].pt()<minPt || fabs(leptons[i].eta())>maxEta) continue;

      //check if this lepton should be vetoed by request
      if(vetoParticles){
        bool skipThisLepton(false);
        for(auto &vetoL : *vetoParticles){
          if(vetoL.originalReference()!=leptons[i].originalReference()) continue;
          skipThisLepton=true;
          break;
        }
        if(skipThisLepton) continue;
      }

      //lepton is selected
      selLeptons.push_back(leptons[i]);
    }

  //all done here
  return selLeptons;
}

//
std::vector<Particle> SelectionTool::getGenLeptons(MiniEvent_t &ev, double minPt, double maxEta){
  std::vector<Particle> leptons;
  
  //loop over leptons from pseudotop producer
  for (int i = 0; i < ev.ng; i++) {
    int absid(abs(ev.g_id[i]));
    if(absid!=11 && absid!=13) continue;

    bool passKin(ev.g_pt[i]>minPt && fabs(ev.g_eta[i])<maxEta);
    if(!passKin) continue;

    TLorentzVector lp4;
    lp4.SetPtEtaPhiM(ev.g_pt[i],ev.g_eta[i],ev.g_phi[i],ev.g_m[i]);
    leptons.push_back( Particle(lp4, -ev.g_id[i]/abs(ev.g_id[i]), ev.g_id[i], 0, 1) );
  }
  
  return leptons;
}

//
std::vector<Jet> SelectionTool::getGoodJets(MiniEvent_t &ev, double minPt, double maxEta, std::vector<Particle> leptons) {
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
    jet.setCSV(ev.j_csv[k]);

    //fill jet constituents
    for (int p = 0; p < ev.npf; p++) {
      if (ev.pf_j[p] == k) {
        TLorentzVector pp4;
        pp4.SetPtEtaPhiM(ev.pf_pt[p],ev.pf_eta[p],ev.pf_phi[p],ev.pf_m[p]);
        jet.addParticle(Particle(pp4, ev.pf_c[p], ev.pf_id[p], 0, p, ev.pf_puppiWgt[p]));
        if (ev.pf_c[p] != 0) jet.addTrack(pp4, ev.pf_id[p]);
      }
    }

    if(debug_) cout << "Jet #" << jets.size() 
		    << " pt=" << jp4.Pt() << " eta=" << jp4.Eta() << " csv=" << ev.j_csv[k] << endl;

    
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

//
std::vector<Jet> SelectionTool::getGenJets(MiniEvent_t &ev, double minPt, double maxEta, std::vector<Particle> *leptons) {
  std::vector<Jet> jets;
  
  for (int i = 0; i < ev.ng; i++) {
    if (abs(ev.g_id[i])>10) continue;
    TLorentzVector jp4;
    jp4.SetPtEtaPhiM(ev.g_pt[i],ev.g_eta[i],ev.g_phi[i],ev.g_m[i]);

    //cross clean with leptons
    bool overlapsWithLepton(false);
    if(leptons){
      for (auto& lepton : *leptons) {
        if(jp4.DeltaR(lepton.p4())<0.4) overlapsWithLepton=true;
      }
      if(overlapsWithLepton) continue;
    }

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
	//if(pp4.DeltaR(jp4)<0.4){
	jet.addParticle(Particle(pp4, ev.gpf_c[p], ev.gpf_id[p], 0, p, 1.));
	if (ev.gpf_c[p] != 0) jet.addTrack(pp4, ev.gpf_id[p]);
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

//
TString SelectionTool::flagFinalState(MiniEvent_t &ev, std::vector<Particle> preselleptons) {

  //clear vectors
  leptons_.clear(); 
  vetoLeptons_.clear();
  jets_.clear();

  //if no set of pre-selected leptons has been passed, use standard top selections
  if(preselleptons.size()==0) preselleptons=getTopFlaggedLeptons(ev);

  //decide the channel based on the lepton multiplicity and set lepton collections
  std::vector<Particle> passllid( getLeptons(preselleptons,PASSLLID) ), passlid( getLeptons(preselleptons,PASSLID) );
  TString chTag("");
  if(passllid.size()>=2){
    int ch( abs(passllid[0].id()*passllid[1].id()) );
    if      (ch==11*13) chTag = "EM";
    else if (ch==13*13) chTag = "MM";
    else if (ch==11*11) chTag = "EE";
    leptons_=passllid;
    vetoLeptons_=getLeptons(preselleptons,SelectionTool::PASSLVETO, 0., 99., &leptons_);
    atleastVetoLeptons_=getLeptons(preselleptons,SelectionTool::PASSLVETO);
  }
  else if(passlid.size()==1){
    int ch(abs(passlid[0].id()) );
    if      (ch==13) chTag = "M";
    else if (ch==11) chTag = "E";
    leptons_=passlid;
    vetoLeptons_=getLeptons(preselleptons,SelectionTool::PASSLVETO, 0., 99., &leptons_);
    atleastVetoLeptons_=getLeptons(preselleptons,SelectionTool::PASSLVETO);
  }

  //select jets based on the leptons
  jets_=getGoodJets(ev,30.,2.4,leptons_);

  //build the met
  met_.SetPtEtaPhiM( ev.met_pt[0], 0, ev.met_phi[0], 0. );


  //check if triggers have fired
  bool hasETrigger(  hasTriggerBit("HLT_Ele32_eta2p1_WPTight_Gsf_v",           ev.triggerBits) );
  bool hasMTrigger(  hasTriggerBit("HLT_IsoMu24_v",                            ev.triggerBits) || 
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
      hasEMTrigger=false;
      if(ev.run<=280385)
	{
	  hasEMTrigger |= hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v", ev.triggerBits);
	  hasEMTrigger |= hasTriggerBit("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",  ev.triggerBits);
	}
      if(ev.run>=278273 && ev.run<=280385)
	{
	  hasEMTrigger |= hasTriggerBit("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v", ev.triggerBits);
	}
      if(ev.run>=278273)
	{
	  hasEMTrigger |= hasTriggerBit("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",  ev.triggerBits);
	  hasEMTrigger |= hasTriggerBit("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v", ev.triggerBits);
	  hasEMTrigger |= hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", ev.triggerBits);
	}
    }


  //check consistency with data
  if(chTag=="EM")
    {
      if(!hasEMTrigger && !hasETrigger && !hasMTrigger)                         chTag="";
      if(isDoubleEGPD_      || isDoubleMuonPD_)                                 chTag="";
      if(isSingleElectronPD_ && (hasEMTrigger || !hasETrigger))                 chTag="";
      if(isSingleMuonPD_    && (hasEMTrigger  || hasETrigger || !hasMTrigger))  chTag="";
      if(isMuonEGPD_        && !hasEMTrigger)                                   chTag="";
    }
  if(chTag=="EE")
    {
      if(!hasEETrigger && !hasETrigger)                          chTag="";
      if(isMuonEGPD_ || isSingleMuonPD_ || isDoubleMuonPD_)      chTag="";
      if(isSingleElectronPD_ && (hasEETrigger || !hasETrigger) ) chTag="";
      if(isDoubleEGPD_      && !hasEETrigger)                    chTag="";
    }
  if(chTag=="MM")
    {
      if(!hasMMTrigger && !hasMTrigger)                         chTag="";
      if(isMuonEGPD_ || isSingleElectronPD_ || isDoubleEGPD_)   chTag="";
      if(isSingleMuonPD_ && (hasMMTrigger || !hasMTrigger) )    chTag="";
      if(isDoubleMuonPD_ && !hasMMTrigger)                      chTag="";
    }
  if(chTag=="M")
    {
      if(!hasMTrigger)                    chTag="";
      if(isMuonEGPD_ || isDoubleMuonPD_ || isDoubleEGPD_ || isSingleElectronPD_)   chTag="";
      if(isSingleMuonPD_ && !hasMTrigger) chTag="";
    }
  if(chTag=="E")
    {
      if(!hasETrigger)                        chTag="";
      if(isMuonEGPD_ || isDoubleMuonPD_ || isDoubleEGPD_ || isSingleMuonPD_)   chTag="";
      if(isSingleElectronPD_ && !hasETrigger) chTag="";
    }
      

  if(debug_) cout << "[flagFinalState] chTag=" << chTag << endl
		  << "\t Pre-selection lepton mult." << preselleptons.size() << endl
		  << "\t 2l cands=" << passllid.size() << " 1l cands=" << passlid.size() << endl
		  << "\t Trig bits. e=" << hasETrigger << " m=" << hasMTrigger << " em=" << hasEMTrigger << " mm=" << hasMMTrigger << " ee=" << hasEETrigger << endl
		  << "\t Sel mult. l=" << leptons_.size() << " vl=" << vetoLeptons_.size() << " j=" << jets_.size() << endl;

  //all done
  return chTag;
}

//
bool SelectionTool::hasTriggerBit(TString triggerName,unsigned int word) 
{ 
  std::map<TString,unsigned int>::iterator it=triggerBits_.find(triggerName);
  if(it==triggerBits_.end()) return false;
  unsigned int bit=it->second;
  return ((word>>bit)&0x1); 
}

//
TString SelectionTool::flagGenFinalState(MiniEvent_t &ev, std::vector<Particle> leptons) 
{
  //update current state
  genLeptons_=leptons;
  if(genLeptons_.size()==0) genLeptons_=getGenLeptons(ev,20.,2.5);
  genJets_=getGenJets(ev,30.,2.4,&genLeptons_);

  //decide the channel
  TString chTag("");
  if(genLeptons_.size()>=2) {
    int chId(abs(genLeptons_[0].id()*genLeptons_[1].id()));
    if      (chId==11*13) chTag = "EM";
    else if (chId==13*13) chTag = "MM";
    else if (chId==11*11) chTag = "EE";
  }
  if(genLeptons_.size()==1) {
    int absid(abs(genLeptons_[0].id()));
    if      (absid==13) chTag = "M";
    else if (absid==11) chTag = "E";
  }
  
  return chTag;
}

