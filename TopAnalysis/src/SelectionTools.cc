#include "TopLJets2015/TopAnalysis/interface/SelectionTools.h"
#include "DataFormats/MuonReco/interface/Muon.h"

#include <iostream>

using namespace std;

//
SelectionTool::SelectionTool(TString dataset_,bool debug,TH1 *triggerList, AnalysisType anType) :
  dataset(dataset_),
  debug_(debug),
  anType_(anType),
  isSingleElectronPD_(dataset.Contains("SingleElectron")), 
  isSingleMuonPD_(dataset.Contains("SingleMuon")), 
  isDoubleEGPD_(dataset.Contains("DoubleEG")), 
  isDoubleMuonPD_(dataset.Contains("DoubleMuon")), 
  isMuonEGPD_(dataset.Contains("MuonEG")),
  isPhotonPD_(dataset.Contains("Photon") || dataset.Contains("EGamma")),
  isJetHTPD_(dataset.Contains("JetHT"))
{
  if(triggerList!=0)
    for(int xbin=0; xbin<triggerList->GetNbinsX(); xbin++)
      triggerBits_[ triggerList->GetXaxis()->GetBinLabel(xbin+1) ] = xbin;  

  setPhotonSelection();
}

//
// RECO LEVEL SELECTORS
//

//
bool SelectionTool::passSingleLeptonTrigger(MiniEvent_t &ev) {
  //check if triggers have fired and are consistent with the offline selection
  bool hasETrigger(  hasTriggerBit("HLT_Ele35_eta2p1_WPTight_Gsf_v",           ev.triggerBits) ||
                     hasTriggerBit("HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v",     ev.triggerBits) ||
                     hasTriggerBit("HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v",     ev.triggerBits) );
  bool hasMTrigger(  //hasTriggerBit("HLT_IsoMu24_2p1_v",                                        ev.triggerBits) || 
		     hasTriggerBit("HLT_IsoMu27_v",                                            ev.triggerBits) );

  if(!hasETrigger && !hasMTrigger) return false;
  if(ev.isData) {
	if(hasETrigger && !hasMTrigger) { if(!isSingleElectronPD_) return false; }
 	if(!hasETrigger && hasMTrigger) { if(!isSingleMuonPD_)     return false; }
  	if(hasETrigger && hasMTrigger)  { if(!isSingleMuonPD_)     return false; }
  }
  return true;
}


//
TString SelectionTool::flagFinalState(MiniEvent_t &ev, std::vector<Particle> preselLeptons,std::vector<Particle> preselPhotons, bool isCR, bool isQCDTemp, bool isSRfake) {

 //clear vectors
  leptons_.clear(); 
  photons_.clear();
  vetoLeptons_.clear();
  jets_.clear();

  //if no set of pre-selected leptons has been passed, use standard top selections
  if(preselLeptons.size()==0) preselLeptons=flaggedLeptons(ev);
  if(preselPhotons.size()==0) preselPhotons=flaggedPhotons(ev);

  //decide the channel based on the lepton multiplicity and set lepton collections
  std::vector<Particle> tightLeptons( selLeptons(preselLeptons,TIGHT,MVA80) );
  std::vector<Particle> tightPhotons( selPhotons(preselPhotons,offlinePhoton_, tightLeptons) );
  std::vector<Particle> inclusivePhotons( selPhotons(preselPhotons,CONTROL, tightLeptons) );
  tmpPhotons          = selPhotons(preselPhotons,QCDTEMP, tightLeptons);
  relaxedTightPhotons = selPhotons(preselPhotons,RELAXEDTIGHT, tightLeptons);
  std::vector<Particle> fakePhotons;
  for(auto a : inclusivePhotons) {
    int idx = a.originalReference();
    if (!this->isFakePhoton(ev,idx)) continue;
    fakePhotons.push_back(a);
  }
  TString chTag("");
  if(anType_==TOP)
    {
      if(tightLeptons.size()>=2){
        int ch( abs(tightLeptons[0].id()*tightLeptons[1].id()) );
        if      (ch==11*13) chTag = "EM";
        else if (ch==13*13) chTag = "MM";
        else if (ch==11*11) chTag = "EE";
        leptons_=tightLeptons;
      }
      else if(tightLeptons.size()==1){
        int ch(abs(tightLeptons[0].id()) );
        if      (ch==13) chTag = "M";
        else if (ch==11) chTag = "E";
        leptons_=tightLeptons;
        vetoLeptons_=selLeptons(preselLeptons,VETO, VETO, 0., 99., leptons_);
      }
    } else if(anType_==VBF){ 
      if (!isCR){
	if(tightLeptons.size()==2){
	  int ch( abs(tightLeptons[0].id()*tightLeptons[1].id()) );
	  float mll( (tightLeptons[0]+tightLeptons[1]).M() );
	  if( ch==13*13 && fabs(mll-91)<15 && (tightLeptons[0].pt()>30 || tightLeptons[1].pt()>30)) chTag="MM";          
	  leptons_=tightLeptons;
	} else {
	  if (!isSRfake) {
	    if( tightPhotons.size()>=1){
	      chTag="A";
	      leptons_   =tightLeptons;
	      photons_   =tightPhotons;
	    }
	  } else {
	    if(fakePhotons.size()>=1){
	      chTag="A";
	      leptons_   =tightLeptons;
	      photons_   =tightPhotons;
	    }
	  }
	}
      } else {
	if(isSRfake) return "";
	bool passPhoton = (!isSRfake && !isQCDTemp && inclusivePhotons.size()>=1) || (!isSRfake && isQCDTemp && tmpPhotons.size()>=1);
	if(passPhoton) {
	  chTag="A";
	  if(!isQCDTemp) photons_   =inclusivePhotons;
	  else                photons_   =tmpPhotons;
	  //cout<< "Number of very loose photons: "<<photons_.size()<<endl;
	  leptons_   =tightLeptons;
	} else 	if(tightLeptons.size()==2){
	  int ch( abs(tightLeptons[0].id()*tightLeptons[1].id()) );
	  float mll( (tightLeptons[0]+tightLeptons[1]).M() );
	  if( ch==13*13 && fabs(mll-91)<15 && (tightLeptons[0].pt()>30 || tightLeptons[1].pt()>30)) chTag="MM";          
	  leptons_=tightLeptons;
	}
      }
    }
  
  //select jets based on the leptons and photon candidates
  float maxJetEta(2.4);
  if(anType_==VBF) maxJetEta=4.7;
  jets_=getGoodJets(ev,30.,maxJetEta,leptons_,photons_);
  //getGoodJets(ev,30.,maxJetEta,leptons_,photons_);

  //build the met
  met_.SetPtEtaPhiM( ev.met_pt, 0, ev.met_phi, 0. );

  //check if triggers have fired and are consistent with the offline selection
  bool hasETrigger(  hasTriggerBit("HLT_Ele35_WPTight_Gsf_v",                               ev.triggerBits) );
  bool hasMTrigger(  hasTriggerBit("HLT_IsoMu24_v",                                         ev.triggerBits) ||
                     hasTriggerBit("HLT_IsoMu24_2p1_v",                                     ev.triggerBits) ||
                     hasTriggerBit("HLT_IsoMu27_v",                                         ev.triggerBits) );
  bool hasEMTrigger( hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",     ev.triggerBits) ||
                     hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",  ev.triggerBits) ||
                     hasTriggerBit("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",  ev.triggerBits) ||
                     hasTriggerBit("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",   ev.triggerBits) );
  bool hasMMTrigger( hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",                   ev.triggerBits) ||
                     hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v",           ev.triggerBits) ||
                     hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v",         ev.triggerBits) );
  bool hasEETrigger( hasTriggerBit("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v",              ev.triggerBits) ||
                     hasTriggerBit("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",           ev.triggerBits) );
  bool hasPhotonTrigger(false);
  for(auto &t:photonTriggers_) hasPhotonTrigger |= hasTriggerBit(t, ev.triggerBits);

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
  if(chTag=="A")
    {
      if(!hasPhotonTrigger) chTag="";
      if(ev.isData && isCR && !isPhotonPD_ && !isJetHTPD_) chTag = "";
      if(ev.isData && !isCR && !isPhotonPD_ ) chTag = "";  
      //if(ev.isData && !isPhotonPD_) chTag="";    
    }
      
  if(debug_) cout << "[flagFinalState] chTag=" << chTag << endl
		  << "\t Pre-selection lepton mult." << preselLeptons.size() << endl
                  << "\t tight lepton cands=" << tightLeptons.size()  << endl
                  << "\t Pre-selection photon mult." << preselPhotons.size()
                  << "\t photon cands=" << tightPhotons.size() << endl               
		  << "\t Trig bits."
                  << " e=" << hasETrigger << " m=" << hasMTrigger 
                  << " em=" << hasEMTrigger << " mm=" << hasMMTrigger << " ee=" << hasEETrigger 
                  << " gamma=" << hasPhotonTrigger << endl
		  << "\t Sel mult. l=" << leptons_.size() << " vl=" << vetoLeptons_.size() 
                  << " photons=" << photons_.size() 
                  << " j=" << jets_.size() << endl;

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
bool SelectionTool::passMETFilters(MiniEvent_t &ev){  
  if(ev.isData) return ev.met_filterBits==0xff;
  else          return ((ev.met_filterBits&0xf)==0xf) && ((ev.met_filterBits>>5)==0x7);
}

//
std::vector<Particle> SelectionTool::flaggedLeptons(MiniEvent_t &ev)
{
  //leptons
  std::vector<Particle> leptons;
  for (int il=0; il<ev.nl; il++) {

    float pt(ev.l_pt[il]);
    float eta(fabs(ev.l_eta[il]));
    unsigned int pid(ev.l_pid[il]);
    float relIso(ev.l_relIso[il]);

    //see bits in plugins/MiniAnalyzer.cc
    int qualityFlagsWord(0);
    Float_t unc(0.);
    if(abs(ev.l_id[il])==11)
      {
	if( pt>20 && eta<2.5 ) {
          if((pid>>1)&0x1)  qualityFlagsWord |= (0x1 << VETO);
          if((pid>>2)&0x1)  qualityFlagsWord |= (0x1 << LOOSEIDONLY);
          if((pid>>3)&0x1)  qualityFlagsWord |= (0x1 << LOOSE);
          if((pid>>4)&0x1)  qualityFlagsWord |= (0x1 << MEDIUMIDONLY);
          if((pid>>5)&0x1)  qualityFlagsWord |= (0x1 << MEDIUM);
          if((pid>>6)&0x1)  qualityFlagsWord |= (0x1 << TIGHTIDONLY);
          if((pid>>7)&0x1)  qualityFlagsWord |= (0x1 << TIGHT);
          if((pid>>9)&0x1)  qualityFlagsWord |= (0x1 << MVA80);
          if((pid>>10)&0x1) qualityFlagsWord |= (0x1 << MVA90);
          if((pid>>10)&0x1) qualityFlagsWord |= (0x1 << MVANONISOWPLOOSE);
        }
        unc = TMath::Sqrt(
                          pow(ev.l_scaleUnc1[il],2)+
                          pow(ev.l_scaleUnc2[il],2)+
                          pow(ev.l_scaleUnc3[il],2)+
                          pow(ev.l_scaleUnc4[il],2)+
                          pow(ev.l_scaleUnc5[il],2)+
                          pow(ev.l_scaleUnc6[il],2)+
                          pow(ev.l_scaleUnc7[il],2)
                          );
      }
    else
      {
        if(pt>20 && eta<2.5) {
          if( (pid&reco::Muon::Selector::CutBasedIdLoose)==reco::Muon::Selector::CutBasedIdLoose ) {
            qualityFlagsWord |= (0x1 << LOOSEIDONLY);
            if( (pid&reco::Muon::Selector::PFIsoLoose)==reco::Muon::Selector::PFIsoLoose ) 
              qualityFlagsWord |= (0x1 << LOOSE);
          }
          if( (pid&reco::Muon::Selector::CutBasedIdMedium)==reco::Muon::Selector::CutBasedIdMedium ) {
            qualityFlagsWord |= (0x1 << MEDIUMIDONLY);
            if( (pid&reco::Muon::Selector::PFIsoMedium)==reco::Muon::Selector::PFIsoMedium ) 
              qualityFlagsWord |= (0x1 << MEDIUM);
          }
          if( (pid&reco::Muon::Selector::CutBasedIdTight)==reco::Muon::Selector::CutBasedIdTight ) {
            qualityFlagsWord |= (0x1 << TIGHTIDONLY);
            if( (pid&reco::Muon::Selector::PFIsoTight)==reco::Muon::Selector::PFIsoTight ) 
              qualityFlagsWord |= (0x1 << TIGHT);
          }
          if( (pid&reco::Muon::Selector::CutBasedIdTrkHighPt)==reco::Muon::Selector::CutBasedIdTrkHighPt ) {
            qualityFlagsWord |= (0x1 << HIGHPTIDONLY);
            if( (pid&reco::Muon::Selector::TkIsoLoose)==reco::Muon::Selector::TkIsoLoose) 
              qualityFlagsWord |= (0x1 << HIGHPT);            
          }
        }
      }

    if(debug_) cout << "Lepton #" << il << " id=" << ev.l_id[il] 
		    << " pt=" << pt << "+/-" << unc << " eta=" << eta << " relIso=" << relIso 
		    << " charge=" << ev.l_charge[il]
                    << " rawId = 0x" << std::hex << pid
		    << " quality flag=0x" << qualityFlagsWord << std::dec << endl;

    if(qualityFlagsWord==0) continue;

    TLorentzVector lp4;
    lp4.SetPtEtaPhiM(ev.l_pt[il],ev.l_eta[il],ev.l_phi[il],ev.l_mass[il]);
    leptons.push_back(Particle(lp4, ev.l_charge[il], ev.l_id[il], qualityFlagsWord, il, 1.0, unc));
  }
  
  return leptons;
}


//
std::vector<Particle> SelectionTool::selLeptons(std::vector<Particle> &leptons,int muQualBit,int eleQualBit,double minPt, double maxEta,std::vector<Particle> veto){
  std::vector<Particle> selLeptons;
  for(size_t i =0; i<leptons.size(); i++)
    {
      //check quality flag
      if(leptons[i].id()==11 && !leptons[i].hasQualityFlag(eleQualBit) ) continue;
      if(leptons[i].id()==13 && !leptons[i].hasQualityFlag(muQualBit) ) continue;

      //check kinematics
      if(leptons[i].pt()<minPt || fabs(leptons[i].eta())>maxEta) continue;

      //check if this lepton should be vetoed by request      
      bool skipThisLepton(false);
      for(auto &vetoL : veto){
        if(vetoL.originalReference()!=leptons[i].originalReference()) continue;
        skipThisLepton=true;
        break;
      }
      if(skipThisLepton) continue;
      
      
      //lepton is selected
      selLeptons.push_back(leptons[i]);
    }

  //all done here
  return selLeptons;
}

//
std::vector<Particle> SelectionTool::flaggedPhotons(MiniEvent_t &ev)
{
  //leptons
  std::vector<Particle> photons;
  for (int ig=0; ig<ev.ngamma; ig++) {
    float pt(ev.gamma_pt[ig]);
    float eta(fabs(ev.gamma_eta[ig]));
    int pid(ev.gamma_pid[ig]);
    int addpid(ev.gamma_idFlags[ig]);

    //see bits in plugins/MiniAnalyzer.cc

    int qualityFlagsWord(0);
    if( pt>30 && eta<2.4)
      {
        if( (pid&0x7f)==0x7f )            qualityFlagsWord |= (0x1 << LOOSE);
        if( ((pid>>10)&0x7f)==0x7f   )    qualityFlagsWord |= (0x1 << MEDIUM);
        if( ((pid>>20)&0x7f)==0x7f   )    qualityFlagsWord |= (0x1 << TIGHT);
        if( ((addpid>>2)&0x1)        )    qualityFlagsWord |= (0x1 << MVA80);
        if( ((addpid>>3)&0x1) )           qualityFlagsWord |= (0x1 << MVA90);
	if( isInclusivePhoton(ev,ig) )    qualityFlagsWord |= (0x1 << CONTROL);
	if( isQCDTemplate(ev,ig))         qualityFlagsWord |= (0x1 << QCDTEMP);
	if( isRelaxedTight(ev,ig)    )    qualityFlagsWord |= (0x1 << RELAXEDTIGHT);
      }
    if(qualityFlagsWord==0) continue;

    float unc = TMath::Sqrt(
                            pow(ev.gamma_scaleUnc1[ig],2)+
                            pow(ev.gamma_scaleUnc2[ig],2)+
                            pow(ev.gamma_scaleUnc3[ig],2)+
                            pow(ev.gamma_scaleUnc4[ig],2)+
                            pow(ev.gamma_scaleUnc5[ig],2)+
                            pow(ev.gamma_scaleUnc6[ig],2)+
                            pow(ev.gamma_scaleUnc7[ig],2)
                            );

    TLorentzVector p4;
    p4.SetPtEtaPhiM(ev.gamma_pt[ig],ev.gamma_eta[ig],ev.gamma_phi[ig],0);
    photons.push_back(Particle(p4, 0, 22, qualityFlagsWord, ig, 1.0, unc));

    if(debug_) std::cout << "Photon #"<< photons.size() 
                         << " pt=" << pt << "+/-" << unc << " eta=" << p4.Eta()
                         << hex << " raw particle id bits=" << pid 
                         << " quality bits=" << qualityFlagsWord 
                         << dec << endl;

  }
  
  return photons;
}

//
std::vector<Particle> SelectionTool::selPhotons(std::vector<Particle> &photons,int qualBit, std::vector<Particle> leptons,double minPt, double maxEta,std::vector<Particle> veto){
  std::vector<Particle> selPhotons;
  for(size_t i =0; i<photons.size(); i++)
    {
      //check quality flag
      if( !photons[i].hasQualityFlag(qualBit) ) continue;
      //      cout<<"Id Passed!"<<endl;
      //check kinematics
      if(photons[i].pt()<minPt || fabs(photons[i].eta())>maxEta) continue;
      //      cout<<"Kinematics Passed!"<<endl;
      //check if this lepton should be vetoed by request      
      bool skipThisPhoton(false);
      for(auto &vetoL : veto){
        if(vetoL.originalReference()!=photons[i].originalReference()) continue;
        skipThisPhoton=true;
        break;
      }
      if(skipThisPhoton) continue;
      //      cout<<"Not-Veto Passed!"<<endl;     
      // cross-cleaning with leptos
      bool overlapsWithPhysicsObject(false);
      for (auto& lepton : leptons) {
	if(photons[i].p4().DeltaR(lepton.p4())<0.4) overlapsWithPhysicsObject=true;
      }
      
      if(overlapsWithPhysicsObject) continue;
      //      cout<<"No overlap Passed!"<<endl;
      //photon is selected
      selPhotons.push_back(photons[i]);
    }

  //all done here
  return selPhotons;
}


//
std::vector<Jet> SelectionTool::getGoodJets(MiniEvent_t &ev, double minPt, double maxEta, std::vector<Particle> leptons,std::vector<Particle> photons) {
  std::vector<Jet> jets;
  
  for (int k=0; k<ev.nj; k++) {
    TLorentzVector jp4;
    jp4.SetPtEtaPhiM(ev.j_pt[k],ev.j_eta[k],ev.j_phi[k],ev.j_mass[k]);

    //cross clean with leptons/photons
    bool overlapsWithPhysicsObject(false);
    for (auto& lepton : leptons) {
      if(jp4.DeltaR(lepton.p4())<0.4) overlapsWithPhysicsObject=true;
    }
    for (auto& photon : photons) {
      if(jp4.DeltaR(photon.p4())<0.4) overlapsWithPhysicsObject=true;
    }
    if(overlapsWithPhysicsObject) continue;
    
    //jet kinematic selection
    if(jp4.Pt() < minPt || abs(jp4.Eta()) > maxEta) continue;

    //ECAL noise cleaning 
    //if(ev.isData && dataset.Contains("2017")){
    // if(dataset.Contains("2017")){
    //   if(abs(jp4.Eta()) > 2.7 && abs(jp4.Eta()) < 3.0)
    // 	if(jp4.Pt() < 60) continue;
    // }
    // if(dataset.Contains("2017")){
    //   if(abs(jp4.Eta()) > 2.7 && abs(jp4.Eta()) < 3.0)
    // 	if(ev.j_emf[k] > 0.5) continue;
    // }

    //flavor based on b tagging
    int flavor = 0;
    if (ev.j_btag[k]) {
      flavor = 5;
    }
    
    Jet jet(jp4, flavor, k);
    jet.setCSV(ev.j_csv[k]);
    jet.setDeepCSV(ev.j_deepcsv[k]);
    jet.setPUMVA(ev.j_pumva[k]);

    //fill jet constituents
    for (int p = 0; p < ev.npf; p++) {
      if (ev.pf_j[p] == k) {
        TLorentzVector pp4;
        pp4.SetPtEtaPhiM(ev.pf_pt[p],ev.pf_eta[p],ev.pf_phi[p],ev.pf_m[p]);
        jet.addParticle(Particle(pp4, ev.pf_c[p], ev.pf_id[p], 0, p, ev.pf_puppiWgt[p]));
        if (ev.pf_c[p] != 0) jet.addTrack(pp4, ev.pf_id[p]);
      }
    }

    //jes/jer uncertainty
    int jflav(abs(ev.j_flav[k]));
    float jecUp(0),jecDn(0);   
    jecUp=pow(1-ev.j_jerUp[k],2);
    jecDn=pow(1-ev.j_jerDn[k],2);
    for(int iunc=0; iunc<30; iunc++){
           
      //see python/miniAnalyzer_cfi.py for these
      if(iunc==6 && jflav!=21) continue; //FlavorPureGluon
      if(iunc==7 && jflav>=4)  continue; //FlavorPureQuark
      if(iunc==8 && jflav!=4)  continue; //FlavorPureCharm
      if(iunc==9 && jflav!=5)  continue; //FlavorPureGluon
      
      if(ev.j_jecUp[iunc][k]!=0) jecUp += pow(1-ev.j_jecUp[iunc][k],2);
      if(ev.j_jecDn[iunc][k]!=0) jecDn += pow(1-ev.j_jecDn[iunc][k],2);
    }
    
    jecUp=TMath::Sqrt(jecUp);
    jecDn=TMath::Sqrt(jecDn);
    jet.setScaleUnc(0.5*(jecUp+jecDn));

    if(debug_)
      cout << "Jet #" << jets.size() 
           << " pt=" << jp4.Pt() << "+/-" << jet.getScaleUnc()*jp4.Pt() << " (jec+jer)"
           << " eta=" << jp4.Eta() << " deepCSV=" << ev.j_deepcsv[k] << " flav=" << jflav << endl;
    
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

/*void SelectionTool::getGoodJets(MiniEvent_t &ev, double minPt, double maxEta, std::vector<Particle> leptons,std::vector<Particle> photons) {
  jets_.clear();
  pujets_.clear();
  for (int k=0; k<ev.nj; k++) {
    TLorentzVector jp4;
    jp4.SetPtEtaPhiM(ev.j_pt[k],ev.j_eta[k],ev.j_phi[k],ev.j_mass[k]);

    //cross clean with leptons/photons
    bool overlapsWithPhysicsObject(false);
    for (auto& lepton : leptons) {
      if(jp4.DeltaR(lepton.p4())<0.4) overlapsWithPhysicsObject=true;
    }
    for (auto& photon : photons) {
      if(jp4.DeltaR(photon.p4())<0.4) overlapsWithPhysicsObject=true;
    }
    if(overlapsWithPhysicsObject) continue;
    
    //jet kinematic selection
    if(jp4.Pt() < minPt || abs(jp4.Eta()) > maxEta) continue;

    //flavor based on b tagging
    int flavor = 0;
    if (ev.j_btag[k]) {
      flavor = 5;
    }
    
    Jet jet(jp4, flavor, k);
    jet.setCSV(ev.j_csv[k]);
    jet.setDeepCSV(ev.j_deepcsv[k]);
    jet.setPUMVA(ev.j_pumva[k]);

    //fill jet constituents
    for (int p = 0; p < ev.npf; p++) {
      if (ev.pf_j[p] == k) {
        TLorentzVector pp4;
        pp4.SetPtEtaPhiM(ev.pf_pt[p],ev.pf_eta[p],ev.pf_phi[p],ev.pf_m[p]);
        jet.addParticle(Particle(pp4, ev.pf_c[p], ev.pf_id[p], 0, p, ev.pf_puppiWgt[p]));
        if (ev.pf_c[p] != 0) jet.addTrack(pp4, ev.pf_id[p]);
      }
    }

    if(debug_) cout << "Jet #" << jets_.size() 
		    << " pt=" << jp4.Pt() << " eta=" << jp4.Eta() << " deepCSV=" << ev.j_deepcsv[k] << endl;
    
    int jid=ev.j_id[k];
    bool passLoosePu((jid>>2)&0x1);
    if(!passLoosePu) 
      pujets_.push_back(jet);
    else
      jets_.push_back(jet);
  }
  
  //additional jet-jet information
  for (unsigned int i = 0; i < jets_.size(); i++) {
    for (unsigned int j = i+1; j < jets_.size(); j++) {
      //flag jet-jet overlaps
      if (jets_[i].p4().DeltaR(jets_[j].p4()) < 0.8) {
        jets_[i].setOverlap(1);
        jets_[j].setOverlap(1);
      }
      //flag non-b jets as part of W boson candidates: flavor 0->1
      if (jets_[i].flavor()==5 or jets_[j].flavor()==5) continue;
      TLorentzVector wCand = jets_[i].p4() + jets_[j].p4();
      if (abs(wCand.M()-80.4) < 15.) {
        jets_[i].setFlavor(1);
        jets_[j].setFlavor(1);
      }
    }
  }
  }*/

//
// PARTICLE LEVEL SELECTORS
// 

//
TString SelectionTool::flagGenFinalState(MiniEvent_t &ev, std::vector<Particle> leptons, std::vector<Particle> photons) 
{
  //update current state
  genLeptons_=leptons;
  genPhotons_=photons;
  if(genLeptons_.size()==0) genLeptons_=getGenLeptons(ev,20.,2.5);
  if(genPhotons_.size()==0) genPhotons_=getGenPhotons(ev,50.,1.442);

  float maxEta(2.4);
  if(anType_==VBF) maxEta=4.7;
  genJets_=getGenJets(ev,30.,maxEta,genLeptons_,genPhotons_);

  //decide the channel
  TString chTag("");
  if(anType_==TOP)
    {
      if(genLeptons_.size()>=2) {
        int chId(abs(genLeptons_[0].id()*genLeptons_[1].id()));
        if      (chId==11*13) chTag = "EM";
        else if (chId==13*13) chTag = "MM";
        else if (chId==11*11) chTag = "EE";
      }
      else if(genLeptons_.size()==1) {
        int absid(abs(genLeptons_[0].id()));
        if      (absid==13) chTag = "M";
        else if (absid==11) chTag = "E";
      }
    }
  if(anType_==VBF)
    {
      if(genLeptons_.size()>=2)
        {
          int chId(abs(genLeptons_[0].id()*genLeptons_[1].id()));
          float mll((genLeptons_[0]+genLeptons_[1]).M());
          if(chId==13*13 && fabs(mll-91)<15) chTag="MM";
        }
      if(genPhotons_.size()>=1) chTag="A";
    }
  
  return chTag;
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
    leptons.push_back( Particle(lp4, -ev.g_id[i]/abs(ev.g_id[i]), ev.g_id[i], 0, 1));
  }
  
  return leptons;
}

//
std::vector<Particle> SelectionTool::getGenPhotons(MiniEvent_t &ev, double minPt, double maxEta){
  std::vector<Particle> photons;
  
  //loop over leptons from pseudotop producer
  for (int i = 0; i < ev.ng; i++) {
    int absid(abs(ev.g_id[i]));
    if(absid!=22) continue;

    bool passKin(ev.g_pt[i]>minPt && fabs(ev.g_eta[i])<maxEta);
    if(!passKin) continue;

    TLorentzVector p4;
    p4.SetPtEtaPhiM(ev.g_pt[i],ev.g_eta[i],ev.g_phi[i],ev.g_m[i]);
	photons.push_back( Particle(p4, 0, 22, 0, 1));
  }
  
  return photons;
}

//
std::vector<Jet> SelectionTool::getGenJets(MiniEvent_t &ev, double minPt, double maxEta, std::vector<Particle> leptons,std::vector<Particle> photons) {
  std::vector<Jet> jets;
  
  for (int i = 0; i < ev.ng; i++) {
    if (abs(ev.g_id[i])>10) continue;
    TLorentzVector jp4;
    jp4.SetPtEtaPhiM(ev.g_pt[i],ev.g_eta[i],ev.g_phi[i],ev.g_m[i]);

    //cross clean with leptons
    bool overlapsWithPhysicsObject(false);
    for (auto& lepton : leptons) {
      if(jp4.DeltaR(lepton.p4())<0.4) overlapsWithPhysicsObject=true;
    }
    for (auto& photon : photons) {
      if(jp4.DeltaR(photon.p4())<0.4) overlapsWithPhysicsObject=true;
    }
    if(overlapsWithPhysicsObject) continue;
    
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



