#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>

#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"
#include "TopLJets2015/TopAnalysis/interface/PhotonAnalyzers.h"

#include <vector>
#include <set>
#include <iostream>
#include <algorithm>
#include <string>

#include "TMath.h"

using namespace std;

//
void RunPhotonTrigEff(TString filename,
                      TString outname,
                      TH1F *normH, 
                      TH1F *genPU,
                      TString era,
                      bool debug) 
{
  /////////////////////
  // INITIALIZATION //
  ///////////////////
  MiniEvent_t ev;

  bool is2016(era.Contains("2016"));
  bool is2017(era.Contains("2017"));
  bool is2018(era.Contains("2018"));

  //READ TREE FROM FILE
  TFile *f = TFile::Open(filename);  
  TH1 *triggerList=(TH1 *)f->Get("analysis/triggerList");
  TTree *t = (TTree*)f->Get("analysis/data");
  attachToMiniEventTree(t,ev,true);
  Int_t nentries(t->GetEntriesFast());
  // if (debug) nentries = 10000; //restrict number of entries for testing
  t->GetEntry(0);

  cout << "...producing " << outname << " from " << nentries << " events" << endl;
  
  //PREPARE OUTPUT
  TString baseName=gSystem->BaseName(outname); 
  TString dirName=gSystem->DirName(outname);
  TFile *fOut=TFile::Open(dirName+"/"+baseName,"RECREATE");
  fOut->cd();
  TTree *outT=new TTree("tree","tree");
  Bool_t passHighPtCtrlTrig(false), passLowPtCtrlTrig(false), passLowPtHighMJJCtrlTrig(false), passHighPtTrig(false), passLowPtHighMJJTrig(false),passLowPtCtrlTrigNoId(false);
  Bool_t lowPtHighMJJCtrTrigActive(true);
  Float_t wgt,vpt,veta,vphi,r9,hoe,chiso,nhiso,phoiso,mjj,detajj,j1pt,j1eta,j2pt,j2eta;
  Int_t VBFreq(0);
  outT->Branch("passHighPtCtrlTrig",        &passHighPtCtrlTrig,        "passHighPtCtrlTrig/O");
  outT->Branch("passLowPtCtrlTrig",         &passLowPtCtrlTrig,         "passLowPtCtrlTrig/O");
  outT->Branch("passLowPtHighMJJCtrlTrig",  &passLowPtHighMJJCtrlTrig,  "passLowPtHighMJJCtrlTrig/O");
  outT->Branch("passHighPtTrig",            &passHighPtTrig,            "passHighPtTrig/O");
  outT->Branch("passLowPtHighMJJTrig",      &passLowPtHighMJJTrig,      "passLowPtHighMJJTrig/O");
  outT->Branch("passLowPtCtrlTrigNoId",     &passLowPtCtrlTrigNoId,     "passLowPtCtrlTrigNoId/O");
  outT->Branch("lowPtHighMJJCtrTrigActive", &lowPtHighMJJCtrTrigActive, "lowPtHighMJJCtrTrigActive/O");
  outT->Branch("wgt",    &wgt,    "wgt/F");
  outT->Branch("vpt",    &vpt,    "vpt/F");
  outT->Branch("veta",   &veta,   "veta/F");
  outT->Branch("r9",     &r9,     "r9/F");
  outT->Branch("hoe",    &hoe,    "hoe/F");
  outT->Branch("chiso",  &chiso,  "chiso/F");
  outT->Branch("phoiso", &phoiso, "phoiso/F");
  outT->Branch("nhiso",  &nhiso,  "nhiso/F");
  outT->Branch("vphi",   &vphi,   "vphi/F");
  outT->Branch("mjj",    &mjj,    "mjj/F");
  outT->Branch("detajj", &detajj, "detajj/F");
  outT->Branch("j1pt",   &j1pt,   "j1pt/F");
  outT->Branch("j1eta",  &j1eta,  "j1eta/F");
  outT->Branch("j2pt",   &j2pt,   "j2pt/F");
  outT->Branch("j2eta",  &j2eta,  "j2eta/F");
  outT->Branch("VBFreq", &VBFreq, "VBFreq/I");
  outT->SetDirectory(fOut);

  std::cout << "init done" << std::endl;
  if (debug){std::cout<<"\n DEBUG MODE"<<std::endl;}

  ///////////////////////
  // LOOP OVER EVENTS //
  /////////////////////
  
  //EVENT SELECTION WRAPPER
  SelectionTool selector(filename, false, triggerList);
  LumiTools lumi(era,genPU);

  for (Int_t iev=0;iev<nentries;iev++)
    {
      t->GetEntry(iev);
      if(iev%1000==0) { printf("\r [%3.0f%%] done", 100.*(float)iev/(float)nentries); fflush(stdout); }
        
      //start weights and pu weight control
      wgt=1.0;
      std::vector<double>plotwgts(1,wgt);
      if(!ev.isData){       
        TString period = lumi.assignRunPeriod();
        double puWgt=(lumi.pileupWeight(ev.g_pu,period)[0]);
        wgt=(normH? normH->GetBinContent(1) : 1.0);
        wgt*=puWgt;
        wgt*=(ev.g_nw>0 ? ev.g_w[0] : 1.0);
      }  
	
      //leptons
      std::vector<Particle> leptons = selector.flaggedLeptons(ev);
      leptons = selector.selLeptons(leptons,SelectionTool::MEDIUM,SelectionTool::MVA80,20,2.5);

      //select offline photons
      std::vector<Particle> photons=selector.flaggedPhotons(ev);
      photons=selector.selPhotons(photons,SelectionTool::TIGHT,leptons,50.,2.4);
      if(photons.size()==0 ) continue;
      vpt=photons[0].Pt();
      veta=photons[0].Eta();
      vphi=photons[0].Phi();
      int pidx = photons[0].originalReference();
      r9       = ev.gamma_r9[pidx];
      hoe      = ev.gamma_hoe[pidx];
      chiso    = ev.gamma_chargedHadronIso[pidx];
      nhiso    = ev.gamma_neutralHadronIso[pidx];
      phoiso   = ev.gamma_photonIso[pidx];


      //jets
      std::vector<Jet> allJets = selector.getGoodJets(ev,50.,4.7,leptons,photons);
      mjj=(allJets.size()>=2 ? (allJets[0]+allJets[1]).M() : -1 );
      detajj=(allJets.size()>=2 ? fabs(allJets[0].eta()-allJets[1].eta()) : -1 );
      j1pt=(allJets.size()>=1 ? allJets[0].Pt() : -999);
      j2pt=(allJets.size()>=2 ? allJets[1].Pt() : -999);
      j1eta=(allJets.size()>=1 ? allJets[0].Eta() : -999);
      j2eta=(allJets.size()>=2 ? allJets[1].Eta() : -999);

      //in 2016,2017 there is no cross cleaning between the photon and jet candidates...
      VBFreq=0;
      if(allJets.size()>=2) {
        TLorentzVector vj1(photons[0]+allJets[0]);
        float detavj1(fabs(photons[0].Eta()-allJets[0].Eta()));
        TLorentzVector vj2(photons[0]+allJets[1]);
        float detavj2(fabs(photons[0].Eta()-allJets[1].Eta()));
        if(mjj>500 && detajj>3)    VBFreq |= 1;
        if(vj1.M()>500 && detavj1) VBFreq |= (1<<1);
        if(vj2.M()>500 && detavj2) VBFreq |= (1<<2);
      }
      
      //online categories
      passHighPtCtrlTrig=false; passLowPtCtrlTrig=false; passLowPtHighMJJCtrlTrig=false;      
      passHighPtTrig=false;     passLowPtHighMJJTrig=false;
      passLowPtCtrlTrigNoId=false;
      lowPtHighMJJCtrTrigActive=true;
      if(is2016) {
        passLowPtCtrlTrigNoId    = selector.hasTriggerBit("HLT_Photon50_v",  ev.triggerBits);                
        passHighPtCtrlTrig       = selector.hasTriggerBit("HLT_Photon90_v",  ev.triggerBits);                
        passHighPtTrig           = selector.hasTriggerBit("HLT_Photon175_v", ev.triggerBits);
        passLowPtCtrlTrig        = selector.hasTriggerBit("HLT_Photon50_R9Id90_HE10_IsoM_v",ev.triggerBits);
        passLowPtHighMJJCtrlTrig = selector.hasTriggerBit("HLT_Photon75_R9Id90_HE10_IsoM_v",ev.triggerBits);
        passLowPtHighMJJTrig     = selector.hasTriggerBit("HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_VBF",ev.triggerBits);
      }else if(is2017){
        if(ev.run<305405 || ev.run>306460) lowPtHighMJJCtrTrigActive=false;
        passHighPtCtrlTrig       = selector.hasTriggerBit("HLT_Photon150_v",ev.triggerBits);
        passHighPtTrig           = selector.hasTriggerBit("HLT_Photon200_v",ev.triggerBits);
        passLowPtCtrlTrig        = selector.hasTriggerBit("HLT_Photon50_R9Id90_HE10_IsoM_v",ev.triggerBits);
        passLowPtHighMJJCtrlTrig = selector.hasTriggerBit("HLT_Photon75_R9Id90_HE10_IsoM_v",ev.triggerBits);
        passLowPtHighMJJTrig     = selector.hasTriggerBit("HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_v",ev.triggerBits);
      }else if(is2018) {
        passLowPtCtrlTrigNoId    = selector.hasTriggerBit("HLT_Photon50_v",  ev.triggerBits);                
        passHighPtCtrlTrig       = selector.hasTriggerBit("HLT_Photon90_v",  ev.triggerBits);                
        passHighPtTrig           = selector.hasTriggerBit("HLT_Photon200_v", ev.triggerBits);
        passLowPtCtrlTrig        = selector.hasTriggerBit("HLT_Photon50_R9Id90_HE10_IsoM_v",ev.triggerBits);
        passLowPtHighMJJCtrlTrig = selector.hasTriggerBit("HLT_Photon75_R9Id90_HE10_IsoM_v",ev.triggerBits);
        passLowPtHighMJJTrig     = selector.hasTriggerBit("HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_v",ev.triggerBits);
      }
      else{
        continue;
      }


      outT->Fill();
    }
      
  //close files
  f->Close();
  fOut->cd();
  outT->Write();
  fOut->Close();
}
