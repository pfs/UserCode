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
void RunPhotonFakeRate(TString filename,
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

  //READ TREE FROM FILE
  TFile *f = TFile::Open(filename);  
  TH1 *triggerList=(TH1 *)f->Get("analysis/triggerList");
  TTree *t = (TTree*)f->Get("analysis/data");
  attachToMiniEventTree(t,ev,true);
  Int_t nentries(t->GetEntriesFast());
  // if (debug) nentries = 10000; //restrict number of entries for testing
  t->GetEntry(0);

  cout << "...producing " << outname << " from " << nentries << " events" << endl;

  TString jetTrigs[] = {
    "HLT_PFJet40_v","HLT_PFJet60_v","HLT_PFJet80_v","HLT_PFJet140_v",
    "HLT_PFJet200_v","HLT_PFJet260_v","HLT_PFJet320_v","HLT_PFJet400_v",
    "HLT_PFJet450_v","HLT_PFJet500_v","HLT_PFJet550_v",
  };
  Int_t nJetTrig(sizeof(jetTrigs)/sizeof(TString));

  //PREPARE OUTPUT
  TString baseName=gSystem->BaseName(outname); 
  TString dirName=gSystem->DirName(outname);
  TFile *fOut=TFile::Open(dirName+"/"+baseName,"RECREATE");
  fOut->cd();
  TTree *outT=new TTree("tree","tree");
  Int_t passJetTrig(0),passid(0);
  Bool_t passHighPtTrig(false), passLowPtHighMJJTrig(false);
  Float_t wgt,vpt,veta,vphi,r9,hoe,chiso,nhiso,phoiso,sieie,njets,mindraj;
  outT->Branch("passJetTrig",          &passJetTrig,          "passJetTrig/I");
  outT->Branch("passHighPtTrig",       &passHighPtTrig,       "passHighPtTrig/O");
  outT->Branch("passLowPtHighMJJTrig", &passLowPtHighMJJTrig, "passLowPtHighMJJTrig/O");
  outT->Branch("wgt",    &wgt,    "wgt/F");
  outT->Branch("vpt",    &vpt,    "vpt/F");
  outT->Branch("veta",   &veta,   "veta/F");
  outT->Branch("r9",     &r9,     "r9/F");
  outT->Branch("hoe",    &hoe,    "hoe/F");
  outT->Branch("chiso",  &chiso,  "chiso/F");
  outT->Branch("phoiso", &phoiso, "phoiso/F");
  outT->Branch("nhiso",  &nhiso,  "nhiso/F");
  outT->Branch("vphi",   &vphi,   "vphi/F");
  outT->Branch("sieie",  &sieie,  "sieie/F");
  outT->Branch("passid", &passid, "passid/F"); 
  outT->SetDirectory(fOut);
  
  //TODO
  //no leptons
  //check first passed jet trigger and passed photon triggers
  //add cut-based ids with/without iso
  //add kinematics

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
      if(leptons.size()>0) continue;

      //select offline photons
      std::vector<Particle> photons=selector.flaggedPhotons(ev);
      if(photons.size()==0 ) continue;
      vpt=photons[0].Pt();
      veta=photons[0].Eta();
      vphi=photons[0].Phi();
      int pidx = photons[0].originalReference();
      r9       = ev.gamma_r9[pidx];
      hoe      = ev.gamma_hoe[pidx];
      //    sieie    = ev.gamma_sihih[pidx];
      chiso    = ev.gamma_chargedHadronIso[pidx];
      nhiso    = ev.gamma_neutralHadronIso[pidx];
      phoiso   = ev.gamma_photonIso[pidx];
      //passid    = (photons[0].hasQualityFlag(LOOSE) )
      //  | (photons[0].hasQualityFlag(TIGHT) << 1 )
      //  | (photons[0].hasQualityFlag(TIGHTIFNOSIHIH) << 2);
      
      //online categories      
      passHighPtTrig=false;     passLowPtHighMJJTrig=false;
      if(is2016) {
        passHighPtTrig           = selector.hasTriggerBit("HLT_Photon175_v", ev.triggerBits);
        passLowPtHighMJJTrig     = selector.hasTriggerBit("HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_VBF",ev.triggerBits);
      }else if(is2017){
        passHighPtTrig           = selector.hasTriggerBit("HLT_Photon200_v",ev.triggerBits);
        passLowPtHighMJJTrig     = selector.hasTriggerBit("HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_v",ev.triggerBits);
      }else {
        //
      }

      outT->Fill();
    }
      
  //close files
  f->Close();
  fOut->cd();
  outT->Write();
  fOut->Close();
}
