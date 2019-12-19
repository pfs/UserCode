#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>

#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"
#include "TopLJets2015/TopAnalysis/interface/TopSummer2019.h"
#include "TopLJets2015/TopAnalysis/interface/L1PrefireEfficiencyWrapper.h"

#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>

#include "TMath.h"

using namespace std;


//
void RunTopSummer2019(const TString in_fname,
                      TString outname,
                      TH1F *normH, 
                      TH1F *genPU, 
                      TString era,
                      Bool_t debug)
{
  /////////////////////
  // INITIALIZATION //
  ///////////////////
  
  bool isLowPUrun(false);
  if(in_fname.Contains("2017H")) {
    isLowPUrun=true;
    cout << "Running with low PU run settings" << endl;
  }

  //preselection cuts to apply
  float minLeptonPt( isLowPUrun ? 20. : 27);
  size_t minJetMultiplicity(4);

  //CORRECTIONS: LUMINOSITY+PILEUP
  LumiTools lumi(era,genPU);
  std::map<Int_t,Float_t> lumiPerRun=lumi.lumiPerRun();
  
  //CORRECTIONS: LEPTON EFFICIENCIES
  EfficiencyScaleFactorsWrapper lepEffH(in_fname.Contains("Data13TeV"),era);

  //CORRECTIONS: L1-prefire 
  L1PrefireEfficiencyWrapper l1PrefireWR(in_fname.Contains("Data13TeV"),era);
  
  //CORRECTIONS: B-TAG CALIBRATION
  BTagSFUtil btvSF(era,BTagEntry::OperatingPoint::OP_MEDIUM,"",0);

  //PREPARE OUTPUT (BOOK SOME HISTOGRAMS)
  TString baseName=gSystem->BaseName(outname); 
  TString dirName=gSystem->DirName(outname);
  TFile *fOut=TFile::Open(dirName+"/"+baseName,"RECREATE");
  fOut->cd();
  HistTool ht;
  ht.setNsyst(0);
  ht.addHist("nvtx",         new TH1F("nvtx",        ";Vertex multiplicity;Events",50,0,100));
  ht.addHist("njets",        new TH1F("njets",       ";Jet multiplicity;Events",10,0,10));
  ht.addHist("mlb",          new TH1F("mlb",         ";m(l,b) [GeV];Events",20,0,250)); 
  ht.addHist("nprotons",     new TH1F("nprotons",    ";Proton multiplicity; Events",6,0,6) );
  ht.addHist("csi",          new TH1F("csi",         ";#xi = #deltap/p; Events",50,0,0.3) );
  ht.addHist("x",            new TH1F("x",           ";x  [cm]; Events",50,0,25) );
  ht.addHist("ratevsrun",    new TH1F("ratevsrun",   ";Run number; #sigma [pb]",int(lumiPerRun.size()),0,float(lumiPerRun.size())));
  int i=0;
  for(auto key : lumiPerRun) {
    i++;
    ht.getPlots()["ratevsrun"]->GetXaxis()->SetBinLabel(i,Form("%d",key.first));
  }


  //INPUT
  MiniEvent_t ev;
  TFile *f = TFile::Open(in_fname);
  if(f==NULL || f->IsZombie()) {
    cout << "Corrupted or missing file " << in_fname << endl;
    return;
  }

  TH1 *triggerList=(TH1 *)f->Get("analysis/triggerList");
  TTree *t = (TTree*)f->Get("analysis/data");
  attachToMiniEventTree(t,ev,true);
  Int_t nentries(t->GetEntriesFast());
  if (debug) nentries = min(100000,nentries); //restrict number of entries for testing
  t->GetEntry(0);
  cout << "...producing " << outname << " from " << nentries << " events" << endl;

  //EVENT SELECTION WRAPPER (GETS LISTS OF PHYSICS OBJECTS FROM THE INPUT)
  SelectionTool selector(in_fname, false, triggerList);
  
  //EVENT LOOP
  //select mu+>=4 jets events triggered by a single muon trigger
  for (Int_t iev=0;iev<nentries;iev++)
    {
      t->GetEntry(iev);
      if(iev%1000==0) { printf("\r [%3.0f%%] done", 100.*(float)iev/(float)nentries); fflush(stdout); }

      //trigger
      bool hasMTrigger(false);
      if(era.Contains("2016")) hasMTrigger=(selector.hasTriggerBit("HLT_IsoMu24_v", ev.triggerBits) );     
      if(era.Contains("2017")) {
        if(isLowPUrun) hasMTrigger=(selector.hasTriggerBit("HLT_HIMu12_v",  ev.addTriggerBits) );   
        else           hasMTrigger=(selector.hasTriggerBit("HLT_IsoMu27_v", ev.triggerBits) );   
      }
      cout << hasMTrigger << endl;
      if(!hasMTrigger) continue;

      //select one offline muon
      std::vector<Particle> leptons = selector.flaggedLeptons(ev);     
      SelectionTool::QualityFlags muId(SelectionTool::TIGHT);
      if(isLowPUrun) muId=SelectionTool::LOOSE;
      leptons = selector.selLeptons(leptons,muId,SelectionTool::MVA90,minLeptonPt,2.1);
      cout << "nl=" << leptons.size() << endl;
      if(leptons.size()!=1) continue;
      if(leptons[0].id()!=13) continue;

      //select jets
      btvSF.addBTagDecisions(ev);
      if(!ev.isData) btvSF.updateBTagDecisions(ev);      
      std::vector<Jet> allJets = selector.getGoodJets(ev,30.,2.4,leptons,{});
      bool passJets(allJets.size()<minJetMultiplicity);

      //met
      TLorentzVector met(0,0,0,0);
      met.SetPtEtaPhiM(ev.met_pt,0,ev.met_phi,0.);

      //event weight
      float evWgt(1.0);
      
      //data specific: check event rates after selection
      if(ev.isData){
        std::map<Int_t,Float_t>::iterator rIt=lumiPerRun.find(ev.run);
        if(rIt!=lumiPerRun.end()){
          int runBin=std::distance(lumiPerRun.begin(),rIt);
          float lumi=1./rIt->second;
          ht.fill("ratevsrun",runBin,lumi,"inc");
        }else{
          cout << "[Warning] Unable to find run=" << ev.run << endl;
        }
      }

      //MC specific: compute event weight
      if (!ev.isData) {

        float normWgt(normH? normH->GetBinContent(1) : 1.0);        
        TString period = lumi.assignRunPeriod();
        double puWgt(lumi.pileupWeight(ev.g_pu,period)[0]);
        EffCorrection_t selSF(1.0,0.0);// = lepEffH.getOfflineCorrection(leptons[0], period);
        EffCorrection_t l1prefireProb=l1PrefireWR.getCorrection(allJets,{});

        evWgt  = normWgt*puWgt*selSF.first*l1prefireProb.first;
        evWgt *= (ev.g_nw>0 ? ev.g_w[0] : 1.0);        
      }
      

      ht.fill("nvtx",       ev.nvtx,        evWgt, "inc");
      ht.fill("njets",      allJets.size(), evWgt, "inc");

      if(!passJets) continue;

      //lepton-b systems
      for(size_t ij=0; ij<allJets.size(); ij++) 
        {
          int idx=allJets[ij].getJetIndex();
          bool passBtag(ev.j_btag[idx]>0);
          if(!passBtag) continue;

          float mlb( (leptons[0]+allJets[ij]).M() );
          std::vector<TString> tags={"inc",leptons[0].charge()>0 ? "plus" : "minus"};
          ht.fill("mlb",mlb,evWgt,tags);
        }

      //roman pots
      int nprotons23(0), nprotons123(0);
      int ntrks( isLowPUrun ? ev.nppstrk : ev.nfwdtrk );
      for (int ift=0; ift<ntrks; ift++) {

        //single pot reconstruction
        if(!isLowPUrun && ev.fwdtrk_method[ift]!=0) continue;

        //only near (pixels) detectors
        const unsigned short pot_raw_id = (isLowPUrun ? ev.ppstrk_pot[ift] : ev.fwdtrk_pot[ift]);
        if (pot_raw_id!=23 && pot_raw_id!=123) continue;            
          
        nprotons23 += (pot_raw_id==23);
        nprotons123 += (pot_raw_id==123);
        std::vector<TString> tags={"inc",Form("rp%d",pot_raw_id)};

        float xi= (isLowPUrun ? 0.               : ev.fwdtrk_xi[ift]);
        float x=  (isLowPUrun ? ev.ppstrk_x[ift] :  0. );

        ht.fill("csi",     xi,                    evWgt,tags);
        ht.fill("x",       x,                     evWgt,tags);
        ht.fill("nprotons",nprotons23+nprotons123,evWgt,"inc");
        ht.fill("nprotons",nprotons23,            evWgt,"rp23");
        ht.fill("nprotons",nprotons123,           evWgt,"rp123");
      }
 
    }
  
  //close input file
  f->Close();
  
  //save histos to file  
  fOut->cd();
  for (auto& it : ht.getPlots())  { 
    if(it.second->GetEntries()==0) continue;
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  for (auto& it : ht.get2dPlots())  { 
    if(it.second->GetEntries()==0) continue;
    it.second->SetDirectory(fOut); it.second->Write(); 
  }  
  fOut->Close();
}
