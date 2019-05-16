#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>

#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"
#include "TopLJets2015/TopAnalysis/interface/PhotonTrigEff.h"

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

  //READ TREE FROM FILE
  TFile *f = TFile::Open(filename);  
  TH1 *triggerList=(TH1 *)f->Get("analysis/triggerList");
  TTree *t = (TTree*)f->Get("analysis/data");
  attachToMiniEventTree(t,ev,true);
  Int_t nentries(t->GetEntriesFast());
  // if (debug) nentries = 10000; //restrict number of entries for testing
  t->GetEntry(0);

  cout << "...producing " << outname << " from " << nentries << " events" << endl;
  
   //BOOK HISTOGRAMS
  HistTool ht;
  ht.setNsyst(0);
  Float_t phoBins[]={25,50,60,65,70,75,80,85,90,100,110,120,140,150,160,170,180,190,200,220,240,260,300,500};
  Int_t nphoBins=sizeof(phoBins)/sizeof(Float_t);
  Float_t mjjBins[]={120,200,400,600,800,1000,1500,2000,5000};
  Int_t nmjjBins=sizeof(mjjBins)/sizeof(Float_t);
  ht.addHist("apt",      new TH1F("apt",      ";Photon transverse momentum [GeV];Events",nphoBins-1,phoBins));
  ht.addHist("mjj",      new TH1F("mjj",      ";Dijet invariant mass [GeV];Events",nmjjBins-1,mjjBins));
  ht.addHist("aptvsmjj", new TH2F("aptvsmjj", ";Photon transverse momentum [GeV];Dijet invariant mass [GeV];Events",nphoBins-1,phoBins,nmjjBins-1,mjjBins));

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
      float wgt(1.0);
      std::vector<double>plotwgts(1,wgt);
      double puWgt(1.0);
      if(!ev.isData){
        ht.fill("puwgtctr",0,plotwgts);
        TString period = lumi.assignRunPeriod();
        puWgt=(lumi.pileupWeight(ev.g_pu,period)[0]);
        std::vector<double>puPlotWgts(1,puWgt);
        ht.fill("puwgtctr",1,puPlotWgts);

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

      //jets
      std::vector<Jet> allJets = selector.getGoodJets(ev,50.,4.7,leptons,photons);
      float mjj(allJets.size()>=2 ? (allJets[0]+allJets[1]).M() : -1 );
      float detajj(allJets.size()>=2 ? fabs(allJets[0].eta()-allJets[1].eta()) : -1 );
	  
      //online categories
      bool passHighPtCtrlTrig(false), passLowPtCtrlTrig(false), passLowPtHighMJJCtrlTrig(false);      
      bool passHighPtTrig(false),     passLowPtHighMJJTrig(false);
      if(is2016) {
        passHighPtCtrlTrig       = selector.hasTriggerBit("HLT_Photon90_v",ev.triggerBits);                
        passHighPtTrig           = selector.hasTriggerBit("HLT_Photon175_v",ev.triggerBits);
        passLowPtCtrlTrig        = selector.hasTriggerBit("HLT_Photon50_R9Id90_HE10_IsoM_v",ev.triggerBits);
        passLowPtHighMJJCtrlTrig = selector.hasTriggerBit("HLT_Photon75_R9Id90_HE10_IsoM_v",ev.triggerBits);
        passLowPtHighMJJTrig     = selector.hasTriggerBit("HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_VBF",ev.triggerBits);
      }else{
        passHighPtCtrlTrig       = selector.hasTriggerBit("HLT_Photon150_v",ev.triggerBits);
        passHighPtTrig       = selector.hasTriggerBit("HLT_Photon200_v",ev.triggerBits);
        passLowPtCtrlTrig        = selector.hasTriggerBit("HLT_Photon50_R9Id90_HE10_IsoM_v",ev.triggerBits);
        passLowPtHighMJJCtrlTrig = selector.hasTriggerBit("HLT_Photon75_R9Id90_HE10_IsoM_v",ev.triggerBits);
        passLowPtHighMJJTrig = selector.hasTriggerBit("HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_v",ev.triggerBits);
      }

      //offline categories
      bool passHighPtOff(false), passLowPtOff(false), passLowPtHighMJJOff(false);
      if(photons[0].Pt()>200 && fabs(photons[0].Eta())<2.4)    passHighPtOff=true;
      if(photons[0].Pt()>75  && fabs(photons[0].Eta())<1.442)  passLowPtOff=true;
      if(passLowPtOff && mjj>300 && detajj>3)                  passLowPtHighMJJOff=true;
 
      std::vector<TString> cats;
      if(passHighPtCtrlTrig){
        if(passHighPtOff)                               cats.push_back("hptoff");
        if(passHighPtTrig && passHighPtOff)             cats.push_back("hpttrig_hptoff");
      }
      if(passLowPtCtrlTrig){
        if(passLowPtOff)                                cats.push_back("lptoff");
        if(passLowPtCtrlTrig && passLowPtOff)           cats.push_back("lpttrig_lptoff");
      }
      if(passLowPtHighMJJCtrlTrig){
        if(passLowPtHighMJJOff)                         cats.push_back("lpthmjjoff");
        if(passLowPtHighMJJTrig && passLowPtHighMJJOff) cats.push_back("lpthmjjtrig_lpthmjjoff");
      }
      if(cats.size()==0) continue;

      plotwgts[0]=wgt;
      ht.fill("apt",        photons[0].pt(),       plotwgts,cats);
      ht.fill("mjj",        mjj,                   plotwgts,cats);
      ht.fill2D("aptvsmjj", photons[0].pt(), mjj,  plotwgts,cats);
    }
      
  //close input file
  f->Close();
  
  //PREPARE OUTPUT
  TString baseName=gSystem->BaseName(outname); 
  TString dirName=gSystem->DirName(outname);
  TFile *fOut=TFile::Open(dirName+"/"+baseName,"RECREATE");
  fOut->cd();
  for (auto& it : ht.getPlots())  { 
    if(!it.second) continue;
    if(it.second->GetEntries()==0) continue;
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  for (auto& it : ht.get2dPlots())  { 
    if(!it.second) continue;
    if(it.second->GetEntries()==0) continue;
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  fOut->Close();
}
