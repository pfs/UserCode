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
  ht.addHist("apt",     new TH1F("apt",      ";Photon transverse momentum [GeV];Events",50,0,500));
  ht.addHist("mjj",     new TH1F("mjj",      ";Dijet invariant mass [GeV];Events",50,0,2000));
  ht.addHist("detajj",  new TH1F("detajj",   ";Dijet rapidity span;Events",50,0,10));
  ht.addHist("j1pt",    new TH1F("j1pt",     ";Jet transverse momentum [GeV];Events",50,50,250));
  ht.addHist("j2pt",    new TH1F("j2pt",     ";Jet transverse momentum [GeV];Events",50,50,250));


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

      //trigger fo data
      if(ev.isData) {        
        bool hasStdMTrigger=(selector.hasTriggerBit("HLT_IsoMu24_v",     ev.triggerBits) ||
                             selector.hasTriggerBit("HLT_IsoMu24_2p1_v", ev.triggerBits) ||
                             selector.hasTriggerBit("HLT_IsoMu27_v",     ev.triggerBits) );     
        if(filename.Contains("2017E") || filename.Contains("2017F")){
          hasStdMTrigger=selector.hasTriggerBit("HLT_IsoMu27_v",ev.triggerBits);
        }      
        if(filename.Contains("2016")) {
          hasStdMTrigger=(selector.hasTriggerBit("HLT_IsoMu24_v",ev.triggerBits) ||
                          selector.hasTriggerBit("HLT_IsoTkMu24_v",ev.triggerBits) );            
        }
        if(!hasStdMTrigger) continue;
      }

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
      float j1pt(allJets.size()>0 ? allJets[0].pt() : -1);
      float j2pt(allJets.size()>1 ? allJets[1].pt() : -1);

      std::vector<TString> cats(1,"offlinephoton");      
      if(selector.hasTriggerBit("HLT_Photon200_v",ev.triggerBits)) {
        cats.push_back("photon200");
      }
      if(selector.hasTriggerBit("HLT_Photon175_v",ev.triggerBits)) {
        cats.push_back("photon175");
      }
      if(mjj>500 && detajj>3 && fabs(photons[0].Eta())<1.442){
        cats.push_back("offlinephotonvbf");
        if(selector.hasTriggerBit("HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_v",ev.triggerBits)) {
          cats.push_back("photon75_vbf2017");
        }
        if(selector.hasTriggerBit("HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_VBF",ev.triggerBits)) {
          cats.push_back("photon75_vbf2016");
        }
      }

      plotwgts[0]=wgt;
      ht.fill("apt",    photons[0].pt(), plotwgts,cats);
      ht.fill("mjj",    mjj,             plotwgts,cats);
      ht.fill("detajj", detajj,          plotwgts,cats);
      ht.fill("j1pt",   j1pt,            plotwgts,cats);
      ht.fill("j2pt",   j2pt,            plotwgts,cats);
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
