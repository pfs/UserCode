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
#include "TopLJets2015/TopAnalysis/interface/GeneratorTools.h"
#include "TopLJets2015/TopAnalysis/interface/EfficiencyScaleFactorsWrapper.h"
#include "TopLJets2015/TopAnalysis/interface/L1PrefireEfficiencyWrapper.h"
#include "TopLJets2015/TopAnalysis/interface/SMP-19-005.h"

#include <vector>
#include <set>
#include <iostream>
#include <algorithm>
#include <string>

#include "TMath.h"

using namespace std;

//
void RunSMP19005(TString filename,
                 TString outname,
                 TH1F *normH, 
                 TH1F *genPU,
                 TString era,
                 float xsec,
                 bool debug,
                 bool skimtree) 
{
  /////////////////////
  // INITIALIZATION //
  ///////////////////
  MiniEvent_t ev;

  bool isSignal(filename.Contains("EWKAJJ"));  
  bool isData(filename.Contains("Data"));
  bool is2016(era.Contains("2016"));
  bool is2017(era.Contains("2017"));
  bool is2018(era.Contains("2018"));

  double normWgt(normH ? normH->GetBinContent(1) : 1.0);

  //efficiency corrections
  EfficiencyScaleFactorsWrapper effSF(isData,era);
  L1PrefireEfficiencyWrapper l1prefire(isData,era);

  //READ TREE FROM FILE
  TFile *f = TFile::Open(filename);  
  TH1 *triggerList=(TH1 *)f->Get("analysis/triggerList");
  TTree *t = (TTree*)f->Get("analysis/data");
  attachToMiniEventTree(t,ev,true);
  Int_t nentries(t->GetEntriesFast());
  if(debug) nentries=2000;
  t->GetEntry(0);

  cout << "...producing " << outname << " from " << nentries << " events" << endl;
  
  //PREPARE OUTPUT
  TString baseName=gSystem->BaseName(outname); 
  TString dirName=gSystem->DirName(outname);
  TFile *fOut=TFile::Open(dirName+"/"+baseName,"RECREATE");
  fOut->cd();

  //reco level 
  VBFDiscriminatorInputs recVBF;
  TTree *recT=NULL;
  Float_t trigSF(1.0),trigSFUnc(0.0),idSF(1.0),idSFUnc(0.0),l1prefireSF(1.0),l1prefireSFUnc(0.0);
  Float_t qgWgt_g(1.0),qgWgt_q(1.0);
  Float_t puWgtUp(1.0), puWgtDn(1.0);
  std::vector<WeightSysts_t> thSysts;
  Float_t thWgt[300];
  for(size_t i=0; i<300; i++) thWgt[i]=0;
  if(!isData) {
    TString proc("");
    if(isSignal)                   proc="EWKAJJ";
    if(filename.Contains("GJets")) proc="GJets";
    thSysts=getWeightSysts(f,proc);
    cout << thSysts.size() << " theory weights assuming " << proc << endl;
    TH1F *hsysts=new TH1F("thsysts",";Theory uncertainty;",thSysts.size(),0,thSysts.size());
    for(size_t i=0; i<thSysts.size(); i++) hsysts->GetXaxis()->SetBinLabel(i+1,thSysts[i].first);
    hsysts->SetDirectory(fOut);
    hsysts->Write();
  }
  if(skimtree) {
    recT=new TTree("rec","rec");
    recVBF.assignTo(recT);
    recT->SetDirectory(fOut);
    recT->Branch("trigSF",          &trigSF,         "trigSF/F");
    recT->Branch("trigSFUnc",       &trigSFUnc,      "trigSFUnc/F");
    recT->Branch("idSF",            &idSF,           "idSF/F");
    recT->Branch("idSFUnc",         &idSFUnc,        "idSFUnc/F");
    recT->Branch("l1prefireSF",     &l1prefireSF,    "l1prefireSF/F");
    recT->Branch("l1prefireSFUnc",  &l1prefireSFUnc, "l1prefireSFUnc/F");
    recT->Branch("qgWgt_g",         &qgWgt_g,        "qgWgt_g/F");
    recT->Branch("qgWgt_q",         &qgWgt_q,        "qgWgt_q/F");
    recT->Branch("puWgtUp",         &puWgtUp,        "puWgtUp/F");
    recT->Branch("puWgtDn",         &puWgtDn,        "puWgtDn/F");
    if(thSysts.size()>0) recT->Branch("thWgt", &thWgt,          Form("thWgt[%d]/F",(int)thSysts.size()));
  }

  VBFDiscriminatorInputs genVBF;
  TTree *genT=0;
  if(!isData && skimtree){
    genT=new TTree("gen","gen");
    genVBF.assignTo(genT);
    genT->SetDirectory(fOut);
  }


  std::cout << "init done" << std::endl;
  if (debug){std::cout<<"\n DEBUG MODE"<<std::endl;}

  ///////////////////////
  // LOOP OVER EVENTS //
  /////////////////////
  SelectionTool selector(filename, false, triggerList);
  LumiTools lumi(era,genPU);
  for (Int_t iev=0;iev<nentries;iev++)
    {
      t->GetEntry(iev);
      if(iev%1000==0) { printf("\r [%3.0f%%] done", 100.*(float)iev/(float)nentries); fflush(stdout); }
        
      //start weights and pu weight control
      std::vector<float> puWgts(3,1.0);
      float genWgt(0);
      if(!ev.isData){       
        TString period = lumi.assignRunPeriod();
        puWgts=lumi.pileupWeight(ev.g_pu,period);
        genWgt=(ev.g_nw>0 ? ev.g_w[0] : 1.0);
      }
      float evWgt(puWgts[0]*genWgt*normWgt);

      //select gen-level objects
      std::vector<Particle> gen_leptons=selector.getGenLeptons(ev,20.,2.5);
      std::vector<Particle> gen_photons=selector.getGenPhotons(ev,50.,2.4);
      std::vector<Jet> gen_jets = selector.getGenJets(ev,20,4.7,gen_leptons,gen_photons);

      //select offline objects
      std::vector<Particle> leptons = selector.flaggedLeptons(ev);
      leptons = selector.selLeptons(leptons,SelectionTool::MEDIUM,SelectionTool::MVA80,20,2.5);
      std::vector<Particle> photons=selector.flaggedPhotons(ev);
      photons=selector.selPhotons(photons,SelectionTool::TIGHT,leptons,75.,2.4);      
      std::vector<Jet> alljets = selector.getGoodJets(ev,20.,4.7,leptons,photons);
      
      //additional cleaning for jets
      std::vector<Jet> jets;
      for(auto &j: alljets) {
 
        int idx=j.getJetIndex();
        float emf=ev.j_emf[idx];

        //ECAL noise (2017)
        if(is2017 && fabs(j.Eta())>2.650 && fabs(j.Eta())<3.139 && emf>0.55) continue;

        //HEM 15/16 issue (2018)
        if(is2018 && j.Eta()>-3.0 && j.Eta()<-1.3 && j.Phi()>-1.57 && j.Phi()<-0.87) continue;

        jets.push_back(j);
      }


      //select trigger bits
      bool passHighPtTrig(false),passLowPtHighMJJTrig(false);
      if(is2016) {
        passHighPtTrig        = selector.hasTriggerBit("HLT_Photon175_v", ev.addTriggerBits); 
        passLowPtHighMJJTrig  = selector.hasTriggerBit("HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_VBF",ev.triggerBits);
      }else if(is2017){
        passHighPtTrig        = selector.hasTriggerBit("HLT_Photon200_v",ev.triggerBits);
        passLowPtHighMJJTrig  = selector.hasTriggerBit("HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_v",ev.triggerBits);
      }else if(is2018) {
        passHighPtTrig        = selector.hasTriggerBit("HLT_Photon200_v", ev.triggerBits);
        passLowPtHighMJJTrig  = selector.hasTriggerBit("HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_v",ev.triggerBits);
      }

      bool fillTree(false);
      if(photons.size()>0) {

        bool passRecVBF = recVBF.fill(photons[0],jets,ev,passLowPtHighMJJTrig,passHighPtTrig,evWgt,false);
        fillTree |= passRecVBF;

        //additional corrections to be stored
        if(passRecVBF && !isData){

          //pileup
          puWgtUp = (puWgts[0]!=0 ? puWgts[1]/puWgts[0] : 0); 
          puWgtDn = (puWgts[0]!=0 ? puWgts[2]/puWgts[0] : 0);

          //trigger
          //TODO

          //id
          EffCorrection_t photonSF=effSF.getPhotonSF(photons[0].pt(),photons[0].eta());
          idSF    = photonSF.first;
          idSFUnc = photonSF.second;
          
          //L1-prefire
          if(is2016 || is2017){
            EffCorrection_t l1prefSF=l1prefire.getCorrection(alljets,photons);
            l1prefireSF    = l1prefSF.first;
            l1prefireSFUnc = l1prefSF.second;
          }

          //q-g discriminator
          qgWgt_g=1.0;
          qgWgt_q=1.0;
          int tagFlav[2]={recVBF.lead_flav, recVBF.sublead_flav};
          float tagQG[2]={recVBF.lead_qg,   recVBF.sublead_qg};
          for(size_t ij=0; ij<2; ij++){
            int jflav=tagFlav[ij];
            float xqg=tagQG[ij];
            if(abs(jflav)==21) {
              double qgsf=-55.7067*pow(xqg,7) + 113.218*pow(xqg,6) -21.1421*pow(xqg,5) -99.927*pow(xqg,4) + 92.8668*pow(xqg,3) -34.3663*pow(xqg,2) + 6.27*xqg + 0.612992;
              qgWgt_g *= min(max(0.5,qgsf),2.);
            }
            else if(jflav!=0) {
              double qgsf=-0.666978*pow(xqg,3) + 0.929524*pow(xqg,2) -0.255505*xqg + 0.981581;
              qgWgt_q *= min(max(0.5,qgsf),2.);
            }

            //theory weights
            if(thSysts.size()>0 && normH) {
              for(size_t is=0; is<thSysts.size(); is++) {
                size_t idx=thSysts[is].second;
                thWgt[is] = (ev.g_w[idx]/ev.g_w[0])*(normH->GetBinContent(idx+1)/normH->GetBinContent(1));
              }
            }
          }
        }
      }
      
      //generator-level analysis
      if(!isData && gen_photons.size()>0){
        fillTree |= genVBF.fill(gen_photons[0],gen_jets,ev,true,true,genWgt*normWgt,true);
      }

      if(fillTree && skimtree){
        recT->Fill();
        if(!isData) genT->Fill();
      }
    }
  
  //close files
  f->Close();
  fOut->cd();
  if(skimtree){
    recT->Write();
    if(!isData) genT->Write();
  }
  fOut->Close();
}
