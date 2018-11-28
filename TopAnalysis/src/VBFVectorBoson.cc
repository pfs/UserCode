#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>

#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"
#include "TopLJets2015/TopAnalysis/interface/VBFVectorBoson.h"
#include "TopLJets2015/TopAnalysis/interface/EfficiencyScaleFactorsWrapper.h"


#include <vector>
#include <set>
#include <iostream>
#include <algorithm>
#include "TMath.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

using namespace std;



//
void VBFVectorBoson::runAnalysis()
{

  //specifics for JetHT
  if(selector_->isJetHTPD()) {
    cout <<"[VBFVectorBoson::RunVBFVectorBoson] this is JetHT data, adapting cuts and triggers" <<endl;
    lowMJJCut_=50;
    highMJJCut_=50;
    lowVPtPhotonTrigs_.clear();
    highVPtPhotonTrigs_.clear();
    TString jetTrigs[] = {"HLT_PFJet40_v","HLT_PFJet60_v","HLT_PFJet80_v","HLT_PFJet140_v",
                          "HLT_PFJet200_v","HLT_PFJet260_v","HLT_PFJet320_v","HLT_PFJet400_v",
                          "HLT_PFJet450_v","HLT_PFJet500_v","HLT_PFJet550_v","HLT_PFJetFwd40_v",
                          "HLT_PFJetFwd60_v","HLT_PFJetFwd80_v","HLT_PFJetFwd140_v","HLT_PFJetFwd200_v",
                          "HLT_PFJetFwd260_v","HLT_PFJetFwd320_v","HLT_PFJetFwd400_v","HLT_PFJetFwd450_v",
                          "HLT_PFJetFwd500_v"};
    size_t nJetTrigs=sizeof(jetTrigs)/sizeof(TString);
    lowVPtPhotonTrigs_.assign(jetTrigs,jetTrigs+nJetTrigs);
    highVPtPhotonTrigs_.assign(jetTrigs,jetTrigs+nJetTrigs);
  }

  //TMVA configuration
  //FIXME: this may need to change for the new training
  std::map<TString,TMVA::Reader *> readers;
  std::map<TString,TGraph *> mvaCDFinv;
  TString method("BDT_VBF0HighMJJ");
  TString weightFile("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/VBF_weights/HighMJJ_BDT_VBF0HighMJJ.weights.xml");
  gSystem->ExpandPathName(weightFile);
  readers[method]=new TMVA::Reader( "!Color:!Silent" );
  readers[method]->AddVariable("D",             &vbfVars_.D);
  readers[method]->AddVariable("C",             &vbfVars_.C);
  readers[method]->AddVariable("circularity",   &vbfVars_.circularity);
  readers[method]->AddVariable("balance",       &vbfVars_.balance);
  readers[method]->AddVariable("jjpt",          &vbfVars_.jjpt);
  readers[method]->AddVariable("ystar",         &vbfVars_.ystar);
  readers[method]->AddVariable("dphijj",        &vbfVars_.dphijj);
  readers[method]->AddVariable("dphivj0",       &vbfVars_.dphivj0);
  readers[method]->AddVariable("j_pt[1]",       &vbfVars_.subleadj_pt);
  readers[method]->AddVariable("j_c2_00[0]",    &vbfVars_.leadj_c2_02);
  readers[method]->AddVariable("j_c2_00[1]",    &vbfVars_.subleadj_c2_02);
  readers[method]->AddVariable("j_gawidth[0]",  &vbfVars_.leadj_gawidth);
  readers[method]->AddVariable("mjj",           &vbfVars_.mjj);
  readers[method]->AddVariable("j_qg[1]",       &vbfVars_.subleadj_qg);
  readers[method]->BookMVA(method,weightFile);
  
  method="BDT_VBF0LowMJJ";
  weightFile="${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/VBF_weights/LowMJJ_BDT_VBF0LowMJJ.weights.xml";
  gSystem->ExpandPathName(weightFile);
  readers[method]=new TMVA::Reader( "!Color:!Silent" );
  readers[method]->AddVariable("D",             &vbfVars_.D);
  readers[method]->AddVariable("C",             &vbfVars_.C);
  readers[method]->AddVariable("balance",       &vbfVars_.balance);
  readers[method]->AddVariable("jjpt",          &vbfVars_.jjpt);
  readers[method]->AddVariable("ystar",         &vbfVars_.ystar);
  readers[method]->AddVariable("dphijj",        &vbfVars_.dphijj);
  readers[method]->AddVariable("dphivj0",       &vbfVars_.dphivj0);
  readers[method]->AddVariable("dphivj1",       &vbfVars_.dphivj1);
  readers[method]->AddVariable("j_pt[1]",       &vbfVars_.subleadj_pt);
  readers[method]->AddVariable("detajj",        &vbfVars_.detajj);
  readers[method]->AddVariable("j_c2_00[0]",    &vbfVars_.leadj_c2_02);
  readers[method]->AddVariable("j_c2_00[1]",    &vbfVars_.subleadj_c2_02);
  readers[method]->AddVariable("mjj",           &vbfVars_.mjj);
  readers[method]->AddVariable("j_qg[0]",       &vbfVars_.leadj_qg);
  readers[method]->AddVariable("j_qg[1]",       &vbfVars_.subleadj_qg);
  readers[method]->BookMVA(method,weightFile);

  //read the transformations based on CDF^{-1}
  weightFile="${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/VBF_weights/inverse_cdfs.root";
  gSystem->ExpandPathName(weightFile);
  TFile *fcdf=TFile::Open(weightFile);
  for(std::map<TString,TMVA::Reader *>::iterator it=readers.begin(); it!=readers.end(); it++) {
    TString key=it->first;
    mvaCDFinv[key]=(TGraph *)fcdf->Get(key+"_cdfinv");
  }

  ///////////////////////
  // LOOP OVER EVENTS //
  /////////////////////
  for (Int_t iev=0;iev<nentries_;iev++)
    {
      t_->GetEntry(iev);
      if(debug_) cout << "Number of event: "<<iev<<endl;
      if(iev%10000==0) printf ("\r [%3.0f%%] done", 100.*(float)iev/(float)nentries_);

      std::vector<double>plotwgts(1,1.0);
      ht_->fill("qscale",ev_.g_qscale,plotwgts);
      
      //assign randomly a run period
      TString period = lumi_->assignRunPeriod();
      
      ///////////////////////////
      // RECO LEVEL SELECTION //
      /////////////////////////

      TString chTag("");
      TLorentzVector boson(0,0,0,0);
      float bosonScaleUnc(0.);
      bool isBosonTrigSafe(false);

      //leptons
      std::vector<Particle> leptons = selector_->flaggedLeptons(ev_);     
      leptons = selector_->selLeptons(leptons,SelectionTool::MEDIUM,SelectionTool::MVA80,20,2.5);
      if(leptons.size()==1) continue;
      if(leptons.size()>=2) {
        if(leptons[0].Pt()<30) continue;
        int selCode=abs(leptons[0].id())*abs(leptons[1].id());
        if(selCode==11*13) continue;
        if( leptons[0].charge()*leptons[1].charge()>0 ) continue;
        boson   += leptons[0];
        boson   += leptons[1];
        bosonScaleUnc= TMath::Sqrt( pow(leptons[0].Pt()*leptons[1].scaleUnc(),2)+
                                    pow(leptons[1].Pt()*leptons[0].scaleUnc(),2) )/boson.Pt();
        float mll=(leptons[0]+leptons[1]).M();
        if( fabs(mll-91)>zMassWindow_ ) continue;

        //check trigger
        if(selCode==11*11) {
          if(ev_.isData && !selector_->isDoubleEGPD() ) continue;
          bool passTrigger=(selector_->hasTriggerBit("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v",    ev_.triggerBits) ||
                            selector_->hasTriggerBit("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", ev_.triggerBits) );
          if(!passTrigger) continue;
          chTag="EE";
          isBosonTrigSafe=true;
        }
        else if(selCode==13*13) {
          if(ev_.isData && !selector_->isSingleMuonPD() ) continue; 
          bool passTrigger=(selector_->hasTriggerBit("HLT_IsoMu24_v",     ev_.triggerBits) ||
                            selector_->hasTriggerBit("HLT_IsoMu24_2p1_v", ev_.triggerBits) ||
                            selector_->hasTriggerBit("HLT_IsoMu27_v",     ev_.triggerBits) );     
          if(filename_.Contains("2017E") || filename_.Contains("2017F")){
            passTrigger=selector_->hasTriggerBit("HLT_IsoMu27_v",ev_.triggerBits);
          }
          if(!passTrigger) continue;
          chTag="MM";
          isBosonTrigSafe=true;
        }
      }
      
      //photons
      photons_.clear();
      relaxedTightPhotons_.clear();
      tmpPhotons_.clear();
      sihih_ = 0, chiso_ = 0, r9_ = 0, hoe_ = 0;
      std::vector<Particle> fakeACR,tightACR,QCDTemplate;
      std::map<TString, int> mults;
      mults["loosefake"]   = 0;
      mults["tightfake"]   = 0;
      mults["looseprompt"] = 0;
      mults["tightprompt"] = 0;
      std::vector<Particle> allPhotons = selector_->flaggedPhotons(ev_);      
      if(chTag=="") {

        //in data photons can only come from these PDs
        if(ev_.isData && !selector_->isPhotonPD() && !selector_->isJetHTPD() ) continue;
        
        //select photons
        photons_              = selector_->selPhotons(allPhotons, SelectionTool::TIGHT, leptons);   
        relaxedTightPhotons_  = selector_->selPhotons(allPhotons,SelectionTool::QCDTEMP, leptons);
        tmpPhotons_           = selector_->selPhotons(allPhotons,SelectionTool::RELAXEDTIGHT, leptons);     

        //if running for control region alter the photon defintion
        if(CR_) {
          if(SRfake_) {
            photons_.clear();
          }
          else {
            if(!QCDTemp_) photons_ = selector_->selPhotons(allPhotons,SelectionTool::CONTROL, leptons);
            else          photons_ = tmpPhotons_;                            
          }
        }
       
        //assign relevant information
        if(photons_.size()>0) {
          boson         = photons_[0];
          bosonScaleUnc = photons_[0].scaleUnc()/photons_[0].Pt();
          int pidx      = photons_[0].originalReference();
          sihih_         = ev_.gamma_sieie[pidx];
          chiso_         = ev_.gamma_chargedHadronIso[pidx];         
          r9_            = ev_.gamma_r9[pidx];
          hoe_           = ev_.gamma_hoe[pidx];
          if(!CR_ && applyTrigSafePhoton_) {
            bool eveto       = ev_.gamma_idFlags[pidx] & 0x1;
            bool pixelseed   = ((ev_.gamma_idFlags[pidx]>>1) & 0x1);
            isBosonTrigSafe = (r9_>0.9 && eveto && !pixelseed && hoe_<0.01); 
          }else{
            isBosonTrigSafe=true;
          }          
          
          //trigger will be verified later when dividing the categories
          chTag="A";
        }
        
        //MC truth for the selected photons
        for(auto a : photons_) {
          int idx=a.originalReference();
          bool isPrompt(ev_.gamma_isPromptFinalState[idx]);
          if (selector_->isFakePhoton(ev_,idx)) {
            fakeACR.push_back(a);
            mults[isPrompt ? "looseprompt" : "loosefake"]++;	 
          } else if (a.hasQualityFlag(SelectionTool::TIGHT)){
            tightACR.push_back(a);
            mults[isPrompt ? "tightprompt" : "tightfake"]++;
          }
        }

      }

      //no interesting channel has been found
      if(chTag=="") continue;
      if(chTag=="A" && vetoPromptPhotons_) {
        if(ev_.gamma_isPromptFinalState[ photons_[0].originalReference() ] ) continue;
      }

      //jet selection
      std::vector<Jet> alljets = selector_->getGoodJets(ev_,30.,4.7,leptons,photons_);
      std::vector<Jet> jets;
      for(auto j : alljets) {

        int idx=j.getJetIndex();
        if(cleanEENoise_ && fabs(j.Eta())>2.7 && fabs(j.Eta())<3 && ev_.j_emf[idx]>0.55) continue;

        int jid=ev_.j_id[idx];
        bool passPu((jid>>jetPuId_)&0x1);
        bool passLoosePu((jid>>2)&0x1);

        //for the control region the second jet is required to fail a loose PU ID
        if(!CR_ && !passPu) continue;
        if(CR_ && jets.size()==1 && passLoosePu) continue;

        if(!passPu) continue;
        jets.push_back(j);
      }


      //variables
      vbfVars_.fillDiscriminatorVariables(boson,jets,ev_);
      mindrl_=9999.;
      for(auto &l: leptons) mindrl_ = min(l.DeltaR(boson),Double_t(mindrl_));
      mindrj_=9999.;
      for(auto &j: jets)    mindrj_ = min(j.DeltaR(boson),Double_t(mindrj_));


      //
      // ASSIGN EVENT CATEGORY
      //
      category_.reset();
      bool passLowVPtTrig(true),passHighVPtTrig(true);
      std::vector<bool> cat(8,false);
      if(chTag == "EE")  cat[0] = true;
      if(chTag == "MM")  cat[1] = true;
      if(chTag == "A")   {
        cat[2] = true;
        passLowVPtTrig=false;
        for(auto tname : lowVPtPhotonTrigs_) passLowVPtTrig |= selector_->hasTriggerBit(tname,ev_.triggerBits);
        passHighVPtTrig=false;
        for(auto tname : highVPtPhotonTrigs_) passHighVPtTrig |= selector_->hasTriggerBit(tname,ev_.triggerBits);
      }
      cat[3]  = ( boson.Pt()>lowVPtCut_  && boson.Pt()<=highVPtCut_ && isBosonTrigSafe);
      cat[4]  = ( boson.Pt()>highVPtCut_ && isBosonTrigSafe);
      if(jets.size()>=2 && jets[0].Pt()>leadJetPt_ && jets[1].Pt()>subLeadJetPt_) {
        cat[5]  =  (vbfVars_.mjj>lowMJJCut_ && vbfVars_.mjj<=highMJJCut_);
        cat[6]  =  (vbfVars_.mjj>highMJJCut_);
      }
      category_.set(cat);
      std::vector<TString> chTags( category_.getChannelTags() );

      //veto untriggerable events for photons
      if(chTag=="A" && cat[3]) {
        if(fabs(boson.Rapidity())>lowVPtMaxRapCut_
           || vbfVars_.mjj<highMJJCut_
           || vbfVars_.detajj<lowVPtDetaJJCut_)
          chTags.clear();
      }
      if(chTags.size()==0) continue;

      TString baseCategory(chTags[chTags.size()-1]);

      //evaluate discriminator MVA for categories of interest
      //FIXME: this probably needs to be modified for the new training
      vbfmva_ = -1000;
      if(cat[4] || cat[5]) {
        TString key(cat[4] ?"BDT_VBF0HighMJJ":"BDT_VBF0LowMJJ");
        vbfmva_ = readers[key]->EvaluateMVA(key);
        if(mvaCDFinv[key]) vbfmva_=max(0.,mvaCDFinv[key]->Eval(vbfmva_));
        if(doBlindAnalysis_ && ev_.isData && vbfmva_>0.8) vbfmva_=-1000;
      }

      ////////////////////
      // EVENT WEIGHTS //
      //////////////////
      float wgt(1.0);
      std::vector<float>puWgts(3,1.0);
      EffCorrection_t trigSF(1.0,0.),selSF(1.0,0.),l1prefireProb(1.0,0.);
      if (!ev_.isData) {

        // norm weight
        wgt  = (normH_? normH_->GetBinContent(1) : 1.0);
            
        // pu weight
        ht_->fill("puwgtctr",0,plotwgts);
        puWgts=lumi_->pileupWeight(ev_.g_pu,period);
        std::vector<double>puPlotWgts(1,puWgts[0]);
        ht_->fill("puwgtctr",1,puPlotWgts);

        //L1 prefire probability
        l1prefireProb=l1PrefireWR_->getCorrection(jets,photons_);
        wgt *= l1prefireProb.first;

        // photon trigger*selection weights        
        if(chTag=="A")
          {
            trigSF = gammaEffWR_->getTriggerCorrection({},photons_,{}, period);
            selSF  = gammaEffWR_->getOfflineCorrection(22,photons_[0].pt(),photons_[0].eta(), period);
          }
        else
          {
            trigSF = gammaEffWR_->getTriggerCorrection(leptons,{},{}, period);
            selSF  = gammaEffWR_->getOfflineCorrection(leptons[0].id(),leptons[0].pt(),leptons[0].eta(), period);
            EffCorrection_t sel2SF=gammaEffWR_->getOfflineCorrection(leptons[1].id(),leptons[1].pt(),leptons[1].eta(), period);
            selSF.first *= sel2SF.first;
            selSF.second = TMath::Sqrt( pow(selSF.second,2)+pow(sel2SF.second,2) );
          }
        wgt *= puWgts[0]*trigSF.first*selSF.first;
        
        // generator level weights
        wgt *= (ev_.g_nw>0 ? ev_.g_w[0] : 1.0);
        
        //update weight for plotter
        plotwgts[0]=wgt;
      }
      
      //fake rate
      if(ev_.isData && chTag=="A" && SRfake_ && (cat[4]||cat[5])) {
        cout << "Fake Rate will be applied! " <<endl;
        TString FRcat(cat[5] ? "LowMJJ" : "HighMJJ");
        plotwgts[0]*=fr_->getWeight(FRcat, vbfVars_.mjj, photons_[0].Eta());
      }

      
      //fill control histograms
      for( auto c : chTags)
	fillControlHistos( boson,  jets,  plotwgts[0], c, mults, fakeACR, tightACR);

      //fill tree
      if(skimtree_) {
        evtWeight_ = plotwgts[0]*xsec_;
        training_ = useForTraining(); 	
        newTree_->Fill();
      }

      //experimental systs cycle: better not to do anything else after this...
      //final category selection is repeated ad nauseam with varied objects/weights and mva is re-evaluated several times
      if(ev_.isData) continue;
      if(!runSysts_) continue;
      std::vector<std::pair<float,float> > mvaWithWeights;
      selector_->setDebug(false);
      vbf::DiscriminatorInputs origVbfVars(vbfVars_);
      for(size_t is=0; is<expSysts_.size(); is++){
        
        //reset to the original values
        vbfVars_.assignValuesFrom(origVbfVars);
        
        //uncertainty
        TString sname=expSysts_[is];
        bool isUpVar(sname.Contains("up"));
        
        //base values and kinematics
        TString icat(baseCategory);
        float imva=vbfmva_;
        float iwgt=(ev_.g_nw>0 ? ev_.g_w[0] : 1.0);
        iwgt *= (normH_? normH_->GetBinContent(1) : 1.0);
        TLorentzVector iBoson(boson);
        std::vector<Jet> ijets(jets);
        bool reSelect(false);        
        
        if(sname=="puup")             iwgt *= puWgts[1]*trigSF.first*selSF.first*l1prefireProb.first;
        else if(sname=="pudn")        iwgt *= puWgts[2]*trigSF.first*selSF.first*l1prefireProb.first;
        else if(sname=="trigup")      iwgt *= puWgts[0]*(trigSF.first+trigSF.second)*selSF.first*l1prefireProb.first;
        else if(sname=="trigdn")      iwgt *= puWgts[0]*(trigSF.first-trigSF.second)*selSF.first*l1prefireProb.first;
        else if(sname=="selup")       iwgt *= puWgts[0]*trigSF.first*(selSF.first+selSF.second)*l1prefireProb.first;
        else if(sname=="seldn")       iwgt *= puWgts[0]*trigSF.first*(selSF.first-selSF.second)*l1prefireProb.first;
        else if(sname=="l1prefireup") iwgt *= puWgts[0]*trigSF.first*selSF.first*(l1prefireProb.first+l1prefireProb.second);
        else if(sname=="l1prefiredn") iwgt *= puWgts[0]*trigSF.first*selSF.first*(l1prefireProb.first-l1prefireProb.second);
        else                          iwgt = wgt;


        if(sname.Contains("aes") && chTag=="A")  {
          reSelect=true;
          iBoson *= (1+(isUpVar?1:-1)*bosonScaleUnc); 
        }
        //technically we should re-select the leptons but given we're looking to high pT Z's
        //assume effect is negligible and all that counts is the Z energy scale?
        if(sname.Contains("mes") && chTag=="MM") {
          reSelect=true;
          iBoson *= (1+(isUpVar?1:-1)*bosonScaleUnc); 
        }
        if(sname.Contains("ees") && chTag=="EE") {
          reSelect=true;
          iBoson *= (1+(isUpVar?1:-1)*bosonScaleUnc); 
        }
        if(sname.Contains("JEC") || sname.Contains("JER") )  {
          reSelect=true;
          int jecIdx=-1;
          if(sname.Contains("AbsoluteStat"))     jecIdx=0;
          if(sname.Contains("AbsoluteScale"))    jecIdx=1; 
          if(sname.Contains("AbsoluteMPFBias"))  jecIdx=2; 
          if(sname.Contains("Fragmentation"))    jecIdx=3; 
          if(sname.Contains("SinglePionECAL"))   jecIdx=4; 
          if(sname.Contains("SinglePionHCAL"))   jecIdx=5; 
          if(sname.Contains("FlavorPureGluon"))  jecIdx=6; 
          if(sname.Contains("FlavorPureQuark"))  jecIdx=7; 
          if(sname.Contains("FlavorPureCharm"))  jecIdx=8; 
          if(sname.Contains("FlavorPureBottom")) jecIdx=9; 
          if(sname.Contains("TimePtEta"))        jecIdx=10; 
          if(sname.Contains("RelativeJEREC1"))   jecIdx=11; 
          if(sname.Contains("RelativeJEREC2"))   jecIdx=12; 
          if(sname.Contains("RelativeJERHF"))    jecIdx=13; 
          if(sname.Contains("RelativePtBB"))     jecIdx=14; 
          if(sname.Contains("RelativePtEC1"))    jecIdx=15; 
          if(sname.Contains("RelativePtEC2"))    jecIdx=16; 
          if(sname.Contains("RelativePtHF"))     jecIdx=17; 
          if(sname.Contains("RelativeBal"))      jecIdx=18; 
          if(sname.Contains("RelativeFSR"))      jecIdx=19; 
          if(sname.Contains("RelativeStatFSR"))  jecIdx=20; 
          if(sname.Contains("RelativeStatEC"))   jecIdx=21; 
          if(sname.Contains("RelativeStatHF"))   jecIdx=22; 
          if(sname.Contains("PileUpDataMC"))     jecIdx=23; 
          if(sname.Contains("PileUpPtRef"))      jecIdx=24; 
          if(sname.Contains("PileUpPtBB"))       jecIdx=25; 
          if(sname.Contains("PileUpPtEC1"))      jecIdx=26; 
          if(sname.Contains("PileUpPtEC2"))      jecIdx=27; 
          if(sname.Contains("PileUpPtHF"))       jecIdx=28;
          
          //re-scale and re-select jets
          std::vector<Jet> newJets = selector_->getGoodJets(ev_,20.,4.7,leptons,photons_);
          ijets.clear();
          for(auto j : alljets) {

            int idx=j.getJetIndex();
            if(cleanEENoise_ && fabs(j.Eta())>2.7 && fabs(j.Eta())<3 && ev_.j_emf[idx]>0.55) continue;

            //shift jet energy
            float scaleVar(1.0);
            if(jecIdx<0) {
              scaleVar=isUpVar ? ev_.j_jerUp[idx] : ev_.j_jerDn[idx];
            } 
            else {
              int jflav(abs(ev_.j_flav[idx]));
              bool flavorMatches(true);
              if(jecIdx==6 && jflav!=21) flavorMatches=false; //FlavorPureGluon
              if(jecIdx==7 && jflav>=4)  flavorMatches=false; //FlavorPureQuark
              if(jecIdx==8 && jflav!=4)  flavorMatches=false; //FlavorPureCharm
              if(jecIdx==9 && jflav!=5)  flavorMatches=false; //FlavorPureGluon
              if(flavorMatches)
                scaleVar=isUpVar ? ev_.j_jecUp[jecIdx][idx] : ev_.j_jecDn[jecIdx][idx];
            }
            j*=scaleVar;
            if(j.Pt()<30.) continue;

            //pileup id
            int jid=ev_.j_id[idx];
            bool passPu((jid>>jetPuId_)&0x1);
            if(!passPu) continue;

            ijets.push_back(j);
          }
        }
        
        //re-select if needed
        if(reSelect) {
          
          if (ijets.size()<2) continue;

          vbfVars_.fillDiscriminatorVariables(iBoson,ijets,ev_);                   
          if(vbfVars_.leadj_pt    < leadJetPt_)    continue;
          if(vbfVars_.subleadj_pt < subLeadJetPt_) continue;

          //set the new category tag
          bool isLowVPt = (passLowVPtTrig                                                      
                           && iBoson.Pt()>lowVPtCut_
                           && iBoson.Pt()<=highVPtCut_
                           && isBosonTrigSafe);
          bool isHighVPt = (passHighVPtTrig
                            && iBoson.Pt()>highVPtCut_
                            && isBosonTrigSafe );
          bool isLowMJJ(vbfVars_.mjj>lowMJJCut_ && vbfVars_.mjj<=highMJJCut_);
          bool isHighMJJ(vbfVars_.mjj>highMJJCut_);

          if(isLowVPt && isHighMJJ 
             && fabs(iBoson.Rapidity())<lowVPtMaxRapCut_ && vbfVars_.mjj>highMJJCut_ && vbfVars_.detajj>lowVPtDetaJJCut_)       
            icat="LowVPtHighMJJ"+chTag;
          else if(isHighVPt && isLowMJJ)  icat="HighVPtLowMJJ"+chTag;
          else if(isHighVPt && isHighMJJ) icat="HighVPtHighMJJ"+chTag;
          
          //re-evaluate MVA         
          //FIXME: this needs probably to be changed for the next version of the training
          TString key(isLowVPt ? "BDT_VBF0HighMJJ" : "BDT_VBF0LowMJJ");
          imva = readers[key]->EvaluateMVA(key);
          if(mvaCDFinv[key]) imva=max(0.,mvaCDFinv[key]->Eval(imva));
        }
        
        if( !icat.Contains("MJJ") ) continue;

        //fill with new values/weights
        std::vector<double> eweights(1,iwgt);
        ht_->fill2D("vbfmva_exp",       imva,                 is,eweights,icat);
        ht_->fill2D("centraleta_exp",   vbfVars_.centraleta,  is,eweights,icat);
        ht_->fill2D("forwardeta_exp",   vbfVars_.forwardeta,  is,eweights,icat);
        ht_->fill2D("leadpt_exp",       vbfVars_.leadj_pt,    is,eweights,icat);
        ht_->fill2D("subleadpt_exp",    vbfVars_.subleadj_pt, is,eweights,icat);
        ht_->fill2D("detajj_exp",       vbfVars_.detajj,      is,eweights,icat);
        ht_->fill2D("dphijj_exp",       vbfVars_.dphijj,      is,eweights,icat);
        ht_->fill2D("mjj_exp", 	        vbfVars_.mjj,         is,eweights,icat);
        ht_->fill2D("vpt_exp", 	        iBoson.Pt(),          is,eweights,icat);
        ht_->fill2D("evcount_exp",      0,                    is,eweights,icat);
        if(imva>0.9)
          ht_->fill2D("evcount_exp",    1,                    is,eweights,icat);
      }
      selector_->setDebug(debug_);
    }
  
  //close input file
  f_->Close();

  //save histos and new tree to file  
  saveHistos();
  if(skimtree_){
    fMVATree_->cd();
    newTree_->Write();
    fMVATree_->Close();
  }
}

void VBFVectorBoson::saveHistos(){
  fOut_->cd();
  for (auto& it : ht_->getPlots())  { 
    if(it.second->GetEntries()==0) continue;
    it.second->SetDirectory(fOut_); 
    it.second->Write(); 
  }
  for (auto& it : ht_->get2dPlots())  { 
    if(it.second->GetEntries()==0) continue;
    it.second->SetDirectory(fOut_); 
    it.second->Write(); 
  }
  fOut_->Close();		
}

//
void VBFVectorBoson::readTree()
{
  f_          = TFile::Open(filename_);
  vetoPromptPhotons_ = filename_.Contains("_QCDEM_") || filename_.Contains("_TTJets");
  triggerList_ = (TH1F *)f_->Get("analysis/triggerList");
  weightSysts_ = getWeightSysts(f_);
  t_           = (TTree*)f_->Get("analysis/data");
  attachToMiniEventTree(t_,ev_,true);
  nentries_   = t_->GetEntriesFast();
  if (debug_) nentries_ = 10000; //restrict number of entries for testing
  t_->GetEntry(0);
}

//
void VBFVectorBoson::prepareOutput()
{
  baseName_=gSystem->BaseName(outname_); 
  TString dirName=gSystem->DirName(outname_);
  
  if(skimtree_) {
    fMVATree_ = TFile::Open(dirName+"/MVATree_"+baseName_,"RECREATE");
    newTree_ = t_->CloneTree(0);
  }
  fOut_=TFile::Open(dirName+"/"+baseName_,"RECREATE");
  fOut_->cd();
}

//
void VBFVectorBoson::bookHistograms() {
  ht_ = new HistTool(0);
  ht_->addHist("puwgtctr",      new TH1F("puwgtctr",         ";Weight sums;Events",                2,0,2));  
  ht_->addHist("qscale",        new TH1F("qscale",           ";Q^{2} scale;Events",                100,0,2000));  
  ht_->addHist("nvtx",          new TH1F("nvtx",             ";Vertex multiplicity;Events",        100,-0.5,99.5));  
  ht_->addHist("vpt", 	       new TH1F("vectorbosonPt",    ";Boson p_{T}[GeV];Events",           25,50,550));  
  ht_->addHist("vy", 	       new TH1F("vectorbosony",     ";Boson rapidity;Events",             25,-3,3));  
  ht_->addHist("mindrl",        new TH1F("mindrl",           ";min #Delta R(boson,lepton);Events", 25,0,6));  
  ht_->addHist("mindrj",        new TH1F("mindrj",           ";min #Delta R(boson,jet);Events",    25,0,6));  
  ht_->addHist("sihih", 	       new TH1F("sihih",            ";#sigma(i#eta,i#eta);Events",        50,0,0.1));  
  ht_->addHist("hoe", 	       new TH1F("hoe",              ";h/e;Events",                        25,0,0.05));  
  ht_->addHist("r9", 	       new TH1F("r9",               ";r9;Events",                         25,0.9,1.0));  
  ht_->addHist("chiso", 	       new TH1F("chiso",            ";Charged isolation [GeV];Events",    50,0,10));  

  ht_->addHist("njets",         new TH1F("njets",            ";Jet multiplicity;Events",            10,-0.5,9.5));  
  ht_->addHist("leadpt",        new TH1F("leadpt",           ";Leading jet p_{T} [GeV];Events",     25,0,500));  
  ht_->addHist("subleadpt",     new TH1F("subleadpt"   ,     ";Sub-leading jet p_{T} [GeV];Events", 25,0,500));  
  ht_->addHist("centraleta",    new TH1F("centraleta",       ";Most central jet |#eta|;Events",     25,0,5));  
  ht_->addHist("forwardeta",    new TH1F("forwardeta",       ";Most forward jet |#eta|;Events",     25,0,5));  
  ht_->addHist("leadeta",       new TH1F("leadeta",          ";Leading jet |#eta|;Events",          25,0,5));  
  ht_->addHist("subleadeta",    new TH1F("subleadeta",       ";Sub-leading jet |#eta|;Events",      25,0,5));  
  ht_->addHist("leadpumva",     new TH1F("leadpumva",        ";Pileup MVA;Events",                  25,-1,1));  
  ht_->addHist("subleadpumva",  new TH1F("subleadpumva"   ,  ";Pileup MVA;Events",                  25,-1,1));  
  ht_->addHist("etaphi",        new TH2F("etaphi",           ";Most central jet |#eta|V|#phi|;Events",30,-4.7,4.7,25,-TMath::Pi(),TMath::Pi()));
  ht_->addHist("jet_raw_pt",    new TH1F("jet_raw_pt",       ";raw PT of jets;Jets",                50,0,200));
  ht_->addHist("jet_raw_empt",  new TH1F("jet_raw_empt",     ";raw e.m. PT of jets;Jets",           50,0,200));
  ht_->addHist("jet_emf",       new TH1F("jet_emf",          ";Electromagnetic fraction;Jets",      50,0,1));
  ht_->addHist("jet_qg",        new TH1F("jet_qg",           ";q/g discriminator;Jets",             50,0,1));
  ht_->addHist("jet_c2_00",     new TH1F("jet_c2_00",        ";Energy correlation C_{2}^{(0)};Jets",50,0,1));  
  ht_->addHist("jet_c2_02",     new TH1F("jet_c2_02",        ";Energy correlation C_{2}^{(2)};Jets",50,0,1));  
  ht_->addHist("jet_c2_05",     new TH1F("jet_c2_05",        ";Energy correlation C_{2}^{(5)};Jets",50,0,1));  
  ht_->addHist("jet_zg",        new TH1F("jet_zg",           ";Jet z_{g};Jets",                     50,0,1));  
  ht_->addHist("jet_gaptd",     new TH1F("jet_gaptd",        ";Jet p_{T}^{D};Jets",                 50,0,1));  
  ht_->addHist("jet_gawidth",   new TH1F("jet_gawidth",      ";Jet width;Jets",                     50,0,1));

  ht_->addHist("mjj", 	       new TH1F("mjj",              ";Dijet invariant mass [GeV];Events", 40,0,4000));  
  ht_->addHist("detajj",        new TH1F("detajj" ,          ";#Delta#eta(J,J);Events",            20,0,8));  
  ht_->addHist("dphijj",        new TH1F("dphijj" ,          ";#Delta#phi(J,J) [rad];Events",      20,-3.15,3.15));  
  ht_->addHist("dijetpt",       new TH1F("dijetpt",          ";Dijet p_{T} [GeV];Events",          20,0,1000));  
  ht_->addHist("ht",            new TH1F("ht",               ";H_{T} [GeV];Events",                20,0,4000));  
  ht_->addHist("mht",           new TH1F("mht",              ";Missing H_{T} [GeV];Events",        20,0,500));  

  ht_->addHist("drj1b",         new TH1F("drj1b",            ";#DeltaR(j_{1},boson);Events",       25,0,8));  
  ht_->addHist("drj2b",         new TH1F("drj2b"   ,         ";#DeltaR(j_{2},boson);Events",       25,0,8));  
  ht_->addHist("vystar",        new TH1F("vectorbosonystar", ";y-(1/2)(y_{j1}+y_{j2});Events",     25,-5,5));  
  ht_->addHist("balance",       new TH1F("balance",          ";System p_{T} balance [GeV];Events", 20,0,300));  
  ht_->addHist("relbpt",        new TH1F("relbpt",           ";#Sigma p_{T}(j)/Boson p_{T};Events",25,0,2));  
  ht_->addHist("dphibjj",       new TH1F("dphibjj",          ";#Delta#phi(JJ,boson) [rad];Events", 20,-3.15,3.15));  
  ht_->addHist("sphericity",    new TH1F("sphericity",       ";Sphericity;Events",                 20,0,1.0));  
  ht_->addHist("aplanarity",    new TH1F("aplanarity",       ";Aplanarity;Events",                 20,0,1.0));  
  ht_->addHist("C",             new TH1F("C",                ";C;Events",                          20,0,1.0));  
  ht_->addHist("D",             new TH1F("D",                ";D;Events",                          20,0,1.0));  
  ht_->addHist("isotropy",      new TH1F("isotropy",         ";Isotropy;Events",                   20,0,1.0));  
  ht_->addHist("circularity",   new TH1F("circularity",      ";Circularity;;Events",               20,0,1.0));
  //additional variables from https://link.springer.com/content/pdf/10.1140/epjc/s10052-017-5315-6.pdf
  ht_->addHist("jjetas",        new TH1F("jjetas",           ";#eta_{j1}#eta_{j2};Events",50,-25,15));  
  ht_->addHist("centjy",	       new TH1F("centjy",           ";Central jet rapidity;Jets",25,0,3));  
  ht_->addHist("ncentj",        new TH1F("ncentjj",          ";Number of central jets;Events",10,-0.5,9.5));  
  ht_->addHist("dphivj0",       new TH1F("dphivj0",          ";#Delta#phi(V,j0);Jets",20,0,TMath::Pi()));  
  ht_->addHist("dphivj1",       new TH1F("dphivj1",          ";#Delta#phi(V,j1);Jets",20,0,TMath::Pi()));  
  ht_->addHist("dphivj2",       new TH1F("dphivj2",          ";#Delta#phi(V,j2);Jets",20,0,TMath::Pi()));  
  ht_->addHist("dphivj3",       new TH1F("dphivj3",          ";#Delta#phi(V,j3);Jets",20,0,TMath::Pi()));
  //final analyses distributions
  ht_->addHist("evcount",         new TH1F("evcount",        ";Pass;Events",2,0,2));  
  ht_->getPlots()["evcount"]->GetXaxis()->SetBinLabel(1,"Inclusive");
  ht_->getPlots()["evcount"]->GetXaxis()->SetBinLabel(2,"MVA>0.9");
  ht_->addHist("vbfmva",          new TH1F("vbfmva",         ";VBF MVA;Events",50,0,1));  

  TString expSystNames[]={"puup","pudn","trigup","trigdn","selup","seldn","l1prefireup","l1prefiredn",
                          "aesup","aesdn",
                          "mesup","mesdn",
                          "JERup","JERdn",
                          "AbsoluteStatJECup","AbsoluteScaleJECup","AbsoluteMPFBiasJECup","FragmentationJECup","SinglePionECALJECup","SinglePionHCALJECup","FlavorPureGluonJECup","FlavorPureQuarkJECup","FlavorPureCharmJECup","FlavorPureBottomJECup","TimePtEtaJECup","RelativeJEREC1JECup","RelativeJEREC2JECup","RelativeJERHFJECup","RelativePtBBJECup","RelativePtEC1JECup","RelativePtEC2JECup","RelativePtHFJECup","RelativeBalJECup","RelativeFSRJECup","RelativeStatFSRJECup","RelativeStatECJECup","RelativeStatHFJECup","PileUpDataMCJECup","PileUpPtRefJECup","PileUpPtBBJECup","PileUpPtEC1JECup","PileUpPtEC2JECup","PileUpPtHFJECup",
                          "AbsoluteStatJECdn","AbsoluteScaleJECdn","AbsoluteMPFBiasJECdn","FragmentationJECdn","SinglePionECALJECdn","SinglePionHCALJECdn","FlavorPureGluonJECdn","FlavorPureQuarkJECdn","FlavorPureCharmJECdn","FlavorPureBottomJECdn","TimePtEtaJECdn","RelativeJEREC1JECdn","RelativeJEREC2JECdn","RelativeJERHFJECdn","RelativePtBBJECdn","RelativePtEC1JECdn","RelativePtEC2JECdn","RelativePtHFJECdn","RelativeBalJECdn","RelativeFSRJECdn","RelativeStatFSRJECdn","RelativeStatECJECdn","RelativeStatHFJECdn","PileUpDataMCJECdn","PileUpPtRefJECdn","PileUpPtBBJECdn","PileUpPtEC1JECdn","PileUpPtEC2JECdn","PileUpPtHFJECdn"};
  
  //instantiate 2D histograms for most relevant variables to trace with systs
  TString hoi[]={"vbfmva","evcount","mjj","detajj","dphijj","leadpt","subleadpt","forwardeta","centraleta","vpt"};
  size_t nexpSysts=sizeof(expSystNames)/sizeof(TString);
  expSysts_=std::vector<TString>(expSystNames,expSystNames+nexpSysts);  
  for(size_t ih=0; ih<sizeof(hoi)/sizeof(TString); ih++)
    {
      TH1 *histo=ht_->getPlots()[hoi[ih]];
      
      //experimental systs
      ht_->addHist(hoi[ih]+"_exp",      
                  new TH2F(hoi[ih]+"_exp", 
                           Form(";%s;Experimental systematic variation;Events",histo->GetName()),
                           histo->GetNbinsX(),histo->GetXaxis()->GetXmin(),histo->GetXaxis()->GetXmax(),
                           nexpSysts,0,nexpSysts));
      for(size_t is=0; is<nexpSysts; is++)
        ht_->get2dPlots()[hoi[ih]+"_exp"]->GetYaxis()->SetBinLabel(is+1,expSystNames[is]);
      
      //theory systs
      size_t nthSysts(weightSysts_.size()); 
      if(nthSysts>0){
        ht_->addHist(hoi[ih]+"_th",      
                    new TH2F(hoi[ih]+"_th", 
                             Form(";%s;Theory systematic variation;Events",histo->GetName()),
                             histo->GetNbinsX(),histo->GetXaxis()->GetXmin(),histo->GetXaxis()->GetXmax(),
                             nthSysts,0,nthSysts));
        for(size_t is=0; is<nthSysts; is++)
          ht_->get2dPlots()[hoi[ih]+"_th"]->GetYaxis()->SetBinLabel(is+1,weightSysts_[is].first);
      }
    }
      
  //Photons in CR
  ht_->addHist("allsihih",         new TH1F("allsihih",         ";All #sigma(i#eta,i#eta);Photons",100,0,0.05));
  ht_->addHist("relaxedTightsihih",new TH1F("relaxedTightsihih",      ";Relaxed tight #sigma(i#eta,i#eta);Photons",100,0,0.05));
  ht_->addHist("fakesihih",        new TH1F("fakesihih",        ";Fake #sigma(i#eta,i#eta);Photons",100,0,0.05));
  ht_->addHist("tightsihih",       new TH1F("tightsihih",       ";Tight #sigma(i#eta,i#eta);Photons",100,0,0.05));
  ht_->addHist("fakechiso",        new TH1F("fakechiso",        ";Fake ch. isolation [GeV];Photons",50,0,10));  
  ht_->addHist("tightchiso",       new TH1F("tightchiso",       ";Tight ch. isolation [GeV];Photons",50,0,10));
  ht_->addHist("fakeneutiso",      new TH1F("fakeneutiso",      ";Fake neut. isolation [GeV];Photons",50,0,10));  
  ht_->addHist("tightneutiso",     new TH1F("tightneutiso",     ";Tight neut. isolation [GeV];Photons",50,0,10));
  ht_->addHist("fakeaiso",         new TH1F("fakeaiso",         ";Fake #gamma isolation [GeV];Photons",50,0,10));  
  ht_->addHist("tightaiso",        new TH1F("tightaiso",        ";Tight #gamma isolation [GeV];Photons",50,0,10));  
  ht_->addHist("nloose",           new TH1F("nloose",           ";Number of loose #gamma; Events",20,-0.5,19.5));
  ht_->addHist("ntight",           new TH1F("ntight",           ";Number of tight #gamma; Events",20,-0.5,19.5));
  ht_->addHist("nloosefake",       new TH1F("nloosefake",       ";Number of fake loose #gamma; Events",20,-0.5,19.5));
  ht_->addHist("ntightfake",       new TH1F("ntightfake",       ";Number of fake tight #gamma; Events",20,-0.5,19.5));
  ht_->addHist("nlooseprompt",     new TH1F("nlooseprompt",     ";Number of prompt loose #gamma; Events",20,-0.5,19.5));
  ht_->addHist("ntightprompt",     new TH1F("ntightprompt",     ";Number of prompt tight #gamma; Events",20,-0.5,19.5));
  //2D's for Mjj-binned FR
  double bins[]={0,500,1000,2000,4000};
  ht_->addHist("relaxedTightMjjEB",  new TH2F("relaxedTightMjjEB",";Relaxed tight #sigma_{i#etai#eta}; m_{jj} (GeV)",100,0,0.05,4,bins)); //80,0,4000
  ht_->addHist("tightMjjEB",         new TH2F("tightMjjEB",";Tight #sigma_{i#etai#eta}; m_{jj} (GeV)",100,0,0.05,4,bins));
  ht_->addHist("looseMjjEB",         new TH2F("looseMjjEB",";Loose #sigma_{i#etai#eta}; m_{jj} (GeV)",100,0,0.05,4,bins));
  ht_->addHist("allMjjEB",           new TH2F("allMjjEB",";All #sigma_{i#etai#eta}; m_{jj} (GeV)",100,0,0.05,4,bins));
  ht_->addHist("tmpQCDMjjEB",        new TH2F("tmpQCDMjjEB",";All #sigma_{i#etai#eta}; m_{jj} (GeV)",100,0,0.05,4,bins));
  ht_->addHist("relaxedTightMjjEE",  new TH2F("relaxedTightMjjEE",";Relaxed tight #sigma_{i#etai#eta}; m_{jj} (GeV)",100,0,0.05,4,bins)); //80,0,4000
  ht_->addHist("tightMjjEE",         new TH2F("tightMjjEE",";Tight #sigma_{i#etai#eta}; m_{jj} (GeV)",100,0,0.05,4,bins));
  ht_->addHist("looseMjjEE",         new TH2F("looseMjjEE",";Loose #sigma_{i#etai#eta}; m_{jj} (GeV)",100,0,0.05,4,bins));
  ht_->addHist("allMjjEE",           new TH2F("allMjjEE",";All #sigma_{i#etai#eta}; m_{jj} (GeV)",100,0,0.05,4,bins));
  ht_->addHist("tmpQCDMjjEE",        new TH2F("tmpQCDMjjEE",";All #sigma_{i#etai#eta}; m_{jj} (GeV)",100,0,0.05,4,bins));
}


void VBFVectorBoson::loadCorrections(){
  fr_ = new FakeRateTool(era_, "fakeRatios.root");
  lumi_ = new LumiTools(era_,genPU_);

  std::map<TString,TString> cfgMap;
  cfgMap["g_id"]     = "Tight";
  cfgMap["m_id"]     = "TightID";
  cfgMap["m_iso"]    = "TightRelIso";
  cfgMap["m_id4iso"] = "TightIDandIPCut";
  cfgMap["e_id"]     = "MVA80";
  gammaEffWR_  = new EfficiencyScaleFactorsWrapper(filename_.Contains("Data13TeV"),era_,cfgMap);
  l1PrefireWR_ = new L1PrefireEfficiencyWrapper(filename_.Contains("Data13TeV"),era_);
}


//
void VBFVectorBoson::addMVAvars(){
  newTree_->Branch("centralEta",       &vbfVars_.centraleta);
  newTree_->Branch("forwardeta",       &vbfVars_.forwardeta);
  newTree_->Branch("leadj_pt",         &vbfVars_.leadj_pt);
  newTree_->Branch("subleadj_pt",      &vbfVars_.subleadj_pt);
  newTree_->Branch("mjj",              &vbfVars_.mjj);
  newTree_->Branch("detajj",           &vbfVars_.detajj);
  newTree_->Branch("jjpt",             &vbfVars_.jjpt);
  newTree_->Branch("dphijj",           &vbfVars_.dphijj);
  newTree_->Branch("ystar",            &vbfVars_.ystar);
  newTree_->Branch("relbpt",           &vbfVars_.relbpt);
  newTree_->Branch("dphibjj",          &vbfVars_.dphibjj);
  newTree_->Branch("balance",          &vbfVars_.balance);
  newTree_->Branch("leadj_gawidth",    &vbfVars_.leadj_gawidth);
  newTree_->Branch("subleadj_gawidth", &vbfVars_.subleadj_gawidth);
  newTree_->Branch("subleadj_c2_02",   &vbfVars_.subleadj_c2_02);
  newTree_->Branch("jjetas",           &vbfVars_.jjetas);
  newTree_->Branch("centjy",           &vbfVars_.centjy);
  newTree_->Branch("ncentjj",          &vbfVars_.ncentj);
  newTree_->Branch("dphivj0",          &vbfVars_.dphivj0);
  newTree_->Branch("dphivj1",          &vbfVars_.dphivj1);
  newTree_->Branch("dphivj2",          &vbfVars_.dphivcentj[0]);
  newTree_->Branch("dphivj3",          &vbfVars_.dphivcentj[0]);
  newTree_->Branch("mht",              &vbfVars_.mht);
  newTree_->Branch("ht",               &vbfVars_.scalarht);
  newTree_->Branch("isotropy",         &vbfVars_.isotropy);
  newTree_->Branch("circularity",      &vbfVars_.circularity);
  newTree_->Branch("sphericity",       &vbfVars_.sphericity);
  newTree_->Branch("aplanarity",       &vbfVars_.aplanarity);
  newTree_->Branch("C",                &vbfVars_.C);
  newTree_->Branch("D",                &vbfVars_.D);
  newTree_->Branch("evtWeight",        &evtWeight_);
  newTree_->Branch("training",         &training_);
  newTree_->Branch("category",         &category_, "EE:MM:A:LowVPt:HighVPt:LowMJJ:HighMJJ");
}

//
void VBFVectorBoson::fillControlHistos(TLorentzVector boson, std::vector<Jet> jets, float wgt, TString c, std::map<TString, int> mults, std::vector<Particle> fakeACR, std::vector<Particle> tightACR){

  //plot weight
  std::vector<double> cplotwgts(1,wgt);

  //mandatory...
  ht_->fill("nvtx",   ev_.nvtx,          cplotwgts,c);        
  
  //boson histos
  ht_->fill("vpt",    boson.Pt(),       cplotwgts,c);
  ht_->fill("vy",     boson.Rapidity(), cplotwgts,c);   
  ht_->fill("mindrl", mindrl_,          cplotwgts,c);   
  ht_->fill("mindrj", mindrj_,          cplotwgts,c);   
  
  //photon specific plots
  if( !(c.Contains("EE") || c.Contains("MM")) ) {
      ht_->fill("sihih",  sihih_,           cplotwgts,c);   
      ht_->fill("r9",     r9_,              cplotwgts,c);   
      ht_->fill("hoe",    hoe_,             cplotwgts,c);   
      ht_->fill("chiso",  chiso_,           cplotwgts,c);   

      for(auto a : photons_) {
        int idx = a.originalReference();
        ht_->fill("allsihih",   ev_.gamma_sieie[idx]             ,cplotwgts,c);
        if (fabs(ev_.gamma_eta[idx]) < 1.442)
          ht_->fill2D("allMjjEB"  ,   ev_.gamma_sieie[idx], vbfVars_.mjj        ,cplotwgts,c); 
        else
          ht_->fill2D("allMjjEE"  ,   ev_.gamma_sieie[idx], vbfVars_.mjj        ,cplotwgts,c); 
      }
      for(auto a : relaxedTightPhotons_) {
        int idx = a.originalReference();
        ht_->fill("relaxedTightsihih",   ev_.gamma_sieie[idx]             ,cplotwgts,c);
        if (fabs(ev_.gamma_eta[idx]) < 1.442)
          ht_->fill2D("relaxedTightMjjEB"  ,   ev_.gamma_sieie[idx], vbfVars_.mjj        ,cplotwgts,c); 
        else
          ht_->fill2D("relaxedTightMjjEE"  ,   ev_.gamma_sieie[idx], vbfVars_.mjj        ,cplotwgts,c); 
      }
      for(auto a : tmpPhotons_) {
        int idx = a.originalReference();
        if (fabs(ev_.gamma_eta[idx]) < 1.442)
          ht_->fill2D("tmpQCDMjjEB"  ,   ev_.gamma_sieie[idx], vbfVars_.mjj        ,cplotwgts,c); 
        else
          ht_->fill2D("tmpQCDMjjEE"  ,   ev_.gamma_sieie[idx], vbfVars_.mjj        ,cplotwgts,c); 
      }
      //bosons in CR and fakes
      for(auto a : fakeACR) {
        int idx = a.originalReference();
        ht_->fill("fakesihih",   ev_.gamma_sieie[idx]             ,cplotwgts,c);
        ht_->fill("fakechiso",   ev_.gamma_chargedHadronIso[idx]  ,cplotwgts,c);
        ht_->fill("fakeneutiso", ev_.gamma_neutralHadronIso[idx]  ,cplotwgts,c);
        ht_->fill("fakeaiso",    ev_.gamma_photonIso[idx]         ,cplotwgts,c);
        if (fabs(ev_.gamma_eta[idx]) < 1.442)
          ht_->fill2D("looseMjjEB"  ,  ev_.gamma_sieie[idx], vbfVars_.mjj        ,cplotwgts,c);
        else
          ht_->fill2D("looseMjjEE"  ,  ev_.gamma_sieie[idx], vbfVars_.mjj        ,cplotwgts,c);
      }
      for(auto a : tightACR) { 
        int idx = a.originalReference();
        ht_->fill("tightsihih",   ev_.gamma_sieie[idx]             ,cplotwgts,c);
        ht_->fill("tightchiso",   ev_.gamma_chargedHadronIso[idx]  ,cplotwgts,c);
        ht_->fill("tightneutiso", ev_.gamma_neutralHadronIso[idx]  ,cplotwgts,c);
        ht_->fill("tightaiso",    ev_.gamma_photonIso[idx]         ,cplotwgts,c);
        if (fabs(ev_.gamma_eta[idx]) < 1.442)
          ht_->fill2D("tightMjjEB"  ,   ev_.gamma_sieie[idx], vbfVars_.mjj        ,cplotwgts,c);
        else
          ht_->fill2D("tightMjjEE"  ,   ev_.gamma_sieie[idx], vbfVars_.mjj        ,cplotwgts,c);
      }
      ht_->fill("nloose",        fakeACR.size(),       cplotwgts,c);
      ht_->fill("ntight",        tightACR.size(),      cplotwgts,c);
      ht_->fill("nloosefake",    mults["loosefake"],   cplotwgts,c);
      ht_->fill("ntightfake",    mults["tightfake"],   cplotwgts,c);
      ht_->fill("nlooseprompt",  mults["looseprompt"], cplotwgts,c);
      ht_->fill("ntightprompt",  mults["tightprompt"], cplotwgts,c);
  }

  //jet histos
  ht_->fill("njets",        jets.size(), cplotwgts,c);
  ht_->fill("ht",           vbfVars_.scalarht,    cplotwgts,c);
  ht_->fill("mht",          vbfVars_.mht,         cplotwgts,c);
  for(size_t ij=0; ij<min(size_t(2),jets.size());ij++) {
    TString jtype(ij==0?"lead":"sublead");
    ht_->fill(jtype+"pt",       jets[ij].Pt(),        cplotwgts,c);          
    ht_->fill(jtype+"eta",      fabs(jets[ij].Eta()), cplotwgts,c);          
    ht_->fill(jtype+"pumva",    jets[ij].getPUMVA(),  cplotwgts,c);
    ht_->fill(Form("drj%db",(int)ij+1),   jets[ij].DeltaR(boson),  cplotwgts,c);
    ht_->fill("jet_c2_00", 	ev_.j_c2_00[jets[ij].getJetIndex()]	  ,  cplotwgts,c);
    ht_->fill("jet_c2_02", 	ev_.j_c2_02[jets[ij].getJetIndex()]	  ,  cplotwgts,c);
    ht_->fill("jet_c2_05",	ev_.j_c2_05[jets[ij].getJetIndex()]	  ,  cplotwgts,c);
    ht_->fill("jet_zg", 		ev_.j_zg[jets[ij].getJetIndex()]	          ,  cplotwgts,c);
    ht_->fill("jet_gaptd", 	ev_.j_gaptd[jets[ij].getJetIndex()]	  ,  cplotwgts,c);
    ht_->fill("jet_gawidth",     ev_.j_gawidth[jets[ij].getJetIndex()]	  ,  cplotwgts,c);
    float j_rawpt   = jets[ij].pt()/ev_.j_rawsf[jets[ij].getJetIndex()];
    float j_emf     = ev_.j_emf[jets[ij].getJetIndex()];
    float j_rawempt = j_rawpt*j_emf;
    ht_->fill("jet_emf"     , j_emf        ,cplotwgts,c);
    ht_->fill("jet_raw_pt" ,j_rawpt        , cplotwgts,c );
    ht_->fill("jet_raw_empt" ,j_rawempt        , cplotwgts,c );
    ht_->fill("jet_qg"   ,ev_.j_qg[jets[ij].getJetIndex()]    , cplotwgts,c);
    ht_->fill2D("etaphi",  jets[ij].Eta(),jets[ij].Phi() ,   cplotwgts,c);    
  }
  ht_->fill("centraleta",   vbfVars_.centraleta,  cplotwgts,c);
  ht_->fill("forwardeta",   vbfVars_.forwardeta,  cplotwgts,c);
  ht_->fill("dijetpt",      vbfVars_.jjpt,        cplotwgts,c);
  ht_->fill("detajj",       vbfVars_.detajj,      cplotwgts,c);
  ht_->fill("dphijj",       vbfVars_.dphijj,      cplotwgts,c);
  ht_->fill("mjj", 	   vbfVars_.mjj,         cplotwgts,c);   
  ht_->fill("jjetas",       vbfVars_.jjetas,   cplotwgts,c);
  ht_->fill("dphivj0",      vbfVars_.dphivj0 ,  cplotwgts,c);
  ht_->fill("dphivj1",      vbfVars_.dphivj1 ,  cplotwgts,c);
  
  //central jet activity
  ht_->fill("ncentj", vbfVars_.ncentj, cplotwgts, c);
  if(vbfVars_.ncentj>0){
    ht_->fill("dphivj2", vbfVars_.dphivcentj[0] ,  cplotwgts,c);
    if(vbfVars_.ncentj>1) 
      ht_->fill("dphivj3", vbfVars_.dphivcentj[1] ,  cplotwgts,c);    
  }
	
  //visible system histos
  ht_->fill("relbpt",       vbfVars_.relbpt,      cplotwgts,c);
  ht_->fill("dphibjj",      vbfVars_.dphibjj,     cplotwgts,c);
  ht_->fill("vystar",       vbfVars_.ystar,              cplotwgts,c);        
  ht_->fill("balance",      vbfVars_.balance,            cplotwgts,c);
  ht_->fill("isotropy",     vbfVars_.isotropy,     cplotwgts,c);
  ht_->fill("circularity",  vbfVars_.circularity,  cplotwgts,c);
  ht_->fill("sphericity",   vbfVars_.sphericity, cplotwgts,c);
  ht_->fill("aplanarity",   vbfVars_.aplanarity, cplotwgts,c);
  ht_->fill("C",            vbfVars_.C,          cplotwgts,c);
  ht_->fill("D",            vbfVars_.D,          cplotwgts,c);

  //final analysis histograms
  ht_->fill("evcount",  0, cplotwgts, c);
  if(vbfmva_>=0) {
    ht_->fill("vbfmva", vbfmva_, cplotwgts,c);
    if(vbfmva_>0.9)
      ht_->fill("evcount",  1, cplotwgts, c);  
  }

  //theory uncertainties are filled only for MC
  if(runSysts_ && !ev_.isData && ev_.g_w[0]!=0 && normH_ && normH_->GetBinContent(1)>0 && c.Contains("MJJ")) {
    
    //replicas for theory systs
    for(size_t is=0; is<weightSysts_.size(); is++){
      std::vector<double> sweights(1,cplotwgts[0]);
      size_t idx=weightSysts_[is].second;
      sweights[0] *= (ev_.g_w[idx]/ev_.g_w[0])*(normH_->GetBinContent(idx+1)/normH_->GetBinContent(1));
      ht_->fill2D("vbfmva_th",       vbfmva_,              is,sweights,c);
      ht_->fill2D("centraleta_th",   vbfVars_.centraleta,  is,sweights,c);
      ht_->fill2D("forwardeta_th",   vbfVars_.forwardeta,  is,sweights,c);
      ht_->fill2D("leadpt_th",       vbfVars_.leadj_pt,    is,sweights,c);
      ht_->fill2D("subleadpt_th",    vbfVars_.subleadj_pt, is,sweights,c);
      ht_->fill2D("detajj_th",       vbfVars_.detajj,      is,sweights,c);
      ht_->fill2D("dphijj_th",       vbfVars_.dphijj,      is,sweights,c);
      ht_->fill2D("mjj_th", 	     vbfVars_.mjj,         is,sweights,c);
      ht_->fill2D("vpt_th", 	     boson.Pt(),           is,sweights,c);
      ht_->fill2D("evcount_th",      0,                    is,sweights,c);
      if(vbfmva_>0.9)
        ht_->fill2D("evcount_th",    1,                    is,sweights,c);
    }  
  }
}
