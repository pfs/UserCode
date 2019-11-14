
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


  TString method("BDT_VBF0LowVPtHighMJJ");
  TString weightFile("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/VBF_weights/LowVPtHighMJJ_BDT_VBF0LowVPtHighMJJ.weights.xml");
  gSystem->ExpandPathName(weightFile);
  readers[method]=new TMVA::Reader( "!Color:!Silent" );
  readers[method]->AddVariable("ystar",                                &vbfVars_.ystar);
  readers[method]->AddVariable("mjj",                                  &vbfVars_.mjj);
  readers[method]->AddVariable("dphijj",                               &vbfVars_.dphijj);
  readers[method]->AddVariable("j_qg[0]",                              &vbfVars_.leadj_qg);
  readers[method]->AddVariable("dphibjj",                              &vbfVars_.dphibjj);
  readers[method]->AddVariable("dphivj0",                              &vbfVars_.dphivj0);
  readers[method]->AddVariable("j_qg[1]",                              &vbfVars_.subleadj_qg);
  readers[method]->AddVariable("C",                                    &vbfVars_.C);
  readers[method]->AddVariable("jjpt",                                 &vbfVars_.jjpt);
  readers[method]->AddVariable("j_pt[0]+j_pt[1]+gamma_pt[0]",          &vbfVars_.syspt);
  readers[method]->BookMVA(method,weightFile);

  method="BDT_VBF0HighVPt";
  weightFile="${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/VBF_weights/HighVPt_BDT_VBF0HighVPt.weights.xml";
  gSystem->ExpandPathName(weightFile);
  readers[method]=new TMVA::Reader( "!Color:!Silent" );
  readers[method]->AddVariable("subleadj_pt",   &vbfVars_.subleadj_pt);
  readers[method]->AddVariable("mjj",           &vbfVars_.mjj);
  readers[method]->AddVariable("jjpt",          &vbfVars_.jjpt);
  readers[method]->AddVariable("detajj",        &vbfVars_.detajj);
  readers[method]->AddVariable("dphijj",        &vbfVars_.dphijj);
  readers[method]->AddVariable("ystar",         &vbfVars_.ystar);
  readers[method]->AddVariable("dphibjj",       &vbfVars_.dphibjj);
  readers[method]->AddVariable("j_qg[0]",       &vbfVars_.leadj_qg);
  readers[method]->AddVariable("j_qg[1]",       &vbfVars_.subleadj_qg);
  readers[method]->AddVariable("dphivj0",       &vbfVars_.dphivj0);
  readers[method]->AddVariable("dphivj1",       &vbfVars_.dphivj1);
  readers[method]->AddVariable("aplanarity",    &vbfVars_.aplanarity);
  readers[method]->AddVariable("C",             &vbfVars_.C);
  readers[method]->BookMVA(method,weightFile);

  method="BDT_VBF0LowVPt";
  weightFile="${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/VBF_weights/LowVPt_BDT_VBF0LowVPt.weights.xml";
  gSystem->ExpandPathName(weightFile);
  readers[method]=new TMVA::Reader( "!Color:!Silent" );
  readers[method]->AddVariable("circularity",   &vbfVars_.circularity);
  readers[method]->AddVariable("mjj",           &vbfVars_.mjj);
  readers[method]->AddVariable("jjpt",          &vbfVars_.jjpt);
  readers[method]->AddVariable("dphijj",        &vbfVars_.dphijj);
  readers[method]->AddVariable("ystar",         &vbfVars_.ystar);
  readers[method]->AddVariable("relbpt",        &vbfVars_.relbpt);
  readers[method]->AddVariable("dphibjj",       &vbfVars_.dphibjj);
  readers[method]->AddVariable("j_qg[0]",       &vbfVars_.leadj_qg);
  readers[method]->AddVariable("j_qg[1]",       &vbfVars_.subleadj_qg);
  readers[method]->AddVariable("dphivj0",       &vbfVars_.dphivj0);
  readers[method]->AddVariable("D",             &vbfVars_.D);
  readers[method]->BookMVA(method,weightFile);


  method="BDT_VBF0HighVPtLowMJJ";
  weightFile="${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/VBF_weights/HighVPtLowMJJ_BDT_VBF0HighVPtLowMJJ.weights.xml";
  gSystem->ExpandPathName(weightFile);
  readers[method]=new TMVA::Reader( "!Color:!Silent" );

  readers[method]->AddVariable("mjj",           &vbfVars_.mjj);
  readers[method]->AddVariable("j_qg[0]",       &vbfVars_.leadj_qg);
  readers[method]->AddVariable("ystar",         &vbfVars_.ystar);
  readers[method]->AddVariable("j_qg[1]",       &vbfVars_.subleadj_qg);
  readers[method]->AddVariable("j_pt[0]",       &vbfVars_.leadj_pt);
  readers[method]->AddVariable("dphivj0",       &vbfVars_.dphivj0);
  readers[method]->AddVariable("C",             &vbfVars_.C);


  readers[method]->BookMVA(method,weightFile);

  method="BDT_VBF0HighVPtHighMJJ";
  weightFile="${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/VBF_weights/HighVPtHighMJJ_BDT_VBF0HighVPtHighMJJ.weights.xml";
  gSystem->ExpandPathName(weightFile);
  readers[method]=new TMVA::Reader( "!Color:!Silent" );
  readers[method]->AddVariable("ystar",         &vbfVars_.ystar);
  readers[method]->AddVariable("mjj",           &vbfVars_.mjj);
  readers[method]->AddVariable("j_qg[0]",       &vbfVars_.leadj_qg);
  readers[method]->AddVariable("circularity",   &vbfVars_.circularity);
  readers[method]->AddVariable("sphericity",    &vbfVars_.sphericity);
  readers[method]->AddVariable("j_qg[1]",       &vbfVars_.subleadj_qg);
  readers[method]->AddVariable("jjpt",          &vbfVars_.jjpt);
  readers[method]->AddVariable("dphijj",        &vbfVars_.dphijj);
  readers[method]->AddVariable("j_pt[0]+j_pt[1]+gamma_pt[0]",          &vbfVars_.syspt);

  readers[method]->BookMVA(method,weightFile);

  
  //read the transformations based on CDF^{-1}
  weightFile="${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/VBF_weights/inverse_cdfs.root";
  gSystem->ExpandPathName(weightFile);
  TFile *fcdf=TFile::Open(weightFile);
  for(std::map<TString,TMVA::Reader *>::iterator it=readers.begin(); it!=readers.end(); it++) {
    TString key=it->first;
    mvaCDFinv[key]=(TGraph *)fcdf->Get(key+"_cdfinv");
     TString key1=it->first;
    mvaCDFinv[key1]=(TGraph *)fcdf->Get(key1+"_cdfinv");
    TString key2=it->first;
    mvaCDFinv[key2]=(TGraph *)fcdf->Get(key2+"_cdfinv");

  }
  

  ///////////////////////
  // LOOP OVER EVENTS //
  /////////////////////
    for (Int_t iev=0;iev<nentries_;iev++)
  // for (Int_t iev=0;iev<1000;iev++)    
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
	  if (era_.Contains("2016"))
	    {
	      bool passTrigger=(selector_->hasTriggerBit("HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v",    ev_.triggerBits) ||
				selector_->hasTriggerBit("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", ev_.triggerBits) );
	    
	      if(!passTrigger) continue;
	      chTag="EE";
	      isBosonTrigSafe=true;
	    }
	  
	  //..........2017
	  else 
	    {

	      bool passTrigger=(selector_->hasTriggerBit("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v",    ev_.triggerBits) ||
				selector_->hasTriggerBit("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", ev_.triggerBits) );
	      
	      if(!passTrigger) continue;
	      chTag="EE";
	      isBosonTrigSafe=true;
	    }
	}
	
	
        else if(selCode==13*13) {
          if(ev_.isData && !selector_->isSingleMuonPD() ) continue;
	  if (era_.Contains("2016"))
	    {
	      //Lowest unprescaled paths from
	      //https://twiki.cern.ch/twiki/bin/view/CMS/TopTriggerYear2016
	      //Run B:    HLT_Iso(Tk)Mu20_v 
	      //Run C--D: HLT_Iso(Tk)Mu22_v
	      //Run E--H: HLT_Iso(Tk)Mu24_v 
	      bool passTrigger=(selector_->hasTriggerBit("HLT_IsoMu24_v",     ev_.triggerBits) ||
				selector_->hasTriggerBit("HLT_IsoTkMu24_v",     ev_.triggerBits)
			       );     
	      if(!passTrigger) continue;
	      chTag="MM";
	      isBosonTrigSafe=true;
	    }
	  
	  else 
	    {
	      //Lowest unprescaled paths from                                                                                                                                                              
              //https://twiki.cern.ch/twiki/bin/view/CMS/TopTriggerYear2017
	      //Run B:    HLT_Iso(Tk)Mu24_v, HLT_IsoMu24_2p1_v, HLT_IsoMu27_v
	      //Run C--D: HLT_IsoMu24_2p1_v, HLT_IsoMu27_v
	      //Run E--F: HLT_IsoMu27_v
	      bool passTrigger=selector_->hasTriggerBit("HLT_IsoMu27_v",ev_.triggerBits);
	      if(!passTrigger) cout<<"trigger failed"<<endl;//continue;
	      chTag="MM";
	      isBosonTrigSafe=true;
	    }
	  
        }
      }
      
      //photons
      photons_.clear();
      relaxedTightPhotons_.clear();
      tmpPhotons_.clear();
      fakePhotons_.clear();
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
        relaxedTightPhotons_  = selector_->selPhotons(allPhotons,SelectionTool::RELAXEDTIGHT, leptons);
        tmpPhotons_           = selector_->selPhotons(allPhotons,SelectionTool::QCDTEMP, leptons);     
	inclusivePhotons_     = selector_->selPhotons(allPhotons,SelectionTool::CONTROL, leptons);
	for(auto a : inclusivePhotons_) {
	  int idx = a.originalReference();
	  if (!selector_->isFakePhoton(ev_,idx)) continue;
	  fakePhotons_.push_back(a);
	}

        //Special treatments for based on 
	//being SR or CR (CR_ option), 
	//aiming to extract the QCD template (QCDTemp option), 
	//aiming to apply the fake rate (SRfake option).
	//Also differences between jet data and photon data
        if(!CR_){
	  if (!SRfake_) {
	    if( photons_.size()>=1) chTag="A";
	  } else {
	    if(fakePhotons_.size()>=1){
	      chTag="A";
	      photons_.clear();
	      photons_   = fakePhotons_;
	    }
	  }
	} else {
	  if(SRfake_) chTag = "";
	  bool passPhoton = (!SRfake_ && !QCDTemp_ && inclusivePhotons_.size()>=1) || (!SRfake_ && QCDTemp_ && tmpPhotons_.size()>=1);
	  if(passPhoton) {
	    chTag="A";
	    photons_.clear();
	    if(!QCDTemp_)      photons_   =inclusivePhotons_;
	    else               photons_   =tmpPhotons_;
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
            isBosonTrigSafe = (r9_>0.9 && eveto && !pixelseed && hoe_<0.1); 
          }else{
            isBosonTrigSafe=true;
          }          
          
          //trigger will be verified later when dividing the categories
        }
        
        //MC truth for the selected photons
	for(auto a : photons_) {
	  int idx=a.originalReference();
	  bool isPrompt(ev_.gamma_isPromptFinalState[idx]);
	  if (selector_->isFakePhoton(ev_,idx)){
	    fakeACR.push_back(a);
	    mults[isPrompt ? "looseprompt" : "loosefake"]++;
	  } else if (selector_->isTight(ev_,idx)){
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
      //qg discriminator re-weighting uncertainty
      //https://twiki.cern.ch/twiki/bin/view/CMS/QuarkGluonLikelihood#Systematics
      std::vector<Jet> alljets = selector_->getGoodJets(ev_,30.,4.7,leptons,photons_);
      std::vector<Jet> jets;

      float qgWgt_g(1.0),qgWgt_q(1.0);
      // bool is_Loose=0;
      // bool is_Mid=0;
      // bool is_tight=0;
      for(auto j : alljets) {
        int idx=j.getJetIndex();
        if(cleanEENoise_ && fabs(j.Eta())>2.7 && fabs(j.Eta())<3 && ev_.j_emf[idx]>0.55) continue;
        int jid=ev_.j_id[idx];
       
        bool passPu((jid>>jetPuId_)&0x1);
        bool passLoosePu((jid>>2)&0x1);

        //for the control region the second jet is required to fail a loose PU ID
	if(!CR_) {
	  if(!passPu) continue;
	} else {
	  if(jets.size()!=1 && !passPu) continue;
	  if(jets.size()==1 && passLoosePu) continue;
	}
        
        if(!ev_.isData && jets.size()<2 && fabs(j.Eta())<2.0){
          float xqg=ev_.j_qg[idx];
          int jflav(abs(ev_.j_flav[idx]));
          if(jflav==21) {
            double qgsf=-55.7067*pow(xqg,7) + 113.218*pow(xqg,6) -21.1421*pow(xqg,5) -99.927*pow(xqg,4) + 92.8668*pow(xqg,3) -34.3663*pow(xqg,2) + 6.27*xqg + 0.612992;
            qgWgt_g *= min(max(0.5,qgsf),2.);
          }
          else if(jflav!=0) {
            double qgsf=-0.666978*pow(xqg,3) + 0.929524*pow(xqg,2) -0.255505*xqg + 0.981581;
            qgWgt_q *= min(max(0.5,qgsf),2.);
          }
        }
        
        jets.push_back(j);
      }


      //variables
      vbfVars_.fillDiscriminatorVariables(boson,jets,ev_);
      //Flexiblity for jet pt cut --- similar to reselect part       
      if(vbfVars_.leadj_pt    < leadJetPt_)    continue;
      if(vbfVars_.subleadj_pt < subLeadJetPt_) continue;

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
      cat[3]  = (passLowVPtTrig && boson.Pt()>lowVPtCut_  && boson.Pt()<=highVPtCut_ && isBosonTrigSafe); 
      cat[4]  = (passHighVPtTrig && boson.Pt()>highVPtCut_ && isBosonTrigSafe);
      if(jets.size()>=2) {

	cat[5]  =  (vbfVars_.mjj>looseMJJCut_  && vbfVars_.mjj<=lowMJJCut_ );//LowMjj200-500
        cat[6]  =  (vbfVars_.mjj>lowMJJCut_ && vbfVars_.mjj<=highMJJCut_  );//HighMjj500-1TeV
	cat[7]  =  (vbfVars_.mjj>lowMJJCut_ );//500GeV
	//	cat[8]  =  (vbfVars_.mjj>200. && vbfVars_.mjj<500.);
      }

      category_.set(cat);


      std::vector<TString> chTags( category_.getChannelTags() );

      //veto untriggerable events for photons
      if(chTag=="A" && cat[3]) {
        if(fabs(boson.Eta())>lowVPtMaxRapCut_
	   // || vbfVars_.mjj<lowVPtMinMJJCut_
           || vbfVars_.detajj<lowVPtDetaJJCut_)
          chTags.clear();
      }

   
     
      if(chTags.size()==0) continue;

      //evaluate discriminator MVA for categories of interest
      //FIXME: this probably needs to be modified for the new training





       for(unsigned int icat = 0; icat<chTags.size(); icat++){
	       
	int pos(chTags[icat].EndsWith("A")? chTags[icat].Sizeof()-1 : chTags[icat].Sizeof()-2);
	std::string s(chTags[icat]);
	//	std::cout << "Debug:=category " << pos << "  "  << s << "  "  << s.substr(0,pos-1) << std::endl;
	TString key("BDT_VBF0"+s.substr(0,pos-1));
	vbfmva_[chTags[icat]]      = (readers[key]?readers[key]->EvaluateMVA(key):-1000);
	flat_vbfmva_[chTags[icat]] = (readers[key]?readers[key]->EvaluateMVA(key):-1000);
	if(mvaCDFinv[key]) flat_vbfmva_[chTags[icat]]=max(0.,mvaCDFinv[key]->Eval(vbfmva_[icat]));


	  if(doBlindAnalysis_ && ev_.isData && chTags[icat].EndsWith("A") && flat_vbfmva_[chTags[icat]]>0.8) flat_vbfmva_[chTags[icat]]=-1000;
	  if(doBlindAnalysis_ && ev_.isData && chTags[icat].EndsWith("A") && vbfmva_[chTags[icat]]>0.8) vbfmva_[chTags[icat]]=-1000;
	  
	 
      }
     

     

      

     
      
      
     
	
      // if (cat[5] || cat[6]) {
      // 	TString key(cat[3] ?"BDT_VBF0LowVPtHighMJJ":(cat[5] ? "BDT_VBF0HighVPtLowMJJ" : "BDT_VBF0HighVPtHighMJJ"));
      // 	vbfmva_      = readers[key]->EvaluateMVA(key);
      // 	flat_vbfmva_ = readers[key]->EvaluateMVA(key);
      // 	if(mvaCDFinv[key]) flat_vbfmva_=max(0.,mvaCDFinv[key]->Eval(vbfmva_));
      // 	if(doBlindAnalysis_ && ev_.isData && flat_vbfmva_>0.8) flat_vbfmva_=-1000;
      // 	if(doBlindAnalysis_ && ev_.isData && vbfmva_>0.8) vbfmva_=-1000;
      // }
            

      ////////////////////
      // EVENT WEIGHTS //
      //////////////////
      float wgt(1.0);
      std::vector<float>puWgts(3,1.0);
      EffCorrection_t trigSF(1.0,0.),selSF(1.0,0.),l1prefireProb(1.0,0.);
      // std::cout << "Debug: Is data?" << ev_.isData << " current wt:" << wgt << std::endl;
      if (!ev_.isData) {

        // norm weight
	//	std::cout << "Norm weight:" << normH_->GetBinContent(1) << std::endl;
	wgt  = (normH_? normH_->GetBinContent(1) : 1.0);
	//	std::cout << "Debug: Norm weight is multiplied, current wt:" << wgt << std::endl;
	
        // pu weight
        ht_->fill("puwgtctr",0,plotwgts);
        puWgts=lumi_->pileupWeight(ev_.g_pu,period);
	
        std::vector<double>puPlotWgts(1,puWgts[0]);
        ht_->fill("puwgtctr",1,puPlotWgts);

        //L1 prefire probability
        l1prefireProb=l1PrefireWR_->getCorrection(jets,photons_);

	wgt *= l1prefireProb.first;
	//	std::cout << "Debug: L1 Prefire weight is multiplied, current wt:" << wgt << std::endl;	

        // photon trigger*selection weights        
        if(chTag=="A")
          {
            trigSF = gammaEffWR_->getPhotonTrigCorrection(photons_[0].pt(),vbfVars_.mjj);
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
	//	std::cout << "Debug: puWt*trigSF*selecionSF are  multiplied, current wt:" << wgt << std::endl;	
        
        // generator level weights
	//	std::cout << "Generator level weight:" << ev_.g_w[0] << std::endl;
	wgt *= (ev_.g_nw>0 ? ev_.g_w[0] : 1.0);//if condition is true ? calculate otherwise do : 

	//	std::cout << "Debug:  gen level weight is multiplied, current wt:" << wgt << std::endl;
	
        //update weight for plotter
        plotwgts[0]=wgt;

	//	std::cout << "Debug:  gen level weight is multiplied with wgt, current wt:" <<  plotwgts[0] << std::endl;
      }

      //gen level
      genMjj_=0;
      genAPt_=0;
      if(!ev_.isData){
        std::vector<Particle> genPhotons=selector_->getGenPhotons(ev_,50.,2.5);
        std::vector<Particle> genLeptons=selector_->getGenPhotons(ev_,20.,2.5);
        std::vector<Jet> genJets=selector_->getGenJets(ev_,30.,4.7,genLeptons,genPhotons);
        genAPt_ = (genPhotons.size() ? genPhotons[0].Pt() : 0.);
        genMjj_ = (genJets.size()>1 ? (genJets[0]+genJets[1]).M() : 0.);        
      }

      
      bool FRmeasured(chTags.size() > 1);
      //fill control histograms
      for( auto c : chTags){
	//	std::cout << "Debug: Filling control histograms for tag:" << c << " with current wt:" << plotwgts[0] << std::endl;
	float myWgt(plotwgts[0]);
	//fake rate      
	if(ev_.isData && chTag=="A" && SRfake_ && FRmeasured) {
	  if(c == "A") continue;
	  TString FRcat(c);
	  FRcat.Remove(FRcat.First(chTag),1);
	  //	  myWgt*=fr_->getWeight(FRcat, vbfVars_.mjj, photons_[0].Eta());	 
	  myWgt*=fr_->getWeight(FRcat, photons_[0].Pt(), photons_[0].Eta());	 
	}
	fillControlHistos( boson,  jets,  myWgt, c, mults, fakeACR, tightACR);

      }

      //fill tree  weight for the training 
      	
      if(skimtree_) {
      
        evtWeight_ = plotwgts[0]*xsec_;
        training_ = useForTraining(); 	
	if(!cat[2]) continue;
	if(!cat[3] && !cat[4]) continue;
	// avoid overtraining
	if (float(iev) > nentries_*0.2) break;
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
	std::map<TString,float> imva(vbfmva_);
	std::map<TString,float> flat_imva(flat_vbfmva_);
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
        else if(sname=="gluonqg")     iwgt *= qgWgt_g*puWgts[0]*trigSF.first*selSF.first*(l1prefireProb.first-l1prefireProb.second);
        else if(sname=="quarkqg")     iwgt *= qgWgt_q*puWgts[0]*trigSF.first*selSF.first*(l1prefireProb.first-l1prefireProb.second);
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
        if(sname.Contains("JEC") || sname.Contains("JER") ) {
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

            //shift jet energy/resolution
            float scaleVar(1.0);
            if(jecIdx<0) {
              float genJet_pt(ev_.j_g[idx]>-1 ? ev_.g_pt[ ev_.j_g[idx] ] : 0);
              TLorentzVector smearP4=jerTool_.getSmearedJet(j,genJet_pt,ev_.rho,isUpVar ? Variation::UP : Variation::DOWN);
              scaleVar=smearP4.Pt()/j.Pt();
            } 
            else {
              int jflav(abs(ev_.j_flav[idx]));
              bool flavorMatches(true);
              if(jecIdx==6 && jflav!=21) flavorMatches=false; //FlavorPureGluon
              if(jecIdx==7 && jflav>=4)  flavorMatches=false; //FlavorPureQuark
              if(jecIdx==8 && jflav!=4)  flavorMatches=false; //FlavorPureCharm
              if(jecIdx==9 && jflav!=5)  flavorMatches=false; //FlavorPureGluon
              if(flavorMatches)
                scaleVar=isUpVar ? ev_.j_jecUp[jecIdx][idx] : 2-ev_.j_jecDn[jecIdx][idx];
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
	std::vector<TString> myCat(chTags);
        if(reSelect) {
	  imva.clear();
	  flat_imva.clear();
	  myCat.clear();
          vbfVars_.fillDiscriminatorVariables(iBoson,ijets,ev_);                   
          if(vbfVars_.leadj_pt    < leadJetPt_)    continue;
          if(vbfVars_.subleadj_pt < subLeadJetPt_) continue;

          //set the new category tag

          bool isLowVPt = (passLowVPtTrig && iBoson.Pt()>lowVPtCut_ && iBoson.Pt()<=highVPtCut_ && isBosonTrigSafe); 
          bool isHighVPt = (passHighVPtTrig  && iBoson.Pt()>highVPtCut_ && isBosonTrigSafe );
	  bool isLowMJJ(false), isHighMJJ(false), isAllMJJ(false); //islooseMJJ(false);
	  if(ijets.size() >= 2){
	    
	   
	    //islooseMJJ=  (vbfVars_.mjj>200. && vbfVars_.mjj<500.);
	    isLowMJJ =  (vbfVars_.mjj>looseMJJCut_ && vbfVars_.mjj<=lowMJJCut_  );//LowMjj200-500
	    isHighMJJ =(vbfVars_.mjj>lowMJJCut_  && vbfVars_.mjj<=highMJJCut_ );//HighMjj500-1TeV
	    isAllMJJ =  (vbfVars_.mjj> lowMJJCut_);//500GeV
	  }

	  if(chTag=="A" && isLowVPt) {
	    if(fabs(iBoson.Eta())>lowVPtMaxRapCut_
	       //  || vbfVars_.mjj<lowVPtMinMJJCut_
	       || vbfVars_.detajj<lowVPtDetaJJCut_)
	      continue;
	  }
	  //  /*-----*/if(vbfVars_.mjj < lowMJJCut_) continue;
          if(isLowVPt) {
	    if(isAllMJJ) myCat.push_back("LowVPt"+chTag);
	    if(isLowMJJ) myCat.push_back("LowVPtLowMJJ"+chTag);
	    if(isHighMJJ) myCat.push_back("LowVPtHighMJJ"+chTag);
	    //  if (isLooseMJJ)myCat.push_back("LowVPtLowMJJ"+chTag);
	  } else if(isHighVPt) {
	    if(isAllMJJ) myCat.push_back("HighVPt"+chTag);
	    if(isLowMJJ)  myCat.push_back("HighVPtLowMJJ"+chTag);
	    if(isHighMJJ) myCat.push_back("HighVPtHighMJJ"+chTag);
	    // if (isLooseMJJ)myCat.push_back("HighVPtLowMJJ"+chTag);
	  }

          
          //re-evaluate MVA         
	  // if (isLowMJJ || isHighMJJ) {
	  //   TString key(isLowVPt ?"BDT_VBF0LowVPtHighMJJ":(isLowMJJ ? "BDT_VBF0HighVPtLowMJJ" : "BDT_VBF0HighVPtHighMJJ"));
	  //   imva = readers[key]->EvaluateMVA(key);	    
	  //   if(mvaCDFinv[key]) {
          //     flat_imva=max(0.,mvaCDFinv[key]->Eval(imva));	    
          //   }
          // }
	  for(unsigned int ic = 0; ic<myCat.size(); ic++){
	    int pos(myCat[ic].EndsWith("A")? myCat[ic].Sizeof()-1 : myCat[ic].Sizeof()-2);
	    std::string s(myCat[ic]);
	    TString key("BDT_VBF0"+s.substr(0,pos-1));
	    imva[myCat[ic]]           = (readers[key]?readers[key]->EvaluateMVA(key):-1000);
	    flat_imva[myCat[ic]]      = (readers[key]?readers[key]->EvaluateMVA(key):-1000);
	    if(mvaCDFinv[key]) flat_imva[myCat[ic]]=max(0.,mvaCDFinv[key]->Eval(imva[myCat[ic]]));
	  }
        }
        
	//        if( !icat.Contains("VPt") ) continue;
	if(myCat.size() == 0) continue;


        std::vector<float> tagJetResol;        
        for(size_t ij=0; ij<min((size_t)2,ijets.size()); ij++) {
          int idx=ijets[ij].getJetIndex();

          if(ev_.j_g[idx]>-1){
            float genJet_pt(ev_.g_pt[ ev_.j_g[idx] ]);
            if(genJet_pt>0)
              tagJetResol.push_back( ijets[ij].Pt()/genJet_pt-1 );
          }
        }
        
	  
        //fill with new values/weights

        double qgwgt(qgWgt_q*qgWgt_g);

        std::vector<double> eweights(1,iwgt*qgwgt);

	for(unsigned int ic = 0; ic < myCat.size(); ic++){
	  for(auto jr : tagJetResol) {
	    ht_->fill2D("tagjetresol_exp",  jr, is,eweights,myCat[ic]);
	  }	  
	  std::vector<double> eweights(1,iwgt);
	  ht_->fill2D("vbfmva_exp",       imva[myCat[ic]],      is,eweights,myCat[ic]);
	  ht_->fill2D("acdfvbfmva_exp",   flat_imva[myCat[ic]], is,eweights,myCat[ic]);
	  ht_->fill2D("centraleta_exp",   vbfVars_.centraleta,  is,eweights,myCat[ic]);
	  ht_->fill2D("forwardeta_exp",   vbfVars_.forwardeta,  is,eweights,myCat[ic]);
	  ht_->fill2D("leadpt_exp",       vbfVars_.leadj_pt,    is,eweights,myCat[ic]);
	  ht_->fill2D("subleadpt_exp",    vbfVars_.subleadj_pt, is,eweights,myCat[ic]);
	  ht_->fill2D("detajj_exp",       vbfVars_.detajj,      is,eweights,myCat[ic]);
	  ht_->fill2D("dphijj_exp",       vbfVars_.dphijj,      is,eweights,myCat[ic]);
	  ht_->fill2D("mjj_exp", 	  vbfVars_.mjj,         is,eweights,myCat[ic]);
	  ht_->fill2D("vpt_exp", 	  iBoson.Pt(),          is,eweights,myCat[ic]);
	  ht_->fill2D("evcount_exp",      0,                    is,eweights,myCat[ic]);
	  if(flat_imva[myCat[ic]]>0.9)
	    ht_->fill2D("evcount_exp",    1,                    is,eweights,myCat[ic]);
	}
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
    // std::cout << "Debug: Hname:" << it.first << "\t" << it.second->GetEntries() << std::endl;
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
  if (debug_)  nentries_ = min(500000,nentries_); //restrict number of entries for testing  
  //  nentries_ = 5000;
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
  ht_->addHist("vpt", 	       new TH1F("vectorbosonPt",    ";Boson p_{T}[GeV];Events",           20,50,550));  
  ht_->addHist("vy", 	       new TH1F("vectorbosony",     ";Boson rapidity;Events",             25,-3,3));  
  ht_->addHist("mindrl",        new TH1F("mindrl",           ";min #Delta R(boson,lepton);Events", 25,0,6));  
  ht_->addHist("mindrj",       new TH1F("mindrj",           ";min #Delta R(boson,jet);Events",    25,0,6));  
  ht_->addHist("sihih",        new TH1F("sihih",            ";#sigma(i#eta,i#eta);Events",        50,0,0.1));  
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

  // Study the jet bump
  ht_->addHist("ptcentj",          new TH2F("ptcentj",          ";Cent Jet p_{T}; #eta_{j}",25,20,280,25,0,5));  
  ht_->addHist("ptfwdj",           new TH2F("ptfwdj",           ";Fwd Jet p_{T}; #eta_{j}",25,30,280,25,0,5)); 
  ht_->addHist("gawidthcentj",     new TH2F("gawidthcentj",     ";gawidth cent jet;#eta_{j}",100,-1,1,25,0,5));
  ht_->addHist("gawidthfwdj",      new TH2F("gawidthfwdj",      ";gawidth fwd jet;#eta_{j}",100,-1,1,25,0,5));
  ht_->addHist("emfcentj",         new TH2F("emfcentj",         ";emf cent jet;#eta_{j}",100,0,2,25,0,5));
  ht_->addHist("emffwdj",          new TH2F("emffwdj",          ";emf fwd jet;#eta_{j}",100,0,2,25,0,5));
  ht_->addHist("etaphi",           new TH2F("etaphi",           ";Most central jet #eta;#phi [rad];Events",30,-4.7,4.7,25,-TMath::Pi(),TMath::Pi()));
  ht_->addHist("centjEtaVphi",      new TH2F("centjEtaVphi",      ";Most central jet #eta;#phi [rad]",30,-4.7,4.7,25,-TMath::Pi(),TMath::Pi()));
  ht_->addHist("centjEtaVpuMva",    new TH2F("centjEtaVpuMva",    ";Most central jet #eta;PUMVA",30,-4.7,4.7,25,-1,1));
  ht_->addHist("fwdjEtaVphi",       new TH2F("fwdjEtaVphi",       ";Most forward jet #eta;#phi [rad]",30,-4.7,4.7,25,-TMath::Pi(),TMath::Pi()));
  ht_->addHist("fwdjEtaVpuMva",     new TH2F("fwdjEtaVpuMva",     ";Most forward jet #eta;PUMVA",30,-4.7,4.7,25,-1,1));

  //other jet variables
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
  ht_->addHist("mjjlog",        new TH1F("mjjlog",              ";log10(Dijet invariant mass [GeV]);Events", 35,2.30,3.70));
  ht_->addHist("mjj", 	        new TH1F("mjj",              ";Dijet invariant mass [GeV];Events", 40,0,4000));   
  ht_->addHist("genapt", 	new TH1F("genapt",           ";Generator level photon p_{T} [GeV];Events", 50,50,550));  
  ht_->addHist("genmjj", 	new TH1F("genmjj",           ";Generator level dijet invariant mass [GeV];Events", 40,0,4000));  
  ht_->addHist("matchedjpt", 	new TH1F("matchedjetpt",     ";Transverse momentum [GeV];Jets", 50,0,200));  
  ht_->addHist("unmatchedjpt", 	new TH1F("unmatchedjetpt",   ";Transverse momentum [GeV];Jets", 50,0,200));  
  ht_->addHist("matchedjeta", 	new TH1F("matchedjeteta",    ";Pseudo-rapidity;Jets", 50,0,5.));  
  ht_->addHist("unmatchedjeta", new TH1F("unmatchedjeteta",  ";Pseudo-rapidity;Jets", 50,0,5.));  
  ht_->addHist("njnmatched",    new TH1F("njnmatched",       ";Matched jets;Jets", 5,0,5.));  
  ht_->addHist("njnmatchedsoft",new TH1F("njnmatchedsoft",   ";Matched jets;Jets", 5,0,5.));  

  ht_->addHist("detajj",        new TH1F("detajj" ,          ";#Delta#eta(J,J);Events",            20,0,8));  
  ht_->addHist("dphijj",        new TH1F("dphijj" ,          ";#Delta#phi(J,J) [rad];Events",      20,-3.15,3.15));
  ht_->addHist("absdphijj",     new TH1F("absdphijj" ,          ";|#Delta#phi(J,J)| [rad];Events",      20,0,3.15)); 
  ht_->addHist("dijetpt",       new TH1F("dijetpt",          ";Dijet p_{T} [GeV];Events",          20,0,1000));  
  ht_->addHist("ht",            new TH1F("ht",               ";H_{T} [GeV];Events",                20,0,4000));
  ht_->addHist("Ght",            new TH1F("Ght",               ";(#gamma+jets) H_{T} [GeV];Events",                20,0,4000)); 
  ht_->addHist("mht",           new TH1F("mht",              ";Missing H_{T} [GeV];Events",        20,0,500));  
  ht_->addHist("dijetht",       new TH1F("dijetht",          ";Dijet H_{T} [GeV];Events",          20,0,4000));  
  ht_->addHist("syspt",         new TH1F("syspt",            ";Dijet H_{T}+boson P_{T} [GeV];Events",         20,0,1000));  

  ht_->addHist("drj1b",         new TH1F("drj1b",            ";#DeltaR(j_{1},boson);Events",       25,0,8));  
  ht_->addHist("drj2b",         new TH1F("drj2b"   ,         ";#DeltaR(j_{2},boson);Events",       25,0,8));  
  ht_->addHist("vystar",        new TH1F("vectorbosonystar", ";y-(1/2)(y_{j1}+y_{j2});Events",     25,-5,5));  
  ht_->addHist("balance",       new TH1F("balance",          ";System p_{T} balance [GeV];Events", 20,0,300));  
  ht_->addHist("relbpt",        new TH1F("relbpt",           ";#Sigma p_{T}(j)/Boson p_{T};Events",25,0,2));  
  ht_->addHist("dphibjj",       new TH1F("dphibjj",          ";#Delta#phi(JJ,boson) [rad];Events", 20,-3.15,3.15));
  ht_->addHist("absdphibjj",    new TH1F("absdphibjj",          ";|#Delta#phi(JJ,boson)| [rad];Events", 20,0,3.15)); 
  ht_->addHist("sphericity",    new TH1F("sphericity",       ";Sphericity;Events",                 20,0,1.0));  
  ht_->addHist("aplanarity",    new TH1F("aplanarity",       ";Aplanarity;Events",                 20,0,1.0));  
  ht_->addHist("C",             new TH1F("C",                ";C;Events",                          20,0,1.0));  
  ht_->addHist("D",             new TH1F("D",                ";D;Events",                          20,0,1.0));  
  ht_->addHist("isotropy",      new TH1F("isotropy",         ";Isotropy;Events",                   20,0,1.0));  
  ht_->addHist("circularity",   new TH1F("circularity",      ";Circularity;Events",               20,0,1.0));
  //additional variables from https://link.springer.com/content/pdf/10.1140/epjc/s10052-017-5315-6.pdf
  ht_->addHist("jjetas",        new TH1F("jjetas",           ";#eta_{j1}#eta_{j2};Events",50,-25,15));  
  ht_->addHist("centjy",	new TH1F("centjy",           ";Central jet rapidity;Jets",25,0,5));  
  ht_->addHist("ncentj",        new TH1F("ncentjj",          ";Number of central jets;Events",5,0,5));
  ht_->addHist("dphivj0",       new TH1F("dphivj0",          ";#Delta#phi(V,lead jet);Jets",20,0,TMath::Pi()));
  ht_->addHist("dphivj1",       new TH1F("dphivj1",          ";#Delta#phi(V,sublead jet);Jets",20,0,TMath::Pi()));  
  ht_->addHist("dphivj2",       new TH1F("dphivj2",          ";#Delta#phi(V,extra jet);Jets",20,0,TMath::Pi()));  
  ht_->addHist("dphivj3",       new TH1F("dphivj3",          ";#Delta#phi(V,next extra jet);Jets",20,0,TMath::Pi()));


  ht_->addHist("cosqj1",       new TH1F("cosqj1",          ";|cos#theta^{*}(V,lead jet)|;Events",20,0,1) );
  ht_->addHist("cosqjj",       new TH1F("cosqjj",          ";|cos#theta^{*}(V,dijet)|;Events",20,0,1) );
  ht_->addHist("betavj2",      new TH1F("betavj2",         ";#beta(V,sublead jet) [rad];Events",20,0,TMath::Pi()) );
  ht_->addHist("betaj1j2",     new TH1F("betaj1j2",        ";#beta(lead jet,sublead jet) [rad];Events",20,0,TMath::Pi()) );
  ht_->addHist("betavj3",      new TH1F("betavj3",         ";#beta(V,extra jet) [rad];Events",20,0,TMath::Pi()) );
  ht_->addHist("betaclosejj3", new TH1F("betaclosejj3",    ";#beta(closest tag jet,extra jet) [rad];Events",20,0,TMath::Pi()) );



  //final analyses distributions
  ht_->addHist("evcount",         new TH1F("evcount",        ";Pass;Events",2,0,2));  
  ht_->getPlots()["evcount"]->GetXaxis()->SetBinLabel(1,"Inclusive");
  ht_->getPlots()["evcount"]->GetXaxis()->SetBinLabel(2,"MVA>0.9");

  ht_->addHist("vbfmva",          new TH1F("vbfmva",         ";VBF MVA;Events",50,-1,1));  
  ht_->addHist("vbfmvaAcc",       new TH1F("vbfmvaAcc",      ";VBF MVA;Events",50,-1,1));  
  ht_->addHist("acdfvbfmva",     new TH1F("acdfvbfmva",    ";CDF^{-1}(VBF MVA);Events",50,0,1));  
  ht_->addHist("tagjetresol", new TH1F("tagjetresol",";#Delta p_{T}/p_{T};Jets",50,-0.5,0.5));
   ht_->addHist("vbfmvamjj",          new TH2F("vbfmvamjj",          ";m_{jj}; vbfmva",20,0,4000,50,-1,1)); 
  //  ht_->addHist("vbfmvaHVPt",          new TH1F("vbfmvaHVPt",         ";VBF MVA;Events",50,-1,1));
  // ht_->addHist("vbfmvaLVPt",          new TH1F("vbfmvaLVPt",         ";VBF MVA;Events",50,-1,1));


  TString expSystNames[]={"puup","pudn","trigup","trigdn","selup","seldn","l1prefireup","l1prefiredn",
                          "gluonqg","quarkqg",
                          "aesup","aesdn",
                          "mesup","mesdn",
                          "JERup","JERdn",
                          "AbsoluteStatJECup","AbsoluteScaleJECup","AbsoluteMPFBiasJECup","FragmentationJECup","SinglePionECALJECup","SinglePionHCALJECup","FlavorPureGluonJECup","FlavorPureQuarkJECup","FlavorPureCharmJECup","FlavorPureBottomJECup","TimePtEtaJECup","RelativeJEREC1JECup","RelativeJEREC2JECup","RelativeJERHFJECup","RelativePtBBJECup","RelativePtEC1JECup","RelativePtEC2JECup","RelativePtHFJECup","RelativeBalJECup","RelativeFSRJECup","RelativeStatFSRJECup","RelativeStatECJECup","RelativeStatHFJECup","PileUpDataMCJECup","PileUpPtRefJECup","PileUpPtBBJECup","PileUpPtEC1JECup","PileUpPtEC2JECup","PileUpPtHFJECup",
                          "AbsoluteStatJECdn","AbsoluteScaleJECdn","AbsoluteMPFBiasJECdn","FragmentationJECdn","SinglePionECALJECdn","SinglePionHCALJECdn","FlavorPureGluonJECdn","FlavorPureQuarkJECdn","FlavorPureCharmJECdn","FlavorPureBottomJECdn","TimePtEtaJECdn","RelativeJEREC1JECdn","RelativeJEREC2JECdn","RelativeJERHFJECdn","RelativePtBBJECdn","RelativePtEC1JECdn","RelativePtEC2JECdn","RelativePtHFJECdn","RelativeBalJECdn","RelativeFSRJECdn","RelativeStatFSRJECdn","RelativeStatECJECdn","RelativeStatHFJECdn","PileUpDataMCJECdn","PileUpPtRefJECdn","PileUpPtBBJECdn","PileUpPtEC1JECdn","PileUpPtEC2JECdn","PileUpPtHFJECdn"};
  
  //instantiate 2D histograms for most relevant variables to trace with systs

  TString hoi[]={"vbfmva","acdfvbfmva","evcount","mjj","detajj","dphijj","leadpt","subleadpt","forwardeta","centraleta","vpt","tagjetresol","vbfmvaAcc"};



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

  double binsPt[]={75,100,150,200,250,300,400,600,1000};
  //2D's for Pt-binned FR
  ht_->addHist("relaxedTightPtEB",  new TH2F("relaxedTightPtEB",";Relaxed tight #sigma_{i#etai#eta}; p_{T,#gamma} (GeV)",100,0,0.05,8,binsPt)); //80,0,4000
  ht_->addHist("tightPtEB",         new TH2F("tightPtEB",";Tight #sigma_{i#etai#eta}; p_{T,#gamma} (GeV)",100,0,0.05,8,binsPt));
  ht_->addHist("loosePtEB",         new TH2F("loosePtEB",";Loose #sigma_{i#etai#eta}; p_{T,#gamma} (GeV)",100,0,0.05,8,binsPt));
  ht_->addHist("tmpQCDPtEB",        new TH2F("tmpQCDPtEB",";All #sigma_{i#etai#eta}; p_{T,#gamma} (GeV)",100,0,0.05,8,binsPt));
  ht_->addHist("relaxedTightPtEE",  new TH2F("relaxedTightPtEE",";Relaxed tight #sigma_{i#etai#eta}; p_{T,#gamma} (GeV)",100,0,0.05,8,binsPt)); //80,0,4000
  ht_->addHist("tightPtEE",         new TH2F("tightPtEE",";Tight #sigma_{i#etai#eta}; p_{T,#gamma} (GeV)",100,0,0.05,8,binsPt));
  ht_->addHist("loosePtEE",         new TH2F("loosePtEE",";Loose #sigma_{i#etai#eta}; p_{T,#gamma} (GeV)",100,0,0.05,8,binsPt));
  ht_->addHist("allPtEE",           new TH2F("allPtEE",";All #sigma_{i#etai#eta}; p_{T,#gamma} (GeV)",100,0,0.05,8,binsPt));
  ht_->addHist("tmpQCDPtEE",        new TH2F("tmpQCDPtEE",";All #sigma_{i#etai#eta}; p_{T,#gamma} (GeV)",100,0,0.05,8,binsPt));
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
  newTree_->Branch("dijetht",          &vbfVars_.dijetht);
  newTree_->Branch("syspt",            &vbfVars_.syspt);
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
        if (fabs(ev_.gamma_eta[idx]) < 1.442){
          ht_->fill2D("allMjjEB"  ,   ev_.gamma_sieie[idx], vbfVars_.mjj        ,cplotwgts,c); 
	  ht_->fill2D("allPtEB"  ,    ev_.gamma_sieie[idx], a.Pt()              ,cplotwgts,c); 
	} else {
          ht_->fill2D("allMjjEE"  ,   ev_.gamma_sieie[idx], vbfVars_.mjj        ,cplotwgts,c); 
	  ht_->fill2D("allPtEE"  ,    ev_.gamma_sieie[idx], a.Pt()              ,cplotwgts,c); 
	}
      }
      for(auto a : relaxedTightPhotons_) {
        int idx = a.originalReference();
        ht_->fill("relaxedTightsihih",   ev_.gamma_sieie[idx]             ,cplotwgts,c);
        if (fabs(ev_.gamma_eta[idx]) < 1.442){
          ht_->fill2D("relaxedTightMjjEB"  ,   ev_.gamma_sieie[idx], vbfVars_.mjj        ,cplotwgts,c); 
	  ht_->fill2D("relaxedTightPtEB"  ,    ev_.gamma_sieie[idx], a.Pt()              ,cplotwgts,c); 
	} else {
          ht_->fill2D("relaxedTightMjjEE"  ,   ev_.gamma_sieie[idx], vbfVars_.mjj        ,cplotwgts,c); 
	  ht_->fill2D("relaxedTightPtEE"  ,    ev_.gamma_sieie[idx], a.Pt()              ,cplotwgts,c); 
	}
      }
      for(auto a : tmpPhotons_) {
        int idx = a.originalReference();
        if (fabs(ev_.gamma_eta[idx]) < 1.442){
          ht_->fill2D("tmpQCDMjjEB"  ,   ev_.gamma_sieie[idx], vbfVars_.mjj        ,cplotwgts,c); 
	  ht_->fill2D("tmpQCDPtEB"  ,    ev_.gamma_sieie[idx], a.Pt()              ,cplotwgts,c); 
	} else {
          ht_->fill2D("tmpQCDMjjEE"  ,   ev_.gamma_sieie[idx], vbfVars_.mjj        ,cplotwgts,c); 
	  ht_->fill2D("tmpQCDPtEE"  ,    ev_.gamma_sieie[idx], a.Pt()              ,cplotwgts,c); 
	}
      }
      //bosons in CR and fakes
      for(auto a : fakeACR) {
        int idx = a.originalReference();
        ht_->fill("fakesihih",   ev_.gamma_sieie[idx]             ,cplotwgts,c);
        ht_->fill("fakechiso",   ev_.gamma_chargedHadronIso[idx]  ,cplotwgts,c);
        ht_->fill("fakeneutiso", ev_.gamma_neutralHadronIso[idx]  ,cplotwgts,c);
        ht_->fill("fakeaiso",    ev_.gamma_photonIso[idx]         ,cplotwgts,c);
        if (fabs(ev_.gamma_eta[idx]) < 1.442){
          ht_->fill2D("looseMjjEB"  ,  ev_.gamma_sieie[idx], vbfVars_.mjj        ,cplotwgts,c);
	  ht_->fill2D("loosePtEB"   ,  ev_.gamma_sieie[idx], a.Pt()              ,cplotwgts,c);
        } else {
          ht_->fill2D("looseMjjEE"  ,  ev_.gamma_sieie[idx], vbfVars_.mjj        ,cplotwgts,c);
	  ht_->fill2D("loosePtEE"   ,  ev_.gamma_sieie[idx], a.Pt()              ,cplotwgts,c);
	}
      }
      for(auto a : tightACR) { 
        int idx = a.originalReference();
        ht_->fill("tightsihih",   ev_.gamma_sieie[idx]             ,cplotwgts,c);
        ht_->fill("tightchiso",   ev_.gamma_chargedHadronIso[idx]  ,cplotwgts,c);
        ht_->fill("tightneutiso", ev_.gamma_neutralHadronIso[idx]  ,cplotwgts,c);
        ht_->fill("tightaiso",    ev_.gamma_photonIso[idx]         ,cplotwgts,c);
        if (fabs(ev_.gamma_eta[idx]) < 1.442){
          ht_->fill2D("tightMjjEB"  ,   ev_.gamma_sieie[idx], vbfVars_.mjj        ,cplotwgts,c);
	  ht_->fill2D("tightPtEB"   ,   ev_.gamma_sieie[idx], a.Pt()              ,cplotwgts,c);
        } else {
          ht_->fill2D("tightMjjEE"  ,   ev_.gamma_sieie[idx], vbfVars_.mjj        ,cplotwgts,c);
	  ht_->fill2D("tightPtEE"   ,   ev_.gamma_sieie[idx], a.Pt()              ,cplotwgts,c);
	}
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
  ht_->fill("dijetht",      vbfVars_.dijetht,     cplotwgts,c);
  ht_->fill("syspt",        vbfVars_.syspt,       cplotwgts,c);

  std::vector<float> tagJetResol;        
  for(size_t ij=0; ij<min((size_t)2,jets.size()); ij++) {
    int idx=jets[ij].getJetIndex();
    if(ev_.j_g[idx]>-1){
      float genJet_pt(ev_.g_pt[ ev_.j_g[idx] ]);
      if(genJet_pt>0)
        tagJetResol.push_back( jets[ij].Pt()/genJet_pt-1 );
    }
  }
  for(auto jr : tagJetResol) {
    ht_->fill("tagjetresol",  jr,cplotwgts,c);
  }

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
     
    //if(is_Loose) 
    TString postfix("centj");
    if(fabs(fabs(jets[ij].Eta())-vbfVars_.forwardeta)<0.05) {
      postfix = "fwdj";
    }
    ht_->fill2D(postfix+"EtaVphi",    jets[ij].Eta(), jets[ij].Phi(),      cplotwgts,c);
    ht_->fill2D(postfix+"EtaVpuMva",  jets[ij].Eta(), jets[ij].getPUMVA(), cplotwgts,c);
    ht_->fill2D("pt"+postfix,      fabs(jets[ij].Eta()), jets[ij].Pt(), cplotwgts,c);
    ht_->fill2D("gawidth"+postfix, ev_.j_gawidth[jets[ij].getJetIndex()],  fabs(jets[ij].Eta()), cplotwgts,c);
    ht_->fill2D("emf"+postfix,     ev_.j_emf[jets[ij].getJetIndex()],      fabs(jets[ij].Eta()), cplotwgts,c);
  }
  ht_->fill("centraleta",   vbfVars_.centraleta,  cplotwgts,c);
  ht_->fill("forwardeta",   vbfVars_.forwardeta,  cplotwgts,c);
  ht_->fill("dijetpt",      vbfVars_.jjpt,        cplotwgts,c);
  ht_->fill("detajj",       vbfVars_.detajj,      cplotwgts,c);
  ht_->fill("dphijj",       vbfVars_.dphijj,      cplotwgts,c);
  ht_->fill("absdphijj",    fabs(vbfVars_.dphijj),      cplotwgts,c);


  bool  isLowMJJ(false);
  if(jets.size()>=2){
   isLowMJJ =  (vbfVars_.mjj>looseMJJCut_ && vbfVars_.mjj<=lowMJJCut_  );
 }
  if ( isLowMJJ ){
  ht_->fill("mjjlog", 	   log10(vbfVars_.mjj),         cplotwgts,c);
  }
 
  ht_->fill("mjj", 	    vbfVars_.mjj,         cplotwgts,c);
 
    
  ht_->fill("genmjj", 	    genMjj_,              cplotwgts,c);
  ht_->fill("genapt", 	    genAPt_,              cplotwgts,c);   
  if(genAPt_>0 && genMjj_<120){
    int nmatched(0),nmatchedsoft(0);
    for(size_t i=0; i<jets.size(); i++) {
      int idx=jets[i].getJetIndex();
      float genJet_pt(ev_.j_g[idx]>-1 ? ev_.g_pt[ ev_.j_g[idx] ] : 0);
      if(genJet_pt>0) {
        if(i<=1) nmatched++;
        else     nmatchedsoft++;
      }
      if(i>1) continue;
      TString pf(genJet_pt>0 ? "matched" : "unmatched");
      ht_->fill(pf+"jpt",  jets[i].Pt(),        cplotwgts,c);
      ht_->fill(pf+"jeta",  fabs(jets[i].Eta()), cplotwgts,c);   
    }
    ht_->fill("njnmatched",     nmatched,     cplotwgts,c);
    ht_->fill("njmatchedsoft",  nmatchedsoft, cplotwgts,c);   
  }

  ht_->fill("jjetas",       vbfVars_.jjetas,      cplotwgts,c);
  ht_->fill("dphivj0",      vbfVars_.dphivj0 ,    cplotwgts,c);
  ht_->fill("dphivj1",      vbfVars_.dphivj1 ,    cplotwgts,c);
  /* ht_->fill("leadeta",      leadeta,              cplotwgts,c);
  ht_->fill("subleadeta",   subleadeta,           cplotwgts,c);
  ht_->fill("leadpt",       leadPt,               cplotwgts,c);
  ht_->fill("subleadpt",    subleadPt,            cplotwgts,c);*/

  //central jet activity
  ht_->fill("ncentj", vbfVars_.ncentj, cplotwgts, c);
  if(vbfVars_.ncentj>0){
    ht_->fill("dphivj2", fabs(vbfVars_.dphivcentj[0]) ,  cplotwgts,c);
    if(vbfVars_.ncentj>1) 
      ht_->fill("dphivj3", fabs(vbfVars_.dphivcentj[1]) ,  cplotwgts,c);    
  }

  //additional colour flow and production
  ht_->fill("cosqj1",   fabs(vbfVars_.cosqj1),     cplotwgts, c);
  ht_->fill("cosqjj",   fabs(vbfVars_.cosqjj),     cplotwgts, c);
  ht_->fill("betavj2",  fabs(vbfVars_.beta_v_j2),  cplotwgts, c);
  ht_->fill("betaj1j2", fabs(vbfVars_.beta_j1_j2), cplotwgts, c);
  if(vbfVars_.ncentj>0){
    ht_->fill("betavj3",      fabs(vbfVars_.beta_v_j3),      cplotwgts, c);
    ht_->fill("betaclosejj3", fabs(vbfVars_.beta_closej_j3), cplotwgts, c);
  }

  //visible system histos
  ht_->fill("relbpt",       vbfVars_.relbpt,      cplotwgts,c);
  ht_->fill("dphibjj",      vbfVars_.dphibjj,     cplotwgts,c);
  ht_->fill("absdphibjj",   fabs(vbfVars_.dphibjj),     cplotwgts,c);
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
 

  if(vbfmva_[c]>-1000)  {
    ht_->fill("vbfmva", vbfmva_[c], cplotwgts,c);
    ht_->fill2D("vbfmvamjj"  ,    vbfVars_.mjj,vbfmva_[c]        ,cplotwgts,c); 
    ht_->fill("acdfvbfmva", flat_vbfmva_[c], cplotwgts,c);
    std::vector<double> AccWeights(1,1);
    AccWeights[0]*=(ev_.g_nw>0 ? ev_.g_w[0] : 1.0);
    ht_->fill("vbfmvaAcc", vbfmva_[c], AccWeights,c);

    if(vbfmva_[c]>0.9)

      ht_->fill("evcount",  1, cplotwgts, c);  
  }
  // std::cout << "Debug: Filling HV and LV pt histograms with values>>>" << vbfmvaHVPt_ << "\t" << vbfmvaLVPt_ << std::endl;
  /* if (vbfmva1_ >-1000){
   ht_->fill("vbfmvaHVPt", vbfmvaHVPt_, cplotwgts,c);
   ht_->fill("vbfmvaLVPt", vbfmvaLVPt_, cplotwgts,c);
   }*/

    
 

  //theory uncertainties are filled only for MC
  if(runSysts_ && !ev_.isData && ev_.g_w[0]!=0 && normH_ && normH_->GetBinContent(1)>0 && c.Contains("VPt")) {
    
    //replicas for theory systs
    for(size_t is=0; is<weightSysts_.size(); is++){
      std::vector<double> sweights(1,cplotwgts[0]);
      std::vector<double> sgenweights(1,1);
      size_t idx=weightSysts_[is].second;
      sweights[0]    *= (ev_.g_w[idx]/ev_.g_w[0])*(normH_->GetBinContent(idx+1)/normH_->GetBinContent(1));
      sgenweights[0] *= ev_.g_w[idx];

      if(vbfmva_[c]>-1000){
	ht_->fill2D("vbfmva_th",       vbfmva_[c],               is,sweights,c);
	ht_->fill2D("vbfmvaAcc_th",    vbfmva_[c],               is,sgenweights,c);
      }
      if(flat_vbfmva_[c]>0.9)
	ht_->fill2D("acdfvbfmva_th",   flat_vbfmva_[c], is,sweights,c);
      ht_->fill2D("centraleta_th",   vbfVars_.centraleta,   is,sweights,c);
      ht_->fill2D("forwardeta_th",   vbfVars_.forwardeta,   is,sweights,c);
      ht_->fill2D("leadpt_th",       vbfVars_.leadj_pt,     is,sweights,c);
      ht_->fill2D("subleadpt_th",    vbfVars_.subleadj_pt,  is,sweights,c);
      ht_->fill2D("detajj_th",       vbfVars_.detajj,       is,sweights,c);
      ht_->fill2D("dphijj_th",       vbfVars_.dphijj,       is,sweights,c);
      ht_->fill2D("mjj_th", 	     vbfVars_.mjj,          is,sweights,c);
      ht_->fill2D("vpt_th", 	     boson.Pt(),            is,sweights,c);
      ht_->fill2D("evcount_th",      0,                     is,sweights,c);
      if(vbfmva_[c]>0.9)
        ht_->fill2D("evcount_th",    1,                     is,sweights,c);
    }  
  }
}
