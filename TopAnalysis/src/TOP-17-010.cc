#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TKey.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>

#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"
#include "TopLJets2015/TopAnalysis/interface/TOP-17-010.h"
#include "TopLJets2015/TopAnalysis/interface/EfficiencyScaleFactorsWrapper.h"

#include <vector>
#include <set>
#include <iostream>
#include <algorithm>
#include "TMath.h"

using namespace std;

//
void TOP17010::init(UInt_t scenario){

  //read input
  f_           = TFile::Open(filename_);
  isSignal_    = filename_.Contains("_TT");
  triggerList_ = (TH1F *)f_->Get("analysis/triggerList");
  if(isSignal_) {
//    weightSysts_ = getWeightSysts(f_,"TTJets2016");
   weightSysts_ = getWeightSysts(f_,filename_.Contains("2016") ? "TTJets2016" : "TTJets2017");   // 2017 weights  
   origGt_=1.31;
    origMt_=172.5;
    if ( filename_.Contains("w0p5") ) {                origGt_*=0.5; }
    if ( filename_.Contains("w4p0") ) {                origGt_*=4.0; }
    if ( filename_.Contains("166v5") ){ origMt_=166.5; origGt_=1.16; }
    if ( filename_.Contains("169v5") ){ origMt_=169.5; origGt_=1.23; }
    if ( filename_.Contains("171v5") ){ origMt_=171.5; origGt_=1.28; }
    if ( filename_.Contains("173v5") ){ origMt_=173.5; origGt_=1.34; }
    if ( filename_.Contains("175v5") ){ origMt_=175.5; origGt_=1.39; }
    if ( filename_.Contains("178v5") ){ origMt_=178.5; origGt_=1.48; }
    rbwigner_=getRBW(origMt_,origGt_);
    
    if(scenario!=0){
      int im=((scenario>>16)&0xffff);
      targetMt_=169+im*0.25;
      int ig=(scenario&0xffff);
      targetGt_=0.7+ig*0.01;
    }else{
      targetGt_=origGt_;
      targetMt_=origMt_;
    }

    cout << "[TOP-17-010] scenario is set to " << scenario << endl
         << "\t (" << origMt_ << ";" << origGt_ << ")"
         << " -> (" << targetMt_ << ";" << targetGt_ << ")" << endl;
    
    //MC 2 MC corrections
    TString mc2mcTag("");
    if( filename_.Contains("TTJets_fsrdn") )   mc2mcTag="MC13TeV_2016_TTJets_fsrdn";
    if( filename_.Contains("TTJets_fsrup") )   mc2mcTag="MC13TeV_2016_TTJets_fsrup";
    if( filename_.Contains("TTJets_hdampup") ) mc2mcTag="MC13TeV_2016_TTJets_hdampup";
    if( filename_.Contains("TTJets_hdampdn") ) mc2mcTag="MC13TeV_2016_TTJets_hdampdn";
    if( filename_.Contains("TTJets_uedn") )    mc2mcTag="MC13TeV_2016_TTJets_uedn";
    if( filename_.Contains("TTJets_ueup") )    mc2mcTag="MC13TeV_2016_TTJets_ueup";
    if( filename_.Contains("TTJets_erdon") )   mc2mcTag="MC13TeV_2016_TTJets_erdon";
    if(mc2mcTag!=""){
      cout << "Reading MC2MC corrections for " << mc2mcTag << endl;
      TFile *mcCorF=TFile::Open("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/top17010/mc2mc_corrections.root");
      TString mc2mcCorNames[]={"Rb","Rc","Rudsg","eb","ec","eudsg"};
      for(size_t icor=0; icor<sizeof(mc2mcCorNames); icor++)
        mc2mcCorr_[mc2mcCorNames[icor]]=(TGraphErrors *)mcCorF->Get(mc2mcTag+"/"+mc2mcCorNames[icor]);
      mcCorF->Close();
    }

    //JER SF breakdown
    TFile *jerSFF=TFile::Open("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/top17010/jer_relunc.root");
    TKey *key;
    TIter next( jerSFF->GetListOfKeys());
    while ((key = (TKey *) next())) {     
      TString kname=key->GetName();
      jerSFBreakdown_[kname]=(TH1 *)key->ReadObj();
      jerSFBreakdown_[kname]->SetDirectory(0);
    }
    jerSFF->Close();
      
  }

  t_ = (TTree*)f_->Get("analysis/data");
  attachToMiniEventTree(t_,ev_,true);
  nentries_ = t_->GetEntriesFast();
  if (debug_) nentries_ = 10000; //restrict number of entries for testing
  t_->GetEntry(0);

  TString baseName=gSystem->BaseName(outname_); 
  TString dirName=gSystem->DirName(outname_);
  fOut_=TFile::Open(dirName+"/"+baseName,"RECREATE");

  //corrections
  lumi_ = new LumiTools(era_,genPU_);

  std::map<TString,TString> cfgMap;
  cfgMap["m_id"]     = "TightID";
  cfgMap["m_iso"]    = "TightRelIso";
  cfgMap["m_id4iso"] = "TightIDandIPCut";
  cfgMap["e_id"]     = "Tight";
  gammaEffWR_  = new EfficiencyScaleFactorsWrapper(filename_.Contains("Data13TeV"),era_,cfgMap);
  l1PrefireWR_ = new L1PrefireEfficiencyWrapper(filename_.Contains("Data13TeV"),era_);  
  btvSF_       = new BTagSFUtil(era_);
  if(mc2mcCorr_.find("eb")!=mc2mcCorr_.end())    btvSF_->setMC2MCCorrection(BTagEntry::FLAV_B,mc2mcCorr_["eb"]);
  if(mc2mcCorr_.find("ec")!=mc2mcCorr_.end())    btvSF_->setMC2MCCorrection(BTagEntry::FLAV_C,mc2mcCorr_["eb"]);
  if(mc2mcCorr_.find("eudsg")!=mc2mcCorr_.end()) btvSF_->setMC2MCCorrection(BTagEntry::FLAV_UDSG,mc2mcCorr_["eb"]);

  //theory uncs
  fragWeights_   = getBFragmentationWeights(era_);
  semilepBRwgts_ = getSemilepBRWeights(era_);

  //start selection tool
  selector_ = new SelectionTool(filename_, debug_, triggerList_);
}


//
void TOP17010::bookHistograms() {

  ht_ = new HistTool(0);
  ht_->addHist("puwgtctr", new TH1F("puwgtctr", ";Weight sums;Events",                        2,0,2));  
  ht_->addHist("genscan",  new TH1F("gennscan", ";Parameter;Value",                           4,0,4));  
  TH1 *gscan=ht_->getPlots()["genscan"];
  gscan->SetBinContent(1,origMt_);    gscan->GetXaxis()->SetBinLabel(1,"m_{t}^{i}");
  gscan->SetBinContent(2,origGt_);    gscan->GetXaxis()->SetBinLabel(2,"#Gamma_{t}^{i}");
  gscan->SetBinContent(3,targetMt_);  gscan->GetXaxis()->SetBinLabel(3,"m_{t}^{f}");
  gscan->SetBinContent(4,targetGt_);  gscan->GetXaxis()->SetBinLabel(4,"#Gamma_{t}^{f}");
  ht_->addHist("genmass",  new TH1F("genmass",     ";Mass [GeV];Events",                         100,169,176));  
  ht_->addHist("nvtx",     new TH1F("nvtx",     ";Vertex multiplicity;Events",                100,-0.5,99.5));
  ht_->addHist("rho",      new TH1F("rho",      ";#rho;Events",                               50,0,30));
  ht_->addHist("mll", 	   new TH1F("mll",      ";Dilepton invariant mass [GeV];Events",      50,0,550));  
  ht_->addHist("ptll", 	   new TH1F("ptll",     ";Dilepton p_{T}[GeV];Events",                50,0,550));  
  ht_->addHist("l1pt",     new TH1F("l1pt",     ";Lepton 1 transverse momentum [GeV];Events", 50,20,200));
  ht_->addHist("l1eta",    new TH1F("l1eta",    ";Lepton 1 pseudo-rapidity;Events",           10,0,2.5));
  ht_->addHist("l2pt",     new TH1F("l2pt",     ";Lepton 2 transverse momentum [GeV];Events", 50,20,200));
  ht_->addHist("l2eta",    new TH1F("l2eta",    ";Lepton 2 pseudo-rapidity;Events",           10,0,2.5));
  ht_->addHist("njets",    new TH1F("njets",    ";Jet multiplicity;Events",                 6,2,8)); 
  ht_->addHist("nbjets",   new TH1F("nbjets",   ";b jet multiplicity;Events",               5,1,6));
  ht_->addHist("j1pt",     new TH1F("j1pt",     ";Jet 1 transverse momentum [GeV];Events",  50,30,200));
  ht_->addHist("j1eta",    new TH1F("j1eta",    ";Jet 1 pseudo-rapidity;Events",            10,0,2.5));
  ht_->addHist("j2pt",     new TH1F("j2pt",     ";Jet 2 transverse momentum [GeV];Events",  50,30,200));
  ht_->addHist("j2eta",    new TH1F("j2eta",    ";Jet 2 pseudo-rapidity;Events",            10,0,2.5));
  ht_->addHist("evcount",  new TH1F("evcount",  ";Pass;Events", 1,0,1));  
  ht_->addHist("drlb",     new TH1F("drlb",     ";#DeltaR(l,b);Events", 15,0,2*TMath::Pi()));  
  TFile *rIn=TFile::Open("$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/top17010/mlbresol.root");  
  std::vector<TString> templates={"mlb","ptlb"};
  for(auto t : templates) {
    TH1 *th=(TH1 *)rIn->Get(t); 
    th->GetYaxis()->SetTitle("Events");
    th->SetTitle("");
    th->SetDirectory(0);
    th->Reset("ICE");
    ht_->addHist(t,th);  
  }
  rIn->Close();

  TString expSystNames[]={"puup",        "pudn",
                          "eetrigup",    "eetrigdn",
                          "emtrigup",    "emtrigdn",
                          "mmtrigup",    "mmtrigdn",
                          "eselup",      "eseldn",
                          "mselup",      "mseldn",
                          "l1prefireup", "l1prefiredn",
                          "ees1up", "ees1dn", "ees2up", "ees2dn", "ees3up", "ees3dn", "ees4up", "ees4dn",  "ees5up", "ees5dn",  "ees6up", "ees6dn",  "ees7up", "ees7dn",
                          "mes1up", "mes1dn", "mes2up", "mes2dn", "mes3up", "mes3dn", "mes4up", "mes4dn",
                          "btagup",  "btagdn",
                          "ltagup",  "ltagdn",
                          "JERup",       "JERdn",
                          "JERstat","JERJEC", "JERPU", "JERPLI", "JERptCut", "JERtrunc", "JERpTdep", "JERSTmFE",
                          "topptup",     "topptdn",
                          "AbsoluteStatJECup","AbsoluteScaleJECup","AbsoluteMPFBiasJECup","FragmentationJECup","SinglePionECALJECup","SinglePionHCALJECup","FlavorPureGluonJECup","FlavorPureQuarkJECup","FlavorPureCharmJECup","FlavorPureBottomJECup","TimePtEtaJECup","RelativeJEREC1JECup","RelativeJEREC2JECup","RelativeJERHFJECup","RelativePtBBJECup","RelativePtEC1JECup","RelativePtEC2JECup","RelativePtHFJECup","RelativeBalJECup","RelativeFSRJECup","RelativeStatFSRJECup","RelativeStatECJECup","RelativeStatHFJECup","PileUpDataMCJECup","PileUpPtRefJECup","PileUpPtBBJECup","PileUpPtEC1JECup","PileUpPtEC2JECup","PileUpPtHFJECup",
                          "AbsoluteStatJECdn","AbsoluteScaleJECdn","AbsoluteMPFBiasJECdn","FragmentationJECdn","SinglePionECALJECdn","SinglePionHCALJECdn","FlavorPureGluonJECdn","FlavorPureQuarkJECdn","FlavorPureCharmJECdn","FlavorPureBottomJECdn","TimePtEtaJECdn","RelativeJEREC1JECdn","RelativeJEREC2JECdn","RelativeJERHFJECdn","RelativePtBBJECdn","RelativePtEC1JECdn","RelativePtEC2JECdn","RelativePtHFJECdn","RelativeBalJECdn","RelativeFSRJECdn","RelativeStatFSRJECdn","RelativeStatECJECdn","RelativeStatHFJECdn","PileUpDataMCJECdn","PileUpPtRefJECdn","PileUpPtBBJECdn","PileUpPtEC1JECdn","PileUpPtEC2JECdn","PileUpPtHFJECdn",
                          "bfragup", "bfragdn",
                          "slbrup",  "slbrdn"};
  
  //instantiate 2D histograms for most relevant variables to trace with systs
  TString hoi[]={"ptlb","drlb","mlb","evcount"};
  size_t nexpSysts=sizeof(expSystNames)/sizeof(TString);
  expSysts_=std::vector<TString>(expSystNames,expSystNames+nexpSysts);  
  for(size_t ih=0; ih<sizeof(hoi)/sizeof(TString); ih++)
    {
      TH1 *histo=ht_->getPlots()[hoi[ih]];
      
      //experimental systs
      ht_->addHist(hoi[ih]+"_exp",      
                  new TH2F(hoi[ih]+"_exp", 
                           Form(";%s;Experimental systematic variation;Events",histo->GetName()),
                           histo->GetNbinsX(),histo->GetXaxis()->GetXbins()->GetArray(),
                           nexpSysts,0,nexpSysts));
      for(size_t is=0; is<nexpSysts; is++)
        ht_->get2dPlots()[hoi[ih]+"_exp"]->GetYaxis()->SetBinLabel(is+1,expSystNames[is]);
      
      //theory systs
      size_t nthSysts(weightSysts_.size()); 
      if(nthSysts>0){
        ht_->addHist(hoi[ih]+"_th",      
                    new TH2F(hoi[ih]+"_th", 
                             Form(";%s;Theory systematic variation;Events",histo->GetName()),
                             histo->GetNbinsX(),histo->GetXaxis()->GetXbins()->GetArray(),
                             nthSysts,0,nthSysts));
        for(size_t is=0; is<nthSysts; is++)
          ht_->get2dPlots()[hoi[ih]+"_th"]->GetYaxis()->SetBinLabel(is+1,weightSysts_[is].first);
      }
    }
}



//
void TOP17010::runAnalysis()
{
  ///////////////////////
  // LOOP OVER EVENTS //
  /////////////////////
  for (Int_t iev=0;iev<nentries_;iev++)
    {
      t_->GetEntry(iev);
      if(debug_) cout << "Event number: "<<iev<<endl;
      if(iev%10000==0) printf ("\r [%3.0f%%] done", 100.*(float)iev/(float)nentries_);
      
      //////////////////
      // CORRECTIONS  //
      //////////////////
      TString period = lumi_->assignRunPeriod();
      btvSF_->addBTagDecisions(ev_);
      if(!ev_.isData) btvSF_->updateBTagDecisions(ev_);      
      
      //TRIGGER
      bool hasETrigger(  selector_->hasTriggerBit("HLT_Ele32_eta2p1_WPTight_Gsf_v",                       ev_.triggerBits) );
      bool hasMTrigger(  selector_->hasTriggerBit("HLT_IsoMu24_v",                                        ev_.triggerBits) || 
                         selector_->hasTriggerBit("HLT_IsoTkMu24_v",                                      ev_.triggerBits) );
      bool hasEMTrigger( selector_->hasTriggerBit("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",     ev_.triggerBits) ||
                         selector_->hasTriggerBit("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",    ev_.triggerBits) ||
                         selector_->hasTriggerBit("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v", ev_.triggerBits) ||
                         selector_->hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",    ev_.triggerBits) ||
                         selector_->hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", ev_.triggerBits) );
      bool hasMMTrigger( selector_->hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",                ev_.triggerBits) ||
                         selector_->hasTriggerBit("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",              ev_.triggerBits) );
      bool hasEETrigger( selector_->hasTriggerBit("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",          ev_.triggerBits) );
      if(era_.Contains("2017")) {
        hasETrigger =(selector_->hasTriggerBit("HLT_Ele35_WPTight_Gsf_v", ev_.triggerBits));
        hasMTrigger =(selector_->hasTriggerBit("HLT_IsoMu24_v",           ev_.triggerBits) ||
                      selector_->hasTriggerBit("HLT_IsoMu24_2p1_v",       ev_.triggerBits) ||
                      selector_->hasTriggerBit("HLT_IsoMu27_v",           ev_.triggerBits) );          
        hasEMTrigger=(selector_->hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",    ev_.triggerBits) ||
                      selector_->hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", ev_.triggerBits) ||
                      selector_->hasTriggerBit("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",    ev_.triggerBits) ||
                      selector_->hasTriggerBit("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v", ev_.triggerBits) ||
                      selector_->hasTriggerBit("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",     ev_.triggerBits) ||
                      selector_->hasTriggerBit("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",  ev_.triggerBits) );
        hasMMTrigger=(selector_->hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",                  ev_.triggerBits) ||
                      selector_->hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v",          ev_.triggerBits) ||
                      selector_->hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v",        ev_.triggerBits) );
        hasEETrigger=(selector_->hasTriggerBit("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v",             ev_.triggerBits) ||
                      selector_->hasTriggerBit("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",          ev_.triggerBits) );
      }
      if(ev_.isData) {
        if(era_.Contains("2016")) {
          hasEMTrigger=false;
          if(ev_.run<=280385)
            {
              hasEMTrigger |= selector_->hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v", ev_.triggerBits);
              hasEMTrigger |= selector_->hasTriggerBit("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",  ev_.triggerBits);
            }
          if(ev_.run>=278273 && ev_.run<=280385)
            {
              hasEMTrigger |= selector_->hasTriggerBit("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v", ev_.triggerBits);
            }
          if(ev_.run>=278273)
            {
              hasEMTrigger |= selector_->hasTriggerBit("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",  ev_.triggerBits);
              hasEMTrigger |= selector_->hasTriggerBit("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v", ev_.triggerBits);
              hasEMTrigger |= selector_->hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", ev_.triggerBits);
            }
        }
        if(filename_.Contains("2017E") || filename_.Contains("2017F")){
          hasMTrigger=selector_->hasTriggerBit("HLT_IsoMu27_v",ev_.triggerBits);
        }
        if(!(filename_.Contains("2017A") || filename_.Contains("2017B"))){
          hasMMTrigger=(selector_->hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v",   ev_.triggerBits) ||
                        selector_->hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v", ev_.triggerBits) );
        }
        
        //trigger consistency with PD
        if( selector_->isDoubleEGPD()       && !hasEETrigger ) continue;
        if( selector_->isDoubleMuonPD()     && !hasMMTrigger ) continue;
        if( selector_->isMuonEGPD()         && !hasEMTrigger ) continue;
        if( selector_->isSingleElectronPD() && !hasETrigger )  continue;
        if( selector_->isSingleMuonPD()     && !hasMTrigger )  continue;          
      }

      ///////////////////////////
      // RECO LEVEL SELECTION //
      /////////////////////////
      std::vector<Particle> flaggedleptons = selector_->flaggedLeptons(ev_);     
      std::vector<Particle> leptons        = selector_->selLeptons(flaggedleptons,SelectionTool::TIGHT,SelectionTool::TIGHT,20,2.5);
      std::vector<Jet> alljets             = selector_->getGoodJets(ev_,30.,2.4,leptons);
      applyMC2MC(alljets);
      TopWidthEvent twe(leptons,alljets);
      if(twe.dilcode==0) {
        continue;
      }
      else if(twe.dilcode==11*11) {
        if(!(hasEETrigger || hasETrigger)) continue;
        if(selector_->isSingleElectronPD() && hasEETrigger) continue;
      }
      else if(twe.dilcode==13*13) {
        if(!(hasMMTrigger || hasMTrigger)) continue;
        if(selector_->isSingleMuonPD() && hasMMTrigger) continue;
      }
      else if(twe.dilcode==11*13) {
        if(!(hasEMTrigger || hasETrigger || hasMTrigger)) continue;
        if(selector_->isSingleMuonPD()     && hasEMTrigger) continue;
        if(selector_->isSingleElectronPD() && (hasEMTrigger || hasMTrigger)) continue;
      }

      ////////////////////
      // EVENT WEIGHTS //
      //////////////////
      float wgt(1.0),widthWgt(1.0);
      std::vector<float>puWgts(3,1.0),topptWgts(2,1.0),bfragWgts(2,1.0),slbrWgts(2,1.0);
      EffCorrection_t trigSF(1.0,0.),l1SF(1.0,0.),l2SF(1.0,0.0),l1trigprefireProb(1.0,0.);      
      if (!ev_.isData) {

        // pu weight
        std::vector<double> plotwgts(1,1);
        ht_->fill("puwgtctr",0,plotwgts); 
        puWgts=lumi_->pileupWeight(ev_.g_pu,period);
        std::vector<double>puPlotWgts(1,puWgts[0]);
        ht_->fill("puwgtctr",1,puPlotWgts);

        //L1 prefire probability
        l1trigprefireProb=l1PrefireWR_->getCorrection(alljets);

        // photon trigger*selection weights        
        TString lperiod("");
        if(era_.Contains("2016")) {
          if(gRandom->Uniform()>16551.4/35874.8) lperiod="GH";
        }
        trigSF = gammaEffWR_->getDileptonTriggerCorrection(leptons);
        l1SF   = gammaEffWR_->getOfflineCorrection(leptons[0].id(),leptons[0].pt(),leptons[0].eta(), lperiod);
        l2SF   = gammaEffWR_->getOfflineCorrection(leptons[1].id(),leptons[1].pt(),leptons[1].eta(), lperiod);

        //for signal top pt weights        
        if(isSignal_) {
          std::vector<float> genmt;
          for(Int_t igen=0; igen<ev_.ngtop; igen++)
            {
              if(abs(ev_.gtop_id[igen])!=6) continue;
              float topsf=TMath::Exp(0.156-0.00137*ev_.gtop_pt[igen]);
              topptWgts[0] *= topsf;
              topptWgts[1] *= 1./topsf;
              genmt.push_back(ev_.gtop_m[igen]);
            }
          if(genmt.size()==2 
             && rbwigner_ 
             && targetGt_>0 && targetMt_>0 && origGt_>0  && origMt_>0
             && (targetGt_!=origGt_ || targetMt_ != origMt_) ) {
            widthWgt=weightBW(rbwigner_,genmt,targetGt_,targetMt_,origGt_,origMt_); 

            //control the re-weighted BW
            for(auto genm:genmt) {
              std::vector<double> bwWgts(1,1.0);
              ht_->fill("genmass",  genm, bwWgts);
              bwWgts[0]=widthWgt;
              ht_->fill("genmass",  genm, bwWgts,"rwgt");
            }
          }
        }

        //b-fragmentation and semi-leptonic branching fractions
        bfragWgts[0] = computeBFragmentationWeight(ev_,fragWeights_["downFrag"]);
        bfragWgts[1] = computeBFragmentationWeight(ev_,fragWeights_["upFrag"]);
        slbrWgts[0]  = computeSemilepBRWeight(ev_,semilepBRwgts_["semilepbrDown"],0,false);
        slbrWgts[1]  = computeSemilepBRWeight(ev_,semilepBRwgts_["semilepbrUp"],0,false);        

        //final nominal weight
        wgt *= (normH_? normH_->GetBinContent(1) : 1.0);
        wgt *= (ev_.g_nw>0 ? ev_.g_w[0] : 1.0);
        wgt *= widthWgt;
        wgt *= puWgts[0]*l1trigprefireProb.first*trigSF.first*l1SF.first*l2SF.first;
      }
      fillControlHistograms(twe,wgt);

      //experimental systs cycle: better not to do anything else after this...
      //final category selection is repeated ad nauseam with varied objects/weights 
      if(ev_.isData) continue;

      selector_->setDebug(false);

      //combined offline efficiencies
      EffCorrection_t combinedESF(1.0,0.0), combinedMSF(1.0,0.0);
      if(abs(leptons[0].id())==11) {
        combinedESF.second+=pow(combinedESF.first*l1SF.second,2)+pow(l1SF.first*combinedESF.second,2);
        combinedESF.first *=l1SF.first;
      }
      if(abs(leptons[0].id())==13) {
        combinedMSF.second+=pow(combinedMSF.first*l1SF.second,2)+pow(l1SF.first*combinedMSF.second,2);
        combinedMSF.first *=l1SF.first;
      }
      if(abs(leptons[1].id())==11) {
        combinedESF.second+=pow(combinedESF.first*l2SF.second,2)+pow(l2SF.first*combinedESF.second,2);
        combinedESF.first *=l2SF.first;
      }
      if(abs(leptons[1].id())==13) {
        combinedMSF.second+=pow(combinedMSF.first*l2SF.second,2)+pow(l2SF.first*combinedMSF.second,2);
        combinedMSF.first *=l2SF.first;
      }
      combinedESF.second=sqrt(combinedESF.second);
      combinedMSF.second=sqrt(combinedMSF.second);
      
      for(size_t is=0; is<expSysts_.size(); is++){
        
        //uncertainty
        TString sname=expSysts_[is];
        bool isUpVar(sname.EndsWith("up"));
        
        //base values and kinematics
        bool reSelect(false);        
        float iwgt=(ev_.g_nw>0 ? ev_.g_w[0] : 1.0);
        iwgt *= (normH_? normH_->GetBinContent(1) : 1.0);
        iwgt *= widthWgt;

        EffCorrection_t selSF(1.0,0.0);
        if(sname=="puup")       iwgt *= puWgts[1]*trigSF.first*selSF.first*l1trigprefireProb.first;
        else if(sname=="pudn")  iwgt *= puWgts[2]*trigSF.first*selSF.first*l1trigprefireProb.first;
        else if( (sname.Contains("eetrig") && twe.dilcode==11*11) ||
                 (sname.Contains("emtrig") && twe.dilcode==11*13) ||
                 (sname.Contains("mmtrig") && twe.dilcode==13*13) ) {
          float newTrigSF( max(float(0.),float(trigSF.first+(isUpVar ? +1 : -1)*trigSF.second)) );
          iwgt *= puWgts[0]*newTrigSF*selSF.first*l1trigprefireProb.first;
        }
        else if(sname.BeginsWith("esel") ) {
          float newESF( max(float(0.),float(combinedESF.first+(isUpVar ? +1 : -1)*combinedESF.second)) );
          iwgt *= puWgts[0]*trigSF.first*newESF*combinedMSF.first*l1trigprefireProb.first;
        }
        else if(sname.BeginsWith("msel") ) {
          float newMSF( max(float(0.),float(combinedMSF.first+(isUpVar ? +1 : -1)*combinedMSF.second)) );
          iwgt *= puWgts[0]*trigSF.first*combinedESF.first*newMSF*l1trigprefireProb.first;
        }
        else if(sname.BeginsWith("l1prefire") ){
          float newL1PrefireProb( max(float(0.),float(l1trigprefireProb.first+(isUpVar ? +1 : -1)*l1trigprefireProb.second)) );
          iwgt *= puWgts[0]*trigSF.first*selSF.first*newL1PrefireProb;
        }
        else if(sname=="topptup")     iwgt = wgt*topptWgts[0];
        else if(sname=="topptdn")     iwgt = wgt*topptWgts[1];
        else if(sname=="bfragup")     iwgt = wgt*bfragWgts[0];
        else if(sname=="bfragdn")     iwgt = wgt*bfragWgts[1];
        else if(sname=="slbrup")      iwgt = wgt*slbrWgts[0];
        else if(sname=="slbrdn")      iwgt = wgt*slbrWgts[1];
        else                          iwgt = wgt;           

        //leptons
        std::vector<Particle> ileptons(leptons);
        if(sname.BeginsWith("ees") || sname.BeginsWith("mes")) {
          reSelect=true;
          ileptons=flaggedleptons;
          for(size_t il=0; il<ileptons.size(); il++) {
            int id=abs(ileptons[il].id());
            int idx=ileptons[il].originalReference();
            float eScale(0.0);
            if( (id==11 && sname.Contains("ees1")) || (id==13 && sname.Contains("mes1")) ) eScale=ev_.l_scaleUnc1[idx];
            if( (id==11 && sname.Contains("ees2")) || (id==13 && sname.Contains("mes2")) ) eScale=ev_.l_scaleUnc2[idx];
            if( (id==11 && sname.Contains("ees3")) || (id==13 && sname.Contains("mes3")) ) eScale=ev_.l_scaleUnc3[idx];
            if( (id==11 && sname.Contains("ees4")) || (id==13 && sname.Contains("mes4")) ) eScale=ev_.l_scaleUnc4[idx];
            if( (id==11 && sname.Contains("ees5")) || (id==13 && sname.Contains("mes5")) ) eScale=ev_.l_scaleUnc5[idx];
            if( (id==11 && sname.Contains("ees6")) || (id==13 && sname.Contains("mes6")) ) eScale=ev_.l_scaleUnc6[idx];
            if( (id==11 && sname.Contains("ees7")) || (id==13 && sname.Contains("mes7")) ) eScale=ev_.l_scaleUnc7[idx];

            if( sname.BeginsWith("ees") && id==11 ) {
              if(isUpVar) eScale=(1+fabs(eScale)/ileptons[il].Pt());
              else        eScale=(1-fabs(eScale)/ileptons[il].Pt());
            } else if (sname.BeginsWith("mes") && id==13 ) {
              if(isUpVar) eScale=(1+fabs(1-eScale));
              else        eScale=(1-fabs(1-eScale));
            } else {
              eScale=1.0;
            }
            ileptons[il] *= eScale;
          }
          ileptons=selector_->selLeptons(ileptons,SelectionTool::TIGHT,SelectionTool::TIGHT,20,2.5);        
        }
        if(ileptons.size()<2) continue;
        
        //jets
        std::vector<Jet> ijets(alljets);
        if(sname.Contains("JEC") || sname.Contains("JER") || sname.BeginsWith("btag") || sname.BeginsWith("ltag") )  {
          reSelect=true;
          btvSF_->addBTagDecisions(ev_);

          if(sname=="btagup") btvSF_->updateBTagDecisions(ev_, "up",      "central"); 
          if(sname=="btagdn") btvSF_->updateBTagDecisions(ev_, "down",    "central"); 
          if(sname=="ltagup") btvSF_->updateBTagDecisions(ev_, "central", "up"); 
          if(sname=="ltagdn") btvSF_->updateBTagDecisions(ev_, "central", "down"); 
            
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
          std::vector<Jet> newJets = selector_->getGoodJets(ev_,20.,2.4,leptons);
          applyMC2MC(newJets);
          ijets.clear();
          for(auto j : newJets) {

            int idx=j.getJetIndex();
            int jflav(abs(ev_.j_flav[idx]));

            //shift jet energy
            float scaleVar(1.0);
            if(jecIdx<0) {
              float jerUnc=1-(isUpVar ? ev_.j_jerUp[idx] : ev_.j_jerDn[idx]);
              float jerUncSgn(jerUnc<0 ? -1 : 1);
              jerUnc=fabs(jerUnc);
              if(sname=="JERstat")   jerUnc *= getJERSFBreakdown("stat", fabs(j.eta()));
              if(sname=="JERJEC")    jerUnc *= max(getJERSFBreakdown("JECup", fabs(j.eta())), getJERSFBreakdown("JECdown", fabs(j.eta())));
              if(sname=="JERPU")     jerUnc *= max(getJERSFBreakdown("PUup",  fabs(j.eta())), getJERSFBreakdown("PUudown", fabs(j.eta())));
              if(sname=="JERPLI")    jerUnc *= max(getJERSFBreakdown("PLIup", fabs(j.eta())), getJERSFBreakdown("PLIdown", fabs(j.eta())));
              if(sname=="JERptCut")  jerUnc *= getJERSFBreakdown("ptCut", fabs(j.eta()));
              if(sname=="JERtrunc")  jerUnc *= getJERSFBreakdown("trunc", fabs(j.eta()));
              if(sname=="JERpTdep")  jerUnc *= getJERSFBreakdown("pTdep", fabs(j.eta()));
              if(sname=="JERSTmFE")  jerUnc *= getJERSFBreakdown("STmFE", fabs(j.eta()));
              scaleVar*=(1+jerUncSgn*jerUnc);
            } 
            else {
              bool flavorMatches(true);
              if(jecIdx==6 && jflav!=21) flavorMatches=false; //FlavorPureGluon
              if(jecIdx==7 && jflav>=4)  flavorMatches=false; //FlavorPureQuark
              if(jecIdx==8 && jflav!=4)  flavorMatches=false; //FlavorPureCharm
              if(jecIdx==9 && jflav!=5)  flavorMatches=false; //FlavorPureGluon
              if(flavorMatches)
                scaleVar=isUpVar ? ev_.j_jecUp[jecIdx][idx] : ev_.j_jecDn[jecIdx][idx];
            }
            j*=scaleVar;           

            ijets.push_back(j);
          }
        }
        if(ijets.size()<2) continue;

        //re-select lb pairs if needed
        TString icat=twe.cat;
        std::vector<LeptonBJetPair> lbPairs=twe.getPairs();
        if(reSelect) {
          TopWidthEvent itwe(ileptons,ijets);
          icat=itwe.cat;
          lbPairs=itwe.getPairs();
        }

        //fill with new values/weights
        std::vector<double> eweights(1,iwgt);

        ht_->fill2D("evcount_exp", 0.,  is, eweights, icat);
        for(auto p4 : lbPairs){
          ht_->fill2D("drlb_exp",  p4.getDR(),  is, eweights, icat);
          ht_->fill2D("ptlb_exp",  p4.M(),      is, eweights, icat);
          ht_->fill2D("mlb_exp",   p4.M(),      is, eweights, icat);

          TString c(icat+(p4.Pt()>100?"highpt":"lowpt"));
          std::vector<TString> cat_vec(2,c);
          cat_vec[1] += (lbPairs.size()==1 ? "1b" : "2b");
          ht_->fill2D("mlb_exp",   p4.M(),      is, eweights, cat_vec);
        }
      }

      selector_->setDebug(debug_);
    }
  
  //close input file
  f_->Close();

  //save histograms to the output
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
void TOP17010::fillControlHistograms(TopWidthEvent &twe,float &wgt) {

  //plot weight
  std::vector<double> cplotwgts(1,wgt);
  ht_->fill("nvtx",    ev_.nvtx,        cplotwgts,twe.cat);        
  ht_->fill("rho",     ev_.rho,         cplotwgts,twe.cat);        
  ht_->fill("mll",     twe.mll,         cplotwgts,twe.cat);
  ht_->fill("ptll",    twe.ptll,        cplotwgts,twe.cat);
  ht_->fill("l1pt",    twe.l1pt,        cplotwgts,twe.cat);
  ht_->fill("l1eta",   fabs(twe.l1eta), cplotwgts,twe.cat);
  ht_->fill("l2pt",    twe.l2pt,        cplotwgts,twe.cat);
  ht_->fill("l2eta",   fabs(twe.l2eta), cplotwgts,twe.cat);
  ht_->fill("njets",   twe.njets,       cplotwgts,twe.cat);
  ht_->fill("nbjets",  twe.nbjets,      cplotwgts,twe.cat);
  ht_->fill("j1pt",    twe.j1pt,        cplotwgts,twe.cat);
  ht_->fill("j1eta",   fabs(twe.j1eta), cplotwgts,twe.cat);
  ht_->fill("j2pt",    twe.j2pt,        cplotwgts,twe.cat);
  ht_->fill("j2eta",   fabs(twe.j2eta), cplotwgts,twe.cat);
  ht_->fill("evcount", 0.,              cplotwgts,twe.cat);
  std::vector<LeptonBJetPair> lbPairs=twe.getPairs();
  for(auto p4 : lbPairs) {
    ht_->fill("drlb", p4.getDR(), cplotwgts, twe.cat);
    ht_->fill("ptlb", p4.Pt(), cplotwgts, twe.cat);
    ht_->fill("mlb",  p4.M(),  cplotwgts, twe.cat);
    TString c(twe.cat+(p4.Pt()>100?"highpt":"lowpt"));
    std::vector<TString> cat_vec(2,c);
    cat_vec[1] += (lbPairs.size()==1 ? "1b" : "2b");
    ht_->fill("mlb",  p4.M(),  cplotwgts, cat_vec);
  }

  //theory uncertainties are filled only for MC
  if(ev_.isData) return;
  if(ev_.g_w[0]==0 || normH_==NULL || normH_->GetBinContent(1)==0) return;
    
  //replicas for theory systs
  for(size_t is=0; is<weightSysts_.size(); is++){
    std::vector<double> sweights(1,cplotwgts[0]);
    size_t idx=weightSysts_[is].second;
    sweights[0] *= (ev_.g_w[idx]/ev_.g_w[0])*(normH_->GetBinContent(idx+1)/normH_->GetBinContent(1));
    ht_->fill2D("evcount_th",  0., is, sweights, twe.cat);
    for(auto p4 : lbPairs) {
      ht_->fill2D("drlb_th",  p4.getDR(), is, sweights, twe.cat);
      ht_->fill2D("ptlb_th",  p4.Pt(), is, sweights, twe.cat);
      ht_->fill2D("mlb_th",   p4.M(),  is, sweights, twe.cat);
      TString c(twe.cat+(p4.Pt()>100?"highpt":"lowpt"));
      std::vector<TString> cat_vec(2,c);
      cat_vec[1] += (lbPairs.size()==1 ? "1b" : "2b");
      ht_->fill2D("mlb_th",   p4.M(),  is, sweights, cat_vec);
    }
  }
}

//
void TOP17010::applyMC2MC(std::vector<Jet> &jetColl) {

  //loop over collection of jets
  for(size_t i=0; i<jetColl.size(); i++) {

    if(jetColl[i].Pt()<20) continue;
    
    //check jet flavour and if MC2MC correction is available
    int idx=jetColl[i].getJetIndex();
    int jflav(abs(ev_.j_flav[idx]));
    TString mcTag("");
    if(jflav==5)      mcTag="Rb";
    else if(jflav==4) mcTag="Rc";
    else              mcTag="Rudsg";
    if(mc2mcCorr_.find(mcTag)==mc2mcCorr_.end()) continue;
 
    //as mc2mc=var MC / nom MC bring back JES to nom MC response
    //do not do it if it exceeds 10% corrections...
    float mc2mcVal=mc2mcCorr_[mcTag]->Eval(jetColl[i].Pt());
    if(mc2mcVal>0.9 && mc2mcVal<1.1) jetColl[i] *= 1./mc2mcVal;
  }

}

//
float TOP17010::getJERSFBreakdown(TString key,float abseta){
  if(jerSFBreakdown_.find(key)==jerSFBreakdown_.end()) return 1.0;
  int xbin=jerSFBreakdown_[key]->GetXaxis()->FindBin(abseta);
  return jerSFBreakdown_[key]->GetBinContent(xbin);
}
