#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>

#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"
#include "TopLJets2015/TopAnalysis/interface/ExclusiveTop.h"
#include "TopQuarkAnalysis/TopTools/interface/MEzCalculator.h"

#include <vector>
#include <set>
#include <iostream>
#include <algorithm>
#include <string>

#include "TMath.h"

#include "TopLJets2015/CTPPSAnalysisTools/interface/LHCConditionsFactory.h"
#include "TopLJets2015/CTPPSAnalysisTools/interface/AlignmentsFactory.h"
#include "TopLJets2015/CTPPSAnalysisTools/interface/XiReconstructor.h"

using namespace std;

#define ADDVAR(x,name,t,tree) tree->Branch(name,x,TString(name)+TString(t))

//TODO
//add jet scale uncertainty
//PPS json
//lumi, puweighting, genweighting
//xsec in the analysis json

//
void RunExclusiveTop(TString filename,
                     TString outname,
                     Int_t channelSelection, 
                     Int_t chargeSelection, 
                     TH1F *normH, 
                     TH1F *genPU, 
                     TString era,
                     Bool_t debug)
{
  /////////////////////
  // INITIALIZATION //
  ///////////////////
  const char* CMSSW_BASE = getenv("CMSSW_BASE");
  // CTPPS reconstruction part
  ctpps::XiReconstructor proton_reco;
  proton_reco.feedDispersions(Form("%s/src/TopLJets2015/CTPPSAnalysisTools/data/2017/dispersions.txt", CMSSW_BASE));

  ctpps::AlignmentsFactory ctpps_aligns;
  ctpps_aligns.feedAlignments(Form("%s/src/TopLJets2015/CTPPSAnalysisTools/data/2017/alignments_21aug2018.txt", CMSSW_BASE));

  ctpps::LHCConditionsFactory lhc_conds;
  lhc_conds.feedConditions(Form("%s/src/TopLJets2015/CTPPSAnalysisTools/data/2017/xangle_tillTS2.csv", CMSSW_BASE));
  lhc_conds.feedConditions(Form("%s/src/TopLJets2015/CTPPSAnalysisTools/data/2017/xangle_afterTS2.csv", CMSSW_BASE));
  
  MiniEvent_t ev;

  //PREPARE OUTPUT
  TString baseName=gSystem->BaseName(outname); 
  TString dirName=gSystem->DirName(outname);
  TFile *fOut=TFile::Open(dirName+"/"+baseName,"RECREATE");
  fOut->cd();
  TTree *outT=new TTree("tree","tree");
  outT->Branch("run",&ev.run,"run/i");
  outT->Branch("event",&ev.event,"event/l");
  outT->Branch("lumi",&ev.lumi,"lumi/i");
  outT->Branch("nvtx",&ev.nvtx,"nvtx/I");
  outT->Branch("metfilters",&ev.met_filterBits,"metfilters/I");  

  bool hasETrigger,hasMTrigger,hasMMTrigger,hasEETrigger,hasEMTrigger;
  outT->Branch("hasETrigger",&hasETrigger,"hasETrigger/O");
  outT->Branch("hasMTrigger",&hasMTrigger,"hasMTrigger/O");
  outT->Branch("hasEMTrigger",&hasEMTrigger,"hasEMTrigger/O");
  outT->Branch("hasMMTrigger",&hasMMTrigger,"hasMMTrigger/O");
  outT->Branch("hasEETrigger",&hasEETrigger,"hasEETrigger/O");

  bool isSS,isSF,isZ;
  outT->Branch("isSS",&isSS,"isSS/O");
  outT->Branch("isSF",&isSF,"isSF/O");
  outT->Branch("isZ",&isZ,"isZ/O");
  
  ADDVAR(&ev.rho,"rho","F",outT);
  ADDVAR(&ev.pf_closestDZnonAssoc,"closestDZnonAssoc","F",outT);
  ADDVAR(&ev.pf_ch_wgtSum,"nch","F",outT);
  ADDVAR(&ev.met_pt[1],"met_pt","F",outT);
  ADDVAR(&ev.met_phi[1],"met_phi","F",outT);
  ADDVAR(&ev.met_sig[1],"met_sig","F",outT);
  TString fvars[]={"evwgt", "dilcode", 
                   "l1pt", "l1eta", "l1phi", "ml1", "l1id", "mt1",
                   "l2pt", "l2eta", "l2phi", "ml2", "l2id", "mt2",
                   "llpt", "lleta", "llphi", "mll", "llht", "llacopl", "llcosthetaCS", "llphistar", "llMR", "llR", "mtll", 
                   "llcsip", "llcsipUnc", "llcsim", "llcsimUnc",
                   "xpt",  "xeta",  "xphi",  "mx",  "xid", "xht", "mtx",
                   "fpt",  "feta",  "fphi",  "mf",  "fht", "facopl", "fcosthetaCS", "fphistar", "fMR","fR", 
                   "fcsip", "fcsipUnc", "fcsim", "fcsimUnc",
                   "nb", "nj", "nl","ng","nch", "ht", "htb", "htj", "closestkdz",
                   "puppirecoil_pt","puppirecoil_phi", "puppirecoil_spher",
  };
  std::map<TString,Float_t> outVars;
  for(size_t i=0; i<sizeof(fvars)/sizeof(TString); i++){
    outVars[fvars[i]]=0.;
    ADDVAR(&(outVars[fvars[i]]),fvars[i],"F",outT);
  }
  float beamXangle(0);
  outT->Branch("beamXangle",&beamXangle,"beamXangle/F");
  int nRPtk(0),RPid[50];
  float RPfarcsi[50],RPnearcsi[50];
  if(filename.Contains("Data13TeV")){
    outT->Branch("nRPtk",&nRPtk,"nRPtk/i");
    outT->Branch("RPid",RPid,"RPid[nRPtk]/i");
    outT->Branch("RPfarcsi",RPfarcsi,"RPfarcsi[nRPtk]/F");
    outT->Branch("RPnearcsi",RPnearcsi,"RPnearcsi[nRPtk]/F");
  }
  outT->SetDirectory(fOut);

  //READ TREE FROM FILE
  TFile *f = TFile::Open(filename);  
  TH1 *triggerList=(TH1 *)f->Get("analysis/triggerList");
  TTree *t = (TTree*)f->Get("analysis/data");
  attachToMiniEventTree(t,ev,true);
  Int_t nentries(t->GetEntriesFast());
  if (debug) nentries = 10000; //restrict number of entries for testing
  t->GetEntry(0);

  cout << "...producing " << outname << " from " << nentries << " events" << endl;
  
  //LUMINOSITY+PILEUP
  LumiTools lumi(era,genPU);
  std::map<Int_t,Float_t> lumiPerRun=lumi.lumiPerRun();

  //auxiliary to solve neutrino pZ using MET
  MEzCalculator neutrinoPzComputer;

  //LEPTON EFFICIENCIES
  EfficiencyScaleFactorsWrapper lepEffH(filename.Contains("Data13TeV"),era);

  //B-TAG CALIBRATION
  BTagSFUtil btvSF(era,BTagEntry::OperatingPoint::OP_MEDIUM,"",0);
  
   //BOOK HISTOGRAMS
  HistTool ht;
  ht.setNsyst(0);
  ht.addHist("puwgtctr",     new TH1F("puwgtctr",    ";Weight sums;Events",2,0,2));
  ht.addHist("nvtx",         new TH1F("nvtx",        ";Vertex multiplicity;Events",50,0,100));
  ht.addHist("nlep",         new TH1F("nlep",        ";Lepton multipliciy;Events",3,2,5));
  ht.addHist("nljets",       new TH1F("nljets",      ";light jet multiplicity;Events",6,0,6)); 
  ht.addHist("nbjets",       new TH1F("nbjets",      ";b jet multiplicity;Events",5,0,5));
  ht.addHist("lmpt",         new TH1F("lmpt",        ";Lepton 1 transverse momentum [GeV];Events",50,20,200));
  ht.addHist("lmeta",        new TH1F("lmeta",       ";Lepton 1 pseudo-rapidity;Events",10,0,2.5));
  ht.addHist("lppt",         new TH1F("lppt",        ";Lepton 2 transverse momentum [GeV];Events",50,20,200));
  ht.addHist("lpeta",        new TH1F("lpeta",       ";Lepton 2 pseudo-rapidity;Events",10,0,2.5));
  Float_t mllbins[]={0,20,50,60,70,75,85,90,95,100,105,115,125,200,300,500};
  ht.addHist("mll",          new TH1F("mll",         ";Dilepton invariant mass [GeV];Events",sizeof(mllbins)/sizeof(Float_t)-1,mllbins));
  ht.addHist("drll",         new TH1F("drll",        ";#DeltaR(l,l');Events",50,0,6));
  Float_t ptllbins[]={0,25,50,75,100,150,200,250,500,750,1000,2000};
  ht.addHist("ptll",         new TH1F("ptll",        ";Dilepton transverse momentum [GeV];Events",sizeof(ptllbins)/sizeof(Float_t)-1,ptllbins));
  ht.addHist("phistar",      new TH1F("phistar",     ";log_{10}(dilepton #phi^{*});Events",50,-3,3));
  ht.addHist("costhetaCS",   new TH1F("costhetaCS",  ";Dilepton cos#theta^{*}_{CS};Events",50,-1,1));
  ht.addHist("met",          new TH1F("met",         ";Missing transverse energy [GeV];Events",50,0,200));
  ht.addHist("mindphijmet",  new TH1F("mindphijmet", ";min#Delta#phi(jet,E_{T}^{miss}) [rad];Events",20,0,TMath::Pi()));
  ht.addHist("acopl",        new TH1F("acopl",       ";Acoplanarity;Events",50,0,1.0));
  ht.addHist("ntkrp",        new TH1F("ntkrp",       ";Track multiplicity; Events",6,0,6) );
  ht.addHist("csirp",        new TH1F("csirp",       ";#csi = #deltap/p; Events",50,0,0.3) );
  ht.addHist("xangle",       new TH1F("xangle",      ";Crossing angle; Events",10,100,200) );
  ht.addHist("evyields",     new TH1F("evyields",    ";Category; Events",6,0,6) );
  ht.getPlots()["evyields"]->GetXaxis()->SetBinLabel(1,"inc");
  ht.getPlots()["evyields"]->GetXaxis()->SetBinLabel(2,"inv");
  ht.getPlots()["evyields"]->GetXaxis()->SetBinLabel(3,"#geq3l");
  ht.getPlots()["evyields"]->GetXaxis()->SetBinLabel(4,"#gamma");
  ht.getPlots()["evyields"]->GetXaxis()->SetBinLabel(5,"jj");
  ht.getPlots()["evyields"]->GetXaxis()->SetBinLabel(6,"bb");
  ht.addHist("ratevsrun",    new TH1F("ratevsrun",   ";Run number; #sigma [pb]",int(lumiPerRun.size()),0,float(lumiPerRun.size())));


  int i=0;
  for(auto key : lumiPerRun) {
    i++;
    ht.getPlots()["ratevsrun"]->GetXaxis()->SetBinLabel(i,Form("%d",key.first));
  }

  std::cout << "init done" << std::endl;
  if (debug){std::cout<<"\n DEBUG MODE"<<std::endl;}

  ///////////////////////
  // LOOP OVER EVENTS //
  /////////////////////
  
  //EVENT SELECTION WRAPPER
  SelectionTool selector(filename, false, triggerList);
  
  for (Int_t iev=0;iev<nentries;iev++)
    {
      t->GetEntry(iev);
      if(iev%1000==0) { printf("\r [%3.0f%%] done", 100.*(float)iev/(float)nentries); fflush(stdout); }

      //particle level dilepton
      std::vector<double> trivialwgts(1,1.0);
      float gen_llpt(-1),gen_mll(-1),gen_drll(-1);
      TString gen_dilCat("");
      if(!ev.isData){
        std::vector<Particle> genLeptons=selector.getGenLeptons(ev,20.,2.5);
        
        if(genLeptons.size()>=2) {
          gen_llpt=(genLeptons[0]+genLeptons[1]).Pt();
          gen_mll=(genLeptons[0]+genLeptons[1]).M();
          gen_drll=genLeptons[0].DeltaR(genLeptons[1]);
          int gen_dilcode=abs(genLeptons[0].id()*genLeptons[1].id());
          if(gen_dilcode==11*11) gen_dilCat="genee";
          if(gen_dilcode==11*13) gen_dilCat="genem";
          if(gen_dilcode==13*13) gen_dilCat="genmm";
          ht.fill("ptll", gen_llpt, trivialwgts, gen_dilCat);
          ht.fill("mll",  gen_mll,  trivialwgts, gen_dilCat);
          ht.fill("drll", gen_drll,  trivialwgts, gen_dilCat);
        }
      }

      //start weights and pu weight control
      float wgt(1.0);
      std::vector<double>plotwgts(1,wgt);
      double puWgt(1.0);
      TString period = lumi.assignRunPeriod();
      if(!ev.isData){
        ht.fill("puwgtctr",0,plotwgts);
        puWgt=(lumi.pileupWeight(ev.g_pu,period)[0]);
        std::vector<double>puPlotWgts(1,puWgt);
        ht.fill("puwgtctr",1,puPlotWgts);
      }  
	


      //////////////////
      // CORRECTIONS  //
      //////////////////
      btvSF.addBTagDecisions(ev);
      if(!ev.isData) btvSF.updateBTagDecisions(ev);      

      ///////////////////////////
      // RECO LEVEL SELECTION  //
      ///////////////////////////

      //trigger
      hasETrigger=(selector.hasTriggerBit("HLT_Ele35_WPTight_Gsf_v", ev.triggerBits));
      hasMTrigger=(selector.hasTriggerBit("HLT_IsoMu24_v",     ev.triggerBits) ||
                   selector.hasTriggerBit("HLT_IsoMu24_2p1_v", ev.triggerBits) ||
                   selector.hasTriggerBit("HLT_IsoMu27_v",     ev.triggerBits) );
      hasMMTrigger=(selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",                  ev.triggerBits) ||
                    selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v",          ev.triggerBits) ||
                    selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v",        ev.triggerBits) );
      hasEETrigger=(selector.hasTriggerBit("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v",             ev.triggerBits) ||
                    selector.hasTriggerBit("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",          ev.triggerBits) );
      hasEMTrigger=(selector.hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",    ev.triggerBits) ||
                    selector.hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", ev.triggerBits) ||
                    selector.hasTriggerBit("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",    ev.triggerBits) ||
                    selector.hasTriggerBit("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v", ev.triggerBits) ||
                    selector.hasTriggerBit("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",     ev.triggerBits) ||
                    selector.hasTriggerBit("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",  ev.triggerBits) );
      if (ev.isData) { 
        //use only these unprescaled triggers for these eras
        if(filename.Contains("2017E") || filename.Contains("2017F")){
          hasMTrigger=selector.hasTriggerBit("HLT_IsoMu27_v",ev.triggerBits);
        }
        if(!(filename.Contains("2017A") || filename.Contains("2017B"))){
          hasMMTrigger=(selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v",   ev.triggerBits) ||
                        selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v", ev.triggerBits) );
        }
      }

      //trigger efficiency
      if( (gen_dilCat=="genee" && (hasEETrigger || hasETrigger)) ||
          (gen_dilCat=="genmm" && (hasMMTrigger || hasMTrigger)) ||
          (gen_dilCat=="genem" && (hasEMTrigger || hasETrigger || hasEMTrigger)) ) { 
        ht.fill("ptll", gen_llpt, trivialwgts, gen_dilCat+"trig");
        ht.fill("mll",  gen_mll,  trivialwgts, gen_dilCat+"trig");
        ht.fill("drll", gen_drll,  trivialwgts, gen_dilCat+"trig");
      }
      if( (gen_dilCat=="genee" && hasETrigger) ||
          (gen_dilCat=="genmm" && hasMTrigger) ||
          (gen_dilCat=="genem" && (hasETrigger || hasEMTrigger)) ) { 
        ht.fill("ptll", gen_llpt, trivialwgts, gen_dilCat+"singletrig");
        ht.fill("mll",  gen_mll,  trivialwgts, gen_dilCat+"singletrig");
        ht.fill("drll", gen_drll,  trivialwgts, gen_dilCat+"singletrig");
      }

      
      //identify the offline final state from the leading leptons
      int l1idx(0),l2idx(1);
      std::vector<Particle> flaggedLeptons = selector.flaggedLeptons(ev);
      std::vector<Particle> leptons  = selector.selLeptons(flaggedLeptons, SelectionTool::MEDIUM, SelectionTool::MVA90,20.,2.4);
      //require the leading lepton to be tight
      bool hasLeadingTightLepton(false), hasSubLeadingTightLepton(false);
      if(leptons.size()){
        if(leptons[0].Pt()>30 && fabs(leptons[0].Eta())<2.1){
          if( (leptons[0].id()==11  && leptons[0].hasQualityFlag(SelectionTool::MVA80) ) ||
              (leptons[0].id()==13 && leptons[0].hasQualityFlag(SelectionTool::TIGHT) ) ) hasLeadingTightLepton=true;
        }
        if(leptons[1].Pt()>30 && fabs(leptons[0].Eta())<2.1){
          if( (leptons[1].id()==11  && leptons[0].hasQualityFlag(SelectionTool::MVA80) ) ||
              (leptons[1].id()==13 && leptons[0].hasQualityFlag(SelectionTool::TIGHT) ) ) hasSubLeadingTightLepton=true;
        }
      }

      //photons
      std::vector<Particle> flaggedPhotons = selector.flaggedPhotons(ev);
      std::vector<Particle> photons        = selector.selPhotons(flaggedPhotons,SelectionTool::MVA80,leptons,30.,2.1);
      
      //pure hadronic sample (no tight leptons, or photons, jet-triggered)
      TString dilCat("");
      if(ev.isData && selector.isJetHTPD()) {
        
        //reset variables
        dilCat="jet";
        for(auto v : outVars) v.second=0.;
       
        //check trigger
        bool hasJetTrigger(false);
        Int_t thr[]={40,60,80,140,200,260,320,400,450,500,550};
        TString trigName("");
        for(size_t ithr=0; ithr<sizeof(thr)/sizeof(Int_t); ithr++) {
          
          TString incTrigName=Form("HLT_PFJet%d_v",thr[ithr]);
          TString fwdTrigName=Form("HLT_PFJetFwd%d_v",thr[ithr]);
          if(selector.hasTriggerBit(trigName, ev.triggerBits)){
            outVars["llpt"]=thr[ithr];            
            trigName=incTrigName;
          }
          else if(selector.hasTriggerBit(fwdTrigName, ev.triggerBits)){
            outVars["llpt"]=thr[ithr];
            outVars["lleta"]=3;  
            trigName=fwdTrigName;
          }
          else {
            continue;
          }
          hasJetTrigger=true;
          break;
        }
        if(!hasJetTrigger) continue;

        outVars["dilcode"]=0;
        if(leptons.size()>=2 && hasLeadingTightLepton){
          outVars["l1pt"]=leptons[0].Pt();
          outVars["l1eta"]=leptons[0].Eta();
          outVars["l1phi"]=leptons[0].Phi(); 
          outVars["ml1"]=leptons[0].M();
          outVars["l1id"]=leptons[0].id();
          outVars["l2pt"]=leptons[1].Pt();
          outVars["l2eta"]=leptons[1].Eta();
          outVars["l2phi"]=leptons[1].Phi();
          outVars["ml2"]=leptons[1].M();
          outVars["l2id"]=leptons[1].id();
        }
      }
      else {
            
        if(leptons.size()<2) continue;
        if(!hasLeadingTightLepton) continue;
        Int_t dilcode=leptons[l1idx].id()*leptons[l2idx].id();
        isSF=( leptons[l1idx].id()==leptons[l2idx].id() );
        isSS=( leptons[l1idx].charge()*leptons[l2idx].charge() > 0 );
        
        //check again origin of the dilepton in data to max. efficiency and avoid double counting
        if(ev.isData) {
          
          if(dilcode==11*11) {
            if( !selector.isDoubleEGPD()      && !selector.isSingleElectronPD()) continue;
            if( selector.isDoubleEGPD()       && !hasEETrigger ) continue;
            if( selector.isSingleElectronPD() && !(hasETrigger && !hasEETrigger) ) continue;
          }
          if(dilcode==13*13) {
            if( !selector.isDoubleMuonPD() && !selector.isSingleMuonPD()) continue;
            if( selector.isDoubleMuonPD()  && !hasMMTrigger ) continue;
            if( selector.isSingleMuonPD()  && !(hasMTrigger && !hasMMTrigger) ) continue;
          }
          if(dilcode==11*13) {
            if( !selector.isMuonEGPD()        && !selector.isSingleElectronPD() && !selector.isSingleMuonPD()) continue;
            if( selector.isMuonEGPD()         && !hasEMTrigger ) continue;
            if( selector.isSingleElectronPD() && !(hasETrigger && !hasEMTrigger) ) continue;
            if( selector.isSingleMuonPD()     && !(hasMTrigger && !hasETrigger && !hasMTrigger) ) continue;
          }

          //check trigger rates and final channel assignment
          std::map<Int_t,Float_t>::iterator rIt=lumiPerRun.find(ev.run);
          if(rIt!=lumiPerRun.end()){
            int runBin=std::distance(lumiPerRun.begin(),rIt);
            float lumi=1./rIt->second;
            if(dilcode==11*11) ht.fill("ratevsrun",runBin,lumi,"ee");
            if(dilcode==13*13) ht.fill("ratevsrun",runBin,lumi,"mm");
            if(dilcode==11*13) ht.fill("ratevsrun",runBin,lumi,"em");
          }else{
            cout << "[Warning] Unable to find run=" << ev.run << endl;
          }
        }

        TLorentzVector lm(leptons[l1idx].charge()>0 ? leptons[l1idx] : leptons[l2idx]);
        float lmScaleUnc(leptons[l1idx].charge()>0 ? leptons[l1idx].scaleUnc() : leptons[l2idx].scaleUnc());
        TLorentzVector lp(leptons[l1idx].charge()>0 ? leptons[l2idx] : leptons[l1idx]);
        float lpScaleUnc(leptons[l1idx].charge()>0 ? leptons[l2idx].scaleUnc() : leptons[l1idx].scaleUnc());
        if(isSS)  { 
          lm=leptons[l1idx]; 
          lmScaleUnc=leptons[l1idx].scaleUnc();
          lp=leptons[l2idx];
          lpScaleUnc=leptons[l2idx].scaleUnc();
        }
        TLorentzVector dil(lm+lp);
        Float_t dilScaleUnc=TMath::Sqrt( pow(lm.Pt()*lmScaleUnc,2)+pow(lp.Pt()*lpScaleUnc,2) )/dil.Pt();
        float mll(dil.M());
        if(mll<20) continue;
        isZ=( isSF && !isSS && fabs(mll-91)<10);
        
        //met
        TLorentzVector met(0,0,0,0);
        met.SetPtEtaPhiM(ev.met_pt[1],0,ev.met_phi[1],0.);
        Float_t metScaleUnc(1./ev.met_sig[1]);
        
        //jets (require PU jet id)
        std::vector<Jet> allJets = selector.getGoodJets(ev,25.,2.5,leptons,photons);
        int njets(0);
        std::vector<Jet> bJets,lightJets,jets;
        float scalarht(0.),scalarhtb(0.),scalarhtj(0.),mindphijmet(99999.);
        for(size_t ij=0; ij<allJets.size(); ij++) 
          {
            int idx=allJets[ij].getJetIndex();
            int jid=ev.j_id[idx];
            bool passLoosePu((jid>>2)&0x1);
            bool passBtag(ev.j_btag[ij]>0);
            if(!passLoosePu) continue;
            if(passBtag) { bJets.push_back(allJets[ij]);     scalarhtb+=allJets[ij].pt();  }
            else         { lightJets.push_back(allJets[ij]); scalarhtj+= allJets[ij].pt(); }
            njets++;
            jets.push_back(allJets[ij]);
            scalarht += jets[ij].pt();

            float dphij2met=fabs(allJets[ij].DeltaPhi(met));
            if(dphij2met>mindphijmet) continue;
            mindphijmet=dphij2met;
          }
        
        
        ///dilepton kinematics
        Float_t llht(lm.Pt()+lp.Pt());
        Float_t llacopl = computeAcoplanarity(lm,lp);
        Float_t llphistar=computePhiStar(lm,lp);
        Float_t llcosthetaCS=computeCosThetaStar(lm,lp);
        Float_t llMR(computeMR(lm,lp));
        Float_t llR(computeRsq(lm,lp,met));
        Float_t mtm(computeMT(lm,met));
        Float_t mtp(computeMT(lp,met));
        
        //baseline categories and additional stuff produced with the dilepton
        std::vector<TString> cats(1,"inc");
        dilCat=(isSS ? "ss" : "os" );
        dilCat+=!isSF ? "em" : dilcode==121 ? "ee" : "mm";        
        cats.push_back(dilCat);       
        
        //selection efficiency
        if(isZ) {
          ht.fill("ptll", gen_llpt, trivialwgts,gen_dilCat+"rec");
          ht.fill("mll",  gen_mll,  trivialwgts,gen_dilCat+"rec");
          ht.fill("drll", gen_drll,  trivialwgts, gen_dilCat+"rec");
          
          if(hasSubLeadingTightLepton){
            ht.fill("ptll", gen_llpt, trivialwgts,gen_dilCat+"tightrec");
            ht.fill("mll",  gen_mll,  trivialwgts,gen_dilCat+"tightrec");
            ht.fill("drll", gen_drll,  trivialwgts, gen_dilCat+"tightrec");
          }
        }
        
        TLorentzVector X(met);
        float xScaleUnc(metScaleUnc);
        float xEtaUnc(0);
        Int_t xid(0);
        Float_t xht(0);
        Float_t mt3l(-1);              
        if(leptons.size()>2){
          int l3idx(2);
          neutrinoPzComputer.SetMET(met);
          neutrinoPzComputer.SetLepton(leptons[l3idx].p4());
          float nupz=neutrinoPzComputer.Calculate();
          float metShiftUp(fabs(1+metScaleUnc));
          TLorentzVector upMET(metShiftUp*met);
          neutrinoPzComputer.SetMET(upMET);
          float nupzUp=neutrinoPzComputer.Calculate();
          float metShiftDn(fabs(1-metScaleUnc));
          TLorentzVector dnMET(metShiftDn*met);
          neutrinoPzComputer.SetMET(dnMET);        
          float nupzDn=neutrinoPzComputer.Calculate();
          
          TLorentzVector neutrinoP4(met.Px(),met.Py(),nupz ,TMath::Sqrt(TMath::Power(met.Pt(),2)+TMath::Power(nupz,2)));
          TLorentzVector neutrinoP4Up(upMET.Px(),upMET.Py(),nupzUp ,TMath::Sqrt(TMath::Power(upMET.Pt(),2)+TMath::Power(nupzUp,2)));
          TLorentzVector neutrinoP4Dn(dnMET.Px(),dnMET.Py(),nupzDn ,TMath::Sqrt(TMath::Power(dnMET.Pt(),2)+TMath::Power(nupzDn,2)));
          
          X=leptons[l3idx]+neutrinoP4;
          xScaleUnc=TMath::Sqrt(
                                pow(leptons[l3idx].scaleUnc()*leptons[l3idx].Pt(),2)+
                                pow(metScaleUnc*neutrinoP4.Pt(),2)
                                )/X.Pt();
          xEtaUnc=max( fabs((leptons[l3idx]+neutrinoP4Up).Eta()-X.Eta()),
                       fabs((leptons[l3idx]+neutrinoP4Dn).Eta()-X.Eta()));
          xht=X.Pt();
          xid=leptons[l3idx].id();
          mt3l=computeMT(leptons[l3idx],met);
        }
        else if(photons.size()>0) {
          X=photons[0].p4();
          xScaleUnc=photons[0].scaleUnc();
          xid=22;
          xht=X.Pt();
        }
        else if (jets.size()>2) {
          X=jets[0]+jets[1];
          xScaleUnc=TMath::Sqrt(pow(jets[0].getScaleUnc()*jets[0].Pt(),2)+
                                pow(jets[1].getScaleUnc()*jets[1].Pt(),2))/X.Pt();
          xid=2121;
          xht=jets[0].Pt()+jets[1].Pt();
          if(bJets.size()>1) {                      
            X=bJets[0]+bJets[1];
            xScaleUnc=TMath::Sqrt(pow(bJets[0].getScaleUnc()*bJets[0].Pt(),2)+
                                  pow(bJets[1].getScaleUnc()*bJets[1].Pt(),2))/X.Pt();
            xid=55;
            xht=bJets[0].Pt()+bJets[1].Pt();
          }
        }

        //final state F=ll+X
        TLorentzVector F(dil+X);
        Float_t fht(llht+xht);
        Float_t facopl(computeAcoplanarity(dil,X));
        Float_t fphiStar(computePhiStar(dil,X));
        Float_t fcosthetaStarCS(computeCosThetaStar(dil,X));
        Float_t fMR(computeMR(dil,X));
        Float_t fR(computeRsq(dil,X,met));
        Float_t mtll(computeMT(dil,met));
        Float_t mtx(computeMT(X,met));
        if(mt3l>0) mtx=mt3l; //for 3-leptons use the MT(3rd lep,MET) instead
        
        //recoil and UE
        int nch(int(ev.pf_ch_wgtSum));
        float closestTrackDZ(ev.pf_closestDZnonAssoc);
        TVector2 puppiRecoil(ev.pf_puppi_px,ev.pf_puppi_py);
        float puppiRecoilHt(ev.pf_puppi_ht);
        TVector2 h2( xid!=0 ? TVector2(F.Px(),F.Py()) : TVector2(dil.Px(),dil.Py()) );
        puppiRecoil   -= h2;
        puppiRecoilHt -= (xid!=0 ? fht : llht);
        

        ////////////////////
        // EVENT WEIGHTS //
        //////////////////
        if (!ev.isData) {
          
          // norm weight
          wgt  = (normH? normH->GetBinContent(1) : 1.0);
          
          // lepton trigger*selection weights
          EffCorrection_t trigSF(1.,0.);
          //EffCorrection_t trigSF = lepEffH.getTriggerCorrection(leptons,{},{},period);
          EffCorrection_t  sel1SF = lepEffH.getOfflineCorrection(leptons[l1idx], period);
          EffCorrection_t  sel2SF = lepEffH.getOfflineCorrection(leptons[l2idx], period);
          
          
          wgt *= puWgt*trigSF.first*sel1SF.first*sel2SF.first;
          
          // generator level weights
          wgt *= (ev.g_nw>0 ? ev.g_w[0] : 1.0);
          
          //update weight for plotter
          plotwgts[0]=wgt;
        }
        
        //control histograms
        ht.fill("nvtx",       ev.nvtx,         plotwgts, cats);
      
        //event yields
        ht.fill("evyields",  0,  plotwgts, cats);
        if(xid==0)         ht.fill("evyields",  1,  plotwgts, cats);
        else if(xid==22)   ht.fill("evyields",  3,  plotwgts, cats);
        else if(xid==2121) ht.fill("evyields",  4,  plotwgts, cats);
        else if(xid==55)   ht.fill("evyields",  5,  plotwgts, cats);
        else               ht.fill("evyields",  2,  plotwgts, cats);

        //dilepton system
        ht.fill("nlep",       leptons.size(),  plotwgts, cats);
        ht.fill("lmpt",       lm.Pt(),         plotwgts, cats);
        ht.fill("lmeta",      fabs(lm.Eta()),  plotwgts, cats);
        ht.fill("lppt",       lp.Pt(),         plotwgts, cats);
        ht.fill("lmeta",      fabs(lp.Eta()),  plotwgts, cats);
        ht.fill("mll",        dil.M(),         plotwgts, cats);        
        ht.fill("drll",       leptons[l1idx].DeltaR(leptons[l2idx]),plotwgts,cats);
        ht.fill("ptll",       dil.Pt(),        plotwgts, cats);
        ht.fill("phistar",    TMath::Log10(llphistar),  plotwgts, cats);
        ht.fill("costhetaCS", llcosthetaCS,      plotwgts, cats);
        ht.fill("acopl",      llacopl,           plotwgts, cats);
        
        //bjets
        ht.fill("nbjets",     bJets.size(),    plotwgts, cats);
        ht.fill("scalarhtb",     scalarhtb,    plotwgts, cats);
        for(auto &b:bJets) {
          ht.fill("bpt",     b.Pt(),           plotwgts, cats);
          ht.fill("beta",    fabs(b.Eta()),    plotwgts, cats);
        }
        
        //light jets
        ht.fill("nljets",     lightJets.size(),    plotwgts, cats);
        ht.fill("scalarhtj",  scalarhtj,    plotwgts, cats);
        for(auto &l:lightJets) {
          ht.fill("jpt",     l.Pt(),           plotwgts, cats);
          ht.fill("jeta",    fabs(l.Eta()),    plotwgts, cats);
        }
        
        //photons
        ht.fill("npho",       photons.size(),  plotwgts, cats);
        for(auto &a:photons) {
          ht.fill("apt",     a.Pt(),           plotwgts, cats);
          ht.fill("aeta",    fabs(a.Eta()),    plotwgts, cats);
        }
        
        ht.fill("met",           met.Pt(),    plotwgts, cats);
        ht.fill("mindphijmet",   mindphijmet, plotwgts, cats);
        ht.fill("nch",           nch,         plotwgts, cats);
            
        //fill tree with central detector information
        outVars["evwgt"]=plotwgts[0];
        outVars["dilcode"]=float(dilcode);

        outVars["l1pt"]=lm.Pt();
        outVars["l1eta"]=lm.Eta();
        outVars["l1phi"]=lm.Phi(); 
        outVars["ml1"]=lm.M();
        outVars["l1id"]=leptons[l1idx].id();
        outVars["mt1"]=mtm;

        outVars["l2pt"]=lp.Pt();
        outVars["l2eta"]=lp.Eta();
        outVars["l2phi"]=lp.Phi();
        outVars["ml2"]=lp.M();
        outVars["l2id"]=leptons[l2idx].id();
        outVars["mt2"]=mtp;
        
        outVars["llpt"]=dil.Pt();
        outVars["lleta"]=dil.Eta();
        outVars["llphi"]=dil.Phi();
        outVars["mll"]=dil.M();
        outVars["llacopl"]=llacopl;
        outVars["llcosthetaCS"]=llcosthetaCS;
        outVars["llphistar"]=llphistar;
        outVars["llMR"]=llMR;
        outVars["llR"]=llR;
        outVars["mtll"]=mtll;
        outVars["llht"]=llht;      

        ValueCollection_t llcsi=calcCsi(lm,lmScaleUnc,lp,lpScaleUnc);
        outVars["llcsip"]=llcsi[0].first;
        outVars["llcsipUnc"]=llcsi[0].second;
        outVars["llcsim"]=llcsi[1].first;
        outVars["llcsimUnc"]=llcsi[1].second;
        outVars["xpt"]=X.Pt();
        outVars["xeta"]=X.Eta();
        outVars["xphi"]=X.Phi(); 
        outVars["mx"]=X.M();
        outVars["xid"]=xid;
        outVars["xht"]=xht;
        outVars["mtx"]=mtx;
        
        outVars["fpt"]=F.Pt();
        outVars["feta"]=F.Eta();
        outVars["fphi"]=F.Phi();
        outVars["mf"]=F.M();
        outVars["fht"]=fht;      
        outVars["facopl"]=facopl;
        outVars["fcosthetaCS"]=fcosthetaStarCS;
        outVars["fphistar"]=fphiStar;
        outVars["fMR"]=fMR;
        outVars["fR"]=fR;
        ValueCollection_t fcsi=calcCsi(dil,dilScaleUnc,X,xScaleUnc,xEtaUnc);
        outVars["fcsip"]=fcsi[0].first;
        outVars["fcsipUnc"]=fcsi[0].second;
        outVars["fcsim"]=fcsi[1].first;
        outVars["fcsimUnc"]=fcsi[1].second;
        outVars["nb"]=bJets.size();
        outVars["nj"]=lightJets.size();
        outVars["nl"]=leptons.size();
        outVars["ng"]=photons.size();
        outVars["nch"]=nch;
        outVars["ht"]=scalarht;
        outVars["htb"]=scalarhtb;
        outVars["htj"]=scalarhtj;
        outVars["closestkdz"]=closestTrackDZ;
        outVars["puppirecoil_pt"]=puppiRecoil.Mod();
        outVars["puppirecoil_phi"]=puppiRecoil.Phi();
        outVars["puppirecoil_spher"]=fabs(puppiRecoil.Mod())/puppiRecoilHt;
      }
      
      //fill data with roman pot information
      nRPtk=0;
      if (ev.isData) {
        
        //reset information
        for(size_t irp=0; irp<50; irp++) { 
          RPid[irp]=0; 
          RPfarcsi[irp]=0; 
          RPnearcsi[irp]=0; 
        }
        
        try{
          const edm::EventID ev_id( ev.run, ev.lumi, ev.event );        
          const ctpps::conditions_t lhc_cond = lhc_conds.get( ev_id );
          beamXangle = lhc_cond.crossing_angle;
          ht.fill("beamXangle", beamXangle, plotwgts, dilCat);
          
          //only dispersions for these angles are available (@ CTPPSAnalysisTools/data/2017/dispersions.txt)
          if(beamXangle==120 || beamXangle==130 || beamXangle==140) {
            
            std::vector< std::pair<int,float> > nearCsis;
            std::map<int,int> ntks;
            ntks[23]=0; ntks[123]=0;
            for (int ift=0; ift<ev.nfwdtrk; ift++) {
            
              //only these roman pots are aligned CTPPSAnalysisTools/data/2017/alignments_21aug2018.txt 
              const unsigned short pot_raw_id = 100*ev.fwdtrk_arm[ift]+10*ev.fwdtrk_station[ift]+ev.fwdtrk_pot[ift];
              if (pot_raw_id!=3 && pot_raw_id!=23 && pot_raw_id!=103 && pot_raw_id!=123) continue;
            
              const ctpps::alignment_t align = ctpps_aligns.get( ev_id, pot_raw_id );
              double xi, xi_error;
              
              proton_reco.reconstruct(beamXangle, pot_raw_id, ev.fwdtrk_x[ift]*100+align.x_align, xi, xi_error);
              
              //information is only saved for the far detectors (pixels)
              if (pot_raw_id==23 || pot_raw_id==123) {
                RPid[nRPtk]=pot_raw_id;
                RPfarcsi[nRPtk]=xi;
                RPnearcsi[nRPtk]=0;              
                nRPtk++;
                
                //monitor track multiplicity and csi values
                if(ntks.find(pot_raw_id)==ntks.end()) ntks[pot_raw_id]=0;
                ntks[pot_raw_id]++;
                ht.fill("csirp",xi,plotwgts, Form("%s_%d",dilCat.Data(),pot_raw_id));
                
              }
              else{
                //save near detector info to match to pixel tracks
                nearCsis.push_back( std::pair<int,float>(pot_raw_id,xi) );
              }
            }
            
            //now try to find the best matches for strip in pixels
            for(auto stk : nearCsis) {
              
              int matchTk(-1);
              float minDcsi=1;
              for(int itk=0; itk<nRPtk; itk++) {
                
                //require on the same side of the beam pipe
                if( !( (RPid[itk]==123 && stk.first==103) || (RPid[itk]==23 && stk.first==3) ) )
                  continue;
              
                float dcsi=fabs(stk.second-RPfarcsi[itk]);
                if(dcsi>minDcsi) continue;
                matchTk=itk;
                minDcsi=dcsi;                
              }
              
              if(matchTk<0) continue;
              RPnearcsi[matchTk]=stk.second;
            }
          
          
            for(auto nit : ntks) 
              ht.fill("ntkrp", nit.second, plotwgts, Form("%s_%d",dilCat.Data(),nit.first));

          }
          
        }catch(...){
        }
      }
      
      outT->Fill();
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
  outT->Write();
  fOut->Close();
}
