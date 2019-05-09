#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>

//#include "TopLJets2015/TopAnalysis/interface/JSONWrapper.h"
#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"
#include "TopLJets2015/TopAnalysis/interface/ExclusiveZX.h"
#include "TopQuarkAnalysis/TopTools/interface/MEzCalculator.h"
#include "TopLJets2015/TopAnalysis/interface/L1PrefireEfficiencyWrapper.h"

#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>

#include "TMath.h"

#include "TopLJets2015/CTPPSAnalysisTools/interface/LHCConditionsFactory.h"


using namespace std;

#define ADDVAR(x,name,t,tree) tree->Branch(name,x,TString(name)+TString(t))

//
void RunExclusiveZX(const TString in_fname,
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
  const TString filename(in_fname);
  bool vetoPromptPhotons = filename.Contains("_QCDEM_") || filename.Contains("_TTJets");
  
  //RP in json
  /*
  std::string RPoutFile(Form("%s/src/TopLJets2015/TopAnalysis/test/analysis/pps/golden_noRP.json", CMSSW_BASE));
  JSONWrapper::Object json(RPoutFile,true);
  for(auto k : json.key){
    cout << k << endl;
    std::vector<JSONWrapper::Object> lumis=json.getObject(k).daughters();
    for(auto ll : lumis) {
      for(auto kk : ll["obj"].obj)
        cout << "\t" << kk.val << endl;
    }
  }
  */

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
  Int_t beamXangle;
  outT->Branch("beamXangle",&beamXangle,"beamXangle/I");
  outT->Branch("nvtx",&ev.nvtx,"nvtx/I");
  outT->Branch("nchPV",&ev.nchPV,"nchPV/I");
  outT->Branch("sumPVChPt",&ev.sumPVChPt,"sumPVChPt/F");
  outT->Branch("sumPVChPz",&ev.sumPVChPz,"sumPVChPz/F");
  outT->Branch("sumPVChHt",&ev.sumPVChHt,"sumPVChHt/F");
  outT->Branch("metfilters",&ev.met_filterBits,"metfilters/I");  

  //variables of doom
  outT->Branch("nrawmu", &ev.nrawmu, "nrawmu/I");
  outT->Branch("rawmu_pt", ev.rawmu_pt, "rawmu_pt[nrawmu]/S");
  outT->Branch("rawmu_eta", ev.rawmu_eta, "rawmu_eta[nrawmu]/S");
  outT->Branch("rawmu_phi", ev.rawmu_phi, "rawmu_phi[nrawmu]/S");
  outT->Branch("rawmu_pid", ev.rawmu_pid, "rawmu_pid[nrawmu]/I");

  bool hasMTrigger,hasMMTrigger,hasEETrigger,hasEMTrigger,hasZBTrigger,hasATrigger;
  outT->Branch("hasMTrigger",  &hasMTrigger," hasMTrigger/O");
  outT->Branch("hasEMTrigger", &hasEMTrigger, "hasEMTrigger/O");
  outT->Branch("hasMMTrigger", &hasMMTrigger, "hasMMTrigger/O");
  outT->Branch("hasEETrigger", &hasEETrigger, "hasEETrigger/O");
  outT->Branch("hasATrigger",  &hasATrigger,  "hasATrigger/O");
  outT->Branch("hasZBTrigger", &hasZBTrigger, "hasZBTrigger/O");

  bool isSS,isSF,isZ,isA;
  outT->Branch("isSS",&isSS,"isSS/O");
  outT->Branch("isSF",&isSF,"isSF/O");
  outT->Branch("isZ",&isZ,"isZ/O");
  outT->Branch("isA",&isA,"isA/O");
  
  ADDVAR(&ev.rho,"rho","F",outT);
  ADDVAR(&ev.met_pt,"met_pt","F",outT);
  ADDVAR(&ev.met_phi,"met_phi","F",outT);
  ADDVAR(&ev.met_sig,"met_sig","F",outT);
  TString fvars[]={"evwgt", "evcat", 
                   "l1pt", "l1eta", "l1phi", "ml1", "l1id", 
                   "l2pt", "l2eta", "l2phi", "ml2", "l2id", 
                   "bosonpt","bosoneta", "bosonphi", "mboson", 
                   "llacopl", "llcosthetaCS", "llphistar", "llMR", "llR", 
                   "llcsip", "llcsim",
                   "j1pt","j1eta","j1phi","j1m",
                   "j2pt","j2eta","j2phi","j2m",
                   "j3pt","j3eta","j3phi","j3m",
                   "nb", "nj", "htb","htj",
                   "lumiDeliv","lumiReco",
                   "PFMultSumEB",    "PFMultSumEE",    "PFMultSumHE",    "PFMultSumHF", 
                   "PFMultDiffEB",   "PFMultDiffEE",   "PFMultDiffHE",   "PFMultDiffHF", 
                   "PFChMultSumEB",  "PFChMultSumEE",  "PFChMultSumHE",  "PFChMultSumHF", 
                   "PFChMultDiffEB", "PFChMultDiffEE", "PFChMultDiffHE", "PFChMultDiffHF", 
                   "PFPzSumEB",      "PFPzSumEE",      "PFPzSumHE",      "PFPzSumHF", 
                   "PFPzDiffEB",     "PFPzDiffEE",     "PFPzDiffHE",     "PFPzDiffHF", 
                   "PFChPzSumEB",    "PFChPzSumEE",    "PFChPzSumHE",    "PFChPzSumHF", 
                   "PFChPzDiffEB",   "PFChPzDiffEE",   "PFChPzDiffHE",   "PFChPzDiffHF", 
                   "PFHtSumEB",      "PFHtSumEE",      "PFHtSumHE",      "PFHtSumHF", 
                   "PFHtDiffEB",     "PFHtDiffEE",     "PFHtDiffHE",     "PFHtDiffHF", 
                   "PFChHtSumEB",    "PFChHtSumEE",    "PFChHtSumHE",    "PFChHtSumHF", 
                   "PFChHtDiffEB",   "PFChHtDiffEE",   "PFChHtDiffHE",   "PFChHtDiffHF",
                   "trainCat"
  };

  std::map<TString,Float_t> outVars;
  for(size_t i=0; i<sizeof(fvars)/sizeof(TString); i++){
    outVars[fvars[i]]=0.;
    ADDVAR(&(outVars[fvars[i]]),fvars[i],"F",outT);
  }
  int nRPtk(0),RPid[50];
  float RPfarcsi[50],RPnearcsi[50];
  if(filename.Contains("Data13TeV")){
    outT->Branch("nRPtk",&nRPtk,"nRPtk/i");
    outT->Branch("RPid",RPid,"RPid[nRPtk]/i");
    outT->Branch("RPfarcsi",RPfarcsi,"RPfarcsi[nRPtk]/F");
    outT->Branch("RPnearcsi",RPnearcsi,"RPnearcsi[nRPtk]/F");
  }
  outT->SetDirectory(fOut);
  
  //LUMINOSITY+PILEUP
  LumiTools lumi(era,genPU);
  std::map<Int_t,Float_t> lumiPerRun=lumi.lumiPerRun();

  //LEPTON EFFICIENCIES
  std::map<TString,TString> cfgMap;
  cfgMap["g_id"]="MVAwp90";
  cfgMap["m_id"]="TightID";
  cfgMap["m_iso"]="TightRelIso";
  cfgMap["m_id4iso"]="TightIDandIPCut";
  cfgMap["e_id"]="MVA90";
  EfficiencyScaleFactorsWrapper lepEffH(filename.Contains("Data13TeV"),era,cfgMap);

  //L1-prefire 
  L1PrefireEfficiencyWrapper l1PrefireWR(filename.Contains("Data13TeV"),era);

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
  Float_t ptbosonbins[]={0,25,50,75,100,150,200,250,500,750,1000,2000};
  ht.addHist("ptboson",      new TH1F("ptboson",        ";Transverse momentum [GeV];Events",sizeof(ptbosonbins)/sizeof(Float_t)-1,ptbosonbins));
  ht.addHist("yboson",       new TH1F("yboson",      ";Rapidity;Events",50,-3,3));
  ht.addHist("phistar",      new TH1F("phistar",     ";Dilepton #phi^{*};Events",50,0,5000));
  ht.addHist("costhetaCS",   new TH1F("costhetaCS",  ";Dilepton cos#theta^{*}_{CS};Events",50,-1,1));
  ht.addHist("met",          new TH1F("met",         ";Missing transverse energy [GeV];Events",50,0,200));
  ht.addHist("mindphijmet",  new TH1F("mindphijmet", ";min#Delta#phi(jet,E_{T}^{miss}) [rad];Events",20,0,TMath::Pi()));
  ht.addHist("acopl",        new TH1F("acopl",       ";Acoplanarity;Events",50,0,1.0));
  ht.addHist("ntkrp",        new TH1F("ntkrp",       ";Track multiplicity; Events",6,0,6) );
  ht.addHist("csirp",        new TH1F("csirp",       ";#xi = #deltap/p; Events",50,0,0.3) );
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
  ofstream debug_out;
  if (debug){
    std::cout<<"\n DEBUG MODE"<<std::endl;
    debug_out.open("selz_info.txt");
    debug_out << "#run lumi event m_trig mm_trig zpt zeta zphi zm nmuons (RPid RPcsi)" << endl;
    debug_out << endl;
  }

  ///////////////////////
  // LOOP OVER EVENTS //
  /////////////////////

  //READ TREE FROM FILE
  TFile *f = TFile::Open(filename);  
  TH1 *triggerList=(TH1 *)f->Get("analysis/triggerList");
  TTree *t = (TTree*)f->Get("analysis/data");
  attachToMiniEventTree(t,ev,true);
  Int_t nentries(t->GetEntriesFast());
  if (debug) nentries = min(100000,nentries); //restrict number of entries for testing
  t->GetEntry(0);
  cout << "...producing " << outname << " from " << nentries << " events" << endl;

  
  //EVENT SELECTION WRAPPER
  SelectionTool selector(filename, false, triggerList);
  
  for (Int_t iev=0;iev<nentries;iev++)
    {
      t->GetEntry(iev);
      if(iev%1000==0) { printf("\r [%3.0f%%] done", 100.*(float)iev/(float)nentries); fflush(stdout); }

      //particle level dilepton or photon
      std::vector<double> trivialwgts(1,1.0);
      float gen_pt(-1),gen_m(-1),gen_dr(-1);
      std::vector<TString> gen_cats;
      if(!ev.isData){
        std::vector<Particle> genLeptons=selector.getGenLeptons(ev,20.,2.4);
        std::vector<Particle> genPhotons=selector.getGenPhotons(ev,50.,1.442);
        
        if(genLeptons.size()>=2 && genLeptons[0].Pt()>30) { 
          gen_pt=(genLeptons[0]+genLeptons[1]).Pt();          
          gen_m=(genLeptons[0]+genLeptons[1]).M();
          bool isZ(fabs(gen_m-91)<10);
          gen_dr=genLeptons[0].DeltaR(genLeptons[1]);
          int gen_dilcode=abs(genLeptons[0].id()*genLeptons[1].id());
          if(gen_dilcode==11*11) { gen_cats.push_back("genee"); if(isZ) gen_cats.push_back("geneez"); };
          if(gen_dilcode==11*13) { gen_cats.push_back("genem"); }
          if(gen_dilcode==13*13) { gen_cats.push_back("genmm"); if(isZ) gen_cats.push_back("genmmz"); };
        }
        else if(genPhotons.size()>=1){
          gen_pt=genPhotons[0].Pt();
          gen_m=0;
          gen_dr=0;
          gen_cats.push_back("gena");
        }
        
        if(gen_cats.size()>0){
          ht.fill("ptboson", gen_pt, trivialwgts, gen_cats);
          ht.fill("mll",  gen_m,  trivialwgts, gen_cats);
          ht.fill("drll", gen_dr,  trivialwgts, gen_cats);
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
      hasZBTrigger=((ev.addTriggerBits>>20)&0x1);
      hasATrigger=selector.hasTriggerBit("HLT_Photon90_R9Id90_HE10_IsoM_v", ev.triggerBits);
      hasMTrigger=(selector.hasTriggerBit("HLT_IsoMu27_v",     ev.triggerBits) );     
      hasMMTrigger=selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v",        ev.triggerBits);
      hasEETrigger=selector.hasTriggerBit("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v",             ev.triggerBits);
      hasEMTrigger=(selector.hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",    ev.triggerBits) ||
                    selector.hasTriggerBit("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", ev.triggerBits) ||
                    selector.hasTriggerBit("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",  ev.triggerBits) );
      
  
      //trigger efficiency
      for(auto gen_cat : gen_cats) {
        if( (gen_cat.Contains("gena") && hasATrigger) || 
            (gen_cat.Contains("genee") && hasEETrigger) ||
            (gen_cat.Contains("genmm") && (hasMMTrigger || hasMTrigger)) ||
            (gen_cat.Contains("genem") && (hasEMTrigger || hasMTrigger)) ) { 
          ht.fill("ptboson", gen_pt, trivialwgts, gen_cat+"trig");
          ht.fill("mll",  gen_m,  trivialwgts, gen_cat+"trig");
          ht.fill("drll", gen_dr, trivialwgts, gen_cat+"trig");
        }
        if( (gen_cat.Contains("gena") && hasATrigger) ||            
            (gen_cat.Contains("genmm") && hasMTrigger) ||
            (gen_cat.Contains("genem") && hasMTrigger) ) { 
          ht.fill("ptboson", gen_pt, trivialwgts, gen_cat+"singletrig");
          ht.fill("mll",  gen_m,  trivialwgts, gen_cat+"singletrig");
          ht.fill("drll", gen_dr, trivialwgts, gen_cat+"singletrig");
        }
      }
      
      //identify the offline final state from the leading leptons or the photon
      int l1idx(0),l2idx(1);
      std::vector<Particle> leptons = selector.flaggedLeptons(ev);     
      leptons = selector.selLeptons(leptons,SelectionTool::LOOSE,SelectionTool::MVA90,20,2.5);
      std::vector<Particle> allPhotons=selector.flaggedPhotons(ev);
      allPhotons=selector.selPhotons(allPhotons,SelectionTool::MVA90,{},50,3);
      std::vector<Particle> photons=selector.selPhotons(allPhotons,SelectionTool::MVA90,leptons,95,1.4442);

      //jets
      std::vector<Jet> allJets = selector.getGoodJets(ev,30.,4.7,leptons,photons);

      //met
      TLorentzVector met(0,0,0,0);
      met.SetPtEtaPhiM(ev.met_pt,0,ev.met_phi,0.);
      
      //
      //OFFLINE SELECTION(S)
      //
      bool passTightSel(false),passMediumSel(false);
      if(leptons.size()>1){
        
        bool isTrigSafe(leptons[0].Pt()>30);

        bool isLeadingTight( (leptons[0].id()==11 && leptons[0].hasQualityFlag(SelectionTool::MVA80)) ||
                             (leptons[0].id()==13 && leptons[0].hasQualityFlag(SelectionTool::TIGHT)) );
        bool isSubLeadingTight( (leptons[1].id()==11 && leptons[1].hasQualityFlag(SelectionTool::MVA80)) ||
                                (leptons[1].id()==13 && leptons[1].hasQualityFlag(SelectionTool::TIGHT)) );
        passTightSel = (isTrigSafe && isLeadingTight && isSubLeadingTight);
        
        bool isLeadingMedium( (leptons[0].id()==11 && leptons[0].hasQualityFlag(SelectionTool::MVA90)) ||
                              (leptons[0].id()==13 && leptons[0].hasQualityFlag(SelectionTool::TIGHT)) );        
        bool isSubLeadingMedium( (leptons[1].id()==11 && leptons[1].hasQualityFlag(SelectionTool::MVA90)) ||
                                 (leptons[1].id()==13 && leptons[1].hasQualityFlag(SelectionTool::TIGHT)) );

        passMediumSel = (isTrigSafe && isLeadingMedium && isSubLeadingMedium);        
      }
      else if(photons.size()>0){
        passTightSel=true;
        passMediumSel=true;
      }

      //apply selection
      TString selCat("");
      int selCode(0);
      if(!selector.isZeroBiasPD()){
        if(leptons.size()<2 && photons.size()==0) continue;
        if(!passMediumSel) continue;
      }

      //set kinematics
      TLorentzVector boson(0,0,0,0);
      TLorentzVector lm(0,0,0,0),lp(0,0,0,0);
      float mass(0);
      ValueCollection_t llcsi;
      float llacopl(0),llphistar(0),llcosthetaCS(0),llMR(0),llR(0),drll(0);
      if(leptons.size()>1){
        selCode=leptons[l1idx].id()*leptons[l2idx].id();
        isSF=( leptons[l1idx].id()==leptons[l2idx].id() );
        isSS=( leptons[l1idx].charge()*leptons[l2idx].charge() > 0 );
        if(selCode==11*11) selCat="ee";
        if(selCode==11*13) selCat="em";
        if(selCode==13*13) selCat="mm";

        lm=leptons[l1idx];
        lp=leptons[l2idx];
        boson=lm+lp;

        //dilepton specific
        mass=boson.M();
        llcsi=calcCsi(lm,lp);
        llacopl = computeAcoplanarity(lm,lp);
        drll=leptons[l1idx].DeltaR(leptons[l2idx]);
        llphistar=computePhiStar(lm,lp);
        llcosthetaCS=computeCosThetaStar(lm,lp);
        llMR=computeMR(lm,lp);
        llR=computeRsq(lm,lp,met);       
      }else if(photons.size()>0) {
        selCode=22;
        isSF=false;
        isSS=false;
        selCat="a";        
        boson=photons[0];
        
        //remove double counting of prompt photons in other samples
        if(vetoPromptPhotons)
          if(ev.gamma_isPromptFinalState[ photons[0].originalReference() ] ) continue;
      }
    
      //further selection for dileptons
      if(!selector.isZeroBiasPD()) {
        if(selCode!=22 && mass<20) continue;
        isZ=( isSF && !isSS && fabs(mass-91)<10);
        isA=(selCode==22);
      }
      else {
        if(!hasZBTrigger) continue;
        selCat="zbias";
      }

      //check again origin of the boson in data to max. efficiency and avoid double counting
      if(ev.isData) {
        if(isA) {
          if( !selector.isPhotonPD() ) continue;
          if( !hasATrigger) continue;
        }
        if(selCode==11*11) {
          if( !selector.isDoubleEGPD()      && !selector.isSingleElectronPD()) continue;
          if( selector.isDoubleEGPD()       && !hasEETrigger ) continue;
          if( selector.isSingleElectronPD()) continue;
        }
        if(selCode==13*13) {
          if( !selector.isDoubleMuonPD() && !selector.isSingleMuonPD()) continue;
          if( selector.isDoubleMuonPD()  && !hasMMTrigger ) continue;
          if( selector.isSingleMuonPD()  && !(hasMTrigger && !hasMMTrigger) ) continue;
        }
        if(selCode==11*13) {
          if( !selector.isMuonEGPD()        && !selector.isSingleElectronPD() && !selector.isSingleMuonPD()) continue;
          if( selector.isMuonEGPD()         && !hasEMTrigger ) continue;
          if( selector.isSingleElectronPD()) continue; 
          if( selector.isSingleMuonPD()     && !(hasMTrigger && !hasEMTrigger) ) continue;
        }
        
        //check trigger rates and final channel assignment
        std::map<Int_t,Float_t>::iterator rIt=lumiPerRun.find(ev.run);
        if(rIt!=lumiPerRun.end()){
          int runBin=std::distance(lumiPerRun.begin(),rIt);
          float lumi=1./rIt->second;
          ht.fill("ratevsrun",runBin,lumi,selCat);
        }else{
          cout << "[Warning] Unable to find run=" << ev.run << endl;
        }
      }

      //jets (require PU jet id)
      int njets(0);
      std::vector<Jet> bJets,lightJets,jets;
      float scalarht(0.),scalarhtb(0.),scalarhtj(0.),mindphijmet(99999.);     
      for(size_t ij=0; ij<allJets.size(); ij++) 
        {
          int idx=allJets[ij].getJetIndex();
          bool passBtag(ev.j_btag[idx]>0);

          int jid=ev.j_id[idx];
          bool passLoosePu((jid>>2)&0x1);          
          if(!passLoosePu) continue;
          jets.push_back(allJets[ij]);
          njets++;

          scalarht += jets[ij].pt();          
          if(passBtag) { bJets.push_back(allJets[ij]);     scalarhtb+=allJets[ij].pt();  }
          else         { lightJets.push_back(allJets[ij]); scalarhtj+= allJets[ij].pt(); }

          float dphij2met=fabs(allJets[ij].DeltaPhi(met));
          if(dphij2met>mindphijmet) continue;
          mindphijmet=dphij2met;
        }
       
      ////////////////////
      // EVENT WEIGHTS //
      //////////////////
      if (!ev.isData) {
        
        // norm weight
        wgt  = (normH? normH->GetBinContent(1) : 1.0);
        
        // lepton trigger*selection weights
        EffCorrection_t trigSF(1.,0.); // = lepEffH.getTriggerCorrection(leptons,{},{},period);
        EffCorrection_t sel1SF(1.,0.),sel2SF(1.,0.);
        if(selCode!=22){
          sel1SF = lepEffH.getOfflineCorrection(leptons[l1idx], period);
          sel2SF = lepEffH.getOfflineCorrection(leptons[l2idx], period);
        }else{
          sel1SF = lepEffH.getOfflineCorrection(photons[0], period);
        }

        //L1 pre-fire
        EffCorrection_t l1prefireProb=l1PrefireWR.getCorrection(allJets,allPhotons);

        wgt *= puWgt*trigSF.first*sel1SF.first*sel2SF.first*l1prefireProb.first;
        
        // generator level weights
        wgt *= (ev.g_nw>0 ? ev.g_w[0] : 1.0);
        
        //update weight for plotter
        plotwgts[0]=wgt;
      }
      
      //baseline categories and additional stuff produced with the dilepton
      std::vector<TString> cats(1,selCat);
        
      //selection efficiency                   
      for(auto gen_cat : gen_cats) {
        ht.fill("mll",  gen_m,  trivialwgts,gen_cat+"rec");
        if(passTightSel) ht.fill("mll",  gen_m,  trivialwgts,gen_cat+"2trec");
        if(isZ || isA) {
          ht.fill("ptboson", gen_pt, trivialwgts,gen_cat+"rec");
          ht.fill("drll", gen_dr,  trivialwgts, gen_cat+"rec");
          if(passTightSel) {
            ht.fill("ptboson", gen_pt, trivialwgts,gen_cat+"2trec");
            ht.fill("drll", gen_dr,  trivialwgts, gen_cat+"2trec");
          }
        }
      }

      //from this point onward require offline object and trigger
      bool hasOffAndTrig( (selector.isZeroBiasPD() && hasZBTrigger) ||
                          (selCat=="a"             && hasATrigger)  || 
                          (selCat=="ee"            && hasEETrigger) ||
                          (selCat=="mm"            && (hasMMTrigger || hasMTrigger)) ||
                          (selCat=="em"            && (hasEMTrigger || hasMTrigger)) );
      if(!hasOffAndTrig) continue;


      //control histograms
      ht.fill("nvtx",       ev.nvtx,         plotwgts, cats);
      
      //event yields
      ht.fill("evyields",  0,  plotwgts, cats);

      //dilepton system
      ht.fill("nlep",       leptons.size(),  plotwgts, cats);
      if(!isA){
        ht.fill("lmpt",       lm.Pt(),         plotwgts, cats);
        ht.fill("lmeta",      fabs(lm.Eta()),  plotwgts, cats);
        ht.fill("lppt",       lp.Pt(),         plotwgts, cats);
        ht.fill("lpeta",      fabs(lp.Eta()),  plotwgts, cats);
        ht.fill("mll",        boson.M(),       plotwgts, cats);        
        ht.fill("drll",       drll,            plotwgts,cats);
        ht.fill("phistar",    llphistar,       plotwgts, cats);
        ht.fill("costhetaCS", llcosthetaCS,    plotwgts, cats);
        ht.fill("acopl",      llacopl,         plotwgts, cats);
      }
      ht.fill("ptboson",       boson.Pt(),     plotwgts, cats);
      ht.fill("yboson",       boson.Rapidity(),     plotwgts, cats);
      
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
      
      //fill tree with central detector information
      outVars["evwgt"]=plotwgts[0];
      if(selector.isZeroBiasPD()) outVars["evwgt"]=float(ev.zeroBiasPS);
      outVars["evcat"]=float(selCode);

      outVars["l1pt"]=lm.Pt();
      outVars["l1eta"]=lm.Eta();
      outVars["l1phi"]=lm.Phi(); 
      outVars["ml1"]=lm.M();
      outVars["l1id"]=leptons.size()>1 ? leptons[l1idx].id() : 0.;
      
      outVars["l2pt"]=lp.Pt();
      outVars["l2eta"]=lp.Eta();
      outVars["l2phi"]=lp.Phi();
      outVars["ml2"]=lp.M();
      outVars["l2id"]=leptons.size()>1 ? leptons[l2idx].id() : 0.;
      
      outVars["bosonpt"]=boson.Pt();
      outVars["bosoneta"]=boson.Eta();
      outVars["bosonphi"]=boson.Phi();
      outVars["mboson"]=boson.M();

      //dilepton specific
      outVars["llacopl"]=llacopl;
      outVars["llcosthetaCS"]=llcosthetaCS;
      outVars["llphistar"]=llphistar;
      outVars["llMR"]=llMR;
      outVars["llR"]=llR;      
      outVars["llcsip"]    = llcsi.size()>0 ? llcsi[0].first  : 0;
      outVars["llcsim"]    = llcsi.size()>1 ? llcsi[1].first  : 1;
      
      outVars["nb"]=bJets.size();
      outVars["nj"]=lightJets.size();
      outVars["nl"]=leptons.size();
      outVars["ng"]=photons.size();

      outVars["ht"]=scalarht;
      outVars["htb"]=scalarhtb;
      outVars["htj"]=scalarhtj;
      outVars["j1pt"]=jets.size()>0 ? jets[0].Pt() : 0.;
      outVars["j1eta"]=jets.size()>0 ? jets[0].Eta() : 0.;
      outVars["j1phi"]=jets.size()>0 ? jets[0].Phi() : 0.;
      outVars["j1m"]=jets.size()>0 ? jets[0].M() : 0.;
      outVars["j2pt"]=jets.size()>1 ? jets[1].Pt() : 0.;
      outVars["j2eta"]=jets.size()>1 ? jets[1].Eta() : 0.;
      outVars["j2phi"]=jets.size()>1 ? jets[1].Phi() : 0.;
      outVars["j2m"]=jets.size()>1 ? jets[1].M() : 0.;
      outVars["j3pt"]=jets.size()>2 ? jets[2].Pt() : 0.;
      outVars["j3eta"]=jets.size()>2 ? jets[2].Eta() : 0.;
      outVars["j3phi"]=jets.size()>2 ? jets[2].Phi() : 0.;
      outVars["j3m"]=jets.size()>2 ? jets[2].M() : 0.;
      
      //flux variables
      if(!isA) {
        ev.nchPV -=2;
        ev.sumPVChHt=max(ev.sumPVChPt-lp.Pt() - lm.Pt(),0.);
        ev.sumPVChPz=max(ev.sumPVChPz- fabs(lp.Pz()) - fabs(lm.Pz()),0.);
        for(size_t i=0; i<2; i++){ 
          TLorentzVector lp4(i==0 ? lm : lp);
          size_t etaidx(0);
          if(lp4.Eta()>-3)   etaidx=1;
          if(lp4.Eta()>-2.5) etaidx=2;
          if(lp4.Eta()>-1.5) etaidx=3;
          if(lp4.Eta()>0)    etaidx=4;
          if(lp4.Eta()>1.5)  etaidx=5;
          if(lp4.Eta()>2.5)  etaidx=6;
          if(lp4.Eta()>3.0)  etaidx=7;
          ev.nPFCands[etaidx]--;
          ev.sumPFHt[etaidx]=max(ev.sumPFHt[etaidx]-lp4.Pt(),0.);
          ev.sumPFEn[etaidx]=max(ev.sumPFEn[etaidx]-lp4.E(),0.);
          ev.sumPFPz[etaidx]=max(ev.sumPFPz[etaidx]-fabs(lp4.Pz()),0.);
          ev.nPFChCands[etaidx]--;
          ev.sumPFChHt[etaidx]=max(ev.sumPFChHt[etaidx]-lp4.Pt(),0.);
          ev.sumPFChEn[etaidx]=max(ev.sumPFChEn[etaidx]-lp4.E(),0.);
          ev.sumPFChPz[etaidx]=max(ev.sumPFChPz[etaidx]-fabs(lp4.Pz()),0.);
        }        
      }else if(isA) {
        TLorentzVector ap4(photons[0]);
        size_t etaidx(0);
        if(ap4.Eta()>-3)   etaidx=1;
        if(ap4.Eta()>-2.5) etaidx=2;
        if(ap4.Eta()>-1.5) etaidx=3;
        if(ap4.Eta()>0)    etaidx=4;
        if(ap4.Eta()>1.5)  etaidx=5;
        if(ap4.Eta()>2.5)  etaidx=6;
        if(ap4.Eta()>3.0)  etaidx=7;
        ev.nPFCands[etaidx]--;
        ev.sumPFHt[etaidx]=max(ev.sumPFHt[etaidx]-ap4.Pt(),0.);
        ev.sumPFEn[etaidx]=max(ev.sumPFEn[etaidx]-ap4.E(),0.);
        ev.sumPFPz[etaidx]=max(ev.sumPFPz[etaidx]-fabs(ap4.Pz()),0.);
      }

      outVars["PFMultSumEB"]    = ev.nPFCands[3]+ev.nPFCands[4];
      outVars["PFMultSumEE"]    = ev.nPFCands[2]+ev.nPFCands[5];   
      outVars["PFMultSumHE"]    = ev.nPFCands[1]+ev.nPFCands[6];    
      outVars["PFMultSumHF"]    = ev.nPFCands[0]+ev.nPFCands[7]; 
      outVars["PFMultDiffEB"]   = fabs(ev.nPFCands[3]-ev.nPFCands[4]);
      outVars["PFMultDiffEE"]   = fabs(ev.nPFCands[2]-ev.nPFCands[5]);   
      outVars["PFMultDiffHE"]   = fabs(ev.nPFCands[1]-ev.nPFCands[6]);    
      outVars["PFMultDiffHF"]   = fabs(ev.nPFCands[0]-ev.nPFCands[7]); 
      outVars["PFChMultSumEB"]  = ev.nPFChCands[3]+ev.nPFChCands[4];
      outVars["PFChMultSumEE"]  = ev.nPFChCands[2]+ev.nPFChCands[5];   
      outVars["PFChMultSumHE"]  = ev.nPFChCands[1]+ev.nPFChCands[6];    
      outVars["PFChMultSumHF"]  = ev.nPFChCands[0]+ev.nPFChCands[7]; 
      outVars["PFChMultDiffEB"] = fabs(ev.nPFChCands[3]-ev.nPFChCands[4]);
      outVars["PFChMultDiffEE"] = fabs(ev.nPFChCands[2]-ev.nPFChCands[5]);   
      outVars["PFChMultDiffHE"] = fabs(ev.nPFChCands[1]-ev.nPFChCands[6]);    
      outVars["PFChMultDiffHF"] = fabs(ev.nPFChCands[0]-ev.nPFChCands[7]); 
      outVars["PFPzSumEB"]      = ev.sumPFPz[3]+ev.sumPFPz[4];
      outVars["PFPzSumEE"]      = ev.sumPFPz[2]+ev.sumPFPz[5];   
      outVars["PFPzSumHE"]      = ev.sumPFPz[1]+ev.sumPFPz[6];    
      outVars["PFPzSumHF"]      = ev.sumPFPz[0]+ev.sumPFPz[7]; 
      outVars["PFPzDiffEB"]     = fabs(ev.sumPFPz[3]-ev.sumPFPz[4]);
      outVars["PFPzDiffEE"]     = fabs(ev.sumPFPz[2]-ev.sumPFPz[5]);   
      outVars["PFPzDiffHE"]     = fabs(ev.sumPFPz[1]-ev.sumPFPz[6]);    
      outVars["PFPzDiffHF"]     = fabs(ev.sumPFPz[0]-ev.sumPFPz[7]); 
      outVars["PFHtSumEB"]      = ev.sumPFHt[3]+ev.sumPFHt[4];
      outVars["PFHtSumEE"]      = ev.sumPFHt[2]+ev.sumPFHt[5];   
      outVars["PFHtSumHE"]      = ev.sumPFHt[1]+ev.sumPFHt[6];    
      outVars["PFHtSumHF"]      = ev.sumPFHt[0]+ev.sumPFHt[7]; 
      outVars["PFHtDiffEB"]     = fabs(ev.sumPFHt[3]-ev.sumPFHt[4]);
      outVars["PFHtDiffEE"]     = fabs(ev.sumPFHt[2]-ev.sumPFHt[5]);   
      outVars["PFHtDiffHE"]     = fabs(ev.sumPFHt[1]-ev.sumPFHt[6]);    
      outVars["PFHtDiffHF"]     = fabs(ev.sumPFHt[0]-ev.sumPFHt[7]); 
      outVars["PFChPzSumEB"]    = ev.sumPFChPz[3]+ev.sumPFChPz[4];
      outVars["PFChPzSumEE"]    = ev.sumPFChPz[2]+ev.sumPFChPz[5];   
      outVars["PFChPzSumHE"]    = ev.sumPFChPz[1]+ev.sumPFChPz[6];    
      outVars["PFChPzSumHF"]    = ev.sumPFChPz[0]+ev.sumPFChPz[7]; 
      outVars["PFChPzDiffEB"]   = fabs(ev.sumPFChPz[3]-ev.sumPFChPz[4]);
      outVars["PFChPzDiffEE"]   = fabs(ev.sumPFChPz[2]-ev.sumPFChPz[5]);   
      outVars["PFChPzDiffHE"]   = fabs(ev.sumPFChPz[1]-ev.sumPFChPz[6]);    
      outVars["PFChPzDiffHF"]   = fabs(ev.sumPFChPz[0]-ev.sumPFChPz[7]); 
      outVars["PFChHtSumEB"]    = ev.sumPFChHt[3]+ev.sumPFChHt[4];
      outVars["PFChHtSumEE"]    = ev.sumPFChHt[2]+ev.sumPFChHt[5];   
      outVars["PFChHtSumHE"]    = ev.sumPFChHt[1]+ev.sumPFChHt[6];    
      outVars["PFChHtSumHF"]    = ev.sumPFChHt[0]+ev.sumPFChHt[7]; 
      outVars["PFChHtDiffEB"]   = fabs(ev.sumPFChHt[3]-ev.sumPFChHt[4]);
      outVars["PFChHtDiffEE"]   = fabs(ev.sumPFChHt[2]-ev.sumPFChHt[5]);   
      outVars["PFChHtDiffHE"]   = fabs(ev.sumPFChHt[1]-ev.sumPFChHt[6]);    
      outVars["PFChHtDiffHF"]   = fabs(ev.sumPFChHt[0]-ev.sumPFChHt[7]); 

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
          beamXangle = std::round(lhc_cond.crossing_angle);
          outVars["lumiDeliv"] = lhc_cond.luminosity.delivered;
          outVars["lumiReco"]  = lhc_cond.luminosity.recorded;
          ht.fill("beamXangle", beamXangle, plotwgts, selCat);
          
          if(beamXangle==120 || beamXangle==130 || beamXangle==140 || beamXangle==150) {
            
            std::vector< std::pair<int,float> > nearCsis;
            std::map<int,int> ntks;
            ntks[23]=0; ntks[123]=0;
            for (int ift=0; ift<ev.nfwdtrk; ift++) {
              if(ev.fwdtrk_method[ift]!=0) continue;
            
              const unsigned short pot_raw_id = ev.fwdtrk_pot[ift];
              float xi=ev.fwdtrk_xi[ift];
              if (pot_raw_id!=3 && pot_raw_id!=23 && pot_raw_id!=103 && pot_raw_id!=123) continue;              
              if (pot_raw_id==23 || pot_raw_id==123) {
                RPid[nRPtk]=pot_raw_id;
                RPfarcsi[nRPtk]=xi;
                RPnearcsi[nRPtk]=0;              
                nRPtk++;
                
                //monitor track multiplicity and csi values
                if(ntks.find(pot_raw_id)==ntks.end()) ntks[pot_raw_id]=0;

                ntks[pot_raw_id]++;
                ht.fill("csirp",xi,plotwgts, Form("%s_%d",selCat.Data(),pot_raw_id));
                
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
              ht.fill("ntkrp", nit.second, plotwgts, Form("%s_%d",selCat.Data(),nit.first));

          }
          
        }catch(...){
        }
      }

      //training category
      float trainCatVal=-1;
      if(nRPtk==0) trainCatVal=0;
      if(nRPtk==2 && RPid[0]*RPid[1]==23*123) trainCatVal=1;
      outVars["trainCat"]=trainCatVal;
      
      if(debug && isZ) {
        debug_out.precision(3);
        debug_out << ev.run << " " << ev.lumi << " " << ev.event << " " 
                  << hasMTrigger << " " << hasMMTrigger << " "
                  << boson.Pt() << " " << boson.Eta() << " " << boson.Phi() << " " << boson.M() << " " 
                  << ev.nrawmu-2 << " ";
        for(int irp=0; irp<nRPtk; irp++)
          debug_out << RPid[irp] << " " << RPfarcsi[irp] << " ";
        debug_out << endl;
      }


      outT->Fill();
    }
      
  if(debug) debug_out.close();

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
