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
  proton_reco.interpolateCrossingAngle( 150., 0.046 );

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

  //variables of doom
  outT->Branch("nrawmu", &ev.nrawmu, "nrawmu/I");
  outT->Branch("rawmu_pt", ev.rawmu_pt, "rawmu_pt[nrawmu]/S");
  outT->Branch("rawmu_eta", ev.rawmu_eta, "rawmu_eta[nrawmu]/S");
  outT->Branch("rawmu_phi", ev.rawmu_phi, "rawmu_phi[nrawmu]/S");
  outT->Branch("rawmu_pid", ev.rawmu_pid, "rawmu_pid[nrawmu]/I");
  outT->Branch("vtxHt", ev.vtxHt, "vtxHt[nvtx]/F");
  outT->Branch("vtxPt", ev.vtxPt, "vtxPt[nvtx]/F");
  outT->Branch("vtxMult", ev.vtxMult, "vtxMult[nvtx]/I");

  bool hasETrigger,hasMTrigger,hasMMTrigger,hasEETrigger,hasEMTrigger,hasATrigger;
  outT->Branch("hasETrigger",&hasETrigger,"hasETrigger/O");
  outT->Branch("hasMTrigger",&hasMTrigger,"hasMTrigger/O");
  outT->Branch("hasEMTrigger",&hasEMTrigger,"hasEMTrigger/O");
  outT->Branch("hasMMTrigger",&hasMMTrigger,"hasMMTrigger/O");
  outT->Branch("hasEETrigger",&hasEETrigger,"hasEETrigger/O");
  outT->Branch("hasATrigger",&hasATrigger,"hasATrigger/O");

  bool isSS,isSF,isZ,isA;
  outT->Branch("isSS",&isSS,"isSS/O");
  outT->Branch("isSF",&isSF,"isSF/O");
  outT->Branch("isZ",&isZ,"isZ/O");
  outT->Branch("isA",&isA,"isA/O");
  
  ADDVAR(&ev.rho,"rho","F",outT);
  ADDVAR(&ev.pf_closestDZnonAssoc,"closestDZnonAssoc","F",outT);
  ADDVAR(&ev.pf_ch_wgtSum,"nch","F",outT);
  ADDVAR(&ev.met_pt,"met_pt","F",outT);
  ADDVAR(&ev.met_phi,"met_phi","F",outT);
  ADDVAR(&ev.met_sig,"met_sig","F",outT);
  TString fvars[]={"evwgt", "evcat", 
                   "l1pt", "l1eta", "l1phi", "ml1", "l1id", "mt1",
                   "l2pt", "l2eta", "l2phi", "ml2", "l2id", "mt2",
                   "bosonpt", "bosonptunc","bosoneta", "bosonphi", "mboson", "bosonht", "mtboson", 
                   "llacopl", "llcosthetaCS", "llphistar", "llMR", "llR", 
                   "llcsip", "llcsipUnc", "llcsim", "llcsimUnc",
                   "nb", "nj", "nl","ng","ht", "htb", "htj",
                   "jsumposen","jsumnegen",
                   "jsumposhfen","jsumneghfen"};
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
  //if (debug) nentries = 10000; //restrict number of entries for testing
  t->GetEntry(0);

  cout << "...producing " << outname << " from " << nentries << " events" << endl;
  
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
  ht.addHist("ptboson",         new TH1F("ptboson",        ";Transverse momentum [GeV];Events",sizeof(ptbosonbins)/sizeof(Float_t)-1,ptbosonbins));
  ht.addHist("phistar",      new TH1F("phistar",     ";Dilepton #phi^{*};Events",50,0,5000));
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

      //particle level dilepton or photon
      std::vector<double> trivialwgts(1,1.0);
      float gen_pt(-1),gen_m(-1),gen_dr(-1);
      std::vector<TString> gen_cats;
      if(!ev.isData){
        std::vector<Particle> genLeptons=selector.getGenLeptons(ev,20.,2.4);
        std::vector<Particle> genPhotons=selector.getGenPhotons(ev,50.,2.5);
        
        if(genLeptons.size()>=2 && genLeptons[0].Pt()>30 && fabs(genLeptons[0].Eta())<2.1) {
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
      hasATrigger=(selector.hasTriggerBit("HLT_Photon200_v", ev.triggerBits));
      hasETrigger=(selector.hasTriggerBit("HLT_Ele35_WPTight_Gsf_v", ev.triggerBits));
      bool hasHighPtMTrigger=(selector.hasTriggerBit("HLT_Mu50_v",     ev.triggerBits));
      bool hasStdMTrigger=(selector.hasTriggerBit("HLT_IsoMu24_v",     ev.triggerBits) ||
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
          hasStdMTrigger=selector.hasTriggerBit("HLT_IsoMu27_v",ev.triggerBits);
        }
        if(!(filename.Contains("2017A") || filename.Contains("2017B"))){
          hasMMTrigger=(selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v",   ev.triggerBits) ||
                        selector.hasTriggerBit("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v", ev.triggerBits) );
        }
      }
      hasMTrigger=(hasHighPtMTrigger || hasStdMTrigger);

      //trigger efficiency
      for(auto gen_cat : gen_cats) {
        if( (gen_cat.Contains("gena") && hasATrigger) || 
            (gen_cat.Contains("genee") && (hasEETrigger || hasETrigger)) ||
            (gen_cat.Contains("genmm") && (hasMMTrigger || hasMTrigger)) ||
            (gen_cat.Contains("genem") && (hasEMTrigger || hasETrigger || hasMTrigger)) ) { 
          ht.fill("ptboson", gen_pt, trivialwgts, gen_cat+"trig");
          ht.fill("mll",  gen_m,  trivialwgts, gen_cat+"trig");
          ht.fill("drll", gen_dr, trivialwgts, gen_cat+"trig");
        }
        if( (gen_cat.Contains("gena") && hasATrigger) || 
            (gen_cat.Contains("genee") && hasETrigger) ||
            (gen_cat.Contains("genmm") && hasStdMTrigger) ||
            (gen_cat.Contains("genem") && (hasETrigger || hasStdMTrigger)) ) { 
          ht.fill("ptboson", gen_pt, trivialwgts, gen_cat+"singletrig");
          ht.fill("mll",  gen_m,  trivialwgts, gen_cat+"singletrig");
          ht.fill("drll", gen_dr, trivialwgts, gen_cat+"singletrig");
        }
      }
      
      //identify the offline final state from the leading leptons or the photon
      int l1idx(0),l2idx(1);
      std::vector<Particle> leptons = selector.flaggedLeptons(ev);     
      leptons = selector.selLeptons(leptons,SelectionTool::LOOSE,SelectionTool::MVA90,20,2.5);
      std::vector<Particle> photons=selector.flaggedPhotons(ev);
      photons=selector.selPhotons(photons,SelectionTool::MVA90,leptons,200.,2.4);

      //jets
      std::vector<Jet> allJets = selector.getGoodJets(ev,15.,4.7,leptons,photons);

      //met
      TLorentzVector met(0,0,0,0);
      met.SetPtEtaPhiM(ev.met_pt,0,ev.met_phi,0.);
      
      //
      //OFFLINE SELECTION(S)
      //
      bool passTightSel(false),passMediumSel(false);
      if(leptons.size()>1){
        
        bool isTrigSafe(leptons[0].Pt()>30 && fabs(leptons[0].Eta())<2.1);

        bool isLeadingTight( (leptons[0].id()==11 && leptons[0].hasQualityFlag(SelectionTool::MVA80)) ||
                             (leptons[0].id()==13 && leptons[0].hasQualityFlag(SelectionTool::TIGHT)) );
        bool isSubLeadingTight( (leptons[1].id()==11 && leptons[1].hasQualityFlag(SelectionTool::MVA80)) ||
                                (leptons[1].id()==13 && leptons[1].hasQualityFlag(SelectionTool::TIGHT)) );
        passTightSel = (isTrigSafe && isLeadingTight && isSubLeadingTight);
        
        //bool isLeadingMedium( (leptons[0].id()==11 && leptons[0].hasQualityFlag(SelectionTool::MVA90)) ||
        //                      (leptons[0].id()==13 && leptons[0].hasQualityFlag(SelectionTool::LOOSE)) );        
        //bool isSubLeadingMedium( (leptons[1].id()==11 && leptons[1].hasQualityFlag(SelectionTool::MVA90)) ||
        //                         (leptons[1].id()==13 && leptons[1].hasQualityFlag(SelectionTool::LOOSE)) );
        bool isLeadingMedium( (leptons[0].id()==11 && leptons[0].hasQualityFlag(SelectionTool::MVA90)) ||
                              (leptons[0].id()==13 && leptons[0].hasQualityFlag(SelectionTool::TIGHT)) );        
        bool isSubLeadingMedium( (leptons[1].id()==11 && leptons[1].hasQualityFlag(SelectionTool::MVA90)) ||
                                 (leptons[1].id()==13 && leptons[1].hasQualityFlag(SelectionTool::TIGHT)) );

        passMediumSel = (isTrigSafe && isLeadingMedium && isSubLeadingMedium);        
      }
      else if(photons.size()>0){
        bool isTrigSafe(photons[0].Pt()>200 && fabs(photons[0].Eta())<2.4);
        bool isTight(photons[0].hasQualityFlag(SelectionTool::MVA80));
        passTightSel=(isTrigSafe && isTight);
        passMediumSel=isTrigSafe;
      }

      //apply selection
      TString selCat("");
      int selCode(0);
      if(leptons.size()<2 && photons.size()==0) continue;
      if(!passMediumSel) continue;

      //set kinematics
      TLorentzVector boson(0,0,0,0);
      TLorentzVector lm(0,0,0,0),lp(0,0,0,0);
      float lmScaleUnc(0),lpScaleUnc(0),bosonScaleUnc(0);
      float mass(0);
      ValueCollection_t llcsi;
      float llacopl(0),llphistar(0),llcosthetaCS(0),llMR(0),llR(0),mtm(0),mtp(0),drll(0),bosonht(0),mtboson(0);
      if(leptons.size()>1){
        selCode=leptons[l1idx].id()*leptons[l2idx].id();
        isSF=( leptons[l1idx].id()==leptons[l2idx].id() );
        isSS=( leptons[l1idx].charge()*leptons[l2idx].charge() > 0 );
        if(selCode==11*11) selCat="ee";
        if(selCode==11*13) selCat="em";
        if(selCode==13*13) selCat="mm";

        lm=TLorentzVector(leptons[l1idx].charge()>0 ? leptons[l1idx] : leptons[l2idx]);
        lmScaleUnc=(leptons[l1idx].charge()>0 ? leptons[l1idx].scaleUnc() : leptons[l2idx].scaleUnc());
        lp=TLorentzVector (leptons[l1idx].charge()>0 ? leptons[l2idx] : leptons[l1idx]);
        lpScaleUnc=(leptons[l1idx].charge()>0 ? leptons[l2idx].scaleUnc() : leptons[l1idx].scaleUnc());
        if(isSS)  { 
          lm=leptons[l1idx]; 
          lmScaleUnc=leptons[l1idx].scaleUnc();
          lp=leptons[l2idx];
          lpScaleUnc=leptons[l2idx].scaleUnc();
        }
        boson=lm+lp;
        bosonScaleUnc=TMath::Sqrt( pow(lm.Pt()*lmScaleUnc,2)+pow(lp.Pt()*lpScaleUnc,2) )/boson.Pt();
        bosonht=lm.Pt()+lp.Pt();

        //dilepton specific
        mass=boson.M();
        llcsi=calcCsi(lm,lmScaleUnc,lp,lpScaleUnc);
        llacopl = computeAcoplanarity(lm,lp);
        drll=leptons[l1idx].DeltaR(leptons[l2idx]);
        llphistar=computePhiStar(lm,lp);
        llcosthetaCS=computeCosThetaStar(lm,lp);
        llMR=computeMR(lm,lp);
        llR=computeRsq(lm,lp,met);
        mtm=computeMT(lm,met);
        mtp=computeMT(lp,met);        
      }else if(photons.size()>0l) {
        selCode=22;
        isSF=false;
        isSS=false;
        selCat="a";
        boson=photons[0];
        bosonScaleUnc=photons[0].scaleUnc();
        bosonht=boson.Pt();
      }
      mtboson=computeMT(boson,met);
    
      //further selection for dileptons
      if(selCode!=22 && mass<20) continue;
      isZ=( isSF && !isSS && fabs(mass-91)<10);
      isA=(selCode==22);

      //check again origin of the boson in data to max. efficiency and avoid double counting
      if(ev.isData) {
        if(isA) {
          if( !selector.isPhotonPD() ) continue;
          if( !hasATrigger) continue;
        }
        if(selCode==11*11) {
          if( !selector.isDoubleEGPD()      && !selector.isSingleElectronPD()) continue;
          if( selector.isDoubleEGPD()       && !hasEETrigger ) continue;
          if( selector.isSingleElectronPD() && !(hasETrigger && !hasEETrigger) ) continue;
        }
        if(selCode==13*13) {
          if( !selector.isDoubleMuonPD() && !selector.isSingleMuonPD()) continue;
          if( selector.isDoubleMuonPD()  && !hasMMTrigger ) continue;
          if( selector.isSingleMuonPD()  && !(hasMTrigger && !hasMMTrigger) ) continue;
        }
        if(selCode==11*13) {
          if( !selector.isMuonEGPD()        && !selector.isSingleElectronPD() && !selector.isSingleMuonPD()) continue;
          if( selector.isMuonEGPD()         && !hasEMTrigger ) continue;
          if( selector.isSingleElectronPD() && !(hasETrigger && !hasEMTrigger) ) continue;
          if( selector.isSingleMuonPD()     && !(hasMTrigger && !hasETrigger && !hasEMTrigger) ) continue;
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
      float jsumposen(0.),jsumnegen(0.),jsumposhfen(0.),jsumneghfen(0.);
      for(size_t ij=0; ij<allJets.size(); ij++) 
        {
          int idx=allJets[ij].getJetIndex();
          bool passBtag(ev.j_btag[idx]>0);

          //int jid=ev.j_id[idx];
          //bool passLoosePu((jid>>2)&0x1);          
          // if(!passLoosePu) continue;

          if(passBtag) { bJets.push_back(allJets[ij]);     scalarhtb+=allJets[ij].pt();  }
          else         { lightJets.push_back(allJets[ij]); scalarhtj+= allJets[ij].pt(); }
          njets++;

          bool isFwd(fabs(allJets[ij].eta())>2.5);
          float jen(allJets[ij].E());
          if(allJets[ij].eta()<0) {
            jsumnegen+= jen; 
            if(isFwd) jsumneghfen += jen;
          } else {
            jsumposen+= jen; 
            if(isFwd) jsumposhfen += jen;
          }
          
          jets.push_back(allJets[ij]);
          scalarht += jets[ij].pt();          
          float dphij2met=fabs(allJets[ij].DeltaPhi(met));
          if(dphij2met>mindphijmet) continue;
          mindphijmet=dphij2met;
        }
      
      //baseline categories and additional stuff produced with the dilepton
      std::vector<TString> cats(1,selCat);
      if(boson.Pt()>50) cats.push_back(selCat+"highpt");
        
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

      //dilepton system
      ht.fill("nlep",       leptons.size(),  plotwgts, cats);
      if(!isA){
        ht.fill("lmpt",       lm.Pt(),         plotwgts, cats);
        ht.fill("lmeta",      fabs(lm.Eta()),  plotwgts, cats);
        ht.fill("lppt",       lp.Pt(),         plotwgts, cats);
        ht.fill("lmeta",      fabs(lp.Eta()),  plotwgts, cats);
        ht.fill("mll",        boson.M(),         plotwgts, cats);        
        ht.fill("drll",       drll,plotwgts,cats);
        ht.fill("phistar",    llphistar,  plotwgts, cats);
        ht.fill("costhetaCS", llcosthetaCS,      plotwgts, cats);
        ht.fill("acopl",      llacopl,           plotwgts, cats);
      }
      ht.fill("ptboson",       boson.Pt(),        plotwgts, cats);
      
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
      ht.fill("nch",           ev.pf_ch_wgtSum,         plotwgts, cats);
      
      //fill tree with central detector information
      outVars["evwgt"]=plotwgts[0];
      outVars["evcat"]=float(selCode);
      
      outVars["l1pt"]=lm.Pt();
      outVars["l1eta"]=lm.Eta();
      outVars["l1phi"]=lm.Phi(); 
      outVars["ml1"]=lm.M();
      outVars["l1id"]=leptons.size()>1 ? leptons[l1idx].id() : 0.;
      outVars["mt1"]=mtm;
      
      outVars["l2pt"]=lp.Pt();
      outVars["l2eta"]=lp.Eta();
      outVars["l2phi"]=lp.Phi();
      outVars["ml2"]=lp.M();
      outVars["l2id"]=leptons.size()>1 ? leptons[l2idx].id() : 0.;
      outVars["mt2"]=mtp;
      
      outVars["bosonpt"]=boson.Pt();
      outVars["bosonptunc"]=bosonScaleUnc;
      outVars["bosoneta"]=boson.Eta();
      outVars["bosonphi"]=boson.Phi();
      outVars["mboson"]=boson.M();
      outVars["mtboson"]=mtboson;
      outVars["bosonht"]=bosonht;      

      //dilepton specific
      outVars["llacopl"]=llacopl;
      outVars["llcosthetaCS"]=llcosthetaCS;
      outVars["llphistar"]=llphistar;
      outVars["llMR"]=llMR;
      outVars["llR"]=llR;      
      outVars["llcsip"]    = llcsi.size()>0 ? llcsi[0].first  : 0;
      outVars["llcsipUnc"] = llcsi.size()>0 ? llcsi[0].second : 0;
      outVars["llcsim"]    = llcsi.size()>1 ? llcsi[1].first  : 1;
      outVars["llcsimUnc"] = llcsi.size()>1 ? llcsi[1].second : 1;
      
      outVars["nb"]=bJets.size();
      outVars["nj"]=lightJets.size();
      outVars["nl"]=leptons.size();
      outVars["ng"]=photons.size();
      outVars["nch"]=float(ev.pf_ch_wgtSum);
      outVars["ht"]=scalarht;
      outVars["htb"]=scalarhtb;
      outVars["htj"]=scalarhtj;
      outVars["jsumnegen"]=jsumnegen;
      outVars["jsumneghfen"]=jsumneghfen;
      outVars["jsumposen"]=jsumposen;
      outVars["jsumposhfen"]=jsumposhfen;

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
          ht.fill("beamXangle", beamXangle, plotwgts, selCat);
          
          //only dispersions for these angles are available (@ CTPPSAnalysisTools/data/2017/dispersions.txt)
          //150 is interpolated
          if(beamXangle==120 || beamXangle==130 || beamXangle==140 || beamXangle==150) {
            
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
