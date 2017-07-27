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
#include "TopLJets2015/TopAnalysis/interface/CorrectionTools.h"
#include "TopLJets2015/TopAnalysis/interface/TOPJetShape.h"
#include "TopLJets2015/TopAnalysis/interface/LeptonEfficiencyWrapper.h"
#include "TopLJets2015/TopAnalysis/interface/BtagUncertaintyComputer.h"

#include <vector>
#include <set>
#include <iostream>
#include <algorithm>

#include "TMath.h"

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>

#include "fastjet/tools/Recluster.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"
using namespace fastjet;
using namespace fastjet::contrib;

#include "Rivet/Math/MatrixN.hh"
#include "Rivet/Math/MatrixDiag.hh"
using Rivet::Matrix;
using Rivet::EigenSystem;

using namespace std;


//
void RunTopJetShape(TString filename,
		    TString outname,
		    Int_t channelSelection, 
		    Int_t chargeSelection, 
		    SelectionTool::FlavourSplitting flavourSplitting,
		    TH1F *normH, 
		    Bool_t runSysts,
		    std::string systVar,
		    TString era,
		    Bool_t debug)
{

  /////////////////////
  // INITIALIZATION //
  ///////////////////

  bool isTTbar( filename.Contains("_TTJets") or (normH and TString(normH->GetTitle()).Contains("_TTJets")));
  bool isData( filename.Contains("Data") );
  
  // explicit systematics
  std::vector<std::string> vSystVar;
  boost::split(vSystVar, systVar, boost::is_any_of("_"));

  //PREPARE OUTPUT
  TopJetShapeEvent_t tjsev;
  TString baseName=gSystem->BaseName(outname); 
  TString dirName=gSystem->DirName(outname);
  TFile *fOut=TFile::Open(dirName+"/"+baseName,"RECREATE");
  fOut->cd();
  TTree *outT=new TTree("tjsev","tjsev");
  createTopJetShapeEventTree(outT,tjsev);
  outT->SetDirectory(fOut);


  //READ TREE FROM FILE
  MiniEvent_t ev;
  TFile *f = TFile::Open(filename);
  TH1 *genPU=(TH1 *)f->Get("analysis/putrue");
  TH1 *triggerList=(TH1 *)f->Get("analysis/triggerList");
  TTree *t = (TTree*)f->Get("analysis/data");
  attachToMiniEventTree(t,ev,true);
  Int_t nentries(t->GetEntriesFast());
  if (debug) nentries = 10000; //restrict number of entries for testing
  t->GetEntry(0);

  cout << "...producing " << outname << " from " << nentries << " events" << endl;

  //lumi
  TH1F *ratevsrunH=0;
  std::map<Int_t,Float_t> lumiMap;
  if( isData )  
    {
      std::pair<std::map<Int_t,Float_t>, TH1F *> result=parseLumiInfo(era);
      lumiMap   = result.first;
      ratevsrunH = result.second;
    }
  
  //PILEUP WEIGHTING
  std::map<TString, std::vector<TGraph *> > puWgtGr;
  if( !isData ) puWgtGr=getPileupWeightsMap(era,genPU);
  
  //LEPTON EFFICIENCIES
  LeptonEfficiencyWrapper lepEffH(filename.Contains("Data13TeV"),era+"GH");


  //B-TAG CALIBRATION
  BTagSFUtil* myBTagSFUtil = new BTagSFUtil();
  std::map<TString, std::map<BTagEntry::JetFlavor, BTagCalibrationReader *> > btvsfReaders = getBTVcalibrationReadersMap(era, BTagEntry::OP_MEDIUM);

  //dummy calls
  btvsfReaders["GH"][BTagEntry::FLAV_B]->eval_auto_bounds("central", BTagEntry::FLAV_B,   0., 30.);
  btvsfReaders["GH"][BTagEntry::FLAV_UDSG]->eval_auto_bounds("central", BTagEntry::FLAV_UDSG,   0., 30.);

  std::map<BTagEntry::JetFlavor, TGraphAsymmErrors *>    expBtagEffPy8 = readExpectedBtagEff(era);
  TString btagExpPostFix("");
  if(isTTbar) {
    if(filename.Contains("_herwig"))    btagExpPostFix="_herwig";
    if(filename.Contains("_scaleup"))   btagExpPostFix="_scaleup";
    if(filename.Contains("_scaledown")) btagExpPostFix="_scaledown";
  }
  std::map<BTagEntry::JetFlavor, TGraphAsymmErrors *> expBtagEff=readExpectedBtagEff(era,btagExpPostFix);
  
  
  //JET ENERGY UNCERTAINTIES
  std::string jecVar = "Total";
  if (vSystVar[0] == "jec") {
    if (vSystVar.size() != 3) {
      std::cout << "Wrong format of JEC uncertainty, expected jec_Source_up/down. Exiting..." << std::endl;
      return;
    }
    jecVar = vSystVar[1];
  }
  TString jecUncUrl(era+"/Summer16_23Sep2016V4_MC_UncertaintySources_AK4PFchs.txt");
  gSystem->ExpandPathName(jecUncUrl);
  JetCorrectorParameters *jecParam = new JetCorrectorParameters(jecUncUrl.Data(), jecVar);
  JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty( *jecParam );
  
  //JER SMEARING
  TString resolutionFile(era+"/Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt");
  gSystem->ExpandPathName(resolutionFile);
  JME::JetResolution *jer = new JME::JetResolution(resolutionFile.Data());
  
  //BFRAG UNCERTAINTIES
  std::map<TString, TGraph*> bfrag = getBFragmentationWeights(era);
  std::map<TString, std::map<int, double> > semilepbr = getSemilepBRWeights(era);
  

  //BOOK HISTOGRAMS
  HistTool ht;
  if (isData) ht.setNsyst(0);
  std::map<TString, TH1 *> allPlots;
  std::map<TString, TH2 *> all2dPlots;
  allPlots["puwgtctr"] = new TH1F("puwgtctr","Weight sums",2,0,2);

  std::vector<TString> stageVec = { "0_pre", "1_1l", "2_1l4j", "3_1l4j2b", "4_1l4j2b2w" };
  std::vector<TString> chTags = { "L", "E", "M" };
  for(auto& stage : stageVec) {
    for(auto& channel : chTags) {  
      TString tag(channel+stage+"_");
      
      if(ratevsrunH)
        allPlots[tag+"ratevsrun"] = (TH1 *)ratevsrunH->Clone(tag+"ratevsrun");
      ht.addHist(tag+"nvtx", new TH1F(tag+"nvtx",";Vertex multiplicity;Events",55,0,55));
      ht.addHist(tag+"nleps", new TH1F(tag+"nleps",";Lepton multiplicity;Events",5,-0.5,4.5));
      ht.addHist(tag+"njets", new TH1F(tag+"njets",";Jet multiplicity;Events",15,-0.5,14.5));
      ht.addHist(tag+"nbjets", new TH1F(tag+"nbjets",";b jet multiplicity;Events",5,-0.5,4.5));
      ht.addHist(tag+"nwjets", new TH1F(tag+"nwjets",";W jet multiplicity;Events",10,-0.5,9.5));
      ht.addHist(tag+"wcandm", new TH1F(tag+"wcandm",";W candidate mass;W candidates",60,0,300));
      ht.addHist(tag+"tcandm", new TH1F(tag+"tcandm",";top candidate mass;top candidates",70,50,400));
      ht.addHist(tag+"tcandwcutm", new TH1F(tag+"tcandwcutm",";top candidate mass;top candidates (W cut)",70,50,400));
      for(int i=0; i<2; i++) {
        TString pf(Form("l%d",i));          
        ht.addHist(tag+pf+"pt", new TH1F(tag+pf+"pt",";Lepton p_{t} [GeV];Events",50,0,250));
        ht.addHist(tag+pf+"eta", new TH1F(tag+pf+"eta",";Lepton pseudo-rapidity;Events",50,-2.5,2.5));
      }
      for(int i=0; i<6; i++) {
        TString pf(Form("j%d",i));
        ht.addHist(tag+pf+"pt", new TH1F(tag+pf+"pt",";Jet transverse momentum [GeV];Events",50,0,250));
        ht.addHist(tag+pf+"eta", new TH1F(tag+pf+"eta",";Jet pseudo-rapidity;Events",50,-5,5));
      }
      ht.addHist(tag+"met", new TH1F(tag+"met",";MET [GeV];Events",50,0,250));
    }
  }
  
  ht.addHist("js_mult_charged", new TH1F("js_mult_charged",";N (charged);Jets",30,0,30));
//  ht.addHist("js_mult_puppi", new TH1F("js_mult_puppi",";N (puppi);Jets",30,0,30));
//  ht.addHist("js_mult_all", new TH1F("js_mult_all",";N (all);Jets",30,0,30));
  ht.addHist("js_width_charged", new TH1F("js_width_charged",";width (charged);Jets",50,0,0.25));
//  ht.addHist("js_width_puppi", new TH1F("js_width_puppi",";width (puppi);Jets",50,0,0.25));
//  ht.addHist("js_width_all", new TH1F("js_width_all",";width (all);Jets",50,0,0.25));
  ht.addHist("js_ptd_charged", new TH1F("js_ptd_charged",";p_{T}D (charged);Jets",50,0,1));
//  ht.addHist("js_ptd_puppi", new TH1F("js_ptd_puppi",";p_{T}D (puppi);Jets",50,0,1));
//  ht.addHist("js_ptd_all", new TH1F("js_ptd_all",";p_{T}D (all);Jets",50,0,1));
  ht.addHist("js_ptds_charged", new TH1F("js_ptds_charged",";scaled p_{T}D (charged);Jets",50,0,1));
//  ht.addHist("js_ptds_puppi", new TH1F("js_ptds_puppi",";scaled p_{T}D (puppi);Jets",50,0,1));
//  ht.addHist("js_ptds_all", new TH1F("js_ptds_all",";scaled p_{T}D (all);Jets",50,0,1));
  ht.addHist("js_ecc_charged", new TH1F("js_ecc_charged",";eccentricity (charged);Jets",50,0,1));
//  ht.addHist("js_ecc_puppi", new TH1F("js_ecc_puppi",";eccentricity (puppi);Jets",50,0,1));
//  ht.addHist("js_ecc_all", new TH1F("js_ecc_all",";eccentricity (all);Jets",50,0,1));
  ht.addHist("js_tau21_charged", new TH1F("js_tau21_charged",";#tau_{21} (charged);Jets",50,0,1));
//  ht.addHist("js_tau21_puppi", new TH1F("js_tau21_puppi",";#tau_{21} (puppi);Jets",50,0,1));
//  ht.addHist("js_tau21_all", new TH1F("js_tau21_all",";#tau_{21} (all);Jets",50,0,1));
  ht.addHist("js_tau32_charged", new TH1F("js_tau32_charged",";#tau_{32} (charged);Jets",50,0,1));
//  ht.addHist("js_tau32_puppi", new TH1F("js_tau32_puppi",";#tau_{32} (puppi);Jets",50,0,1));
//  ht.addHist("js_tau32_all", new TH1F("js_tau32_all",";#tau_{32} (all);Jets",50,0,1));
  ht.addHist("js_tau43_charged", new TH1F("js_tau43_charged",";#tau_{43} (charged);Jets",50,0,1));
//  ht.addHist("js_tau43_puppi", new TH1F("js_tau43_puppi",";#tau_{43} (puppi);Jets",50,0,1));
//  ht.addHist("js_tau43_all", new TH1F("js_tau43_all",";#tau_{43} (all);Jets",50,0,1));
  ht.addHist("js_zg_charged", new TH1F("js_zg_charged",";z_{g} (charged);Jets",40,0.1,0.5));
//  ht.addHist("js_zg_puppi", new TH1F("js_zg_puppi",";z_{g} (puppi);Jets",40,0.1,0.5));
//  ht.addHist("js_zg_all", new TH1F("js_zg_all",";z_{g} (all);Jets",40,0.1,0.5));
  ht.addHist("js_zgxdr_charged", new TH1F("js_zgxdr_charged",";z_{g} #times #DeltaR (charged);Jets",50,0,0.25));
//  ht.addHist("js_zgxdr_puppi", new TH1F("js_zgxdr_puppi",";z_{g} #times #DeltaR (puppi);Jets",50,0,0.25));
//  ht.addHist("js_zgxdr_all", new TH1F("js_zgxdr_all",";z_{g} #times #DeltaR (all);Jets",50,0,0.25));
  ht.addHist("js_zgdr_charged", new TH1F("js_zgdr_charged",";z_{g} #DeltaR (charged);Jets",50,0,0.5));
//  ht.addHist("js_zgdr_puppi", new TH1F("js_zgdr_puppi",";z_{g} #DeltaR (puppi);Jets",50,0,0.5));
//  ht.addHist("js_zgdr_all", new TH1F("js_zgdr_all",";z_{g} #DeltaR (all);Jets",50,0,0.5));
  ht.addHist("js_ga_width_charged", new TH1F("js_ga_width_charged",";#lambda_{ 1}^{1} (width) (charged);Jets",50,0,1));
//  ht.addHist("js_ga_width_puppi", new TH1F("js_ga_width_puppi",";#lambda_{ 1}^{1} (width) (puppi);Jets",50,0,1));
//  ht.addHist("js_ga_width_all", new TH1F("js_ga_width_all",";#lambda_{ 1}^{1} (width) (all);Jets",50,0,1));
  ht.addHist("js_ga_lha_charged", new TH1F("js_ga_lha_charged",";#lambda_{ 0.5}^{1} (LHA) (charged);Jets",50,0,1));
//  ht.addHist("js_ga_lha_puppi", new TH1F("js_ga_lha_puppi",";#lambda_{ 0.5}^{1} (LHA) (puppi);Jets",50,0,1));
//  ht.addHist("js_ga_lha_all", new TH1F("js_ga_lha_all",";#lambda_{ 0.5}^{1} (LHA) (all);Jets",50,0,1));
  ht.addHist("js_ga_thrust_charged", new TH1F("js_ga_thrust_charged",";#lambda_{ 2}^{1} (thrust) (charged);Jets",50,0,0.5));
//  ht.addHist("js_ga_thrust_puppi", new TH1F("js_ga_thrust_puppi",";#lambda_{ 2}^{1} (thrust) (puppi);Jets",50,0,0.5));
//  ht.addHist("js_ga_thrust_all", new TH1F("js_ga_thrust_all",";#lambda_{ 2}^{1} (thrust) (all);Jets",50,0,0.5));
  ht.addHist("js_c1_02_charged", new TH1F("js_c1_02_charged",";C_{ 1}^{ (0.2)} (charged);Jets",50,0,0.5));
//  ht.addHist("js_c1_02_puppi", new TH1F("js_c1_02_puppi",";C_{ 1}^{ (0.2)} (puppi);Jets",50,0,0.5));
//  ht.addHist("js_c1_02_all", new TH1F("js_c1_02_all",";C_{ 1}^{ (0.2)} (all);Jets",50,0,0.5));
  ht.addHist("js_c1_05_charged", new TH1F("js_c1_05_charged",";C_{ 1}^{ (0.5)} (charged);Jets",60,0,0.3));
//  ht.addHist("js_c1_05_puppi", new TH1F("js_c1_05_puppi",";C_{ 1}^{ (0.5)} (puppi);Jets",60,0,0.3));
//  ht.addHist("js_c1_05_all", new TH1F("js_c1_05_all",";C_{ 1}^{ (0.5)} (all);Jets",60,0,0.3));
  ht.addHist("js_c1_10_charged", new TH1F("js_c1_10_charged",";C_{ 1}^{ (1.0)} (charged);Jets",40,0,0.2));
//  ht.addHist("js_c1_10_puppi", new TH1F("js_c1_10_puppi",";C_{ 1}^{ (1.0)} (puppi);Jets",40,0,0.2));
//  ht.addHist("js_c1_10_all", new TH1F("js_c1_10_all",";C_{ 1}^{ (1.0)} (all);Jets",40,0,0.2));
  ht.addHist("js_c1_20_charged", new TH1F("js_c1_20_charged",";C_{ 1}^{ (2.0)} (charged);Jets",50,0,0.1));
//  ht.addHist("js_c1_20_puppi", new TH1F("js_c1_20_puppi",";C_{ 1}^{ (2.0)} (puppi);Jets",50,0,0.1));
//  ht.addHist("js_c1_20_all", new TH1F("js_c1_20_all",";C_{ 1}^{ (2.0)} (all);Jets",50,0,0.1));
  ht.addHist("js_c2_02_charged", new TH1F("js_c2_02_charged",";C_{ 2}^{ (0.2)} (charged);Jets",35,0,0.7));
//  ht.addHist("js_c2_02_puppi", new TH1F("js_c2_02_puppi",";C_{ 2}^{ (0.2)} (puppi);Jets",35,0,0.7));
//  ht.addHist("js_c2_02_all", new TH1F("js_c2_02_all",";C_{ 2}^{ (0.2)} (all);Jets",35,0,0.7));
  ht.addHist("js_c2_05_charged", new TH1F("js_c2_05_charged",";C_{ 2}^{ (0.5)} (charged);Jets",40,0,0.4));
//  ht.addHist("js_c2_05_puppi", new TH1F("js_c2_05_puppi",";C_{ 2}^{ (0.5)} (puppi);Jets",40,0,0.4));
//  ht.addHist("js_c2_05_all", new TH1F("js_c2_05_all",";C_{ 2}^{ (0.5)} (all);Jets",40,0,0.4));
  ht.addHist("js_c2_10_charged", new TH1F("js_c2_10_charged",";C_{ 2}^{ (1.0)} (charged);Jets",50,0,0.25));
//  ht.addHist("js_c2_10_puppi", new TH1F("js_c2_10_puppi",";C_{ 2}^{ (1.0)} (puppi);Jets",50,0,0.25));
//  ht.addHist("js_c2_10_all", new TH1F("js_c2_10_all",";C_{ 2}^{ (1.0)} (all);Jets",50,0,0.25));
  ht.addHist("js_c2_20_charged", new TH1F("js_c2_20_charged",";C_{ 2}^{ (2.0)} (charged);Jets",30,0,0.15));
//  ht.addHist("js_c2_20_puppi", new TH1F("js_c2_20_puppi",";C_{ 2}^{ (2.0)} (puppi);Jets",30,0,0.15));
//  ht.addHist("js_c2_20_all", new TH1F("js_c2_20_all",";C_{ 2}^{ (2.0)} (all);Jets",30,0,0.15));
  ht.addHist("js_c3_02_charged", new TH1F("js_c3_02_charged",";C_{ 3}^{ (0.2)} (charged);Jets",35,0,0.7));
//  ht.addHist("js_c3_02_puppi", new TH1F("js_c3_02_puppi",";C_{ 3}^{ (0.2)} (puppi);Jets",35,0,0.7));
//  ht.addHist("js_c3_02_all", new TH1F("js_c3_02_all",";C_{ 3}^{ (0.2)} (all);Jets",35,0,0.7));
  ht.addHist("js_c3_05_charged", new TH1F("js_c3_05_charged",";C_{ 3}^{ (0.5)} (charged);Jets",40,0,0.4));
//  ht.addHist("js_c3_05_puppi", new TH1F("js_c3_05_puppi",";C_{ 3}^{ (0.5)} (puppi);Jets",40,0,0.4));
//  ht.addHist("js_c3_05_all", new TH1F("js_c3_05_all",";C_{ 3}^{ (0.5)} (all);Jets",40,0,0.4));
  ht.addHist("js_c3_10_charged", new TH1F("js_c3_10_charged",";C_{ 3}^{ (1.0)} (charged);Jets",50,0,0.25));
//  ht.addHist("js_c3_10_puppi", new TH1F("js_c3_10_puppi",";C_{ 3}^{ (1.0)} (puppi);Jets",50,0,0.25));
//  ht.addHist("js_c3_10_all", new TH1F("js_c3_10_all",";C_{ 3}^{ (1.0)} (all);Jets",50,0,0.25));
  ht.addHist("js_c3_20_charged", new TH1F("js_c3_20_charged",";C_{ 3}^{ (2.0)} (charged);Jets",30,0,0.15));
//  ht.addHist("js_c3_20_puppi", new TH1F("js_c3_20_puppi",";C_{ 3}^{ (2.0)} (puppi);Jets",30,0,0.15));
//  ht.addHist("js_c3_20_all", new TH1F("js_c3_20_all",";C_{ 3}^{ (2.0)} (all);Jets",30,0,0.15));
  
  for (auto& it : allPlots)   { it.second->Sumw2(); it.second->SetDirectory(0); }
  for (auto& it : all2dPlots) { it.second->Sumw2(); it.second->SetDirectory(0); }

  std::cout << "init done" << std::endl;

  ///////////////////////
  // LOOP OVER EVENTS //
  /////////////////////
  
  //EVENT SELECTION WRAPPER
  SelectionTool selector(filename, false, triggerList);
  
  for (Int_t iev=0;iev<nentries;iev++)
    {
      t->GetEntry(iev);
      resetTopJetShapeEvent(tjsev);
      if(iev%int(nentries/100)==0) printf ("[%3.0f%%] done\n", 100.*(float)iev/(float)nentries);
      
      //assign run period GH
      TString period = "GH";
      
      //////////////////
      // CORRECTIONS //
      ////////////////
      
      double csvm = 0.8484;
      if (vSystVar[0] == "csv") {
          if (vSystVar[1] == "heavy") {
              //heavy flavor uncertainty +/-3.5%
              if (vSystVar[2] == "up")   addBTagDecisions(ev, 0.8726, csvm);
              if (vSystVar[2] == "down") addBTagDecisions(ev, 0.8190, csvm);
          }
          if (vSystVar[1] == "light") {
              //light flavor uncertainty +/-10%
              if (vSystVar[2] == "up")   addBTagDecisions(ev, csvm, 0.8557);
              if (vSystVar[2] == "down") addBTagDecisions(ev, csvm, 0.8415);
          }
      }
      else addBTagDecisions(ev, csvm, csvm);
      
      if(!ev.isData) {
        //jec
        if (vSystVar[0] == "jec") {
          applyJetCorrectionUncertainty(ev, jecUnc, jecVar, vSystVar[2]);
        }
        //jer
        if (vSystVar[0] == "jer") {
          smearJetEnergies(ev, jer, vSystVar[1]);
        }
        else smearJetEnergies(ev, jer);
        //b tagging
        if (vSystVar[0] == "btag") {
          if (vSystVar[1] == "heavy") {
            updateBTagDecisions(ev, btvsfReaders[period],expBtagEff,expBtagEffPy8,myBTagSFUtil,vSystVar[2],"central");
          }
          if (vSystVar[1] == "light") {
            updateBTagDecisions(ev, btvsfReaders[period],expBtagEff,expBtagEffPy8,myBTagSFUtil,"central",vSystVar[2]);
          }
        }
        else updateBTagDecisions(ev, btvsfReaders[period],expBtagEff,expBtagEffPy8,myBTagSFUtil);
        //tracking efficiency
        if (vSystVar[0] == "tracking") {
          // "up": no correction
          // "down": apply twice
          if (vSystVar[1] == "down") applyTrackingEfficiencySF(ev, pow(lepEffH.getTrackingCorrection(ev.nvtx, period).first,2));
        }
        else applyTrackingEfficiencySF(ev, lepEffH.getTrackingCorrection(ev.nvtx, period).first);
      }
      
      ///////////////////////////
      // RECO LEVEL SELECTION //
      /////////////////////////
      
      //decide the lepton channel and get selected objects
      TString chTag = selector.flagFinalState(ev);
      std::vector<Particle> &leptons     = selector.getSelLeptons(); 
      std::vector<Jet>      &jets        = selector.getJets();  
      
      //count b and W candidates
      int sel_nbjets = 0;
      int sel_nwjets = 0;

      for (auto& jet : jets) {
        if (jet.flavor() == 5) ++sel_nbjets;
        if (jet.flavor() == 1) ++sel_nwjets;
      }
      
      //event selected on reco level?
      bool preselected          (true);
      bool singleLepton         ((chTag=="E" or chTag=="M") and
                                 (selector.getVetoLeptons().size() == 0));
      bool singleLepton4Jets    (singleLepton and jets.size()>=4);
      bool singleLepton4Jets2b  (singleLepton4Jets and sel_nbjets==2);
      bool singleLepton4Jets2b2W(singleLepton4Jets2b and sel_nwjets>0);
      
      std::vector<bool> recoPass; recoPass.push_back(preselected); recoPass.push_back(singleLepton); recoPass.push_back(singleLepton4Jets); recoPass.push_back(singleLepton4Jets2b); recoPass.push_back(singleLepton4Jets2b2W); 
      
      if (singleLepton4Jets2b2W) tjsev.reco_sel = 1;
      
      tjsev.nj=jets.size();
      
      ////////////////////
      // EVENT WEIGHTS //
      //////////////////
      
      float wgt(1.0);
      // Pairs for systematic uncertainty weights
      // double: weight value (divided by central weight)
      // bool: use weight for plotting, otherwise just save to tree
      std::vector<std::pair<double, bool> > varweights;
      std::vector<double> plotwgts;
      
      allPlots["puwgtctr"]->Fill(0.,1.0);
      if (!ev.isData) {
        // norm weight
        wgt  = (normH? normH->GetBinContent(1) : 1.0);
        
        // pu weight
        double puWgt(puWgtGr[period][0]->Eval(ev.g_pu));
        if (std::isnan(puWgt)) puWgt = 1.;
        if (puWgt == 0.) puWgt = 1.;
        allPlots["puwgtctr"]->Fill(1,puWgt);
        wgt *= puWgt;
        double puWgt1(puWgtGr[period][1]->Eval(ev.g_pu));
        if (std::isnan(puWgt1)) puWgt1 = 1.;
        varweights.push_back(std::make_pair(puWgt1/puWgt, true)); // 1
        double puWgt2(puWgtGr[period][2]->Eval(ev.g_pu));
        if (std::isnan(puWgt2)) puWgt2 = 1.;
        varweights.push_back(std::make_pair(puWgt2/puWgt, true)); // 2
        
        // lepton trigger*selection weights
        if (singleLepton) {
          EffCorrection_t trigSF = lepEffH.getTriggerCorrection(leptons, period);
          varweights.push_back(std::make_pair(1.+trigSF.second, true)); // 3
          varweights.push_back(std::make_pair(1.-trigSF.second, true)); // 4
          EffCorrection_t selSF = lepEffH.getOfflineCorrection(leptons[0], period);
          varweights.push_back(std::make_pair(1.+selSF.second, true));  // 5
          varweights.push_back(std::make_pair(1.-selSF.second, true));  // 6
          wgt *= trigSF.first*selSF.first;
        }
        else varweights.insert(varweights.end(), 4, std::make_pair(1., true));
        
        // bfrag weights
        varweights.push_back(std::make_pair(computeBFragmentationWeight(ev, bfrag["upFrag"]), true));       // 7
        varweights.push_back(std::make_pair(computeBFragmentationWeight(ev, bfrag["downFrag"]), true));     // 8
        varweights.push_back(std::make_pair(computeBFragmentationWeight(ev, bfrag["PetersonFrag"]), true)); // 9
        // weights for semilep BR
        // simultaneous variation for all hadrons, average over particle and antiparticle
        // divide by mean weight from 100k events to avoid change in cross section
        varweights.push_back(std::make_pair(computeSemilepBRWeight(ev, semilepbr["semilepbrUp"], 0, true)/1.00480, true)); // 10
        varweights.push_back(std::make_pair(computeSemilepBRWeight(ev, semilepbr["semilepbrDown"], 0, true)/0.992632, true)); // 11
        
    	  //top pt weighting
        double topptsf = 1.0;
        if(isTTbar) {
          for (int igen=0; igen<ev.ngtop; igen++) {
            if(abs(ev.gtop_id[igen])!=6) continue;
            topptsf *= TMath::Exp(0.0615-0.0005*ev.gtop_pt[igen]);
          }
        }
        varweights.push_back(std::make_pair(topptsf, true)); // TODO: should have been 12 but prevented by false `else`
        
        // lhe weights
        wgt *= (ev.g_nw>0 ? ev.g_w[0] : 1.0);
        
        std::set<std::string> scalesForPlotter = {
          "id1002muR1muF2hdampmt272.7225",    // 13 (12=1,1)
          "id1003muR1muF0.5hdampmt272.7225",  // 14
          "id1004muR2muF1hdampmt272.7225",    // 15
          "id1005muR2muF2hdampmt272.7225",    // 16
          "id1007muR0.5muF1hdampmt272.7225",  // 18
          "id1009muR0.5muF0.5hdampmt272.7225" // 20
        };
        for (int i = 1; i < ev.g_nw; ++i) {
          bool forPlotter = (normH and scalesForPlotter.count(normH->GetXaxis()->GetBinLabel(i)) != 0);
          varweights.push_back(std::make_pair(ev.g_w[i]/ev.g_w[0], forPlotter));
        }
        
        tjsev.nw = 1 + varweights.size();
        tjsev.weight[0]=wgt;
        for (unsigned int i = 0; i < varweights.size(); ++i) {
          tjsev.weight[i+1] = varweights[i].first;
        }
        plotwgts.push_back(wgt);
        for (auto vw : varweights)
          if (vw.second)
            plotwgts.push_back(vw.first);
      }
      else {
        tjsev.nw=1;
        tjsev.weight[0]=wgt;
        plotwgts.push_back(wgt);
      }
      
      //////////////////////////
      // RECO LEVEL ANALYSIS //
      ////////////////////////
      
      //W and top masses
      std::vector<TLorentzVector> wCands;
      for (unsigned int i = 0; i < jets.size(); i++) {
        for (unsigned int j = i+1; j < jets.size(); j++) {
          if (jets[i].flavor()==5 or jets[j].flavor()==5) continue;
          TLorentzVector wCand = jets[i].p4() + jets[j].p4();
          wCands.push_back(wCand);
        }
      }
      std::vector<TLorentzVector> tCands;
      for (unsigned int i = 0; i < jets.size(); i++) {
        if (jets[i].flavor()!=5) continue;
        for (auto& wCand : wCands) {
          TLorentzVector tCand = jets[i].p4() + wCand;
          tCands.push_back(tCand);
        }
      }
      std::vector<TLorentzVector> tCandsWcut;
      for (unsigned int i = 0; i < jets.size(); i++) {
        if (jets[i].flavor()!=5) continue;
        for (auto& wCand : wCands) {
          if (abs(wCand.M()-80.4) > 15.) continue;
          TLorentzVector tCand = jets[i].p4() + wCand;
          tCandsWcut.push_back(tCand);
        }
      }
      
      //control histograms
      for(size_t istage=0; istage<stageVec.size(); istage++) { 
        for(auto& channel : chTags) { 
          if (not recoPass[istage]) continue;
          if (channel == "E" and chTag != "E") continue;
          if (channel == "M" and chTag != "M") continue;
          TString tag(channel+stageVec[istage]+"_");
          
          std::map<Int_t,Float_t>::iterator rIt=lumiMap.find(ev.run);
          if(rIt!=lumiMap.end() && ratevsrunH) allPlots[tag+"ratevsrun"]->Fill(std::distance(lumiMap.begin(),rIt),1./rIt->second);
          
          ht.fill(tag+"nvtx", ev.nvtx, plotwgts);
          ht.fill(tag+"nleps", leptons.size(), plotwgts);
          ht.fill(tag+"njets", jets.size(), plotwgts);
          ht.fill(tag+"nbjets", sel_nbjets, plotwgts);
          ht.fill(tag+"nwjets", sel_nwjets, plotwgts);
          for (auto& wCand : wCands) ht.fill(tag+"wcandm", wCand.M(), plotwgts);
          for (auto& tCand : tCands) ht.fill(tag+"tcandm", tCand.M(), plotwgts);
          for (auto& tCand : tCandsWcut) ht.fill(tag+"tcandwcutm", tCand.M(), plotwgts);
          for(unsigned int i=0; i<leptons.size(); i++) {
            if (i>1) break;
            TString pf(Form("l%d",i));          
            ht.fill(tag+pf+"pt", leptons[i].pt(), plotwgts);
            ht.fill(tag+pf+"eta", leptons[i].eta(), plotwgts);
          }
          for(unsigned int i=0; i<jets.size(); i++) {
            if (i>5) break;
            TString pf(Form("j%d",i));
            ht.fill(tag+pf+"pt", jets[i].pt(), plotwgts);
            ht.fill(tag+pf+"eta", jets[i].eta(), plotwgts);
          }
          ht.fill(tag+"met", ev.met_pt[0], plotwgts);
        }
      }

      //fill leptons
      tjsev.nl=leptons.size();
      int il = 0;
      for(auto& lepton : leptons) {
        tjsev.l_pt[il]  = lepton.pt();
        tjsev.l_eta[il] = lepton.eta();
        tjsev.l_phi[il] = lepton.phi();
        tjsev.l_m[il]   = lepton.m();
        tjsev.l_id[il]  = lepton.id();
        il++;
      }
      
      //fill MET
      tjsev.met_pt=ev.met_pt[0];
      tjsev.met_phi=ev.met_phi[0];
      
      //fill jets (with jet shapes)
      for(int ij=0; ij<(int)jets.size(); ij++) {
        tjsev.j_pt[ij]      = jets[ij].p4().Pt();
        tjsev.j_eta[ij]     = jets[ij].p4().Eta();
        tjsev.j_phi[ij]     = jets[ij].p4().Phi();
        tjsev.j_m[ij]       = jets[ij].p4().M(); 
        tjsev.j_flavor[ij]  = jets[ij].flavor();
        tjsev.j_overlap[ij] = jets[ij].overlap();
        
        if (tjsev.reco_sel != 1) continue;
        if (abs(tjsev.j_eta[ij]) > 2. or tjsev.j_overlap[ij]) continue;
        
        tjsev.j_mult_charged[ij] = getMult(jets[ij]);
//        tjsev.j_mult_all[ij]     = getMult(jets[ij], true);
//        tjsev.j_mult_puppi[ij]   = getMult(jets[ij], true, true);
        
        tjsev.j_ptd_charged[ij] = getPtD(jets[ij]);
//        tjsev.j_ptd_all[ij]     = getPtD(jets[ij], true);
//        tjsev.j_ptd_puppi[ij]   = getPtD(jets[ij], true, true);
        
        tjsev.j_ptds_charged[ij] = getPtDs(jets[ij]);
//        tjsev.j_ptds_all[ij]     = getPtDs(jets[ij], true);
//        tjsev.j_ptds_puppi[ij]   = getPtDs(jets[ij], true, true);
        
        tjsev.j_width_charged[ij] = getWidth(jets[ij]);
//        tjsev.j_width_all[ij]     = getWidth(jets[ij], true);
//        tjsev.j_width_puppi[ij]   = getWidth(jets[ij], true, true);
        
        tjsev.j_ecc_charged[ij] = getEcc(jets[ij]);
//        tjsev.j_ecc_all[ij]     = getEcc(jets[ij], true);
//        tjsev.j_ecc_puppi[ij]   = getEcc(jets[ij], true, true);
        
        tjsev.j_tau21_charged[ij] = getTau(2, 1, jets[ij]);
//        tjsev.j_tau21_all[ij]     = getTau(2, 1, jets[ij], true);
//        tjsev.j_tau21_puppi[ij]   = getTau(2, 1, jets[ij], true, true);
        
        tjsev.j_tau32_charged[ij] = getTau(3, 2, jets[ij]);
//        tjsev.j_tau32_all[ij]     = getTau(3, 2, jets[ij], true);
//        tjsev.j_tau32_puppi[ij]   = getTau(3, 2, jets[ij], true, true);
        
        tjsev.j_tau43_charged[ij] = getTau(4, 3, jets[ij]);
//        tjsev.j_tau43_all[ij]     = getTau(4, 3, jets[ij], true);
//        tjsev.j_tau43_puppi[ij]   = getTau(4, 3, jets[ij], true, true);
        
        std::vector<double> zgResult_charged = getZg(jets[ij]);
//        std::vector<double> zgResult_all     = getZg(jets[ij], true);
//        std::vector<double> zgResult_puppi   = getZg(jets[ij], true, true);
        
        tjsev.j_zg_charged[ij]    = zgResult_charged[0];
//        tjsev.j_zg_all[ij]        = zgResult_all[0];
//        tjsev.j_zg_puppi[ij]      = zgResult_puppi[0];
        
        tjsev.j_zgxdr_charged[ij] = zgResult_charged[1];
//        tjsev.j_zgxdr_all[ij]     = zgResult_all[1];
//        tjsev.j_zgxdr_puppi[ij]   = zgResult_puppi[1];
        
        tjsev.j_zgdr_charged[ij]  = zgResult_charged[2];
//        tjsev.j_zgdr_all[ij]      = zgResult_all[2];
//        tjsev.j_zgdr_puppi[ij]    = zgResult_puppi[2];
        
        tjsev.j_ga_width_charged[ij]  = calcGA(1., 1., jets[ij]);
//        tjsev.j_ga_width_all[ij]      = calcGA(1., 1., jets[ij], true);
//        tjsev.j_ga_width_puppi[ij]    = calcGA(1., 1., jets[ij], true, true);
        
        tjsev.j_ga_lha_charged[ij]    = calcGA(0.5, 1., jets[ij]);
//        tjsev.j_ga_lha_all[ij]        = calcGA(0.5, 1., jets[ij], true);
//        tjsev.j_ga_lha_puppi[ij]      = calcGA(0.5, 1., jets[ij], true, true);
        
        tjsev.j_ga_thrust_charged[ij] = calcGA(2., 1., jets[ij]);
//        tjsev.j_ga_thrust_all[ij]     = calcGA(2., 1., jets[ij], true);
//        tjsev.j_ga_thrust_puppi[ij]   = calcGA(2., 1., jets[ij], true, true);
        
        tjsev.j_c1_02_charged[ij] = getC(1, 0.2, jets[ij]);
//        tjsev.j_c1_02_all[ij]     = getC(1, 0.2, jets[ij], true);
//        tjsev.j_c1_02_puppi[ij]   = getC(1, 0.2, jets[ij], true, true);
        tjsev.j_c1_05_charged[ij] = getC(1, 0.5, jets[ij]);
//        tjsev.j_c1_05_all[ij]     = getC(1, 0.5, jets[ij], true);
//        tjsev.j_c1_05_puppi[ij]   = getC(1, 0.5, jets[ij], true, true);
        tjsev.j_c1_10_charged[ij] = getC(1, 1.0, jets[ij]);
//        tjsev.j_c1_10_all[ij]     = getC(1, 1.0, jets[ij], true);
//        tjsev.j_c1_10_puppi[ij]   = getC(1, 1.0, jets[ij], true, true);
        tjsev.j_c1_20_charged[ij] = getC(1, 2.0, jets[ij]);
//        tjsev.j_c1_20_all[ij]     = getC(1, 2.0, jets[ij], true);
//        tjsev.j_c1_20_puppi[ij]   = getC(1, 2.0, jets[ij], true, true);
        tjsev.j_c2_02_charged[ij] = getC(2, 0.2, jets[ij]);
//        tjsev.j_c2_02_all[ij]     = getC(2, 0.2, jets[ij], true);
//        tjsev.j_c2_02_puppi[ij]   = getC(2, 0.2, jets[ij], true, true);
        tjsev.j_c2_05_charged[ij] = getC(2, 0.5, jets[ij]);
//        tjsev.j_c2_05_all[ij]     = getC(2, 0.5, jets[ij], true);
//        tjsev.j_c2_05_puppi[ij]   = getC(2, 0.5, jets[ij], true, true);
        tjsev.j_c2_10_charged[ij] = getC(2, 1.0, jets[ij]);
//        tjsev.j_c2_10_all[ij]     = getC(2, 1.0, jets[ij], true);
//        tjsev.j_c2_10_puppi[ij]   = getC(2, 1.0, jets[ij], true, true);
        tjsev.j_c2_20_charged[ij] = getC(2, 2.0, jets[ij]);
//        tjsev.j_c2_20_all[ij]     = getC(2, 2.0, jets[ij], true);
//        tjsev.j_c2_20_puppi[ij]   = getC(2, 2.0, jets[ij], true, true);
        //39 min without C3
        tjsev.j_c3_02_charged[ij] = getC(3, 0.2, jets[ij]);
//        tjsev.j_c3_02_all[ij]     = getC(3, 0.2, jets[ij], true);
//        tjsev.j_c3_02_puppi[ij]   = getC(3, 0.2, jets[ij], true, true);
        tjsev.j_c3_05_charged[ij] = getC(3, 0.5, jets[ij]);
//        tjsev.j_c3_05_all[ij]     = getC(3, 0.5, jets[ij], true);
//        tjsev.j_c3_05_puppi[ij]   = getC(3, 0.5, jets[ij], true, true);
        tjsev.j_c3_10_charged[ij] = getC(3, 1.0, jets[ij]);
//        tjsev.j_c3_10_all[ij]     = getC(3, 1.0, jets[ij], true);
//        tjsev.j_c3_10_puppi[ij]   = getC(3, 1.0, jets[ij], true, true);
        tjsev.j_c3_20_charged[ij] = getC(3, 2.0, jets[ij]);
//        tjsev.j_c3_20_all[ij]     = getC(3, 2.0, jets[ij], true);
//        tjsev.j_c3_20_puppi[ij]   = getC(3, 2.0, jets[ij], true, true);
        
        ht.fill("js_mult_charged", tjsev.j_mult_charged[ij], plotwgts);
//        ht.fill("js_mult_puppi", tjsev.j_mult_puppi[ij], plotwgts);
//        ht.fill("js_mult_all", tjsev.j_mult_all[ij], plotwgts);
        ht.fill("js_width_charged", tjsev.j_width_charged[ij], plotwgts);
//        ht.fill("js_width_puppi", tjsev.j_width_puppi[ij], plotwgts);
//        ht.fill("js_width_all", tjsev.j_width_all[ij], plotwgts);
        ht.fill("js_ptd_charged", tjsev.j_ptd_charged[ij], plotwgts);
//        ht.fill("js_ptd_puppi", tjsev.j_ptd_puppi[ij], plotwgts);
//        ht.fill("js_ptd_all", tjsev.j_ptd_all[ij], plotwgts);
        ht.fill("js_ptds_charged", tjsev.j_ptds_charged[ij], plotwgts);
//        ht.fill("js_ptds_puppi", tjsev.j_ptds_puppi[ij], plotwgts);
//        ht.fill("js_ptds_all", tjsev.j_ptds_all[ij], plotwgts);
        ht.fill("js_ecc_charged", tjsev.j_ecc_charged[ij], plotwgts);
//        ht.fill("js_ecc_puppi", tjsev.j_ecc_puppi[ij], plotwgts);
//        ht.fill("js_ecc_all", tjsev.j_ecc_all[ij], plotwgts);
        ht.fill("js_tau21_charged", tjsev.j_tau21_charged[ij], plotwgts);
//        ht.fill("js_tau21_puppi", tjsev.j_tau21_puppi[ij], plotwgts);
//        ht.fill("js_tau21_all", tjsev.j_tau21_all[ij], plotwgts);
        ht.fill("js_tau32_charged", tjsev.j_tau32_charged[ij], plotwgts);
//        ht.fill("js_tau32_puppi", tjsev.j_tau32_puppi[ij], plotwgts);
//        ht.fill("js_tau32_all", tjsev.j_tau32_all[ij], plotwgts);
        ht.fill("js_tau43_charged", tjsev.j_tau43_charged[ij], plotwgts);
//        ht.fill("js_tau43_puppi", tjsev.j_tau43_puppi[ij], plotwgts);
//        ht.fill("js_tau43_all", tjsev.j_tau43_all[ij], plotwgts);
        ht.fill("js_zg_charged", tjsev.j_zg_charged[ij], plotwgts);
//        ht.fill("js_zg_puppi", tjsev.j_zg_puppi[ij], plotwgts);
//        ht.fill("js_zg_all", tjsev.j_zg_all[ij], plotwgts);
        ht.fill("js_zgxdr_charged", tjsev.j_zgxdr_charged[ij], plotwgts);
//        ht.fill("js_zgxdr_puppi", tjsev.j_zgxdr_puppi[ij], plotwgts);
//        ht.fill("js_zgxdr_all", tjsev.j_zgxdr_all[ij], plotwgts);
        ht.fill("js_zgdr_charged", tjsev.j_zgdr_charged[ij], plotwgts);
//        ht.fill("js_zgdr_puppi", tjsev.j_zgdr_puppi[ij], plotwgts);
//        ht.fill("js_zgdr_all", tjsev.j_zgdr_all[ij], plotwgts);
        ht.fill("js_ga_width_charged", tjsev.j_ga_width_charged[ij], plotwgts);
//        ht.fill("js_ga_width_puppi", tjsev.j_ga_width_puppi[ij], plotwgts);
//        ht.fill("js_ga_width_all", tjsev.j_ga_width_all[ij], plotwgts);
        ht.fill("js_ga_lha_charged", tjsev.j_ga_lha_charged[ij], plotwgts);
//        ht.fill("js_ga_lha_puppi", tjsev.j_ga_lha_puppi[ij], plotwgts);
//        ht.fill("js_ga_lha_all", tjsev.j_ga_lha_all[ij], plotwgts);
        ht.fill("js_ga_thrust_charged", tjsev.j_ga_thrust_charged[ij], plotwgts);
//        ht.fill("js_ga_thrust_puppi", tjsev.j_ga_thrust_puppi[ij], plotwgts);
//        ht.fill("js_ga_thrust_all", tjsev.j_ga_thrust_all[ij], plotwgts);
        ht.fill("js_c1_02_charged", tjsev.j_c1_02_charged[ij], plotwgts);
//        ht.fill("js_c1_02_puppi", tjsev.j_c1_02_puppi[ij], plotwgts);
//        ht.fill("js_c1_02_all", tjsev.j_c1_02_all[ij], plotwgts);
        ht.fill("js_c1_05_charged", tjsev.j_c1_05_charged[ij], plotwgts);
//        ht.fill("js_c1_05_puppi", tjsev.j_c1_05_puppi[ij], plotwgts);
//        ht.fill("js_c1_05_all", tjsev.j_c1_05_all[ij], plotwgts);
        ht.fill("js_c1_10_charged", tjsev.j_c1_10_charged[ij], plotwgts);
//        ht.fill("js_c1_10_puppi", tjsev.j_c1_10_puppi[ij], plotwgts);
//        ht.fill("js_c1_10_all", tjsev.j_c1_10_all[ij], plotwgts);
        ht.fill("js_c1_20_charged", tjsev.j_c1_20_charged[ij], plotwgts);
//        ht.fill("js_c1_20_puppi", tjsev.j_c1_20_puppi[ij], plotwgts);
//        ht.fill("js_c1_20_all", tjsev.j_c1_20_all[ij], plotwgts);
        ht.fill("js_c2_02_charged", tjsev.j_c2_02_charged[ij], plotwgts);
//        ht.fill("js_c2_02_puppi", tjsev.j_c2_02_puppi[ij], plotwgts);
//        ht.fill("js_c2_02_all", tjsev.j_c2_02_all[ij], plotwgts);
        ht.fill("js_c2_05_charged", tjsev.j_c2_05_charged[ij], plotwgts);
//        ht.fill("js_c2_05_puppi", tjsev.j_c2_05_puppi[ij], plotwgts);
//        ht.fill("js_c2_05_all", tjsev.j_c2_05_all[ij], plotwgts);
        ht.fill("js_c2_10_charged", tjsev.j_c2_10_charged[ij], plotwgts);
//        ht.fill("js_c2_10_puppi", tjsev.j_c2_10_puppi[ij], plotwgts);
//        ht.fill("js_c2_10_all", tjsev.j_c2_10_all[ij], plotwgts);
        ht.fill("js_c2_20_charged", tjsev.j_c2_20_charged[ij], plotwgts);
//        ht.fill("js_c2_20_puppi", tjsev.j_c2_20_puppi[ij], plotwgts);
//        ht.fill("js_c2_20_all", tjsev.j_c2_20_all[ij], plotwgts);
        ht.fill("js_c3_02_charged", tjsev.j_c3_02_charged[ij], plotwgts);
//        ht.fill("js_c3_02_puppi", tjsev.j_c3_02_puppi[ij], plotwgts);
//        ht.fill("js_c3_02_all", tjsev.j_c3_02_all[ij], plotwgts);
        ht.fill("js_c3_05_charged", tjsev.j_c3_05_charged[ij], plotwgts);
//        ht.fill("js_c3_05_puppi", tjsev.j_c3_05_puppi[ij], plotwgts);
//        ht.fill("js_c3_05_all", tjsev.j_c3_05_all[ij], plotwgts);
        ht.fill("js_c3_10_charged", tjsev.j_c3_10_charged[ij], plotwgts);
//        ht.fill("js_c3_10_puppi", tjsev.j_c3_10_puppi[ij], plotwgts);
//        ht.fill("js_c3_10_all", tjsev.j_c3_10_all[ij], plotwgts);
        ht.fill("js_c3_20_charged", tjsev.j_c3_20_charged[ij], plotwgts);
//        ht.fill("js_c3_20_puppi", tjsev.j_c3_20_puppi[ij], plotwgts);
//        ht.fill("js_c3_20_all", tjsev.j_c3_20_all[ij], plotwgts);
      }
      
      
      ///////////////////////
      // GENERATOR LEVEL  //
      /////////////////////
      
      if (isTTbar) {
        //////////////////////////
        // GEN LEVEL SELECTION //
        ////////////////////////
       
        
        //decide the lepton channel at particle level
        std::vector<Particle> genVetoLeptons = selector.getGenLeptons(ev,15.,2.4);
        std::vector<Particle> genLeptons     = selector.getGenLeptons(ev,30.,2.1);
        TString genChTag = selector.flagGenFinalState(ev, genLeptons);
        std::vector<Jet> genJets = selector.getGenJets();
        
        //count b and W candidates
        int sel_ngbjets = 0;
        int sel_ngwcand = 0;

        for (auto& jet : genJets) {
          if (jet.flavor() == 5) ++sel_ngbjets;
          if (jet.flavor() == 1) ++sel_ngwcand;
        }
        
        //event selected on gen level?
        bool genSingleLepton((genChTag=="E" or genChTag=="M") and
                             (genVetoLeptons.size() == 1)); // only selected lepton in veto collection
        if (sel_ngbjets==2 && sel_ngwcand>0 && genSingleLepton) tjsev.gen_sel = 1;
        
        tjsev.ngj = genJets.size();
            
        /////////////////////////
        // GEN LEVEL ANALYSIS //
        ///////////////////////

        //store jets to tree
        for (int i = 0; i < tjsev.ngj; i++) {
          tjsev.gj_pt     [i] = genJets[i].p4().Pt();
          tjsev.gj_eta    [i] = genJets[i].p4().Eta();
          tjsev.gj_phi    [i] = genJets[i].p4().Phi();
          tjsev.gj_m      [i] = genJets[i].p4().M();
          tjsev.gj_flavor [i] = genJets[i].flavor();
          tjsev.gj_overlap[i] = genJets[i].overlap();
          
          //matching to reco jet
          for(unsigned int ij = 0; ij< jets.size(); ij++) {
            int ig = i;
            if(jets[ij].p4().DeltaR(genJets[ig].p4())>0.4) continue;
            tjsev.j_gj[ij] = ig;
            tjsev.gj_j[ig] = ij;
            break;
          }
          
          if (tjsev.gen_sel != 1) continue;
          if (abs(tjsev.gj_eta[i]) > 2. or tjsev.gj_overlap[i]) continue;
          
          //calculate jet properties            
          tjsev.gj_mult_charged[i] = getMult(genJets[i]);
//          tjsev.gj_mult_all[i]     = getMult(genJets[i], true);
//          tjsev.gj_mult_puppi[i]   = getMult(genJets[i], true, true);
          
          tjsev.gj_ptd_charged[i] = getPtD(genJets[i]);
//          tjsev.gj_ptd_all[i]     = getPtD(genJets[i], true);
//          tjsev.gj_ptd_puppi[i]   = getPtD(genJets[i], true, true);
          
          tjsev.gj_ptds_charged[i] = getPtDs(genJets[i]);
//          tjsev.gj_ptds_all[i]     = getPtDs(genJets[i], true);
//          tjsev.gj_ptds_puppi[i]   = getPtDs(genJets[i], true, true);
          
          tjsev.gj_width_charged[i] = getWidth(genJets[i]);
//          tjsev.gj_width_all[i]     = getWidth(genJets[i], true);
//          tjsev.gj_width_puppi[i]   = getWidth(genJets[i], true, true);
          
          tjsev.gj_ecc_charged[i] = getEcc(genJets[i]);
//          tjsev.gj_ecc_all[i]     = getEcc(genJets[i], true);
//          tjsev.gj_ecc_puppi[i]   = getEcc(genJets[i], true, true);
          
          tjsev.gj_tau21_charged[i] = getTau(2, 1, genJets[i]);
//          tjsev.gj_tau21_all[i]     = getTau(2, 1, genJets[i], true);
//          tjsev.gj_tau21_puppi[i]   = getTau(2, 1, genJets[i], true, true);
          
          tjsev.gj_tau32_charged[i] = getTau(3, 2, genJets[i]);
//          tjsev.gj_tau32_all[i]     = getTau(3, 2, genJets[i], true);
//          tjsev.gj_tau32_puppi[i]   = getTau(3, 2, genJets[i], true, true);
          
          tjsev.gj_tau43_charged[i] = getTau(4, 3, genJets[i]);
//          tjsev.gj_tau43_all[i]     = getTau(4, 3, genJets[i], true);
//          tjsev.gj_tau43_puppi[i]   = getTau(4, 3, genJets[i], true, true);
          
          std::vector<double> zgResult_charged = getZg(genJets[i]);
//          std::vector<double> zgResult_all     = getZg(genJets[i], true);
//          std::vector<double> zgResult_puppi   = getZg(genJets[i], true, true);
          
          tjsev.gj_zg_charged[i]    = zgResult_charged[0];
//          tjsev.gj_zg_all[i]        = zgResult_all[0];
//          tjsev.gj_zg_puppi[i]      = zgResult_puppi[0];
          
          tjsev.gj_zgxdr_charged[i] = zgResult_charged[1];
//          tjsev.gj_zgxdr_all[i]     = zgResult_all[1];
//          tjsev.gj_zgxdr_puppi[i]   = zgResult_puppi[1];
          
          tjsev.gj_zgdr_charged[i]  = zgResult_charged[2];
//          tjsev.gj_zgdr_all[i]      = zgResult_all[2];
//          tjsev.gj_zgdr_puppi[i]    = zgResult_puppi[2];
          
          tjsev.gj_ga_width_charged[i]  = calcGA(1., 1., genJets[i]);
//          tjsev.gj_ga_width_all[i]      = calcGA(1., 1., genJets[i], true);
//          tjsev.gj_ga_width_puppi[i]    = calcGA(1., 1., genJets[i], true, true);
          
          tjsev.gj_ga_lha_charged[i]    = calcGA(0.5, 1., genJets[i]);
//          tjsev.gj_ga_lha_all[i]        = calcGA(0.5, 1., genJets[i], true);
//          tjsev.gj_ga_lha_puppi[i]      = calcGA(0.5, 1., genJets[i], true, true);
          
          tjsev.gj_ga_thrust_charged[i] = calcGA(2., 1., genJets[i]);
//          tjsev.gj_ga_thrust_all[i]     = calcGA(2., 1., genJets[i], true);
//          tjsev.gj_ga_thrust_puppi[i]   = calcGA(2., 1., genJets[i], true, true);
          
          tjsev.gj_c1_02_charged[i] = getC(1, 0.2, genJets[i]);
//          tjsev.gj_c1_02_all[i]     = getC(1, 0.2, genJets[i], true);
//          tjsev.gj_c1_02_puppi[i]   = getC(1, 0.2, genJets[i], true, true);
          tjsev.gj_c1_05_charged[i] = getC(1, 0.5, genJets[i]);
//          tjsev.gj_c1_05_all[i]     = getC(1, 0.5, genJets[i], true);
//          tjsev.gj_c1_05_puppi[i]   = getC(1, 0.5, genJets[i], true, true);
          tjsev.gj_c1_10_charged[i] = getC(1, 1.0, genJets[i]);
//          tjsev.gj_c1_10_all[i]     = getC(1, 1.0, genJets[i], true);
//          tjsev.gj_c1_10_puppi[i]   = getC(1, 1.0, genJets[i], true, true);
          tjsev.gj_c1_20_charged[i] = getC(1, 2.0, genJets[i]);
//          tjsev.gj_c1_20_all[i]     = getC(1, 2.0, genJets[i], true);
//          tjsev.gj_c1_20_puppi[i]   = getC(1, 2.0, genJets[i], true, true);
          tjsev.gj_c2_02_charged[i] = getC(2, 0.2, genJets[i]);
//          tjsev.gj_c2_02_all[i]     = getC(2, 0.2, genJets[i], true);
//          tjsev.gj_c2_02_puppi[i]   = getC(2, 0.2, genJets[i], true, true);
          tjsev.gj_c2_05_charged[i] = getC(2, 0.5, genJets[i]);
//          tjsev.gj_c2_05_all[i]     = getC(2, 0.5, genJets[i], true);
//          tjsev.gj_c2_05_puppi[i]   = getC(2, 0.5, genJets[i], true, true);
          tjsev.gj_c2_10_charged[i] = getC(2, 1.0, genJets[i]);
//          tjsev.gj_c2_10_all[i]     = getC(2, 1.0, genJets[i], true);
//          tjsev.gj_c2_10_puppi[i]   = getC(2, 1.0, genJets[i], true, true);
          tjsev.gj_c2_20_charged[i] = getC(2, 2.0, genJets[i]);
//          tjsev.gj_c2_20_all[i]     = getC(2, 2.0, genJets[i], true);
//          tjsev.gj_c2_20_puppi[i]   = getC(2, 2.0, genJets[i], true, true);
          tjsev.gj_c3_02_charged[i] = getC(3, 0.2, genJets[i]);
//          tjsev.gj_c3_02_all[i]     = getC(3, 0.2, genJets[i], true);
//          tjsev.gj_c3_02_puppi[i]   = getC(3, 0.2, genJets[i], true, true);
          tjsev.gj_c3_05_charged[i] = getC(3, 0.5, genJets[i]);
//          tjsev.gj_c3_05_all[i]     = getC(3, 0.5, genJets[i], true);
//          tjsev.gj_c3_05_puppi[i]   = getC(3, 0.5, genJets[i], true, true);
          tjsev.gj_c3_10_charged[i] = getC(3, 1.0, genJets[i]);
//          tjsev.gj_c3_10_all[i]     = getC(3, 1.0, genJets[i], true);
//          tjsev.gj_c3_10_puppi[i]   = getC(3, 1.0, genJets[i], true, true);
          tjsev.gj_c3_20_charged[i] = getC(3, 2.0, genJets[i]);
//          tjsev.gj_c3_20_all[i]     = getC(3, 2.0, genJets[i], true);
//          tjsev.gj_c3_20_puppi[i]   = getC(3, 2.0, genJets[i], true, true);
        }
      }
      else tjsev.gen_sel = -1;
      
      //proceed only if event is selected on gen or reco level
      if (tjsev.gen_sel + tjsev.reco_sel == -2) continue;
      
      outT->Fill();
    }
  
  //close input file
  f->Close();

  //save histos to file  
  fOut->cd();
  outT->Write();
  for (auto& it : ht.getPlots())  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  for (auto& it : ht.get2dPlots())  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  for (auto& it : allPlots)  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  for (auto& it : all2dPlots)  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  fOut->Close();
}

//Calculate generalized angularities
/*
  Definition from arXiv:1408.3122
  (beta,kappa) =
    (0,0) -> multiplicity
    (0,2) -> pt dispersion
    (1,1) -> broadening/width/girth
    (2,1) -> thrust (m^2/e)
*/
double calcGA(double beta, double kappa, Jet jet, bool includeNeutrals, bool usePuppi, double ptcut) {
  int mult = 0;
  double sumpt = 0.;
  TLorentzVector axis(0., 0., 0., 0.);
  for (auto p : jet.particles()) {
    if (not includeNeutrals and p.charge() == 0) continue;
    if (ptcut > p.pt()) continue;
    double weight = usePuppi ? p.puppi() : 1.;
    if (weight > 0.) ++mult;
    sumpt  += p.pt()*weight;
    axis += p.momentum()*weight;
  }
  if (mult < 2) return -1.;
  //std::cout << "sumpt" << beta << kappa << iptcut << icharge << ": " << sumpt << std::endl;
  
  double ga = 0.;
  for (auto p : jet.particles()) {
    if (not includeNeutrals and p.charge() == 0) continue;
    if (ptcut > p.pt()) continue;
    double weight = usePuppi ? p.puppi() : 1.;
    ga += weight * pow(p.p4().Pt()/sumpt, kappa) * pow(deltaR(axis, p.momentum())/0.4, beta);
  }
  
  //std::cout << "ga" << beta << kappa << iptcut << icharge << ": " << ga << std::endl;
  
  /*
  if (ga > 1. && beta == 1 && kappa == 1 && iptcut == 0 && icharge == 0) {
    std::cout << "ga(1,1) = " << ga << std::endl;
    std::cout << "pt/sum(pt) deltaR" << std::endl;
    double sumpt_test = 0;
    for (auto p : jet.particles()) {
      if (p.p4().Pt() < (iptcut+1)*0.500) continue;
      if(p.charge()==0) continue;
      std::cout << p.p4().Pt()/sumpt << " " << jet.p4().DeltaR(p.p4()) << std::endl;
      sumpt_test += p.p4().Pt()/sumpt;
    }
    std::cout << "sumpt_test = " << sumpt_test << "\n" << std::endl;
  }
  */
  
  return ga;
}

double getMult(Jet jet, bool includeNeutrals, bool usePuppi, double ptcut) {
  //std::cout << "getMult()" << std::endl;
  double sumWeight = 0.;
  for (auto p : jet.particles()) {
    if (not includeNeutrals and p.charge() == 0) continue;
    if (ptcut > p.pt()) continue;
    double weight = usePuppi ? p.puppi() : 1.;
    sumWeight += weight;
  }
  return sumWeight;
}

double getPtD(Jet jet, bool includeNeutrals, bool usePuppi, double ptcut) {
  //std::cout << "getPtD()" << std::endl;
  int mult = 0;
  double sumpt  = 0.;
  double sumpt2 = 0.;
  for (auto p : jet.particles()) {
    if (not includeNeutrals and p.charge() == 0) continue;
    if (ptcut > p.pt()) continue;
    double weight = usePuppi ? p.puppi() : 1.;
    if (weight > 0.) ++mult;
    sumpt  += p.pt()*weight;
    sumpt2 += pow(p.pt()*weight, 2);
  }
  if (mult < 2) return -1.;
  double ptd = sumpt2/pow(sumpt,2);
  return ptd;
}

double getPtDs(Jet jet, bool includeNeutrals, bool usePuppi, double ptcut) {
  //std::cout << "getPtDs()" << std::endl;
  double mult   = 0.;
  double sumpt  = 0.;
  double sumpt2 = 0.;
  for (auto p : jet.particles()) {
    if (not includeNeutrals and p.charge() == 0) continue;
    if (ptcut > p.pt()) continue;
    double weight = usePuppi ? p.puppi() : 1.;
    if (weight > 0.) mult+=1.;
    sumpt  += p.pt()*weight;
    sumpt2 += pow(p.pt()*weight, 2);
  }
  if (mult < 2.) return -1.;
  double ptd = sumpt2/pow(sumpt,2);
  return max(0., sqrt((ptd-1./mult) * mult/(mult-1.)));
}

double getWidth(Jet jet, bool includeNeutrals, bool usePuppi, double ptcut) {
  //std::cout << "getWidth()" << std::endl;
  int mult = 0;
  double sumpt   = 0.;
  double sumptdr = 0.;
  TLorentzVector axis(0., 0., 0., 0.);
  for (auto p : jet.particles()) {
    if (not includeNeutrals and p.charge() == 0) continue;
    if (ptcut > p.pt()) continue;
    double weight = usePuppi ? p.puppi() : 1.;
    if (weight > 0.) ++mult;
    axis += p.momentum()*weight;
  }
  if (mult < 2) return -1.;
  for (auto p : jet.particles()) {
    if (not includeNeutrals and p.charge() == 0) continue;
    if (ptcut > p.pt()) continue;
    double weight = usePuppi ? p.puppi() : 1.;
    sumpt   += p.pt()*weight;
    sumptdr += p.pt()*weight * deltaR(axis, p.momentum());
  }
  double width = sumptdr/sumpt;
  return width;
}

double getEcc(Jet jet, bool includeNeutrals, bool usePuppi, double ptcut) {
  //std::cout << "getEcc()" << std::endl;
  // Get mean axis
  int mult = 0;
  TLorentzVector axis(0., 0., 0., 0.);
  for (auto p : jet.particles()) {
    if (not includeNeutrals and p.charge() == 0) continue;
    if (ptcut > p.pt()) continue;
    double weight = usePuppi ? p.puppi() : 1.;
    if (weight > 0.) ++mult;
    axis += p.momentum()*weight;
  }
  if (mult < 4) return -1.;
  // Covariance matrix
  Matrix<2> M;
  for (auto p : jet.particles()) {
    if (not includeNeutrals and p.charge() == 0) continue;
    if (ptcut > p.pt()) continue;
    double weight = usePuppi ? p.puppi() : 1.;
    Matrix<2> MPart;
    MPart.set(0, 0, (p.eta() - axis.Eta()) * (p.eta() - axis.Eta()));
    MPart.set(0, 1, (p.eta() - axis.Eta()) * mapAngleMPiToPi(p.phi() - axis.Phi()));
    MPart.set(1, 0, mapAngleMPiToPi(p.phi() - axis.Phi()) * (p.eta() - axis.Eta()));
    MPart.set(1, 1, mapAngleMPiToPi(p.phi() - axis.Phi()) * mapAngleMPiToPi(p.phi() - axis.Phi()));
    M += MPart * p.energy() * weight;
  }
  // Calculate eccentricity from eigenvalues
  const EigenSystem<2> eigen = diagonalize(M);
  double ecc = 1. - eigen.getEigenValues()[1]/eigen.getEigenValues()[0];
  
  return ecc;
}

double getTau(int N, int M, Jet jet, bool includeNeutrals, bool usePuppi, double ptcut) {
  // Recluster constituents with CA
  int mult = 0.;
  vector<PseudoJet> particles;
  for (auto p : jet.particles()) {
    if (not includeNeutrals and p.charge() == 0) continue;
    if (ptcut > p.pt()) continue;
    double weight = usePuppi ? p.puppi() : 1.;
    if (weight > 0.) ++mult;
    particles.push_back( PseudoJet(p.px(), p.py(), p.pz(), p.e())*weight );
  }
  if (mult < N+2) return -1.;
  
  JetDefinition jet_def(fastjet::cambridge_algorithm, fastjet::JetDefinition::max_allowable_R);
  
  ClusterSequence cs(particles, jet_def);
  vector<PseudoJet> jets = sorted_by_pt(cs.exclusive_jets(1));
  
  PseudoJet cajet = jets[0];
  
  NsubjettinessRatio tau_ratio(N, M, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(1.0));
  
  return tau_ratio(cajet);
}

double getC(int N, double beta, Jet jet, bool includeNeutrals, bool usePuppi, double ptcut) {
  // Recluster constituents with CA
  int mult = 0.;
  vector<PseudoJet> particles;
  for (auto p : jet.particles()) {
    if (not includeNeutrals and p.charge() == 0) continue;
    if (ptcut > p.pt()) continue;
    double weight = usePuppi ? p.puppi() : 1.;
    if (weight > 0.) ++mult;
    particles.push_back( PseudoJet(p.px(), p.py(), p.pz(), p.e())*weight );
  }
  if (mult < N+2) return -1.;
  
  JetDefinition jet_def(fastjet::cambridge_algorithm, fastjet::JetDefinition::max_allowable_R);
  
  ClusterSequence cs(particles, jet_def);
  vector<PseudoJet> jets = sorted_by_pt(cs.exclusive_jets(1));
  
  PseudoJet cajet = jets[0];
  
  EnergyCorrelatorDoubleRatio C(N, beta);
  
  return C(cajet);
}

std::vector<double> getZg(Jet jet, bool includeNeutrals, bool usePuppi, double ptcut) {
  //std::cout << "getZg()" << std::endl;
  // Recluster constituents with CA
  int mult = 0.;
  vector<PseudoJet> particles;
  for (auto p : jet.particles()) {
    if (not includeNeutrals and p.charge() == 0) continue;
    if (ptcut > p.pt()) continue;
    double weight = usePuppi ? p.puppi() : 1.;
    if (weight > 0.) ++mult;
    particles.push_back( PseudoJet(p.px(), p.py(), p.pz(), p.e())*weight );
  }
  if (mult < 2) return {-1., -1., -1.};
  
  JetDefinition jet_def(fastjet::cambridge_algorithm, fastjet::JetDefinition::max_allowable_R);
  
  ClusterSequence cs(particles, jet_def);
  vector<PseudoJet> jets = sorted_by_pt(cs.exclusive_jets(1));
  
  PseudoJet cajet = jets[0];
  PseudoJet cajetp1, cajetp2;
  double zg = 0.;
  while (zg < 0.1 and cajet.has_parents(cajetp1, cajetp2)) {
    zg    = cajetp2.pt()/cajet.pt();
    cajet = cajetp1;
  }
  if (zg < 0.1) return {-1., -1., -1.};
  //std::cout << "zg = " << zg << std::endl;
  std::vector<double> results;
  results.push_back(zg);
  results.push_back(zg*cajetp1.delta_R(cajetp2));
  results.push_back(cajetp1.delta_R(cajetp2));
  return results;
}

double deltaR(TLorentzVector v1, TLorentzVector v2) { return v1.DeltaR(v2); }
double mapAngleMPiToPi(double phi) { return TVector2::Phi_mpi_pi(phi); }

//
void createTopJetShapeEventTree(TTree *t,TopJetShapeEvent_t &tjsev)
{
  //event weights
  t->Branch("nw",  &tjsev.nw, "nw/I");
  t->Branch("weight",  tjsev.weight, "weight[nw]/F");

  //met
  t->Branch("met_pt",  &tjsev.met_pt, "met_pt/F");
  t->Branch("met_phi",  &tjsev.met_phi, "met_phi/F");

  //leptons
  t->Branch("nl",  &tjsev.nl, "nl/I");
  t->Branch("l_pt",  tjsev.l_pt ,  "l_pt[nl]/F");
  t->Branch("l_eta", tjsev.l_eta , "l_eta[nl]/F");
  t->Branch("l_phi", tjsev.l_phi , "l_phi[nl]/F");
  t->Branch("l_m",   tjsev.l_m ,   "l_m[nl]/F");
  t->Branch("l_id",   tjsev.l_id ,   "l_id[nl]/I");
  
  //gen leptons
  t->Branch("ngl",  &tjsev.ngl, "ngl/I");
  t->Branch("gl_pt",  tjsev.gl_pt ,  "gl_pt[ngl]/F");
  t->Branch("gl_eta", tjsev.gl_eta , "gl_eta[ngl]/F");
  t->Branch("gl_phi", tjsev.gl_phi , "gl_phi[ngl]/F");
  t->Branch("gl_m",   tjsev.gl_m ,   "gl_m[ngl]/F");
  t->Branch("gl_id",  tjsev.gl_id ,  "gl_id[ngl]/I");

  //jets
  t->Branch("nj",  &tjsev.nj, "nj/I");
  t->Branch("j_pt",  tjsev.j_pt ,  "j_pt[nj]/F");
  t->Branch("j_eta", tjsev.j_eta , "j_eta[nj]/F");
  t->Branch("j_phi", tjsev.j_phi , "j_phi[nj]/F");
  t->Branch("j_m",   tjsev.j_m ,   "j_m[nj]/F");
  t->Branch("j_flavor",  tjsev.j_flavor ,  "j_flavor[nj]/I");
  t->Branch("j_overlap",  tjsev.j_overlap ,  "j_overlap[nj]/I");
  t->Branch("j_gj",  tjsev.j_gj ,  "j_gj[nj]/I");
  
  //gen jets
  t->Branch("ngj",  &tjsev.ngj, "ngj/I");
  t->Branch("gj_pt",  tjsev.gj_pt ,  "gj_pt[ngj]/F");
  t->Branch("gj_eta", tjsev.gj_eta , "gj_eta[ngj]/F");
  t->Branch("gj_phi", tjsev.gj_phi , "gj_phi[ngj]/F");
  t->Branch("gj_m",   tjsev.gj_m ,   "gj_m[ngj]/F");
  t->Branch("gj_flavor",  tjsev.gj_flavor ,  "gj_flavor[ngj]/I");
  t->Branch("gj_overlap",  tjsev.gj_overlap ,  "gj_overlap[ngj]/I");
  t->Branch("gj_j",  tjsev.gj_j ,  "gj_j[ngj]/I");
  
  t->Branch("gen_sel", &tjsev.gen_sel ,  "gen_sel/I");
  t->Branch("reco_sel", &tjsev.reco_sel ,  "reco_sel/I");
  
  //observables
  t->Branch("j_mult_charged",   tjsev.j_mult_charged,   "j_mult_charged[nj]/F");
  t->Branch("j_mult_puppi",   tjsev.j_mult_puppi,   "j_mult_puppi[nj]/F");
  t->Branch("j_mult_all",   tjsev.j_mult_all,   "j_mult_all[nj]/F");
  t->Branch("j_width_charged",   tjsev.j_width_charged,   "j_width_charged[nj]/F");
  t->Branch("j_width_puppi",   tjsev.j_width_puppi,   "j_width_puppi[nj]/F");
  t->Branch("j_width_all",   tjsev.j_width_all,   "j_width_all[nj]/F");
  t->Branch("j_ptd_charged",   tjsev.j_ptd_charged,   "j_ptd_charged[nj]/F");
  t->Branch("j_ptd_puppi",   tjsev.j_ptd_puppi,   "j_ptd_puppi[nj]/F");
  t->Branch("j_ptd_all",   tjsev.j_ptd_all,   "j_ptd_all[nj]/F");
  t->Branch("j_ptds_charged",   tjsev.j_ptds_charged,   "j_ptds_charged[nj]/F");
  t->Branch("j_ptds_puppi",   tjsev.j_ptds_puppi,   "j_ptds_puppi[nj]/F");
  t->Branch("j_ptds_all",   tjsev.j_ptds_all,   "j_ptds_all[nj]/F");
  t->Branch("j_ecc_charged",   tjsev.j_ecc_charged,   "j_ecc_charged[nj]/F");
  t->Branch("j_ecc_puppi",   tjsev.j_ecc_puppi,   "j_ecc_puppi[nj]/F");
  t->Branch("j_ecc_all",   tjsev.j_ecc_all,   "j_ecc_all[nj]/F");
  t->Branch("j_tau21_charged",   tjsev.j_tau21_charged,   "j_tau21_charged[nj]/F");
  t->Branch("j_tau21_puppi",   tjsev.j_tau21_puppi,   "j_tau21_puppi[nj]/F");
  t->Branch("j_tau21_all",   tjsev.j_tau21_all,   "j_tau21_all[nj]/F");
  t->Branch("j_tau32_charged",   tjsev.j_tau32_charged,   "j_tau32_charged[nj]/F");
  t->Branch("j_tau32_puppi",   tjsev.j_tau32_puppi,   "j_tau32_puppi[nj]/F");
  t->Branch("j_tau32_all",   tjsev.j_tau32_all,   "j_tau32_all[nj]/F");
  t->Branch("j_tau43_charged",   tjsev.j_tau43_charged,   "j_tau43_charged[nj]/F");
  t->Branch("j_tau43_puppi",   tjsev.j_tau43_puppi,   "j_tau43_puppi[nj]/F");
  t->Branch("j_tau43_all",   tjsev.j_tau43_all,   "j_tau43_all[nj]/F");
  t->Branch("j_zg_charged",   tjsev.j_zg_charged,   "j_zg_charged[nj]/F");
  t->Branch("j_zg_puppi",   tjsev.j_zg_puppi,   "j_zg_puppi[nj]/F");
  t->Branch("j_zg_all",   tjsev.j_zg_all,   "j_zg_all[nj]/F");
  t->Branch("j_zgxdr_charged",   tjsev.j_zgxdr_charged,   "j_zgxdr_charged[nj]/F");
  t->Branch("j_zgxdr_puppi",   tjsev.j_zgxdr_puppi,   "j_zgxdr_puppi[nj]/F");
  t->Branch("j_zgxdr_all",   tjsev.j_zgxdr_all,   "j_zgxdr_all[nj]/F");
  t->Branch("j_zgdr_charged",   tjsev.j_zgdr_charged,   "j_zgdr_charged[nj]/F");
  t->Branch("j_zgdr_puppi",   tjsev.j_zgdr_puppi,   "j_zgdr_puppi[nj]/F");
  t->Branch("j_zgdr_all",   tjsev.j_zgdr_all,   "j_zgdr_all[nj]/F");
  t->Branch("j_ga_width_charged",   tjsev.j_ga_width_charged,   "j_ga_width_charged[nj]/F");
  t->Branch("j_ga_width_puppi",   tjsev.j_ga_width_puppi,   "j_ga_width_puppi[nj]/F");
  t->Branch("j_ga_width_all",   tjsev.j_ga_width_all,   "j_ga_width_all[nj]/F");
  t->Branch("j_ga_lha_charged",   tjsev.j_ga_lha_charged,   "j_ga_lha_charged[nj]/F");
  t->Branch("j_ga_lha_puppi",   tjsev.j_ga_lha_puppi,   "j_ga_lha_puppi[nj]/F");
  t->Branch("j_ga_lha_all",   tjsev.j_ga_lha_all,   "j_ga_lha_all[nj]/F");
  t->Branch("j_ga_thrust_charged",   tjsev.j_ga_thrust_charged,   "j_ga_thrust_charged[nj]/F");
  t->Branch("j_ga_thrust_puppi",   tjsev.j_ga_thrust_puppi,   "j_ga_thrust_puppi[nj]/F");
  t->Branch("j_ga_thrust_all",   tjsev.j_ga_thrust_all,   "j_ga_thrust_all[nj]/F");
  t->Branch("j_c1_02_charged",   tjsev.j_c1_02_charged,   "j_c1_02_charged[nj]/F");
  t->Branch("j_c1_02_puppi",   tjsev.j_c1_02_puppi,   "j_c1_02_puppi[nj]/F");
  t->Branch("j_c1_02_all",   tjsev.j_c1_02_all,   "j_c1_02_all[nj]/F");
  t->Branch("j_c1_05_charged",   tjsev.j_c1_05_charged,   "j_c1_05_charged[nj]/F");
  t->Branch("j_c1_05_puppi",   tjsev.j_c1_05_puppi,   "j_c1_05_puppi[nj]/F");
  t->Branch("j_c1_05_all",   tjsev.j_c1_05_all,   "j_c1_05_all[nj]/F");
  t->Branch("j_c1_10_charged",   tjsev.j_c1_10_charged,   "j_c1_10_charged[nj]/F");
  t->Branch("j_c1_10_puppi",   tjsev.j_c1_10_puppi,   "j_c1_10_puppi[nj]/F");
  t->Branch("j_c1_10_all",   tjsev.j_c1_10_all,   "j_c1_10_all[nj]/F");
  t->Branch("j_c1_20_charged",   tjsev.j_c1_20_charged,   "j_c1_20_charged[nj]/F");
  t->Branch("j_c1_20_puppi",   tjsev.j_c1_20_puppi,   "j_c1_20_puppi[nj]/F");
  t->Branch("j_c1_20_all",   tjsev.j_c1_20_all,   "j_c1_20_all[nj]/F");
  t->Branch("j_c2_02_charged",   tjsev.j_c2_02_charged,   "j_c2_02_charged[nj]/F");
  t->Branch("j_c2_02_puppi",   tjsev.j_c2_02_puppi,   "j_c2_02_puppi[nj]/F");
  t->Branch("j_c2_02_all",   tjsev.j_c2_02_all,   "j_c2_02_all[nj]/F");
  t->Branch("j_c2_05_charged",   tjsev.j_c2_05_charged,   "j_c2_05_charged[nj]/F");
  t->Branch("j_c2_05_puppi",   tjsev.j_c2_05_puppi,   "j_c2_05_puppi[nj]/F");
  t->Branch("j_c2_05_all",   tjsev.j_c2_05_all,   "j_c2_05_all[nj]/F");
  t->Branch("j_c2_10_charged",   tjsev.j_c2_10_charged,   "j_c2_10_charged[nj]/F");
  t->Branch("j_c2_10_puppi",   tjsev.j_c2_10_puppi,   "j_c2_10_puppi[nj]/F");
  t->Branch("j_c2_10_all",   tjsev.j_c2_10_all,   "j_c2_10_all[nj]/F");
  t->Branch("j_c2_20_charged",   tjsev.j_c2_20_charged,   "j_c2_20_charged[nj]/F");
  t->Branch("j_c2_20_puppi",   tjsev.j_c2_20_puppi,   "j_c2_20_puppi[nj]/F");
  t->Branch("j_c2_20_all",   tjsev.j_c2_20_all,   "j_c2_20_all[nj]/F");
  t->Branch("j_c3_02_charged",   tjsev.j_c3_02_charged,   "j_c3_02_charged[nj]/F");
  t->Branch("j_c3_02_puppi",   tjsev.j_c3_02_puppi,   "j_c3_02_puppi[nj]/F");
  t->Branch("j_c3_02_all",   tjsev.j_c3_02_all,   "j_c3_02_all[nj]/F");
  t->Branch("j_c3_05_charged",   tjsev.j_c3_05_charged,   "j_c3_05_charged[nj]/F");
  t->Branch("j_c3_05_puppi",   tjsev.j_c3_05_puppi,   "j_c3_05_puppi[nj]/F");
  t->Branch("j_c3_05_all",   tjsev.j_c3_05_all,   "j_c3_05_all[nj]/F");
  t->Branch("j_c3_10_charged",   tjsev.j_c3_10_charged,   "j_c3_10_charged[nj]/F");
  t->Branch("j_c3_10_puppi",   tjsev.j_c3_10_puppi,   "j_c3_10_puppi[nj]/F");
  t->Branch("j_c3_10_all",   tjsev.j_c3_10_all,   "j_c3_10_all[nj]/F");
  t->Branch("j_c3_20_charged",   tjsev.j_c3_20_charged,   "j_c3_20_charged[nj]/F");
  t->Branch("j_c3_20_puppi",   tjsev.j_c3_20_puppi,   "j_c3_20_puppi[nj]/F");
  t->Branch("j_c3_20_all",   tjsev.j_c3_20_all,   "j_c3_20_all[nj]/F");


  t->Branch("gj_mult_charged",   tjsev.gj_mult_charged,   "gj_mult_charged[ngj]/F");
  t->Branch("gj_mult_puppi",   tjsev.gj_mult_puppi,   "gj_mult_puppi[ngj]/F");
  t->Branch("gj_mult_all",   tjsev.gj_mult_all,   "gj_mult_all[ngj]/F");
  t->Branch("gj_width_charged",   tjsev.gj_width_charged,   "gj_width_charged[ngj]/F");
  t->Branch("gj_width_puppi",   tjsev.gj_width_puppi,   "gj_width_puppi[ngj]/F");
  t->Branch("gj_width_all",   tjsev.gj_width_all,   "gj_width_all[ngj]/F");
  t->Branch("gj_ptd_charged",   tjsev.gj_ptd_charged,   "gj_ptd_charged[ngj]/F");
  t->Branch("gj_ptd_puppi",   tjsev.gj_ptd_puppi,   "gj_ptd_puppi[ngj]/F");
  t->Branch("gj_ptd_all",   tjsev.gj_ptd_all,   "gj_ptd_all[ngj]/F");
  t->Branch("gj_ptds_charged",   tjsev.gj_ptds_charged,   "gj_ptds_charged[ngj]/F");
  t->Branch("gj_ptds_puppi",   tjsev.gj_ptds_puppi,   "gj_ptds_puppi[ngj]/F");
  t->Branch("gj_ptds_all",   tjsev.gj_ptds_all,   "gj_ptds_all[ngj]/F");
  t->Branch("gj_ecc_charged",   tjsev.gj_ecc_charged,   "gj_ecc_charged[ngj]/F");
  t->Branch("gj_ecc_puppi",   tjsev.gj_ecc_puppi,   "gj_ecc_puppi[ngj]/F");
  t->Branch("gj_ecc_all",   tjsev.gj_ecc_all,   "gj_ecc_all[ngj]/F");
  t->Branch("gj_tau21_charged",   tjsev.gj_tau21_charged,   "gj_tau21_charged[ngj]/F");
  t->Branch("gj_tau21_puppi",   tjsev.gj_tau21_puppi,   "gj_tau21_puppi[ngj]/F");
  t->Branch("gj_tau21_all",   tjsev.gj_tau21_all,   "gj_tau21_all[ngj]/F");
  t->Branch("gj_tau32_charged",   tjsev.gj_tau32_charged,   "gj_tau32_charged[ngj]/F");
  t->Branch("gj_tau32_puppi",   tjsev.gj_tau32_puppi,   "gj_tau32_puppi[ngj]/F");
  t->Branch("gj_tau32_all",   tjsev.gj_tau32_all,   "gj_tau32_all[ngj]/F");
  t->Branch("gj_tau43_charged",   tjsev.gj_tau43_charged,   "gj_tau43_charged[ngj]/F");
  t->Branch("gj_tau43_puppi",   tjsev.gj_tau43_puppi,   "gj_tau43_puppi[ngj]/F");
  t->Branch("gj_tau43_all",   tjsev.gj_tau43_all,   "gj_tau43_all[ngj]/F");
  t->Branch("gj_zg_charged",   tjsev.gj_zg_charged,   "gj_zg_charged[ngj]/F");
  t->Branch("gj_zg_puppi",   tjsev.gj_zg_puppi,   "gj_zg_puppi[ngj]/F");
  t->Branch("gj_zg_all",   tjsev.gj_zg_all,   "gj_zg_all[ngj]/F");
  t->Branch("gj_zgxdr_charged",   tjsev.gj_zgxdr_charged,   "gj_zgxdr_charged[ngj]/F");
  t->Branch("gj_zgxdr_puppi",   tjsev.gj_zgxdr_puppi,   "gj_zgxdr_puppi[ngj]/F");
  t->Branch("gj_zgxdr_all",   tjsev.gj_zgxdr_all,   "gj_zgxdr_all[ngj]/F");
  t->Branch("gj_zgdr_charged",   tjsev.gj_zgdr_charged,   "gj_zgdr_charged[ngj]/F");
  t->Branch("gj_zgdr_puppi",   tjsev.gj_zgdr_puppi,   "gj_zgdr_puppi[ngj]/F");
  t->Branch("gj_zgdr_all",   tjsev.gj_zgdr_all,   "gj_zgdr_all[ngj]/F");
  t->Branch("gj_ga_width_charged",   tjsev.gj_ga_width_charged,   "gj_ga_width_charged[ngj]/F");
  t->Branch("gj_ga_width_puppi",   tjsev.gj_ga_width_puppi,   "gj_ga_width_puppi[ngj]/F");
  t->Branch("gj_ga_width_all",   tjsev.gj_ga_width_all,   "gj_ga_width_all[ngj]/F");
  t->Branch("gj_ga_lha_charged",   tjsev.gj_ga_lha_charged,   "gj_ga_lha_charged[ngj]/F");
  t->Branch("gj_ga_lha_puppi",   tjsev.gj_ga_lha_puppi,   "gj_ga_lha_puppi[ngj]/F");
  t->Branch("gj_ga_lha_all",   tjsev.gj_ga_lha_all,   "gj_ga_lha_all[ngj]/F");
  t->Branch("gj_ga_thrust_charged",   tjsev.gj_ga_thrust_charged,   "gj_ga_thrust_charged[ngj]/F");
  t->Branch("gj_ga_thrust_puppi",   tjsev.gj_ga_thrust_puppi,   "gj_ga_thrust_puppi[ngj]/F");
  t->Branch("gj_ga_thrust_all",   tjsev.gj_ga_thrust_all,   "gj_ga_thrust_all[ngj]/F");
  t->Branch("gj_c1_02_charged",   tjsev.gj_c1_02_charged,   "gj_c1_02_charged[ngj]/F");
  t->Branch("gj_c1_02_puppi",   tjsev.gj_c1_02_puppi,   "gj_c1_02_puppi[ngj]/F");
  t->Branch("gj_c1_02_all",   tjsev.gj_c1_02_all,   "gj_c1_02_all[ngj]/F");
  t->Branch("gj_c1_05_charged",   tjsev.gj_c1_05_charged,   "gj_c1_05_charged[ngj]/F");
  t->Branch("gj_c1_05_puppi",   tjsev.gj_c1_05_puppi,   "gj_c1_05_puppi[ngj]/F");
  t->Branch("gj_c1_05_all",   tjsev.gj_c1_05_all,   "gj_c1_05_all[ngj]/F");
  t->Branch("gj_c1_10_charged",   tjsev.gj_c1_10_charged,   "gj_c1_10_charged[ngj]/F");
  t->Branch("gj_c1_10_puppi",   tjsev.gj_c1_10_puppi,   "gj_c1_10_puppi[ngj]/F");
  t->Branch("gj_c1_10_all",   tjsev.gj_c1_10_all,   "gj_c1_10_all[ngj]/F");
  t->Branch("gj_c1_20_charged",   tjsev.gj_c1_20_charged,   "gj_c1_20_charged[ngj]/F");
  t->Branch("gj_c1_20_puppi",   tjsev.gj_c1_20_puppi,   "gj_c1_20_puppi[ngj]/F");
  t->Branch("gj_c1_20_all",   tjsev.gj_c1_20_all,   "gj_c1_20_all[ngj]/F");
  t->Branch("gj_c2_02_charged",   tjsev.gj_c2_02_charged,   "gj_c2_02_charged[ngj]/F");
  t->Branch("gj_c2_02_puppi",   tjsev.gj_c2_02_puppi,   "gj_c2_02_puppi[ngj]/F");
  t->Branch("gj_c2_02_all",   tjsev.gj_c2_02_all,   "gj_c2_02_all[ngj]/F");
  t->Branch("gj_c2_05_charged",   tjsev.gj_c2_05_charged,   "gj_c2_05_charged[ngj]/F");
  t->Branch("gj_c2_05_puppi",   tjsev.gj_c2_05_puppi,   "gj_c2_05_puppi[ngj]/F");
  t->Branch("gj_c2_05_all",   tjsev.gj_c2_05_all,   "gj_c2_05_all[ngj]/F");
  t->Branch("gj_c2_10_charged",   tjsev.gj_c2_10_charged,   "gj_c2_10_charged[ngj]/F");
  t->Branch("gj_c2_10_puppi",   tjsev.gj_c2_10_puppi,   "gj_c2_10_puppi[ngj]/F");
  t->Branch("gj_c2_10_all",   tjsev.gj_c2_10_all,   "gj_c2_10_all[ngj]/F");
  t->Branch("gj_c2_20_charged",   tjsev.gj_c2_20_charged,   "gj_c2_20_charged[ngj]/F");
  t->Branch("gj_c2_20_puppi",   tjsev.gj_c2_20_puppi,   "gj_c2_20_puppi[ngj]/F");
  t->Branch("gj_c2_20_all",   tjsev.gj_c2_20_all,   "gj_c2_20_all[ngj]/F");
  t->Branch("gj_c3_02_charged",   tjsev.gj_c3_02_charged,   "gj_c3_02_charged[ngj]/F");
  t->Branch("gj_c3_02_puppi",   tjsev.gj_c3_02_puppi,   "gj_c3_02_puppi[ngj]/F");
  t->Branch("gj_c3_02_all",   tjsev.gj_c3_02_all,   "gj_c3_02_all[ngj]/F");
  t->Branch("gj_c3_05_charged",   tjsev.gj_c3_05_charged,   "gj_c3_05_charged[ngj]/F");
  t->Branch("gj_c3_05_puppi",   tjsev.gj_c3_05_puppi,   "gj_c3_05_puppi[ngj]/F");
  t->Branch("gj_c3_05_all",   tjsev.gj_c3_05_all,   "gj_c3_05_all[ngj]/F");
  t->Branch("gj_c3_10_charged",   tjsev.gj_c3_10_charged,   "gj_c3_10_charged[ngj]/F");
  t->Branch("gj_c3_10_puppi",   tjsev.gj_c3_10_puppi,   "gj_c3_10_puppi[ngj]/F");
  t->Branch("gj_c3_10_all",   tjsev.gj_c3_10_all,   "gj_c3_10_all[ngj]/F");
  t->Branch("gj_c3_20_charged",   tjsev.gj_c3_20_charged,   "gj_c3_20_charged[ngj]/F");
  t->Branch("gj_c3_20_puppi",   tjsev.gj_c3_20_puppi,   "gj_c3_20_puppi[ngj]/F");
  t->Branch("gj_c3_20_all",   tjsev.gj_c3_20_all,   "gj_c3_20_all[ngj]/F");

}

//
void resetTopJetShapeEvent(TopJetShapeEvent_t &tjsev)
{
  tjsev.nw=0;   tjsev.nl=0;   tjsev.nj=0;   tjsev.ngj=0;   tjsev.ngl=0;   tjsev.met_pt=0; tjsev.met_phi=0;
  for(int i=0; i<1000; i++) tjsev.weight[i]=0;
  for(int i=0; i<5; i++) { tjsev.l_pt[i]=0;   tjsev.l_eta[i]=0;   tjsev.l_phi[i]=0;   tjsev.l_m[i]=0; tjsev.l_id[i]=0; tjsev.gl_pt[i]=0;   tjsev.gl_eta[i]=0;   tjsev.gl_phi[i]=0;   tjsev.gl_m[i]=0; tjsev.gl_id[i]=0; }
  for(int i=0; i<50; i++) {
    tjsev.j_pt[i]=0;   tjsev.j_eta[i]=0;   tjsev.j_phi[i]=0;   tjsev.j_m[i]=0; tjsev.j_flavor[i]=0; tjsev.j_overlap[i]=0; tjsev.j_gj[i]=-1;
    tjsev.gj_pt[i]=0;   tjsev.gj_eta[i]=0;   tjsev.gj_phi[i]=0;   tjsev.gj_m[i]=0; tjsev.gj_flavor[i]=0; tjsev.gj_overlap[i]=0; tjsev.gj_j[i]=-1;
    
    tjsev.j_mult_charged[i]=-99;
    tjsev.j_mult_puppi[i]=-99;
    tjsev.j_mult_all[i]=-99;
    tjsev.j_width_charged[i]=-99;
    tjsev.j_width_puppi[i]=-99;
    tjsev.j_width_all[i]=-99;
    tjsev.j_ptd_charged[i]=-99;
    tjsev.j_ptd_puppi[i]=-99;
    tjsev.j_ptd_all[i]=-99;
    tjsev.j_ptds_charged[i]=-99;
    tjsev.j_ptds_puppi[i]=-99;
    tjsev.j_ptds_all[i]=-99;
    tjsev.j_ecc_charged[i]=-99;
    tjsev.j_ecc_puppi[i]=-99;
    tjsev.j_ecc_all[i]=-99;
    tjsev.j_tau21_charged[i]=-99;
    tjsev.j_tau21_puppi[i]=-99;
    tjsev.j_tau21_all[i]=-99;
    tjsev.j_tau32_charged[i]=-99;
    tjsev.j_tau32_puppi[i]=-99;
    tjsev.j_tau32_all[i]=-99;
    tjsev.j_tau43_charged[i]=-99;
    tjsev.j_tau43_puppi[i]=-99;
    tjsev.j_tau43_all[i]=-99;
    tjsev.j_zg_charged[i]=-99;
    tjsev.j_zg_puppi[i]=-99;
    tjsev.j_zg_all[i]=-99;
    tjsev.j_zgxdr_charged[i]=-99;
    tjsev.j_zgxdr_puppi[i]=-99;
    tjsev.j_zgxdr_all[i]=-99;
    tjsev.j_zgdr_charged[i]=-99;
    tjsev.j_zgdr_puppi[i]=-99;
    tjsev.j_zgdr_all[i]=-99;
    tjsev.j_ga_width_charged[i]=-99;
    tjsev.j_ga_width_puppi[i]=-99;
    tjsev.j_ga_width_all[i]=-99;
    tjsev.j_ga_lha_charged[i]=-99;
    tjsev.j_ga_lha_puppi[i]=-99;
    tjsev.j_ga_lha_all[i]=-99;
    tjsev.j_ga_thrust_charged[i]=-99;
    tjsev.j_ga_thrust_puppi[i]=-99;
    tjsev.j_ga_thrust_all[i]=-99;
    tjsev.j_c1_02_charged[i]=-99;
    tjsev.j_c1_02_puppi[i]=-99;
    tjsev.j_c1_02_all[i]=-99;
    tjsev.j_c1_05_charged[i]=-99;
    tjsev.j_c1_05_puppi[i]=-99;
    tjsev.j_c1_05_all[i]=-99;
    tjsev.j_c1_10_charged[i]=-99;
    tjsev.j_c1_10_puppi[i]=-99;
    tjsev.j_c1_10_all[i]=-99;
    tjsev.j_c1_20_charged[i]=-99;
    tjsev.j_c1_20_puppi[i]=-99;
    tjsev.j_c1_20_all[i]=-99;
    tjsev.j_c2_02_charged[i]=-99;
    tjsev.j_c2_02_puppi[i]=-99;
    tjsev.j_c2_02_all[i]=-99;
    tjsev.j_c2_05_charged[i]=-99;
    tjsev.j_c2_05_puppi[i]=-99;
    tjsev.j_c2_05_all[i]=-99;
    tjsev.j_c2_10_charged[i]=-99;
    tjsev.j_c2_10_puppi[i]=-99;
    tjsev.j_c2_10_all[i]=-99;
    tjsev.j_c2_20_charged[i]=-99;
    tjsev.j_c2_20_puppi[i]=-99;
    tjsev.j_c2_20_all[i]=-99;
    tjsev.j_c3_02_charged[i]=-99;
    tjsev.j_c3_02_puppi[i]=-99;
    tjsev.j_c3_02_all[i]=-99;
    tjsev.j_c3_05_charged[i]=-99;
    tjsev.j_c3_05_puppi[i]=-99;
    tjsev.j_c3_05_all[i]=-99;
    tjsev.j_c3_10_charged[i]=-99;
    tjsev.j_c3_10_puppi[i]=-99;
    tjsev.j_c3_10_all[i]=-99;
    tjsev.j_c3_20_charged[i]=-99;
    tjsev.j_c3_20_puppi[i]=-99;
    tjsev.j_c3_20_all[i]=-99;

    tjsev.gj_mult_charged[i]=-99;
    tjsev.gj_mult_puppi[i]=-99;
    tjsev.gj_mult_all[i]=-99;
    tjsev.gj_width_charged[i]=-99;
    tjsev.gj_width_puppi[i]=-99;
    tjsev.gj_width_all[i]=-99;
    tjsev.gj_ptd_charged[i]=-99;
    tjsev.gj_ptd_puppi[i]=-99;
    tjsev.gj_ptd_all[i]=-99;
    tjsev.gj_ptds_charged[i]=-99;
    tjsev.gj_ptds_puppi[i]=-99;
    tjsev.gj_ptds_all[i]=-99;
    tjsev.gj_ecc_charged[i]=-99;
    tjsev.gj_ecc_puppi[i]=-99;
    tjsev.gj_ecc_all[i]=-99;
    tjsev.gj_tau21_charged[i]=-99;
    tjsev.gj_tau21_puppi[i]=-99;
    tjsev.gj_tau21_all[i]=-99;
    tjsev.gj_tau32_charged[i]=-99;
    tjsev.gj_tau32_puppi[i]=-99;
    tjsev.gj_tau32_all[i]=-99;
    tjsev.gj_tau43_charged[i]=-99;
    tjsev.gj_tau43_puppi[i]=-99;
    tjsev.gj_tau43_all[i]=-99;
    tjsev.gj_zg_charged[i]=-99;
    tjsev.gj_zg_puppi[i]=-99;
    tjsev.gj_zg_all[i]=-99;
    tjsev.gj_zgxdr_charged[i]=-99;
    tjsev.gj_zgxdr_puppi[i]=-99;
    tjsev.gj_zgxdr_all[i]=-99;
    tjsev.gj_zgdr_charged[i]=-99;
    tjsev.gj_zgdr_puppi[i]=-99;
    tjsev.gj_zgdr_all[i]=-99;
    tjsev.gj_ga_width_charged[i]=-99;
    tjsev.gj_ga_width_puppi[i]=-99;
    tjsev.gj_ga_width_all[i]=-99;
    tjsev.gj_ga_lha_charged[i]=-99;
    tjsev.gj_ga_lha_puppi[i]=-99;
    tjsev.gj_ga_lha_all[i]=-99;
    tjsev.gj_ga_thrust_charged[i]=-99;
    tjsev.gj_ga_thrust_puppi[i]=-99;
    tjsev.gj_ga_thrust_all[i]=-99;
    tjsev.gj_c1_02_charged[i]=-99;
    tjsev.gj_c1_02_puppi[i]=-99;
    tjsev.gj_c1_02_all[i]=-99;
    tjsev.gj_c1_05_charged[i]=-99;
    tjsev.gj_c1_05_puppi[i]=-99;
    tjsev.gj_c1_05_all[i]=-99;
    tjsev.gj_c1_10_charged[i]=-99;
    tjsev.gj_c1_10_puppi[i]=-99;
    tjsev.gj_c1_10_all[i]=-99;
    tjsev.gj_c1_20_charged[i]=-99;
    tjsev.gj_c1_20_puppi[i]=-99;
    tjsev.gj_c1_20_all[i]=-99;
    tjsev.gj_c2_02_charged[i]=-99;
    tjsev.gj_c2_02_puppi[i]=-99;
    tjsev.gj_c2_02_all[i]=-99;
    tjsev.gj_c2_05_charged[i]=-99;
    tjsev.gj_c2_05_puppi[i]=-99;
    tjsev.gj_c2_05_all[i]=-99;
    tjsev.gj_c2_10_charged[i]=-99;
    tjsev.gj_c2_10_puppi[i]=-99;
    tjsev.gj_c2_10_all[i]=-99;
    tjsev.gj_c2_20_charged[i]=-99;
    tjsev.gj_c2_20_puppi[i]=-99;
    tjsev.gj_c2_20_all[i]=-99;
    tjsev.gj_c3_02_charged[i]=-99;
    tjsev.gj_c3_02_puppi[i]=-99;
    tjsev.gj_c3_02_all[i]=-99;
    tjsev.gj_c3_05_charged[i]=-99;
    tjsev.gj_c3_05_puppi[i]=-99;
    tjsev.gj_c3_05_all[i]=-99;
    tjsev.gj_c3_10_charged[i]=-99;
    tjsev.gj_c3_10_puppi[i]=-99;
    tjsev.gj_c3_10_all[i]=-99;
    tjsev.gj_c3_20_charged[i]=-99;
    tjsev.gj_c3_20_puppi[i]=-99;
    tjsev.gj_c3_20_all[i]=-99;

  } 
  
  tjsev.gen_sel = -1; tjsev.reco_sel = -1;
}
