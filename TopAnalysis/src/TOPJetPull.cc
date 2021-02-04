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
#include "TopLJets2015/TopAnalysis/interface/CorrectionTools.h"
#include "TopLJets2015/TopAnalysis/interface/BtagUncertaintyComputer_old.h"
#include "TopLJets2015/TopAnalysis/interface/LeptonEfficiencyWrapper.h"
//#include "TopLJets2015/TopAnalysis/interface/GeneratorTools.h"
#include "TopLJets2015/TopAnalysis/interface/JECTools.h"
#include "TopLJets2015/TopAnalysis/interface/EfficiencyScaleFactorsWrapper.h"
#include "TopQuarkAnalysis/TopTools/interface/MEzCalculator.h"
#include "TopLJets2015/TopAnalysis/interface/L1PrefireEfficiencyWrapper.h"


#include <vector>
#include <set>
#include <iostream>
#include <algorithm>

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>

#include "TMath.h"
#include "TProfile.h"
#include "TopLJets2015/TopAnalysis/interface/CFAT_cmssw.hh"

#include "TopLJets2015/TopAnalysis/interface/JetConstituentAnalysisTool.hh"

#include "CFAT_Event.hh"
#include "TopLJets2015/TopAnalysis/interface/CFAT_AssignHistograms.hh"
#include "TopLJets2015/TopAnalysis/interface/CFAT_Core_cmssw.hh"
#include "TopLJets2015/TopAnalysis/interface/Definitions_cmssw.hh"

using namespace std;

void RunTopJetPull( TString                         inputfilename,
		    TString                         outputfilename,
		    TH1F *                          normH, 
		    Bool_t                          runSysts,
		    std::string                     systVar,
		    TString                         era,
		    float                           xsec)
{
  ////////////////////
  // INITIALIZATION //
  ////////////////////

  printf("Starting TOPJetPull CMSSW_9_4_11\n"                        );
  printf("input file name  [%s]\n",   inputfilename.Data()          );
  printf("output file name [%s]\n",   outputfilename.Data()         );
  printf("normH            [%p]\n",   normH                         );
  printf("ruSysts          [%s]\n",   runSysts ? "True" : "False"   );
  printf("systVar          [%s]\n",   systVar.c_str()               );
  printf("era              [%s]\n",   era.Data()                    );
  printf("xsec             [%f]\n",   xsec                          ); 
  //  dbcut = 1;
  TFile * const inputfile             = TFile::Open(inputfilename);
  TH1   * const genPU                 = (TH1 *)  inputfile -> Get("analysis/putrue");
  TH1   * const triggerList           = (TH1 *)  inputfile -> Get("analysis/triggerList");
  TTree * const inputtree             = (TTree*) inputfile -> Get("analysis/data");
  const unsigned long nentries        = inputtree-> GetEntriesFast();
  
  printf("inputfile   [%p]\n",  inputfile     );
  printf("genPu       [%p]\n",  genPU         );
  printf("triggerList [%p]\n",  triggerList   );
  printf("inputtree   [%p]\n",  inputtree     );
  printf("nentries    [%lu]\n", nentries      );

  const bool isTTbar ( inputfilename.Contains("_TTJets") );
  const bool isData  ( inputfilename.Contains("Data")    );

  const TString baseName = gSystem -> BaseName(outputfilename); 
  const TString dirName  = gSystem -> DirName(outputfilename);

  MiniEvent_t ev;
  attachToMiniEventTree(inputtree, ev, true);

  //SYSTEMATICS

  std::vector<std::string> vSystVar;
  printf("systVar.length() %lu\n", systVar.length());
  boost::split(vSystVar, systVar, boost::is_any_of("_"));
  printf("vSystVar.size() %lu\n", vSystVar.size());
  for (std::vector<std::string>::iterator it = vSystVar.begin(); it != vSystVar.end(); it ++)
    {
      printf("[%s]\n", (*it).c_str()); 
    }

  TH1 * ratevsrunH = nullptr;
  std::map<Int_t, Float_t> lumiMap;
  if (isData)  
    {
      std::pair<std::map<Int_t, Float_t>, TH1F *> result = parseLumiInfo(era);
      lumiMap    = result.first;
      ratevsrunH = result.second;
    }


  // Variation option = Variation::NOMINAL;
  // if (runSysts)
  //   {
  //     if (vSystVar.back() == "up")
  // 	option = Variation::UP;
  //     if (vSystVar.back() == "down")
  // 	option = Variation::DOWN;

  //   }

  //JET ENERGY UNCERTAINTIES                                                                                                                                                                                
  std::string jecVar = "Total";
  if (vSystVar[0] == "jec") 
    {
      if (vSystVar.size() != 3) 
	{
	  std::cout << "Wrong format of JEC uncertainty, expected jec_Source_up/down. Exiting..." << std::endl;
	  return;
	}
      jecVar = vSystVar[1];
    }

  //LEPTON EFFICIENCIES
  LeptonEfficiencyWrapper lepEffH(inputfilename.Contains("Data13TeV"), era);

  const TString jecUncUrl(era + "/Summer16_23Sep2016V4_MC_UncertaintySources_AK4PFchs.txt");
  gSystem -> ExpandPathName(jecUncUrl);
  JetCorrectorParameters * jecParam = new JetCorrectorParameters(jecUncUrl.Data(), jecVar);
  JetCorrectionUncertainty * jecUnc = new JetCorrectionUncertainty( *jecParam );

  MEzCalculator neutrinoPzComputer;


  //  JECTools jectool(era, jecVar);


  //BFRAG UNCERTAINTIES                                                                                                                                                                                     
  std::map<TString, TGraph * > bfrag                   = getBFragmentationWeights_old(era);
  std::map<TString, std::map<int, double> > semilepbr  = getSemilepBRWeights_old(era);

  //lumi
  // LumiTools lumi(era, genPU);
  // std::map<Int_t, Float_t> lumiPerRun = lumi.lumiPerRun();


  TRandom *  random = new TRandom(0); // random seed for period selection

  //MEzCalculator neutrinoPzComputer;

  //PILEUP WEIGHTING

  std::map<TString, std::vector<TGraph *> > puWgtGr;
  if( !isData ) 
    puWgtGr = getPileupWeightsMap(era, genPU);
  
  
  //LEPTON EFFICIENCIES
  // std::map<TString, TString> cfgMap;
  // //  cfgMap["g_id"]     = "MVAwp90";
  // cfgMap["m_id"]     = "TightID";
  // cfgMap["m_iso"]    = "TightRelIso";
  // cfgMap["m_id4iso"] = "TightIDandIPCut";
  // // cfgMap["e_id"]     = "MVA90";
  // cfgMap["e_id"]     = "Tight";
  //EfficiencyScaleFactorsWrapper lepEffH(inputfilename.Contains("Data13TeV"), era, cfgMap);

  //L1-prefire 
  L1PrefireEfficiencyWrapper l1PrefireWR(inputfilename.Contains("Data13TeV"), era);

  //B-TAG CALIBRATION
  //  BTagSFUtil btvSF(era, BTagEntry::OperatingPoint::OP_MEDIUM, "", 0);

  BTagSFUtil_old * myBTagSFUtil = new BTagSFUtil_old();
  std::map<TString, std::map<BTagEntry::JetFlavor, BTagCalibrationReader *> > btvsfReaders = getBTVcalibrationReadersMap(era, BTagEntry::OP_MEDIUM);
  const std::vector<RunPeriod_t> runPeriods = getRunPeriods(era);

  btvsfReaders[runPeriods[0].first][BTagEntry::FLAV_B]    -> eval_auto_bounds("central", BTagEntry::FLAV_B,   0.0, 30.0);
  btvsfReaders[runPeriods[0].first][BTagEntry::FLAV_UDSG] -> eval_auto_bounds("central", BTagEntry::FLAV_UDSG,  0.0, 30.0);

  std::map<BTagEntry::JetFlavor, TGraphAsymmErrors *>    expBtagEffPy8 = readExpectedBtagEff(era);
  TString btagExpPostFix("");
  if(isTTbar) 
    {
      if(inputfilename.Contains("_herwig"))    btagExpPostFix = "_herwig";
      if(inputfilename.Contains("_scaleup"))   btagExpPostFix = "_scaleup";
      if(inputfilename.Contains("_scaledown")) btagExpPostFix = "_scaledown";
    }
  std::map<BTagEntry::JetFlavor, TGraphAsymmErrors *> expBtagEff = readExpectedBtagEff(era, btagExpPostFix);


  //BOOK HISTOGRAMS
  HistTool ht;
  
  if (not isTTbar) 
    {
      ht.setNsyst(0);
      Definitions::nsyst_ = 0;
    }
  // std::map<TString, TH1 *> allPlots;
  // std::map<TString, TH2 *> all2dPlots;
  ht.addHist("puwgtctr", new TH1F("puwgtctr", "Weight sums", 2, 0.0, 2.0));

  const std::vector<TString> stageVec = { "0_pre", "1_1l", "2_1l4j", "3_1l4j2b", "4_1l4j2b2w" };
  const std::vector<TString> chTags   = { "L", "E", "M" };

  const unsigned short Nbin_labels = 4;
  for (unsigned char channel_ind = 0; channel_ind < N_channels_types_; channel_ind ++)
    {
      for (unsigned char level_ind = 0; level_ind < N_levels_types_; level_ind ++)
 	{
 	  const TString key(TString(tag_channels_types_[channel_ind]) + "_" + tag_levels_types_[level_ind] + "_selection");
 	  TH1 * h = new TH1F(key, TString(tag_channels_types_[channel_ind]) + "_" + tag_levels_types_[level_ind] + "_selection; Selection stage; Events", Nbin_labels, -0.5, Nbin_labels - 0.5);
 	  ht.addHist(key, h);
 	  for (unsigned short bin_ind = 1; bin_ind <= Nbin_labels; bin_ind ++)
 	    {
 	      h -> GetXaxis() -> SetBinLabel(bin_ind, title_selection_stages_[bin_ind -1]); 
 	    }

 	}
    }
  for (unsigned short bin_ind = 0; bin_ind < Nbin_labels; bin_ind ++)
    {
      const unsigned char Nbin_labels_2D = 3;
      const TString key = TString("selection2D_") + Form("%u", bin_ind);
      TH2 * h2 = new TH2F(key, TString(title_selection_stages_[bin_ind]) + "; reco; gen", 
 			  Nbin_labels_2D, 0.0, Nbin_labels_2D, 
 			  Nbin_labels_2D, 0.0, Nbin_labels_2D);
      h2 -> GetXaxis() -> SetBinLabel(1, "no reco"); 
      h2 -> GetYaxis() -> SetBinLabel(1, "no gen"); 
      h2 -> GetXaxis() -> SetBinLabel(2, "E"); 
      h2 -> GetYaxis() -> SetBinLabel(2, "E"); 
      h2 -> GetXaxis() -> SetBinLabel(3, "M"); 
      h2 -> GetYaxis() -> SetBinLabel(3, "M"); 
      ht.addHist(key, h2);

    }

 
  for(auto & stage: stageVec) 
    {
      for(auto & channel: chTags) 
	{  
	  const TString tag(channel + stage + "_");
	  ht.addHist(tag + "nvtx",            new TH1F(tag + "nvtx",             "; Vertex multiplicity; Events",                          55,   0.0,   55.0));
	  ht.addHist(tag + "rho",             new TH1F(tag + "rho",              "; #rho; Events",                                         40,   0.0,   40.0));
	  ht.addHist(tag + "nleps",           new TH1F(tag + "nleps",            "; Lepton multiplicity; Events",                           5,  -0.5,    4.5));
	  ht.addHist(tag + "njets",           new TH1F(tag + "njets",            "; Jet multiplicity;    Events",                           15, -0.5,   14.5));
	  ht.addHist(tag + "nbjets",          new TH1F(tag + "nbjets",           "; b jet multiplicity;  Events",                           5,  -0.5,    4.5));
	  ht.addHist(tag + "nwjets",          new TH1F(tag + "nwjets",           "; W jet multiplicity;  Events",                           10, -0.5,    9.5));
	  ht.addHist(tag + "wcandm",          new TH1F(tag + "wcandm",           "; W candidate mass;    W candidates",                     60,  0.0,  300.0));
	  ht.addHist(tag + "tcandm",          new TH1F(tag + "tcandm",           "; top candidate mass;  top candidates",                   70, 50.0,  400.0));
	  ht.addHist(tag + "tcandwcutm",      new TH1F(tag + "tcandwcutm",       "; top candidate mass;  top candidates (W cut)",           70, 50.0,  400.0));
	  ht.addHist(tag + "invmleadingjets", new TH1F(tag + "invmleadingjets",  "; W candidate mass;    W candidates (W cut)",             60,  0.0,  300.0));
	  ht.addHist(tag + "invmwjets",       new TH1F(tag + "invmwjets",        "; W candidate mass;    W candidates (W cut)",             60,  0.0,  300.0));

	  for(unsigned short ind = 0; ind < 2; ind ++) 
	    {
	      const TString pf(Form("l%d", ind));          
	      ht.addHist(tag + pf + "pt",       new TH1F(tag + pf + "pt",   "; Lepton p_{t} [GeV];             Events",            50,  0.0, 250.0));
	      ht.addHist(tag + pf + "eta",      new TH1F(tag + pf + "eta",  "; Lepton pseudo-rapidity;         Events",            50, -2.5,   2.5));
	    }
	  for(unsigned short ind = 0; ind < 6; ind ++)
	    {
	      const TString pf(Form("j%d", ind));
	      ht.addHist(tag + pf + "pt",        new TH1F(tag + pf + "pt",   "; Jet transverse momentum [GeV]; Events",            50,  0.0, 250.0));
	      ht.addHist(tag + pf + "eta",       new TH1F(tag + pf + "eta",  "; Jet pseudo-rapidity;           Events",            50, -5.5,   5.5));
	    }
	  ht.addHist(tag + "met",               new TH1F(tag + "met",              "; MET [GeV];  Events",       50,                    0.0, 250.0));
	  if(ratevsrunH)
	    ht.addHist(tag + "ratevsrun", (TH1 *) ratevsrunH -> Clone(tag + "ratevsrun"));
	  // ht.addHist(tag + "ratevsrun",         new TH1F(tag + "ratevsrun",        "; Run number; #sigma [pb]", int(lumiPerRun.size()), 0.0, float(lumiPerRun.size())));
	  int ind = 0;
	  // for(auto key : lumiPerRun) 
	  //   {
	  //     ind ++;
	  //     ht.getPlots()[tag + "ratevsrun"] -> GetXaxis() -> SetBinLabel(ind, Form("%d", key.first));
	  //   }

	}
    }
  CFAT_cmssw colour_flow_analysis_tool;
  map<TString, pair<TH1 *, TH2 *> *> allPlotsCFAT;
  AssignHistograms(allPlotsCFAT);
  for (map<TString, pair<TH1 *, TH2 *> *>::iterator it = allPlotsCFAT.begin(); it != allPlotsCFAT.end(); it ++)
    {
      if (it -> second -> second) 
	ht.addHist(it -> first, it -> second -> second); 
      ht.addHist(it -> first, it -> second -> first);
    }
  AssignSpecificHistograms2D(ht.get2dPlots());

  colour_flow_analysis_tool.plots_ptr_ = & ht.getPlots();
  colour_flow_analysis_tool.plots2D_ptr_ = & ht.get2dPlots();
  colour_flow_analysis_tool.SetMigrationOutput(dirName + "/migration_" + baseName);

  for (auto & it : ht.getPlots())   
    { 
      it.second -> Sumw2(); 
      it.second -> SetDirectory(nullptr); 
    }
  for (auto & it : ht.get2dPlots()) 
    { 
      it.second -> Sumw2(); 
      it.second -> SetDirectory(nullptr); 
    }

  std::cout << "init done" << std::endl;

  unsigned long nGEN_events  = 0;  
  unsigned long nRECO_events = 0; 
  unsigned long nBOTH_events = 0;
  bool file_completed        = false;

  //EVENT SELECTION WRAPPER
  SelectionTool selector(inputfilename, false, triggerList);

  //////////////////////
  // LOOP OVER EVENTS //
  //////////////////////

  unsigned int dbcut = 1;
  for (unsigned long iev = 0; iev < nentries/dbcut; iev ++)
    {
      //      const TString period                   = lumi.assignRunPeriod();
      const TString period = assignRunPeriod(runPeriods, random);

      //printf("-> period [%s]\n", period.Data());
      if (iev % 1000 == 0) 
	printf ("\r [%3.0f%%] done", 100.0 * (float) iev / (float) nentries);
      inputtree -> GetEntry(iev);
      //      ev.ls();
      bool gen_singleLepton4Jets2b2W = false;
      bool singleLepton4Jets2b2W     = false;
      bool GEN_selected              = false;
      bool RECO_selected             = false;
      map<TString, bool*> selection_reco_map;
      map<TString, bool*> selection_gen_map;
      bool selection[2][2][4];
      selection_reco_map["E"] = selection[0][0];
      selection_reco_map["M"] = selection[0][1];
      selection_gen_map["E"]  = selection[1][0];
      selection_gen_map["M"]  = selection[1][1];
      for (unsigned char ind = 0; ind < 4; ind ++)
	{
	  selection[0][0][ind] = false;
	  selection[0][1][ind] = false;
	  selection[1][0][ind] = false;
	  selection[1][1][ind] = false;
	}
      
      if (not ev.isData)
	{
	  //jec
	  if (vSystVar[0] == "jec") 
	    {
	      applyJetCorrectionUncertainty(ev, jecUnc, jecVar, vSystVar[2]);
	      // jectool.applyJetCorrectionUncertainty(ev, jecVar, option);
	    }

	  //jer
	  if (vSystVar[0] == "jer") 
	    {
	      smearJetEnergies(ev, vSystVar[1]);
	      // jectool.smearJetEnergies(ev, option);
	      // jectool.smearJetEnergies(ev, option);
	    }
	  else 
	    {
	      smearJetEnergies(ev);
	      // jectool.smearJetEnergies(ev);
	    }

	  //csb
	  const float csvm = 0.8484;
	  if (vSystVar[0] == "csv") 
	    {
	      if (vSystVar[1] == "heavy") 
		{
		  //heavy flavor uncertainty +/-3.5%
		  if (vSystVar[2] == "up")   addBTagDecisions(ev, 0.8726, csvm);
		  if (vSystVar[2] == "down") addBTagDecisions(ev, 0.8190, csvm);
		}
	      if (vSystVar[1] == "light") 
		{
		  //light flavor uncertainty +/-10%
		  if (vSystVar[2] == "up")   addBTagDecisions(ev, csvm, 0.8557);
		  if (vSystVar[2] == "down") addBTagDecisions(ev, csvm, 0.8415);
		}
	    }
	  else addBTagDecisions(ev, csvm, csvm);


	  // if (vSystVar[0] == "csv") 
	  //   {
	  //     if (vSystVar[1] == "heavy") 
	  // 	{
	  // 	  //heavy flavor uncertainty +/-3.5%
	  // 	  if (vSystVar[2] == "up")   
	  // 	    btvSF.addBTagDecisions(ev, 0.8726, csvm);
	  // 	  if (vSystVar[2] == "down") 
	  // 	    btvSF.addBTagDecisions(ev, 0.8190, csvm);
	  // 	}
	  //     if (vSystVar[1] == "light") 
	  // 	{
	  // 	  //light flavor uncertainty +/-10%
	  // 	  if (vSystVar[2] == "up")   
	  // 	    btvSF.addBTagDecisions(ev, csvm, 0.8557);
	  // 	  if (vSystVar[2] == "down") 
	  // 	    btvSF.addBTagDecisions(ev, csvm, 0.8415);
	  // 	}
	  //   }
	  // else 
	  //   btvSF.addBTagDecisions(ev, csvm, csvm);

	  //b tagging


	if (vSystVar[0] == "btag")
	  {
	    if (vSystVar[1] == "heavy") 
	      {
		updateBTagDecisions_old(ev, btvsfReaders[period], expBtagEff, expBtagEffPy8, myBTagSFUtil, vSystVar[2], "central");
	      }
	    if (vSystVar[1] == "light") 
	      {
		updateBTagDecisions_old(ev, btvsfReaders[period], expBtagEff, expBtagEffPy8, myBTagSFUtil, "central", vSystVar[2]);
	      }
	  }
	else
	  {
	    //	    BTagSFUtil_old * myBTagSFUtil_check = new BTagSFUtil_old();
	    //    printf("period [%s]\n", period.Data());
	    updateBTagDecisions_old(ev, btvsfReaders[period], expBtagEff, expBtagEffPy8,  myBTagSFUtil);
	  }

	  // if (vSystVar[0] == "btag")
	  //   {
	  //     if (vSystVar[1] == "heavy") 
	  // 	{
	  // 	  btvSF.updateBTagDecisions(ev, vSystVar[2], "central");
	  // 	}
	  //     if (vSystVar[1] == "light") 
	  // 	{
	  // 	  btvSF.updateBTagDecisions(ev, "central", vSystVar[2]);
	  // 	}
	  //   }
	  // else
	  //   {
	  //     btvSF.updateBTagDecisions(ev);
	  //   }
	}

      ///////////////////////////
      // RECO LEVEL SELECTION //
      /////////////////////////
      
      // decide the lepton channel and get selected objects
      const TString chTag                     = selector.flagFinalState_old(ev);
      // printf("new chTag %s\n", chTag.Data());
      std::vector<Particle>     & leptons     = selector.getSelLeptons(); 
      std::vector<Jet>          & jets        = selector.getJets();  
      std::vector<unsigned int> & jet_indices = selector.getJetIndices();
      //count b and W candidates
      //      printf("%lu \n", jets.size());
      int sel_nbjets = 0;
      int sel_nwjets = 0;

      vector<TLorentzVector> bJets, lightJets, nonbJets;

      vector<unsigned short> bJets_index, lightJets_index;
      unsigned char jet_index = 0;
      
      for (auto & jet : jets) 
 	{
 	  const TLorentzVector jp4 = jet.p4();
 	  if (jet.flavor() == 5) 
 	    {
 	      bJets.push_back(jp4);
 	      bJets_index.push_back(jet_indices[jet_index]);
 	      ++ sel_nbjets;
 	    }
 	  else 
 	    {
 	      nonbJets.push_back(jp4);
 	    }
 	  if (jet.flavor() == 1) 
 	    {
 	      lightJets.push_back(jp4);
 	      lightJets_index.push_back(jet_indices[jet_index]);
 	      ++ sel_nwjets;
 	    }
 	  jet_index ++;
 	}
      //printf("%u %u\n", sel_nbjets, sel_nwjets);

      TLorentzVector met(0.0, 0.0, 0.0, 0.0);
      met.SetPtEtaPhiM(ev.met_pt, 0.0, ev.met_phi, 0.0);
      met.SetPz(0.0); 
      met.SetE(met.Pt());

      const bool preselected             (true);
      const bool singleLepton            ((chTag == "E" or chTag == "M") and
                                         (selector.getVetoLeptons().size() == 0));
      const bool singleLepton4Jets       (singleLepton and jets.size() >= 4);
      const bool singleLepton4Jets2b     (singleLepton4Jets and sel_nbjets == 2);
                 singleLepton4Jets2b2W = (singleLepton4Jets2b and sel_nwjets == 2);
      std::vector<bool> recoPass; 
      recoPass.push_back(preselected); 
      recoPass.push_back(singleLepton); 
      recoPass.push_back(singleLepton4Jets); 
      recoPass.push_back(singleLepton4Jets2b); 
      recoPass.push_back(singleLepton4Jets2b2W); 
      
      ///////////////////
      // EVENT WEIGHTS //
      ///////////////////
      //start weights and pu weight control
      std::vector<std::pair<double, bool> >    varweights;
      std::vector<double>                      plotwgts;
      double wgt                               (1.0);
      ht.fill("puwgtctr", 0.0, 1.0);
      if (!ev.isData)
	{
	  // Pairs for systematic uncertainty weights
	  // double: weight value (divided by central weight)
	  // bool: use weight for plotting, otherwise just save to tree
      
	  // norm weight
	  if (normH)
	    {
	      const double normwgt  = normH -> GetBinContent(1);
	      // if (singleLepton) printf("normwgt %.9f\n", normwgt);
	      wgt *= normwgt;
	    }
	  // pu weight
	  const double puWgt(puWgtGr[period][0] -> Eval(ev.g_pu));

	  // double puWgt = lumi.pileupWeight(ev.g_pu, period)[0];
	  ht.fill("puwgtctr", 1.0, puWgt);
	  if (not inputfilename.Contains("dire2002") and not inputfilename.Contains("sherpa") and not inputfilename.Contains("herwig7") ) 
	    {
	      wgt *= puWgt;
	    }
	  
	  varweights.push_back(std::make_pair(puWgtGr[period][1] -> Eval(ev.g_pu), true));                                         // 1
	  varweights.push_back(std::make_pair(puWgtGr[period][2] -> Eval(ev.g_pu), true));                                         // 2
        
	  // lepton trigger*selection weights
	  if (singleLepton) 
	    {
	      const EffCorrection_t trigSF = lepEffH.getTriggerCorrection(leptons, period);
	      varweights.push_back(std::make_pair(1.0 + trigSF.second, true));                                                      //3
	      varweights.push_back(std::make_pair(1.0 - trigSF.second, true));                                                      //4
	      const EffCorrection_t selSF = lepEffH.getOfflineCorrection(leptons[0], period);
	      varweights.push_back(std::make_pair(1.0 + selSF.second, true));                                                       //5
	      varweights.push_back(std::make_pair(1.0 - selSF.second, true));                                                       //6
	      // const EffCorrection_t trigSF = lepEffH.getTriggerCorrection(leptons, {}, {}, period);
	      // varweights.push_back(std::make_pair(1.0 + trigSF.second, true));                                                  
	      // varweights.push_back(std::make_pair(1.0 - trigSF.second, true));                                                  
	      // const EffCorrection_t selSF = lepEffH.getOfflineCorrection(leptons[0], period);
	      // varweights.push_back(std::make_pair(1.0 + selSF.second, true));                                                   
	      // varweights.push_back(std::make_pair(1.0 - selSF.second, true));                                                   
	      const EffCorrection_t l1prefireProb = l1PrefireWR.getCorrection(jets);
	      wgt *= trigSF.first * selSF.first * l1prefireProb.first;
	      // printf("puWgt               %.9f\n", puWgt);
	      // printf("trigSF.first        %.9f\n", trigSF.first);
	      // printf("selSF.first         %.9f\n", selSF.first);
	      // printf("l1prefireProb.first %.9f\n", l1prefireProb.first); 
	      //	      getchar();
	    }

	  else varweights.insert(varweights.end(), 4, std::make_pair(1.0, true));                                                    // 7
	
	  // bfrag weights
	  varweights.push_back(std::make_pair(computeBFragmentationWeight_old(ev, bfrag["upFrag"]), true));                          // 8
	  varweights.push_back(std::make_pair(computeBFragmentationWeight_old(ev, bfrag["downFrag"]), true));                        // 9
	  varweights.push_back(std::make_pair(computeBFragmentationWeight_old(ev, bfrag["PetersonFrag"]), true));                    // 10

	  // weights for semilep BR
	  // simultaneous variation for all hadrons, average over particle and antiparticle
	  // divide by mean weight from 100k events to avoid change in cross section
	  varweights.push_back(std::make_pair(computeSemilepBRWeight_old(ev, semilepbr["semilepbrUp"], 0, true)/1.00480, true));     // 11
	  varweights.push_back(std::make_pair(computeSemilepBRWeight_old(ev, semilepbr["semilepbrDown"], 0, true)/0.992632, true));  // 12

	  double topptsf = 1.0;
	  if (isTTbar) 
	    {
	      for (int igen = 0; igen < ev.ngtop; igen ++) 
		{
		  if (abs(ev.gtop_id[igen]) != 6) 
		    continue;
		  topptsf *= TMath::Exp(0.0615 - 0.0005 * ev.gtop_pt[igen]);
		}
	    }
	  varweights.push_back(std::make_pair(topptsf, true));                                                                        // 13

	  // lhe weights
	  const double gw = (ev.g_nw > 0 ? ev.g_w[0] : 1.0);
	  wgt *= gw;
	  // if (singleLepton)
	  //   {
	  //     printf("gw                  %.9f\n", gw);
	  //     printf("wgt                 %.9f\n", wgt);
	  //   }
	  const std::set<std::string> scalesForPlotter = 
	    {
	      "id1002muR1muF2hdampmt272.7225",                                                                                         // 14
	      "id1003muR1muF0.5hdampmt272.7225",                                                                                       // 15
	      "id1004muR2muF1hdampmt272.7225",                                                                                         // 16
	      "id1005muR2muF2hdampmt272.7225",                                                                                         // 17
	      "id1007muR0.5muF1hdampmt272.7225",                                                                                       // 18
	      "id1009muR0.5muF0.5hdampmt272.7225",                                                                                     // 19
	      "id3001PDFset13100",                                                                                                     // 20   
	      "id4001PDFset25200"                                                                                                      //21
	    };
	  for (int g_nwind = 1; g_nwind < ev.g_nw; ++ g_nwind) 
	    {
	      // printf("%s\n", normH -> GetXaxis() -> GetBinLabel(g_nwind));
	      const bool forPlotter = (normH and scalesForPlotter.count(normH -> GetXaxis() -> GetBinLabel(g_nwind)) != 0);
	      varweights.push_back(std::make_pair(ev.g_w[g_nwind]/ev.g_w[0], forPlotter));
	    }
	  float wgtfull = wgt*xsec*35874.8;
	  if (wgtfull > 10000.0)
	    {
	      printf("wgtfull %.9f. Continuing\n", wgtfull);
	      continue;
	    }
	  plotwgts.push_back(wgt);
	  for (auto vw : varweights)
	    if (vw.second)
	      plotwgts.push_back(vw.first);
	  //printf("plotwgts.size() %lu\n", plotwgts.size());
	}
      else 
	{
	  plotwgts.push_back(wgt);
	}
      //////////////////////////
      // RECO LEVEL ANALYSIS //
      ////////////////////////
      //W and top masses
      // printf("probe A\n");
      std::vector<TLorentzVector> wCands;
      for (unsigned int ijind = 0; ijind < jets.size(); ijind ++) 
	{
	  for (unsigned int jjind = ijind + 1; jjind < jets.size(); jjind ++) 
	    {
	      if (jets[ijind].flavor() == 5 or jets[jjind].flavor() == 5) 
		continue;
	      const TLorentzVector wCand = jets[ijind].p4() + jets[jjind].p4();
	      wCands.push_back(wCand);
	    }
	}
      std::vector<TLorentzVector> tCands;
      for (unsigned int jind = 0; jind < jets.size(); jind ++) 
	{
	  if (jets[jind].flavor() != 5) 
	    continue;
	  for (auto & wCand : wCands) 
	    {
	      const TLorentzVector tCand = jets[jind].p4() + wCand;
	      tCands.push_back(tCand);
	    }
	}
      std::vector<TLorentzVector> tCandsWcut;
      for (unsigned int jind = 0; jind < jets.size(); jind ++) 
	{
	  if (jets[jind].flavor() != 5) 
	    continue;
	  for (auto & wCand : wCands) 
	    {
	      if (abs(wCand.M() - 80.4) > 15.0) 
		continue;
	      const TLorentzVector tCand = jets[jind].p4() + wCand;
	      tCandsWcut.push_back(tCand);
	    }
	}
      //control histograms
      for(size_t istage = 0; istage < stageVec.size(); istage ++ ) 
 	{ 
	  for(auto & channel : chTags) 
	    { 
	      if (not recoPass[istage]) 
		{
		  continue;
		}
	      if (channel == "E" and chTag != "E") 
		{
		  continue;
		}
	      if (channel == "M" and chTag != "M") 
		{
		  continue;
		}
	      const TString tag(channel + stageVec[istage] + "_");

	      // std::map<Int_t, Float_t>::iterator rIt = lumiPerRun.find(ev.run);
	      // //check trigger rates and final channel assignment
	      // // if (rIt != lumiPerRun.end()) 
	      // // 	allPlots[tag + "ratevsrun"] -> Fill(std::distance(lumiPerRun.begin(), rIt), 1.0/rIt -> second);
	      // if (rIt != lumiPerRun.end())
	      // 	{
	      // 	  const int runBin = std::distance(lumiPerRun.begin(), rIt);
	      // 	  const float lumi = 1.0/rIt->second;
	      // 	  ht.fill(tag + "ratevsrun", runBin, lumi);
	      // 	  //allPlots[tag + "ratevsrun"] -> Fill(std::distance(lumiPerRun.begin(), rIt), 1.0/rIt -> second);
	      // 	}
	      // else
	      // 	{
	      // 	  //		  cout << "[Warning] Unable to find run=" << ev.run << endl;
	      // 	}   

	      std::map<Int_t, Float_t>::iterator rIt = lumiMap.find(ev.run);
	      if (rIt != lumiMap.end() && ratevsrunH) 
		ht.fill(tag + "ratevsrun", std::distance(lumiMap.begin(), rIt), 1.0/rIt->second);


	      // std::map<Int_t, Float_t>::iterator rIt = lumiPerRun.find(ev.run);
	      // if (rIt != lumiPerRun.end())
	      // 	{
	      // 	  const int runBin = std::distance(lumiPerRun.begin(), rIt);
	      // 	  const float lumi = 1.0/rIt -> second;
	      // 	  ht.fill(tag + "ratevsrun", runBin, lumi);
	      // 	}
	      // else
	      // 	{
	      // 	  cout << "[Warning] Unable to find run=" << ev.run << endl;
	      // 	}   

	      ht.fill(tag + "nvtx",   ev.nvtx,        plotwgts);
	      ht.fill(tag + "rho",    ev.rho,         plotwgts);
	      ht.fill(tag + "nleps",  leptons.size(), plotwgts);
	      ht.fill(tag + "njets",  jets.size(),    plotwgts);
	      ht.fill(tag + "nbjets", sel_nbjets,     plotwgts);
	      ht.fill(tag + "nwjets", sel_nwjets,     plotwgts);
	      for (auto& wCand : wCands) 
		ht.fill(tag + "wcandm",     wCand.M(), plotwgts);
	      for (auto& tCand : tCands) 
		ht.fill(tag + "tcandm",     tCand.M(), plotwgts);
	      for (auto& tCand : tCandsWcut) 
		ht.fill(tag + "tcandwcutm", tCand.M(), plotwgts);
	      for(unsigned int ind = 0; ind < leptons.size(); ind ++) 
		{
		  if (ind > 1) 
		    break;
		  const TString pf(Form("l%d", ind));          
		  ht.fill(tag + pf + "pt",  leptons[ind].pt(),  plotwgts);
		  ht.fill(tag + pf + "eta", leptons[ind].eta(), plotwgts);
		}
	      if (nonbJets.size() > 1)
		{
		  const float mw =  (nonbJets[0] + nonbJets[1]).M();
		  ht.fill(tag + "invmleadingjets", mw, plotwgts);  
		}
	      if (lightJets.size() == 2)
		{
		  const float mw =  (lightJets[0] + lightJets[1]).M();
		  ht.fill(tag + "invmwjets", mw, plotwgts);  
		}	    
	      for(unsigned int ind = 0; ind < jets.size(); ind ++) 
		{
		  if (ind > 5) 
		    break;
		  const TString pf(Form("j%d", ind));
		  ht.fill(tag + pf + "pt", jets[ind].pt(),   plotwgts);
		  ht.fill(tag + pf + "eta", jets[ind].eta(), plotwgts);
		}
	      ht.fill(tag + "met", ev.met_pt, plotwgts);

	      if (istage > 0 and channel != "L")
		{
		  fill_selection_histo(ht.getPlots(), chTag, tag_levels_types_[RECO], title_selection_stages_[istage - 1], plotwgts[0]); //CommonTools
		  selection_reco_map[TString(chTag)][istage - 1] = true;
		}

	    }

	}
      colour_flow_analysis_tool.ResetMigrationValues();
      if(singleLepton4Jets2b2W)
	{
	  colour_flow_analysis_tool.SetPeriod(period);
	  TLorentzVector lp4(leptons[0].p4());
	  //fill MET
	  const float mt(computeMT(lp4, met));

	  //compute neutrino kinematics
	  neutrinoPzComputer.SetMET(met);
	  neutrinoPzComputer.SetLepton(lp4);
      
	  const float nupz = neutrinoPzComputer.Calculate();
	  const TLorentzVector neutrinoHypP4(met.Px(), met.Py(), nupz, TMath::Sqrt(TMath::Power(met.Pt(), 2) + TMath::Power(nupz, 2)));
	  RECO_selected = true;
	  
	  CFAT_Core_cmssw core_reco;
	  CFAT_Event event_reco;
	  
	  core_reco.SetEvent(ev);
	  event_reco.SetCore(core_reco);
	  core_reco.AddLightJets(lightJets, lightJets_index);
 
	  core_reco.AddVector(Definitions::LEPTON, &lp4);
	  core_reco.AddVector(Definitions::NEUTRINO, &neutrinoHypP4);
	  core_reco.AddBJets(bJets, bJets_index);
	  core_reco.SetEventDisplayMode(kFALSE);
	  event_reco.CompleteVectors();
	  event_reco.SetWeights(plotwgts);
	  event_reco.SetEventNumber(iev);

	  colour_flow_analysis_tool.SetEvent(event_reco);

	  colour_flow_analysis_tool.SetWorkMode(Definitions::RECO);
	  core_reco.RecomputeJetsFromParticles();
	  try
	    {
	      if (chTag == "E")
		colour_flow_analysis_tool.SetChannel(E);
	      else if (chTag == "M")
		colour_flow_analysis_tool.SetChannel(M);
	      else
		throw "No RECO channel";
	      
	      colour_flow_analysis_tool.Work();
	    }
	  catch (const char *e)
	    {
	      	      printf("Exception %s\n", e);
	    }
	}

      //////////////////////
      // GENERATOR LEVEL  //
      //////////////////////
      if (isTTbar) 
 	{
 	  /////////////////////////
 	  // GEN LEVEL SELECTION //
 	  /////////////////////////
       
         
 	  //decide the lepton channel at particle level
 	  std::vector<Particle> genVetoLeptons = selector.getGenLeptons(ev, 15.0, 2.5);
 	  std::vector<Particle> genLeptons     = selector.getGenLeptons(ev, 30.0, 2.1);
	

	  const TString genChTag                 = selector.flagGenFinalState(ev, genLeptons);
	  vector<Jet> & genJets                  = selector.getGenJets();
	  vector<unsigned int> & genJets_indices = selector.getGenJetIndices();
	  //count b and W candidates
	  int sel_ngbjets = 0;
	  int sel_ngwcand = 0;

	  vector<TLorentzVector> gen_bJets, gen_lightJets;
 
	  vector<unsigned short> gen_bJets_index, gen_lightJets_index;
	  unsigned char gen_jet_index = 0;

	  for (auto & genjet : genJets) 
	    {
	      const TLorentzVector gen_jp4 = genjet.p4();
	      if (genjet.flavor() == 5) 
		{
		  gen_bJets.push_back(gen_jp4);
		  gen_bJets_index.push_back(genJets_indices[gen_jet_index]);
		  ++ sel_ngbjets;
		}
	      if (genjet.flavor() == 1) 
		{
		  gen_lightJets.push_back(gen_jp4);
		  gen_lightJets_index.push_back(genJets_indices[gen_jet_index]);
		  ++ sel_ngwcand;
		}
	      gen_jet_index ++;
	    }

	  static const unsigned char gtop_size = 25;
	  TLorentzVector  * gen_nu = nullptr; 
	  TLorentzVector    gen_nu_store;
	  bool gen_nu_found = false;
	  for (unsigned char gtop_ind = 0; gtop_ind < gtop_size; gtop_ind ++)
	    {
	      if (abs(ev.gtop_id[gtop_ind]) == 12000 or abs(ev.gtop_id[gtop_ind]) == 14000)
		{
		  gen_nu_store.SetPtEtaPhiM(ev.gtop_pt[gtop_ind], ev.gtop_eta[gtop_ind], ev.gtop_phi[gtop_ind], ev.gtop_m[gtop_ind]);
		  if (gen_nu_found)
		    {
		      gen_nu_found = false;
		      break;
		    }
		  gen_nu_found = true;
		  gen_nu = & gen_nu_store;
		}
	    
	    }

	  const bool gen_preselected             (true);
	  const bool genSingleLepton             ((genChTag=="E" or genChTag=="M") and
						  (genVetoLeptons.size() == 1)); // only selected lepton in veto collection
	  const bool gen_singleLepton4Jets       (genSingleLepton and genJets.size() >= 4);
	  const bool gen_singleLepton4Jets2b     (gen_singleLepton4Jets and sel_ngbjets == 2/* and gen_nu_found*/);
	             gen_singleLepton4Jets2b2W = (gen_singleLepton4Jets2b and sel_ngwcand == 2);
	  const bool gen_pass[4] = {genSingleLepton, gen_singleLepton4Jets, gen_singleLepton4Jets2b, gen_singleLepton4Jets2b2W};
	  for (unsigned char gen_stage_ind = 0; gen_stage_ind < 4 and gen_pass[gen_stage_ind]; gen_stage_ind ++)
	    {
	    
	      fill_selection_histo(ht.getPlots(), genChTag, tag_levels_types_[GEN], title_selection_stages_[gen_stage_ind], wgt);
	      selection_gen_map[genChTag][gen_stage_ind] = true;

	    }
	  for (unsigned char stageind = 0; stageind < 4; stageind ++)
	    {
	      const TString key = TString("selection2D_") + Form("%u", stageind);
	      TH2 * h2 = ht.get2dPlots()[key];
	      TString reco_ch;
	      if (selection_reco_map["E"][stageind]) 
		{
		  reco_ch = "E";
		}	      
	      else if (selection_reco_map["M"][stageind]) 
		{
		  reco_ch = "M";
		}
	      else
		{
		  reco_ch = "no reco";
		}
	      TString gen_ch;
	      if (selection_gen_map["E"][stageind]) 
		{
		  gen_ch = "E";
		}	      
	      else if (selection_gen_map["M"][stageind]) 
		{
		  gen_ch = "M";
		}
	      else
		{
		  gen_ch = "no gen";
		}
	      h2 -> Fill(reco_ch, gen_ch, wgt);
	    }
	  if (gen_singleLepton4Jets2b2W)
	    {
	      colour_flow_analysis_tool.SetPeriod(period);
	      const TLorentzVector gen_lp4(genVetoLeptons[0].p4()); 
	      GEN_selected = true;
	      CFAT_Core_cmssw core_gen;
	      CFAT_Event event_gen;
	      core_gen.SetEvent(ev);
	      event_gen.SetCore(core_gen);
	      core_gen.AddLightJets(gen_lightJets, gen_lightJets_index);
	      core_gen.AddVector(Definitions::LEPTON, & gen_lp4);
	      core_gen.AddVector(Definitions::NEUTRINO, gen_nu);
	      core_gen.AddBJets(gen_bJets, gen_bJets_index);
	      core_gen.SetEventDisplayMode(0);
	      event_gen.CompleteVectors();
	      event_gen.SetWeights(plotwgts);
	      event_gen.SetEventNumber(iev);

	      colour_flow_analysis_tool.SetEvent(event_gen);
	  
	      colour_flow_analysis_tool.SetWorkMode(Definitions::GEN);
	  
	      colour_flow_analysis_tool.SetChannel(genChTag == "E" ? E : M);
	      core_gen.RecomputeJetsFromParticles();
	      colour_flow_analysis_tool.Work();
	    }
	}
      //      printf("plotting values %lu\n", plotwgts.size());
      colour_flow_analysis_tool.PlotMigrationValues();
      if (GEN_selected and not RECO_selected)
	{
	  nGEN_events ++;
	  
	}
      if (not GEN_selected and RECO_selected)
	{
	  nRECO_events ++;
	  
	}
      if (GEN_selected and RECO_selected)
	{
	  nBOTH_events ++;
	}
    }
  
  //close input file
  inputfile -> Close();

  //save histos to file  
  TFile * fOut = TFile::Open(dirName + "/" + baseName, "RECREATE");

  fOut -> cd();
  for (auto & it : ht.getPlots())  
    { 
      it.second->SetDirectory(fOut); 
      it.second->Write(); 
    }
  for (auto & it : ht.get2dPlots())  
    { 
      it.second->SetDirectory(fOut); 
      it.second->Write(); 
    }
  // for (auto & it : allPlots)  
  //   { 
  //     it.second->SetDirectory(fOut); 
  //     it.second->Write(); 
  //   }
  // for (auto & it : all2dPlots)  
  //   { 
  //     it.second->SetDirectory(fOut); it.second->Write(); 
  //   }
  colour_flow_analysis_tool.WriteMigrationTree();
  fOut -> Close();
}
