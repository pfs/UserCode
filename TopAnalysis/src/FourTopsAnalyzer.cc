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
#include "TopLJets2015/TopAnalysis/interface/FourTopsAnalyzer.h"
#include "TopLJets2015/TopAnalysis/interface/LeptonEfficiencyWrapper.h"
#include "TopLJets2015/TopAnalysis/interface/TMVAUtils.h"

#include <vector>
#include <set>
#include <iostream>
#include <algorithm>

#include "TMath.h"

using namespace std;


//
void RunFourTopsAnalyzer(TString filename,
                          TString outname,
                          Int_t channelSelection, 
                          Int_t chargeSelection, 
                          TH1F *normH, 
                          TString era,
                          Bool_t debug)
{
  /////////////////////
  // INITIALIZATION //
  ///////////////////
  TRandom* random = new TRandom(0); // random seed for period selection
  std::vector<RunPeriod_t> runPeriods=getRunPeriods(era);
  bool isTTbar( filename.Contains("_TTJets") or (normH and TString(normH->GetTitle()).Contains("_TTJets")));
  bool isData( filename.Contains("Data") );
  
  TMVA::Reader *tmvaReader = new TMVA::Reader( "!Color:!Silent" );
  //FIXME TOP: add all the variables you need to use here and update the xml file location
  Float_t lepton_ht;
  Float_t lepton_3_chIso;
  Float_t lepton_4_chIso;
  Float_t lepton_3_minIso;
  Float_t lepton_4_minIso;
  Float_t jet_1_highest_csv;
  Float_t jet_2_highest_csv;
  Float_t jet_3_highest_csv;
  Float_t jet_4_highest_csv;
  Float_t jet_ht;
  Float_t jet_btag_ht;
  Float_t jet_non_btag_ht;
  Float_t jet_smallest_angle;
  Float_t jet_smallest_angle_2b;
  Float_t met_pt;
  Float_t n_jets;
  Float_t n_bjets;
  Float_t n_leptons;
  Float_t n_mu_p;
  Float_t n_mu_m;
  Float_t n_ele_p;
  Float_t n_ele_m;

  Float_t event_weight;

  tmvaReader->AddVariable("jet_1_highest_csv",     &jet_1_highest_csv);
  tmvaReader->AddVariable("jet_2_highest_csv",     &jet_2_highest_csv);
  tmvaReader->AddVariable("jet_3_highest_csv",     &jet_3_highest_csv);
  tmvaReader->AddVariable("jet_4_highest_csv",     &jet_4_highest_csv);
  tmvaReader->AddVariable("jet_ht",                &jet_ht);
  tmvaReader->AddVariable("jet_btag_ht",           &jet_btag_ht);
  tmvaReader->AddVariable("jet_non_btag_ht",       &jet_non_btag_ht);
  tmvaReader->AddVariable("jet_smallest_angle",    &jet_smallest_angle);
  tmvaReader->AddVariable("jet_smallest_angle_2b", &jet_smallest_angle_2b);
  tmvaReader->AddVariable("lepton_ht",             &lepton_ht);
  tmvaReader->AddVariable("lepton_3_chIso",        &lepton_3_chIso);
  tmvaReader->AddVariable("lepton_4_chIso",        &lepton_4_chIso);
  tmvaReader->AddVariable("lepton_3_minIso",       &lepton_3_minIso);
  tmvaReader->AddVariable("lepton_4_minIso",       &lepton_4_minIso);
  tmvaReader->AddVariable("met_pt",                &met_pt);
  tmvaReader->AddVariable("n_jets",                &n_jets);
  tmvaReader->AddVariable("n_bjets",               &n_bjets);
  tmvaReader->AddVariable("n_leptons",             &n_leptons);
  tmvaReader->AddVariable("n_mu_p",                &n_mu_p);
  tmvaReader->AddVariable("n_mu_m",                &n_mu_m);
  tmvaReader->AddVariable("n_ele_p",               &n_ele_p);
  tmvaReader->AddVariable("n_ele_m",               &n_ele_m);
  tmvaReader->BookMVA( "BDT", "/afs/cern.ch/user/v/vwachira/CMSSW_8_0_28/src/TopLJets2015/TopAnalysis/test/summer2017/FourTopsML_BDT_606.weights.xml");
  tmvaReader->BookMVA("MLP_ANN","/afs/cern.ch/user/v/vwachira/CMSSW_8_0_28/src/TopLJets2015/TopAnalysis/test/summer2017/FourTopsML_MLP_ANN_606.weights.xml");

  //CREATE OUTPUT TREE
  TTree outTree("TMVAanalysis","NTuple tree for TMVA analysis");
  outTree.Branch("lepton_ht", &lepton_ht,"lepton_ht/F");
  outTree.Branch("lepton_3_chIso", &lepton_3_chIso,"lepton_3_chIso/F");
  outTree.Branch("lepton_4_chIso", &lepton_4_chIso,"lepton_4_chIso/F");
  outTree.Branch("lepton_3_minIso", &lepton_3_minIso,"lepton_3_minIso/F");
  outTree.Branch("lepton_4_minIso", &lepton_4_minIso,"lepton_4_minIso/F");
  outTree.Branch("jet_1_highest_csv", &jet_1_highest_csv,"jet_1_highest_csv/F");
  outTree.Branch("jet_2_highest_csv", &jet_2_highest_csv,"jet_2_highest_csv/F");
  outTree.Branch("jet_3_highest_csv", &jet_3_highest_csv,"jet_3_highest_csv/F");
  outTree.Branch("jet_4_highest_csv", &jet_4_highest_csv,"jet_4_highest_csv/F");
  outTree.Branch("jet_ht", &jet_ht,"jet_ht/F");
  outTree.Branch("jet_btag_ht", &jet_btag_ht,"jet_btag_ht/F");
  outTree.Branch("jet_non_btag_ht", &jet_non_btag_ht,"jet_non_btag_ht/F");
  outTree.Branch("jet_smallest_angle", &jet_smallest_angle,"jet_smallest_angle/F");
  outTree.Branch("jet_smallest_angle_2b", &jet_smallest_angle_2b,"jet_smallest_angle_2b/F");
  outTree.Branch("met_pt",&met_pt,"met_pt/F");
  outTree.Branch("event_weight",&event_weight,"event_weight/F");
  outTree.Branch("n_jets",&n_jets,"n_jets/F");
  outTree.Branch("n_bjets",&n_bjets,"n_bjets/F");
  outTree.Branch("n_leptons",&n_leptons,"n_leptons/F");

  outTree.Branch("n_mu_p",&n_mu_p,"n_mu_p/F");
  outTree.Branch("n_mu_m",&n_mu_m,"n_mu_m/F");
  outTree.Branch("n_ele_p",&n_ele_p,"n_ele_p/F");
  outTree.Branch("n_ele_m",&n_ele_m,"n_ele_m/F");

  //PREPARE OUTPUT
  TString baseName=gSystem->BaseName(outname); 
  TString dirName=gSystem->DirName(outname);
  TFile *fOut=TFile::Open(dirName+"/"+baseName,"RECREATE");
  fOut->cd();

  //READ TREE FROM FILE
  MiniEvent_t ev;
  TFile *f = TFile::Open(filename);
  TH1 *genPU=(TH1 *)f->Get("analysis/putrue");
  TH1 *triggerList=(TH1 *)f->Get("analysis/triggerList");
  TTree *t = (TTree*)f->Get("analysis/data");
  attachToMiniEventTree(t,ev,true,false);
  Int_t nentries(t->GetEntriesFast());
  if (debug) nentries = 10000; //restrict number of entries for testing
  t->GetEntry(0);

  cout << "...producing " << outname << " from " << nentries << " events" << endl;
  
  //PILEUP WEIGHTING
  std::map<TString, std::vector<TGraph *> > puWgtGr;
  if( !isData ) puWgtGr=getPileupWeightsMap(era,genPU);
    
  //LEPTON EFFICIENCIES
  LeptonEfficiencyWrapper lepEffH(filename.Contains("Data13TeV"),era);

  //B-TAG CALIBRATION
  std::map<TString, std::map<BTagEntry::JetFlavor, BTagCalibrationReader *> > btvsfReaders = getBTVcalibrationReadersMap(era, BTagEntry::OP_MEDIUM);
  std::map<BTagEntry::JetFlavor, TGraphAsymmErrors *>    expBtagEffPy8 = readExpectedBtagEff(era);
  
   //BOOK HISTOGRAMS
  HistTool ht;
  ht.setNsyst(0);
  ht.addHist("puwgtctr", new TH1F("puwgtctr",";Weight sums;Events",2,0,2));
  ht.addHist("nvtx",     new TH1F("nvtx",";Vertex multiplicity;Events",50,-0.5,49.5));
  ht.addHist("njets",    new TH1F("njets",";Jet multiplicity;Events",15,-0.5,14.5));
  ht.addHist("nbjets",   new TH1F("nbjets",";b jet multiplicity;Events",10,-0.5,9.5));
  ht.addHist("nleptons", new TH1F("nleptons",";Lepton multiplicity;Events",6,-0.5,5.5));
  ht.addHist("bdt",      new TH1F("bdt",";BDT;Events",50,-1,1));
  ht.addHist("mlp_ann",  new TH1F("mlp_ann",";MLP_ANN; Events",50,0.,1.));

  std::cout << "init done" << std::endl;

  ///////////////////////
  // LOOP OVER EVENTS //
  /////////////////////
  
  //EVENT SELECTION WRAPPER
  SelectionTool selector(filename, false, triggerList);
  
  for (Int_t iev=0;iev<nentries;iev++)
    {
      t->GetEntry(iev);
      if(iev%10==0) printf ("\r [%3.0f%%] done", 100.*(float)iev/(float)nentries);
      
      //assign randomly a run period
      TString period = assignRunPeriod(runPeriods,random);
      
      //////////////////
      // CORRECTIONS //
      ////////////////
      double csvl = 0.5426;
      double csvm = 0.8484;
      double csvt = 0.9535;
      addBTagDecisions(ev, csvm, csvm);
      if(!ev.isData) smearJetEnergies(ev);
           
      ///////////////////////////
      // RECO LEVEL SELECTION //
      /////////////////////////
      TString chTag = selector.flagFinalState(ev);
      if(chTag=="") continue;
      std::vector<Particle> &leptons     = selector.getAtLeastVetoLeptons();
      std::vector<Particle> &LLIDleptons = selector.getSelLeptons();
      std::vector<Particle> &vetoLeptons = selector.getVetoLeptons(); 
      std::vector<Jet>      &jets        = selector.getJets();  
      //if(leptons.size()<2) continue;
      //if (leptons.size() < 1) continue;

      if (leptons.size() < 4 || LLIDleptons.size() < 2) continue;
      
      //count n b-jets
      int sel_nbjets(0);
      for(size_t ij=0; ij<jets.size(); ij++) sel_nbjets += (jets[ij].flavor()==5);
      if (jets.size() < 4 || sel_nbjets < 2) continue;

      ////////////////////
      // EVENT WEIGHTS //
      //////////////////
      float wgt(1.0);
      std::vector<double>plotwgts(1,wgt);
      ht.fill("puwgtctr",0,plotwgts);
      if (!ev.isData) {

        // norm weight
        wgt  = (normH? normH->GetBinContent(1) : 1.0);
        
        // pu weight
        //if (puWgtGr[period][0]? false : true) continue;
        double puWgt(puWgtGr[period][0]->Eval(ev.g_pu));
        std::vector<double>puPlotWgts(1,puWgt);
        ht.fill("puwgtctr",1,puPlotWgts);
        wgt *= puWgt;
        
        // lepton trigger*selection weights
        EffCorrection_t trigSF = lepEffH.getTriggerCorrection(leptons, period);
        EffCorrection_t selSF(1.,0.);
        for(size_t il=0; il<leptons.size(); il++)
          { 
            EffCorrection_t  ilSelSF= lepEffH.getOfflineCorrection(leptons[il], period);
            selSF.first *= ilSelSF.first;
            selSF.second += sqrt(pow(ilSelSF.second,2)+pow(selSF.second,2));
          }
        wgt *= trigSF.first*selSF.first;
        
        //top pt weighting
        double topptsf = 1.0;
        if(isTTbar) {
          for (int igen=0; igen<ev.ngtop; igen++) {
            if(abs(ev.gtop_id[igen])!=6) continue;
            topptsf *= TMath::Exp(0.0615-0.0005*ev.gtop_pt[igen]);
          }
        }
        
        // lhe weights
        wgt *= (ev.g_nw>0 ? ev.g_w[0] : 1.0);

        //update weight for plotter
        plotwgts[0]=wgt;
      }

      //control histograms
      ht.fill("nvtx",     ev.nvtx,        plotwgts);
      ht.fill("njets",    jets.size(),    plotwgts);
      ht.fill("nbjets",   sel_nbjets,     plotwgts);
      ht.fill("nleptons", leptons.size(), plotwgts);

      //FIXME TOP: assign all the values to the variables here before evaluating the bdt discriminator
      //control histograms

      n_jets = (float) jets.size();
      n_bjets = (float) sel_nbjets;
      n_leptons = (float) leptons.size();

      n_mu_p = 0.;
      n_mu_m = 0.;
      n_ele_p = 0.;
      n_ele_m = 0.;

      for (size_t i=0;i<leptons.size();i++)
      {
        if (leptons[i].id() == 13 && leptons[i].charge() > 0.) n_mu_p++;
        else if (leptons[i].id() == 13 && leptons[i].charge() < 0.) n_mu_m++;
        else if (leptons[i].id() == 11 && leptons[i].charge() > 0.) n_ele_p++;
        else if (leptons[i].id() == 11 && leptons[i].charge() < 0.) n_ele_m++;
      }

      // Fill only lepton with highest pt
      
      lepton_3_chIso = 0.;
      lepton_4_chIso = 0.;
      lepton_3_minIso = 0.;
      lepton_4_minIso = 0.;

      if (leptons.size() > 0)
      {

        if (leptons.size() >= 4)
        {
          auto chIsoSort = [ev](Particle &p1, Particle &p2)
          { 
            return ev.l_chargedHadronIso[p1.originalReference()] < ev.l_chargedHadronIso[p2.originalReference()];
          };
          auto minIsoSort = [ev](Particle &p1, Particle &p2)
          { 
            return ev.l_miniIso[p1.originalReference()] < ev.l_miniIso[p2.originalReference()];
          };
          std::vector<Particle> leptons_chIsoSorted  = leptons;
          std::vector<Particle> leptons_minIsoSorted = leptons;
          std::sort(leptons_chIsoSorted.begin(), leptons_chIsoSorted.end(), chIsoSort);
          std::sort(leptons_minIsoSorted.begin(), leptons_minIsoSorted.end(), minIsoSort);

          lepton_3_chIso = ev.l_chargedHadronIso[leptons_chIsoSorted[2].originalReference()];
          lepton_4_chIso = ev.l_chargedHadronIso[leptons_chIsoSorted[3].originalReference()];
          lepton_3_minIso = ev.l_miniIso[leptons_minIsoSorted[2].originalReference()];
          lepton_4_minIso = ev.l_miniIso[leptons_minIsoSorted[3].originalReference()];
        }
      }

      lepton_ht = 0.;
      for (size_t i=0;i<leptons.size();i++) lepton_ht += leptons[i].pt();

      jet_1_highest_csv = 0;
      jet_2_highest_csv = 0;
      jet_3_highest_csv = 0;
      jet_4_highest_csv = 0;

      auto greater_than = [](double i, double j) { return i > j; };

      auto sort_jet_csv = [greater_than](auto ev, auto jets)
      {
        std::vector<float> csv_sorted;
        for (size_t i=0;i<jets.size();i++) csv_sorted.push_back(ev.j_csv[jets[i].getJetIndex()]);
        std::sort(csv_sorted.begin(),csv_sorted.end(),greater_than);
        return csv_sorted;
      };

      std::vector<float> csv_sorted = sort_jet_csv(ev,jets);
      if (csv_sorted.size() > 0) jet_1_highest_csv = csv_sorted[0];
      if (csv_sorted.size() > 1) jet_2_highest_csv = csv_sorted[1];
      if (csv_sorted.size() > 2) jet_3_highest_csv = csv_sorted[2];
      if (csv_sorted.size() > 3) jet_4_highest_csv = csv_sorted[3];

      jet_smallest_angle = TMath::Pi();
      for (size_t i=0;i<jets.size();i++) 
      {
        if (jets[i].pt() == 0) continue;
        for (size_t j=0;j<leptons.size();j++)
        {
          if (leptons[j].pt() == 0) continue;
          if (jets[i].p4().Angle(leptons[j].p4().Vect()) < jet_smallest_angle)
            jet_smallest_angle = jets[i].p4().Angle(leptons[j].p4().Vect());
        }
      }

      jet_smallest_angle_2b = TMath::Pi();
      if (sel_nbjets >= 2) for (size_t i=0;i<jets.size();i++)
      {
        if (jets[i].flavor() != 5) continue;
        for (size_t j=i+1;j<jets.size();j++)
        {
          if (jets[j].flavor() != 5) continue;
          if (jets[i].p4().Angle(jets[j].p4().Vect()) < jet_smallest_angle_2b)
            jet_smallest_angle_2b = jets[i].p4().Angle(jets[j].p4().Vect());
        }
      }

      jet_ht = 0.;
      for (size_t i=0;i<jets.size();i++) jet_ht += jets[i].pt();

      jet_btag_ht = 0.;
      jet_non_btag_ht = 0.;
      for (size_t i=0;i<jets.size();i++) 
      {
        if (ev.j_btag[jets[i].getJetIndex()]) jet_btag_ht += jets[i].pt();
        else jet_non_btag_ht += jets[i].pt();
      }

      met_pt = ev.met_pt[0];

      event_weight = plotwgts[0];
      
      float bdt = tmvaReader->EvaluateMVA("BDT");
      ht.fill("bdt", bdt, plotwgts);
      float mlp_ann = tmvaReader->EvaluateMVA("MLP_ANN");
      ht.fill("mlp_ann", mlp_ann, plotwgts);

      outTree.Fill();
    }
  
  //close input file
  f->Close();
  
  //save histos to file  
  fOut->cd();
  for (auto& it : ht.getPlots())  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  for (auto& it : ht.get2dPlots())  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  outTree.Write();
  fOut->Close();
}
