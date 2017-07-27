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
#include "TopLJets2015/TopAnalysis/interface/ExclusiveTop.h"
#include "TopLJets2015/TopAnalysis/interface/LeptonEfficiencyWrapper.h"

#include "TopLJets2015/TopAnalysis/interface/FillNumberLUTHandler.h"
#include "TopLJets2015/TopAnalysis/interface/AlignmentLUTHandler.h"
#include "TopLJets2015/TopAnalysis/interface/ProtonReconstruction.h"

#include <vector>
#include <set>
#include <iostream>
#include <algorithm>

#include "TMath.h"
#include "TopQuarkAnalysis/TopTools/interface/MEzCalculator.h"


using namespace std;


//
void RunExclusiveTop(TString filename,
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

  const char* CMSSW_BASE = getenv("CMSSW_BASE");
  CTPPSAlCa::AlignmentLUTHandler pots_align(Form("%s/src/TopLJets2015/TopAnalysis/data/era2016/alignment_collection_v2.out", CMSSW_BASE));
  CTPPSAlCa::FillNumberLUTHandler run_to_fill(Form("%s/src/TopLJets2015/TopAnalysis/data/era2016/fill_run_lut_v2.dat", CMSSW_BASE));
  XiInterpolator proton_reco(Form("%s/src/TopLJets2015/TopAnalysis/data/era2016/ctpps_optics_9mar2017.root", CMSSW_BASE));

  bool isTTbar( filename.Contains("_TTJets") or (normH and TString(normH->GetTitle()).Contains("_TTJets")));
  bool isData( filename.Contains("Data") );
  
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
  attachToMiniEventTree(t,ev,true);
  Int_t nentries(t->GetEntriesFast());
  if (debug) nentries = 10000; //restrict number of entries for testing
  t->GetEntry(0);

  cout << "...producing " << outname << " from " << nentries << " events" << endl;
  
  //auxiliary to solve neutrino pZ using MET
  MEzCalculator neutrinoPzComputer;

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
  ht.addHist("puwgtctr",     new TH1F("puwgtctr",    ";Weight sums;Events",2,0,2));
  ht.addHist("nvtx",         new TH1F("nvtx",        ";Vertex multiplicity;Events",55,-0.5,49.5));
  ht.addHist("njets",        new TH1F("njets",       ";Jet multiplicity;Events",15,-0.5,14.5));
  ht.addHist("nbjets",       new TH1F("nbjets",      ";b jet multiplicity;Events",10,-0.5,9.5));
  ht.addHist("ht",           new TH1F("ht",          ";H_{T} [GeV];Events",50,0,250));
  ht.addHist("mttbar_cen",   new TH1F("mttbar_cen",  ";M_{ttbar} [GeV];Events",50,300,500));

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

      int fill_number = run_to_fill.getFillNumber(ev.run);
      proton_reco.setAlignmentConstants(pots_align.getAlignmentConstants(fill_number));
      
      //assign randomly a run period
      TString period = assignRunPeriod(runPeriods,random);
      
      //////////////////
      // CORRECTIONS //
      ////////////////
      double csvm = 0.8484;
      addBTagDecisions(ev, csvm, csvm);
      if(!ev.isData) smearJetEnergies(ev);
           
      ///////////////////////////
      // RECO LEVEL SELECTION //
      /////////////////////////
      TString chTag = selector.flagFinalState(ev);
      if(chTag=="") continue;
      std::vector<Particle> &leptons     = selector.getSelLeptons(); 
      std::vector<Jet>      &jets        = selector.getJets();  

      //count n b-jets
      std::vector<Jet> bJets,lightJets;
      float scalarht(0.);
      for(size_t ij=0; ij<jets.size(); ij++) 
        {
          if(jets[ij].flavor()==5) bJets.push_back(jets[ij]);
          else                     lightJets.push_back(jets[ij]);
          scalarht += jets[ij].pt();
        }

      //require one good lepton
      if(leptons.size()!=1) continue;
      bool passJets(jets.size()>=4);
      bool passBJets(bJets.size()>=2);
      
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
        double puWgt(puWgtGr[period][0]->Eval(ev.g_pu));
        std::vector<double>puPlotWgts(1,puWgt);
        ht.fill("puwgtctr",1,puPlotWgts);
        
        // lepton trigger*selection weights
        EffCorrection_t trigSF = lepEffH.getTriggerCorrection(leptons, period);
        EffCorrection_t  selSF= lepEffH.getOfflineCorrection(leptons[0], period);

        wgt *= puWgt*trigSF.first*selSF.first;
        
        //top pt weighting
        double topptsf = 1.0;
        if(isTTbar) {
          for (int igen=0; igen<ev.ngtop; igen++) {
            if(abs(ev.gtop_id[igen])!=6) continue;
            topptsf *= TMath::Exp(0.0615-0.0005*ev.gtop_pt[igen]);
          }
        }
        
        // generator level weights
        wgt *= (ev.g_nw>0 ? ev.g_w[0] : 1.0);

        //update weight for plotter
        plotwgts[0]=wgt;
      }

      if (ev.isData) {
        for (int ift=0; ift<ev.nfwdtrk; ift++) {
          float xi, xi_error;
          proton_reco.computeXiSpline(ev.fwdtrk_arm[ift], ev.fwdtrk_pot[ift], ev.fwdtrk_x[ift], xi, xi_error);
        }
      }

      //control histograms
      ht.fill("nvtx",     ev.nvtx,        plotwgts);
      if(passJets)   ht.fill("nbjets", bJets.size(), plotwgts);
      if(passBJets)  ht.fill("njets",  jets.size(),  plotwgts);
      if(bJets.size()>=2 && lightJets.size()>=2)
        {
          //visible system
          TLorentzVector visSystem(leptons[0].p4()+bJets[0].p4()+bJets[1].p4()+lightJets[0].p4()+lightJets[1].p4());
          
          //determine the neutrino kinematics
          TLorentzVector met(0,0,0,0);
          met.SetPtEtaPhiM(ev.met_pt[0],0,ev.met_phi[0],0.);
          neutrinoPzComputer.SetMET(met);
          neutrinoPzComputer.SetLepton(leptons[0].p4());
          float nupz=neutrinoPzComputer.Calculate();
          TLorentzVector neutrinoP4(met.Px(),met.Py(),nupz ,TMath::Sqrt(TMath::Power(met.Pt(),2)+TMath::Power(nupz,2)));
          
          //ttbar system
          TLorentzVector ttbarSystem(visSystem+neutrinoP4);

          ht.fill("ht",         scalarht, plotwgts);          
          ht.fill("mttbar_cen", ttbarSystem.M(),   plotwgts);
        }
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
  fOut->Close();
}
