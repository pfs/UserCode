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

#include <vector>
#include <set>
#include <iostream>
#include <algorithm>

#include "TMath.h"
#include "TopQuarkAnalysis/TopTools/interface/MEzCalculator.h"

#include "TopLJets2015/CTPPSAnalysisTools/interface/LHCConditionsFactory.h"
#include "TopLJets2015/CTPPSAnalysisTools/interface/AlignmentsFactory.h"
#include "TopLJets2015/CTPPSAnalysisTools/interface/XiReconstructor.h"

#include "TopLJets2015/TopAnalysis/CODESamplesExclusiveTop/header/Impacts.h"

using namespace std;


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
  ctpps_aligns.feedAlignments(Form("%s/src/TopLJets2015/CTPPSAnalysisTools/data/2017/alignments_30jan2017.txt", CMSSW_BASE));

  ctpps::LHCConditionsFactory lhc_conds;
  lhc_conds.feedConditions(Form("%s/src/TopLJets2015/CTPPSAnalysisTools/data/2017/xangle_tillTS2.csv", CMSSW_BASE));
  lhc_conds.feedConditions(Form("%s/src/TopLJets2015/CTPPSAnalysisTools/data/2017/xangle_afterTS2.csv", CMSSW_BASE));

  bool isTTbar( filename.Contains("_TTJets") or (normH and TString(normH->GetTitle()).Contains("_TTJets")));
  
  //PREPARE OUTPUT
  TString baseName=gSystem->BaseName(outname); 
  TString dirName=gSystem->DirName(outname);
  TFile *fOut=TFile::Open(dirName+"/"+baseName,"RECREATE");
  fOut->cd();

  //READ TREE FROM FILE
  MiniEvent_t ev;
  TFile *f = TFile::Open(filename);  
  TH1 *triggerList=(TH1 *)f->Get("analysis/triggerList");
  TTree *t = (TTree*)f->Get("analysis/data");
  attachToMiniEventTree(t,ev,true);
  Int_t nentries(t->GetEntriesFast());
  if (debug) nentries = 10000; //restrict number of entries for testing
  t->GetEntry(0);

  cout << "...producing " << outname << " from " << nentries << " events" << endl;
  
  //auxiliary to solve neutrino pZ using MET
  MEzCalculator neutrinoPzComputer;

  //LUMINOSITY+PILEUP
  LumiTools lumi(era,genPU);
    
  //LEPTON EFFICIENCIES
  EfficiencyScaleFactorsWrapper lepEffH(filename.Contains("Data13TeV"),era);

  //B-TAG CALIBRATION
  BTagSFUtil btvSF(era,"DeepCSV",BTagEntry::OperatingPoint::OP_MEDIUM,"",0);
  
  //JEC/JER
  JECTools jec(era);

   //BOOK HISTOGRAMS
  HistTool ht;
  ht.setNsyst(0);
  ht.addHist("puwgtctr",     new TH1F("puwgtctr",    ";Weight sums;Events",2,0,2));
  ht.addHist("nvtx",         new TH1F("nvtx",        ";Vertex multiplicity;Events",100,-0.5,99.5));
  ht.addHist("njets",        new TH1F("njets",       ";Jet multiplicity;Events",15,-0.5,14.5));
  ht.addHist("nbjets",       new TH1F("nbjets",      ";b jet multiplicity;Events",10,-0.5,9.5));
  ht.addHist("ht",           new TH1F("ht",          ";H_{T} [GeV];Events",100,-5,1995));
  ht.addHist("mttbar_cen",   new TH1F("mttbar_cen",  ";M_{ttbar} [GeV];Events",200,-5,1995));
  ht.addHist("nfwdtrk",	     new TH1F("nfwdtrk",     ";nfwdtrk;Events",16,-0.5,15.5));
  ht.addHist("j_mass",  new TH1F("j_mass",	";jet mass;Events",60,-0.5,59.5));
  ht.addHist("j_pt",	new TH1F("j_pt",	";jet pt;Events",200,-1,400));
  ht.addHist("rho",	new TH1F("rho",		";rho;Events",100,-0.5,99.5));


  ht.addHist("RP003: impacts",	new TH2F("RP003: impacts",	";x;y;Events",100,-0.0005,0.0995,100,-0.0505,0.0495));
  ht.addHist("RP103: impacts",	new TH2F("RP103: impacts",	";x;y;Events",100,-0.0005,0.0995,100,-0.0505,0.0495));
  ht.addHist("RP023: impacts",	new TH2F("RP023: impacts",	";x;y;Events",100,-0.0005,0.0995,100,-0.0505,0.0495));
  ht.addHist("RP123: impacts",	new TH2F("RP123: impacts",	";x;y;Events",100,-0.0005,0.0995,100,-0.0505,0.0495));
  ht.addHist("RP003: Xi",	new TH1F("RP003: Xi",		";Fractional momentum loss Xi;Events",100,-0.005,0.995));
  ht.addHist("RP103: Xi",	new TH1F("RP103: Xi",		";Fractional momentum loss Xi;Events",100,-0.005,0.995));
  ht.addHist("RP023: Xi",	new TH1F("RP023: Xi",		";Fractional momentum loss Xi;Events",100,-0.005,0.995));
  ht.addHist("RP123: Xi",	new TH1F("RP123: Xi",		";Fractional momentum loss Xi;Events",100,-0.005,0.995));

  ht.addHist("ngdP_low",	new TH1F("ngdP_low",		";Number of good protons;Events",20,-0.5,19.5));
  ht.addHist("ngdP_mid",	new TH1F("ngdP_mid",		";Number of good protons;Events",20,-0.5,19.5));
  ht.addHist("ngdP_high",	new TH1F("ngdP_high",		";Number of good protons;Events",20,-0.5,19.5));

  ht.addHist("run",		new TH1F("run",			";Run number;Events",10000,290000,310000));
  ht.addHist("run_nosel",	new TH1F("run_nosel",		";Run number;Events",10000,290000,310000));

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
      if(debug){if(iev%10==0) printf ("\r [%3.0f%%] done", 100.*(float)iev/(float)nentries);}

      //int fill_number = run_to_fill.getFillNumber(ev.run);
      //if(debug){std::cout<<"fill number = "<<fill_number<<std::endl;}
      //proton_reco.setAlignmentConstants(pots_align.getAlignmentConstants(fill_number));
      
      //assign randomly a run period
      TString period = lumi.assignRunPeriod();

      //////////////////
      // CORRECTIONS //
      ////////////////      
      btvSF.addBTagDecisions(ev);
      btvSF.updateBTagDecisions(ev);
      jec.smearJetEnergies(ev);
           
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
        double puWgt(lumi.pileupWeight(ev.g_pu,period)[0]);
        std::vector<double>puPlotWgts(1,puWgt);
        ht.fill("puwgtctr",1,puPlotWgts);

      	if (jets.size()>=10){
		std::cout<<"\nFound an entry with >10 jets:\t"<<iev+1<<"\tpuWgt = "<<puWgt<<"\tnvtx = "<<ev.nvtx<<"\tev.g_pu = "<<ev.g_pu<<std::endl;}
	

        // lepton trigger*selection weights
        EffCorrection_t trigSF = lepEffH.getTriggerCorrection(leptons,{},{},period);
        EffCorrection_t  selSF = lepEffH.getOfflineCorrection(leptons[0], period);

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


	if(bJets.size()>=2 && lightJets.size()>=2){
		
		ImpactsRP(ev, ht, plotwgts);
		
		/*const edm::EventID ev_id( ev.run, ev.lumi, ev.event );
		// LHC information retrieval from LUT
		const ctpps::conditions_t lhc_cond = lhc_conds.get( ev_id );
		const double xangle = lhc_cond.crossing_angle;
		std::cout<<"xangle = "<<xangle<<std::endl;
		XiPerRP(ev, ht, xangle, ctpps_aligns, plotwgts, proton_reco);*/

		int nvtxlow = 15;
		int nvtxhigh = 25;
		NGDPvsPU(ev, ht, nvtxlow, nvtxhigh, plotwgts);
		
	}


      }//end of if (ev.isData)


      //control histograms

	if(ev.isData){	ht.fill("run_nosel",ev.run,plotwgts);}

	if(bJets.size()>=2 && lightJets.size()>=2){

	  	ht.fill("nbjets", bJets.size(), plotwgts);
	  	ht.fill("njets",  jets.size(),  plotwgts);
          	ht.fill("nvtx",     ev.nvtx,    plotwgts);
		ht.fill("rho",	ev.rho,		plotwgts);

		if(ev.isData){
			ht.fill("run",ev.run,plotwgts);}

		for (size_t i=0;i<jets.size();i++){
      	     		ht.fill("j_mass",ev.j_mass[i],plotwgts);
       			ht.fill("j_pt",ev.j_pt[i],plotwgts);
			std::cout<<"\tj_id = "<<ev.j_id[i]<<"\tj_btag = "<<ev.j_btag[i]<<std::endl;
		}
	
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

    }//end of loop over events
  
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
