#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>

#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TopLJets2015/TopAnalysis/interface/TOP-16-006.h"
#include "TopLJets2015/TopAnalysis/interface/BtagUncertaintyComputer.h"

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"


#include <vector>
#include <iostream>
#include <algorithm>

#include "TMath.h"
#include "TopQuarkAnalysis/TopTools/interface/MEzCalculator.h"

using namespace std;


//
void RunTop16006(TString filename,
		 TString outname,
		 Int_t channelSelection, 
		 Int_t chargeSelection, 
		 FlavourSplitting flavourSplitting,
		 TH1F *normH, 
		 Bool_t runSysts)
{

  bool isTTbar( filename.Contains("_TTJets") );

  //READ TREE FROM FILE
  MiniEvent_t ev;
  TFile *f = TFile::Open(filename);
  TH1 *puTrue=(TH1 *)f->Get("analysis/putrue");
  puTrue->SetDirectory(0);
  puTrue->Scale(1./puTrue->Integral());
  TTree *t = (TTree*)f->Get("analysis/data");
  attachToMiniEventTree(t,ev);
  Int_t nentries(t->GetEntriesFast());
  t->GetEntry(0);
  if(ev.isData) runSysts=false;
  cout << "...producing " << outname << " from " << nentries << " events" << (runSysts ? " syst variations will be considered" : "") << endl;

  //auxiliary to solve neutrino pZ using MET
  MEzCalculator neutrinoPzComputer;

  //for data only get the lumi per run map
  std::map<Int_t,Float_t> lumiMap;
  if(ev.isData) lumiMap=lumiPerRun();

  //PILEUP WEIGHTING
  std::vector<TGraph *>puWgtGr;
  if(!ev.isData)
    {
      TString puWgtUrl("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/pileupWgts.root");
      gSystem->ExpandPathName(puWgtUrl);
      TFile *fIn=TFile::Open(puWgtUrl);
      for(size_t i=0; i<3; i++)
	{
	  TString grName("pu_nom");
	  if(i==1) grName="pu_down";
	  if(i==2) grName="pu_up";
	  TGraph *puData=(TGraph *)fIn->Get(grName);
	  Float_t totalData=puData->Integral();
	  TH1 *tmp=(TH1 *)puTrue->Clone("tmp");
	  for(Int_t xbin=1; xbin<=tmp->GetXaxis()->GetNbins(); xbin++)
	    {
	      Float_t yexp=puTrue->GetBinContent(xbin);
	      Double_t xobs,yobs;
	      puData->GetPoint(xbin-1,xobs,yobs);
	      tmp->SetBinContent(xbin, yexp>0 ? yobs/(totalData*yexp) : 0. );
	    }
	  TGraph *gr=new TGraph(tmp);
	  grName.ReplaceAll("pu","puwgts");
	  gr->SetName(grName);
	  puWgtGr.push_back( gr );
	  tmp->Delete();
	}
      
      /*
	if(fIn)
	{
	puWgtGr.push_back( (TGraph *)fIn->Get("puwgts_nom") );
	puWgtGr.push_back( (TGraph *)fIn->Get("puwgts_down") );
	puWgtGr.push_back( (TGraph *)fIn->Get("puwgts_up") );
	fIn->Close();
	}
      */
    }

  //LEPTON EFFICIENCIES
  TString lepEffUrl("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/leptonEfficiencies.root");
  gSystem->ExpandPathName(lepEffUrl);
  std::map<TString,TH2 *> lepEffH;
  if(!ev.isData)
    {
      TFile *fIn=TFile::Open(lepEffUrl);
      lepEffH["m_sel"]=(TH2 *)fIn->Get("m_sel");
      lepEffH["m_trig"]=(TH2 *)fIn->Get("m_trig");      
      for(auto& it : lepEffH) it.second->SetDirectory(0);
      fIn->Close();
    }

  lepEffUrl="${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/CutBasedID_TightWP_76X_18Feb.txt_SF2D.root";
  gSystem->ExpandPathName(lepEffUrl);
  if(!ev.isData)
    {
      TFile *fIn=TFile::Open(lepEffUrl);
      lepEffH["e_sel"]=(TH2 *)fIn->Get("EGamma_SF2D");
      for(auto& it : lepEffH) it.second->SetDirectory(0);
      fIn->Close();
    }

  //B-TAG CALIBRATION
  TString btagUncUrl("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/CSVv2.csv");
  gSystem->ExpandPathName(btagUncUrl);
  std::vector<BTagCalibrationReader *> sfbReaders, sflReaders;
  TString btagEffExpUrl("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/expTageff.root");
  gSystem->ExpandPathName(btagEffExpUrl);
  std::map<TString, TGraphAsymmErrors *> expBtagEff, expBtagEffPy8;
  BTagSFUtil myBTagSFUtil;
  if(!ev.isData)
    {
      BTagCalibration btvcalib("csvv2", btagUncUrl.Data());
      sfbReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "mujets", "central") );
      sfbReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "mujets", "down") ); 
      sfbReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "mujets", "up") );

      sflReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "incl", "central") );
      sflReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "incl", "down") ); 
      sflReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "incl", "up") );
      
      TFile *beffIn=TFile::Open(btagEffExpUrl);
      expBtagEffPy8["b"]=(TGraphAsymmErrors *)beffIn->Get("b");
      expBtagEffPy8["c"]=(TGraphAsymmErrors *)beffIn->Get("c");
      expBtagEffPy8["udsg"]=(TGraphAsymmErrors *)beffIn->Get("udsg");
      beffIn->Close();
      
      TString btagExpPostFix("");
      if(isTTbar)
	{
	  if(filename.Contains("_herwig")) btagExpPostFix="_herwig";
	  if(filename.Contains("_scaleup")) btagExpPostFix="_scaleup";
	  if(filename.Contains("_scaledown")) btagExpPostFix="_scaledown";
	}
      btagEffExpUrl.ReplaceAll(".root",btagExpPostFix+".root");
      beffIn=TFile::Open(btagEffExpUrl);
      expBtagEff["b"]=(TGraphAsymmErrors *)beffIn->Get("b");
      expBtagEff["c"]=(TGraphAsymmErrors *)beffIn->Get("c");
      expBtagEff["udsg"]=(TGraphAsymmErrors *)beffIn->Get("udsg");
      beffIn->Close();
    }

  //JET ENERGY SCALE: https://twiki.cern.ch/twiki/bin/view/CMS/JECUncertaintySources#Summer15_uncertainties
  TString jecUncUrl("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/Fall15_25nsV2_DATA_UncertaintySources_AK4PFchs.txt");
  gSystem->ExpandPathName(jecUncUrl);
  //FactorizedJetCorrector *jetCorr=getFactorizedJetEnergyCorrector("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/jectxt",!ev.isData);
  std::vector<TString> jecUncSrcs;
  std::vector<JetCorrectionUncertainty*> jecUncs;
  if(runSysts)
    {
      jecUncSrcs.push_back("AbsoluteStat");        
      jecUncSrcs.push_back("AbsoluteScale");        
      jecUncSrcs.push_back("AbsoluteFlavMap");        
      jecUncSrcs.push_back("AbsoluteMPFBias");        
      jecUncSrcs.push_back("Fragmentation");        
      jecUncSrcs.push_back("SinglePionECAL");  
      jecUncSrcs.push_back("SinglePionHCAL"); 
      jecUncSrcs.push_back("TimeEta");
      jecUncSrcs.push_back("TimePt");
      jecUncSrcs.push_back("RelativeJEREC1");  
      jecUncSrcs.push_back("RelativeJEREC2"); 
      jecUncSrcs.push_back("RelativeJERHF");
      jecUncSrcs.push_back("RelativePtBB");    
      jecUncSrcs.push_back("RelativePtEC1");  
      jecUncSrcs.push_back("RelativePtEC2");   
      jecUncSrcs.push_back("RelativePtHF");  
      jecUncSrcs.push_back("RelativeFSR");
      jecUncSrcs.push_back("RelativeStatEC"); 
      jecUncSrcs.push_back("RelativeStatHF");
      jecUncSrcs.push_back("PileUpDataMC");    
      jecUncSrcs.push_back("PileUpPtRef");    
      jecUncSrcs.push_back("PileUpPtBB");     
      jecUncSrcs.push_back("PileUpPtEC1");      
      jecUncSrcs.push_back("PileUpPtEC2");      
      jecUncSrcs.push_back("PileUpPtHF");    
      jecUncSrcs.push_back("FlavorPureGluon"); 
      jecUncSrcs.push_back("FlavorPureQuark");
      jecUncSrcs.push_back("FlavorPureCharm"); 
      jecUncSrcs.push_back("FlavorPureBottom");
      for(size_t i=0; i<jecUncSrcs.size(); i++)
	{
	  JetCorrectorParameters *p = new JetCorrectorParameters(jecUncUrl.Data(), jecUncSrcs[i].Data());
	  jecUncs.push_back( new JetCorrectionUncertainty(*p) );
	}
    }

  //LIST OF SYSTEMATICS
  Int_t nGenSysts(0);
  std::vector<TString> expSysts;
  if(runSysts)
    {
      if(normH)
	for(Int_t xbin=1; xbin<=normH->GetNbinsX(); xbin++) 
	  nGenSysts += (normH->GetBinContent(xbin)!=0);	
      
      expSysts=jecUncSrcs;
      expSysts.push_back("JER");
      expSysts.push_back("Pileup");
      expSysts.push_back("MuTrigger");
      expSysts.push_back("MuEfficiency");
      expSysts.push_back("MuScale");
      expSysts.push_back("EleTrigger");
      expSysts.push_back("EleEfficiency");
      expSysts.push_back("EleScale");
      expSysts.push_back("BtagEff");
      expSysts.push_back("CtagEff");
      expSysts.push_back("LtagEff");
      expSysts.push_back("UncMET");
      expSysts.push_back("topPt");
      
      cout << "\t..." << nGenSysts << "/" << expSysts.size() << " generator level/experimental systematics will be considered" << endl;
    }


  //BOOK HISTOGRAMS
  std::map<TString, TH1 *> allPlots;
  std::map<TString, TH2 *> all2dPlots;
  allPlots["puwgtctr"] = new TH1F("puwgtctr","Weight sums",2,0,2);
  for(Int_t ij=0; ij<=4; ij++)
    {
      for(Int_t itag=-1; itag<=2; itag++)
	{
	  if(itag>ij) continue;
	  TString tag(itag<0 ? Form("%dj",ij) : Form("%dj%dt",ij,itag));
	  if(lumiMap.size()) allPlots["ratevsrun_"+tag] = new TH1F("ratevsrun_"+tag,";Run number; Events/pb",lumiMap.size(),0,lumiMap.size());
	  Int_t runCtr(0);
	  for(std::map<Int_t,Float_t>::iterator it=lumiMap.begin(); it!=lumiMap.end(); it++,runCtr++)
	    allPlots["ratevsrun_"+tag]->GetXaxis()->SetBinLabel(runCtr+1,Form("%d",it->first));
	  allPlots["lpt_"+tag]        = new TH1F("lpt_"+tag,";Transverse momentum [GeV];Events" ,20,0.,300.);
	  allPlots["lsip3d_"+tag]     = new TH1F("lsip3d_"+tag,";3d impact parameter significance;Events" ,40,0.,4.);
	  allPlots["lreliso_"+tag]     = new TH1F("lreliso_"+tag,";Relative isolation;Events" ,25,0.,0.5);
	  allPlots["leta_"+tag]       = new TH1F("leta_"+tag,";Pseudo-rapidity;Events" ,12,0.,3.);
	  allPlots["jpt_"+tag]        = new TH1F("jpt_"+tag,";Transverse momentum [GeV];Events" ,20,0.,300.);
	  allPlots["jeta_"+tag]       = new TH1F("jeta_"+tag,";Pseudo-rapidity;Events" ,12,0.,3.);
	  allPlots["ht_"+tag]         = new TH1F("ht_"+tag,";H_{T} [GeV];Events",40,0,800);
	  allPlots["csv_"+tag]        = new TH1F("csv_"+tag,";CSV discriminator;Events",100,0,1.0);
	  allPlots["rho_"+tag]        = new TH1F("rho_"+tag,";#rho [GeV];Events" ,20,0.,50.);
	  allPlots["nvtx_"+tag]       = new TH1F("nvtx_"+tag,";Vertex multiplicity;Events" ,80,0.,80.);
	  allPlots["nvtxup_"+tag]     = new TH1F("nvtxup_"+tag,";Vertex multiplicity;Events" ,80,0.,80.);
	  allPlots["nvtxdn_"+tag]     = new TH1F("nvtxdn_"+tag,";Vertex multiplicity;Events" ,80,0.,80.);
	  allPlots["metpt_"+tag]      = new TH1F("metpt_"+tag,";Missing transverse energy [GeV];Events" ,10,0.,200.);
	  allPlots["metphi_"+tag]     = new TH1F("metphi_" + tag,";MET #phi [rad];Events" ,50,-3.2,3.2);
	  allPlots["mttbar_"+tag]     = new TH1F("mttbar_"+tag,";#sqrt{#hat{s}} [GeV];Events" ,50,0.,1000.);
	  allPlots["mt_"+tag]         = new TH1F("mt_"+tag,";Transverse Mass [GeV];Events" ,20,0.,200.);
	  allPlots["minmlb_"+tag]     = new TH1F("minmlb_"+tag,";Mass(lepton,b) [GeV];Events" ,25,0.,250.);
	  allPlots["drlb_"+tag]       = new TH1F("drlb_"+tag,";#DeltaR(lepton,b);Events" ,25,0.,6.);
	  allPlots["passSIP3d_"+tag]  = new TH1F("passSIP3d_"+tag,";Pass SIP3d requirement;Events" ,2,0.,2.);
	  allPlots["RMPF_"+tag]       = new TH1F("RMPF_"+tag,";R_{MPF};Events" ,20,0.,2.);
	  allPlots["alpha_"+tag]      = new TH1F("alpha_"+tag,";#alpha;Events" ,20,0.,2.);
	  if(itag==-1)
	    {
	      allPlots["nbtags_"+tag]     = new TH1F("nbtags_"+tag,";Category;Events" ,3,0.,3.);
	      allPlots["nbtags_"+tag]->GetXaxis()->SetBinLabel(1, "=0b");
	      allPlots["nbtags_"+tag]->GetXaxis()->SetBinLabel(2, "=1b");
	      allPlots["nbtags_"+tag]->GetXaxis()->SetBinLabel(3, "#geq2b");
	    }

	  //shape analysis
	  if(runSysts)
	    {
	      //gen level systematics
	      if(normH)
		{
		  all2dPlots["metptshapes_"+tag+"_gen"]                   
		    = new TH2F("metptshapes_"+tag+"_gen", ";Missing transverse energy [GeV];Events" ,    10,0.,200., nGenSysts,0,nGenSysts);
		  all2dPlots["minmlbshapes_"+tag+"_gen"]               
		    = new TH2F("minmlbshapes_"+tag+"_gen", ";Mass(lepton,b) [GeV];Events" , 25,0.,250., nGenSysts,0,nGenSysts);
		  all2dPlots["RMPFshapes_"+tag+"_gen"]               
		    = new TH2F("RMPFshapes_"+tag+"_gen", ";R_{MPF};Events" , 20,0.,2., nGenSysts,0,nGenSysts);
		  if(itag==-1) 
		    all2dPlots["nbtagsshapes_"+tag+"_gen"] 
		      = new TH2F("nbtagsshapes_"+tag+"_gen", ";Category;Events" , 3, 0.,3.,   nGenSysts,0,nGenSysts);		  
		  for(Int_t igen=0; igen<nGenSysts; igen++)
		    {
		      TString label(normH->GetXaxis()->GetBinLabel(igen+3));
		      all2dPlots["metptshapes_"+tag+"_gen"]  ->GetYaxis()->SetBinLabel(igen+1,label);
		      all2dPlots["minmlbshapes_"+tag+"_gen"]->GetYaxis()->SetBinLabel(igen+1,label);
		      all2dPlots["RMPFshapes_"+tag+"_gen"]->GetYaxis()->SetBinLabel(igen+1,label);
		      if(itag!=-1) continue;
		      all2dPlots["nbtagsshapes_"+tag+"_gen"]->GetYaxis()->SetBinLabel(igen+1,label);
		    }
		}
	      
	      //experimental systematics
	      Int_t nExpSysts=expSysts.size();
	      if(nExpSysts>0)
		{
		  all2dPlots["metptshapes_"+tag+"_exp"]                  
		    = new TH2F("metptshapes_"+tag+"_exp",  ";Missing transverse energy [GeV];Events" ,   10,0.,200., 2*nExpSysts,0,2*nExpSysts);
		  all2dPlots["minmlbshapes_"+tag+"_exp"]              
		    = new TH2F("minmlbshapes_"+tag+"_exp", ";Mass(lepton,b) [GeV];Events" ,25,0.,250., 2*nExpSysts,0,2*nExpSysts);
		  all2dPlots["RMPFshapes_"+tag+"_exp"]               
		    = new TH2F("RMPFshapes_"+tag+"_exp", ";R_{MPF};Events" , 20,0.,2., 2*nExpSysts,0,2*nExpSysts);
		  if(itag==-1) 
		    all2dPlots["nbtagsshapes_"+tag+"_exp"] 
		      = new TH2F("nbtagsshapes_"+tag+"_exp", ";Category;Events" ,  3,0.,3.,    2*nExpSysts,0,2*nExpSysts);
		  for(Int_t isyst=0; isyst<nExpSysts; isyst++)
		    {
		      for(Int_t ivar=0; ivar<2; ivar++)
			{
			  TString label(expSysts[isyst] + (ivar==0 ? "Down" : "Up"));
			  all2dPlots["metptshapes_"+tag+"_exp"]->GetYaxis()->SetBinLabel(2*isyst+ivar+1, label);
			  all2dPlots["minmlbshapes_"+tag+"_exp"]->GetYaxis()->SetBinLabel(2*isyst+ivar+1,label);
			  all2dPlots["RMPFshapes_"+tag+"_exp"]->GetYaxis()->SetBinLabel(2*isyst+ivar+1,label);
			  if(itag!=-1) continue;
			  all2dPlots["nbtagsshapes_"+tag+"_exp"]->GetYaxis()->SetBinLabel(2*isyst+ivar+1,label);
			}
		    }
		}
	    }
	}
    }

  for (auto& it : allPlots)   { it.second->Sumw2(); it.second->SetDirectory(0); }
  for (auto& it : all2dPlots) { it.second->Sumw2(); it.second->SetDirectory(0); }

  //LOOP OVER EVENTS
  for (Int_t iev=0;iev<nentries;iev++)
    {
      t->GetEntry(iev);
      if(iev%10000==0) printf ("\r [%3.0f/100] done",100.*(float)(iev)/(float)(nentries));

      //select 1 good lepton
      std::vector<int> tightLeptonsIso, tightLeptonsNonIso, vetoLeptons;
      for(int il=0; il<ev.nl; il++)
	{
	  bool passTightKin(ev.l_pt[il]>30 && fabs(ev.l_eta[il])<2.1);
	  bool passVetoKin(ev.l_pt[il]>10 && fabs(ev.l_eta[il])<2.5);
	  bool passTightId(ev.l_id[il]==13 ? (ev.l_pid[il]>>1)&0x1  : (ev.l_pid[il]>>2)&0x1);
	  float relIso(ev.l_relIso[il]);
	  bool passIso( ev.l_id[il]==13 ? relIso<0.15 : (ev.l_pid[il]>>1)&0x1 );
	  bool passNonIso(relIso>0.4);
	  if( ev.l_id[il]==11 && (passIso || relIso<0.4) ) passNonIso=false;
	  bool passVetoIso(  ev.l_id[il]==13 ? relIso<0.25 : true); 
	  bool passSIP3d(ev.l_ip3dsig[il]<4);
	  if(channelSelection==21) passSIP3d=true;

	  if(passTightKin && passTightId && passSIP3d)
	    {
	      if(passIso)         tightLeptonsIso.push_back(il);
	      else if(passNonIso) tightLeptonsNonIso.push_back(il);
	    }
	  else if(passVetoKin && passVetoIso) vetoLeptons.push_back(il);
	}

      //one good lepton either isolated or in the non-isolated sideband or a Z candidate
      Int_t lepIdx=-1;
      Bool_t isZ(false),isZPassingSIP3d(false);
      TLorentzVector l1p4,l2p4,dilp4;
      if(tightLeptonsIso.size()==1)                                       lepIdx=tightLeptonsIso[0];
      else if (tightLeptonsIso.size()==0 && tightLeptonsNonIso.size()==1) lepIdx=tightLeptonsNonIso[0];
      else if(tightLeptonsIso.size()==2)
	{	  
	  l1p4.SetPtEtaPhiM(ev.l_pt[tightLeptonsIso[0]],ev.l_eta[tightLeptonsIso[0]],ev.l_phi[tightLeptonsIso[0]],ev.l_mass[tightLeptonsIso[0]]);
	  l2p4.SetPtEtaPhiM(ev.l_pt[tightLeptonsIso[1]],ev.l_eta[tightLeptonsIso[1]],ev.l_phi[tightLeptonsIso[1]],ev.l_mass[tightLeptonsIso[1]]);
	  dilp4=l1p4+l2p4;
	  if(ev.l_id[tightLeptonsIso[0]]==ev.l_id[tightLeptonsIso[1]]          && 
	     ev.l_charge[tightLeptonsIso[0]]*ev.l_charge[tightLeptonsIso[1]]<0 && 
	     fabs(dilp4.M()-91)<10 && 
	     dilp4.Pt()>30)
	    { 
	      isZ=true; 
	      isZPassingSIP3d=(ev.l_ip3dsig[0]<4 && ev.l_ip3dsig[1]<4);
	    }
	  lepIdx=tightLeptonsIso[0];
	}

      if(lepIdx<0) continue;
      
      //no extra isolated leptons
      if(vetoLeptons.size()>0) continue;
            
      //apply trigger requirement
      if(ev.l_id[lepIdx]==13)
	{
	  if(ev.isData  && (ev.muTrigger & 0x3)==0) continue;
	  if(!ev.isData && (ev.muTrigger & 0x3)==0) continue;
	}
      if(ev.l_id[lepIdx]==11)
	{ 
	  if( ((ev.elTrigger>>0)&0x1)==0 ) continue;
	}

      //select according to the lepton id/charge
      Int_t lid=ev.l_id[lepIdx];
      if(isZ) lid=2100+ev.l_id[lepIdx];
      else if(tightLeptonsNonIso.size()==1) lid=100*ev.l_id[lepIdx];
      if(channelSelection!=0)
	{
	  if(channelSelection==21) { if(!isZ) continue; }
	  else                     { if(lid!=channelSelection) continue; }
	}
      if(chargeSelection!=0  && ev.l_charge[lepIdx]!=chargeSelection) continue;

      //lepton kinematics
      TLorentzVector lp4;
      lp4.SetPtEtaPhiM(ev.l_pt[lepIdx],ev.l_eta[lepIdx],ev.l_phi[lepIdx],ev.l_mass[lepIdx]);

      //select jets
      Float_t htsum(0);
      TLorentzVector jetDiff(0,0,0,0);
      std::vector<TLorentzVector> bJets,lightJets;
      TLorentzVector visSystem(isZ ? dilp4 : lp4);
      int nbjets(0),ncjets(0),nljets(0),leadingJetIdx(-1);
      std::vector<int> resolvedJetIdx;
      std::vector<TLorentzVector> resolvedJetP4;
      for (int k=0; k<ev.nj;k++)
	{
	  //check kinematics
	  TLorentzVector jp4;
	  jp4.SetPtEtaPhiM(ev.j_pt[k],ev.j_eta[k],ev.j_phi[k],ev.j_mass[k]);
	  //jp4=updateJES(jp4,ev.j_rawsf[k],ev.j_area[k],ev.rho,ev.nvtx,jetCorr);

	  //cross clean with respect to leptons 
	  if(jp4.DeltaR(lp4)<0.5) continue;
	  if(isZ && jp4.DeltaR(l2p4)<0.5)continue;

	  //smear jet energy resolution for MC
	  //jetDiff -= jp4;
	  float genJet_pt(0);
	  if(ev.j_g[k]>-1) genJet_pt=ev.g_pt[ ev.j_g[k] ];
	  if(!ev.isData && genJet_pt>0) 
	    {
	      float jerSmear=getJetResolutionScales(jp4.Pt(),jp4.Eta(),genJet_pt)[0];
	      jp4 *= jerSmear;
	    }
	  //jetDiff += jp4;
	  resolvedJetIdx.push_back(k);
	  resolvedJetP4.push_back(jp4);

	  //require back-to-back configuration with Z
	  if(isZ && jp4.DeltaPhi(dilp4)<2.7) continue;

	  // re-inforce kinematics cuts
	  if(jp4.Pt()<30) continue;
	  if(fabs(jp4.Eta()) > 2.4) continue;
	  
	  if(leadingJetIdx<0) leadingJetIdx=k;
	  htsum += jp4.Pt();
	  if(bJets.size()+lightJets.size()<4) visSystem += jp4;

	  //b-tag
	  float csv = ev.j_csv[k];	  
	  bool isBTagged(csv>0.800);
	  if(!ev.isData)
	    {
	      float jptForBtag(jp4.Pt()>1000. ? 999. : jp4.Pt()), jetaForBtag(fabs(jp4.Eta()));
	      float expEff(1.0), jetBtagSF(1.0);
	      if(abs(ev.j_hadflav[k])==4) 
		{ 
		  ncjets++;
		  expEff    = expBtagEff["c"]->Eval(jptForBtag); 
		  jetBtagSF = sfbReaders[0]->eval( BTagEntry::FLAV_C, jetaForBtag, jptForBtag);
		  jetBtagSF *= expEff>0 ? expBtagEffPy8["c"]->Eval(jptForBtag)/expBtagEff["c"]->Eval(jptForBtag) : 0.;
		}
	      else if(abs(ev.j_hadflav[k])==5) 
		{ 
		  nbjets++;
		  expEff    = expBtagEff["b"]->Eval(jptForBtag); 
		  jetBtagSF = sfbReaders[0]->eval( BTagEntry::FLAV_B, jetaForBtag, jptForBtag);
		  jetBtagSF *= expEff>0 ? expBtagEffPy8["b"]->Eval(jptForBtag)/expBtagEff["b"]->Eval(jptForBtag) : 0.;
		}
	      else
		{
		  nljets++;
		  expEff    = expBtagEff["udsg"]->Eval(jptForBtag);
                  jetBtagSF = sflReaders[0]->eval( BTagEntry::FLAV_UDSG, jetaForBtag, jptForBtag);
		  jetBtagSF *= expEff> 0 ? expBtagEffPy8["udsg"]->Eval(jptForBtag)/expBtagEff["udsg"]->Eval(jptForBtag) : 0.;
		}
	      
	      //updated b-tagging decision with the data/MC scale factor
	      myBTagSFUtil.modifyBTagsWithSF(isBTagged,    jetBtagSF,     expEff);
	    }

	  //save jet
	  if(isBTagged) bJets.push_back(jp4);
	  else          lightJets.push_back(jp4);
	}

      //check if flavour splitting was required
      if(!ev.isData)
	{
	  if(flavourSplitting!=FlavourSplitting::NOFLAVOURSPLITTING)
	    {
	      if(flavourSplitting==FlavourSplitting::BSPLITTING)         { if(nbjets==0)    continue; }
	      else if(flavourSplitting==FlavourSplitting::CSPLITTING)    { if(ncjets==0 || nbjets!=0)    continue; }
	      else if(flavourSplitting==FlavourSplitting::UDSGSPLITTING) { if(nljets==0 || ncjets!=0 || nbjets!=0) continue; }
	    }
	}

      //MET and transverse mass
      TLorentzVector met(0,0,0,0);
      met.SetPtEtaPhiM(ev.met_pt[0],0,ev.met_phi[0],0.);
      met+=jetDiff;
      met.SetPz(0.); met.SetE(met.Pt());
      float mt( computeMT(isZ ? dilp4: lp4,met) );

      //compute neutrino kinematics
      neutrinoPzComputer.SetMET(met);
      neutrinoPzComputer.SetLepton(isZ ? dilp4 : lp4);
      
      float nupz=neutrinoPzComputer.Calculate();
      TLorentzVector neutrinoHypP4(met.Px(),met.Py(),nupz ,TMath::Sqrt(TMath::Power(met.Pt(),2)+TMath::Power(nupz,2)));
      visSystem+=neutrinoHypP4;

      //balancing variable and ISR control variable
      float RMPF(0.0),alpha(0.0);
      TVector2 metT(met.Px(),met.Py());
      TVector2 visT(isZ ? dilp4.Px() : lp4.Px(), isZ ? dilp4.Py() : lp4.Py());
      RMPF=1.0+(metT*visT)/visT.Mod2();
      for(size_t ij=0; ij<resolvedJetIdx.size(); ij++)
	{
	  Int_t k=resolvedJetIdx[ij];
	  if(k==leadingJetIdx) continue;
	  if(ev.j_pt[k]<15 || fabs(ev.j_eta[k])>3.0) continue;
	  alpha=ev.j_pt[k]/(isZ ? dilp4.Pt() : lp4.Pt());
	}

      //event weight
      float wgt(1.0);
      std::vector<float> puWeight(3,1.0),lepTriggerSF(3,1.0),lepSelSF(3,1.0), topPtWgt(3,1.0);
      if(!ev.isData)
	{
	  //update lepton selection scale factors, if found
	  TString prefix("m");
	  if(lid==11 || lid==1100) prefix="e";
	  if(lepEffH.find(prefix+"_sel")!=lepEffH.end() && !isZ)
	    {
	      for(UInt_t il=0; il<TMath::Min((UInt_t)1,(UInt_t)tightLeptonsIso.size()); il++)
		{
		  Int_t ilIdx=tightLeptonsIso[il];
		  float minEtaForEff( lepEffH[prefix+"_sel"]->GetXaxis()->GetXmin() ), maxEtaForEff( lepEffH[prefix+"_sel"]->GetXaxis()->GetXmax()-0.01 );
		  float etaForEff=TMath::Max(TMath::Min(float(fabs(ev.l_eta[ilIdx])),maxEtaForEff),minEtaForEff);
		  Int_t etaBinForEff=lepEffH[prefix+"_sel"]->GetXaxis()->FindBin(etaForEff);
		  		  
		  float minPtForEff( lepEffH[prefix+"_sel"]->GetYaxis()->GetXmin() ), maxPtForEff( lepEffH[prefix+"_sel"]->GetYaxis()->GetXmax()-0.01 );
		  float ptForEff=TMath::Max(TMath::Min(float(ev.l_pt[ilIdx]),maxPtForEff),minPtForEff);
		  Int_t ptBinForEff=lepEffH[prefix+"_sel"]->GetYaxis()->FindBin(ptForEff);
		  		  
		  float selSF(lepEffH[prefix+"_sel"]->GetBinContent(etaBinForEff,ptBinForEff));
		  float selSFUnc(lepEffH[prefix+"_sel"]->GetBinError(etaBinForEff,ptBinForEff));
		  lepSelSF[0]*=selSF;      lepSelSF[1]*=(selSF-selSFUnc);       lepSelSF[2]*=(selSF+selSFUnc);
		  
		  float trigSF(1.0), trigSFUnc(0.03);
		  if(prefix=="m")
		    {
		      trigSF=(lepEffH[prefix+"_trig"]->GetBinContent(etaBinForEff,ptBinForEff));
		      trigSFUnc=(lepEffH[prefix+"_trig"]->GetBinError(etaBinForEff,ptBinForEff));
		    }

		  lepTriggerSF[0]*=trigSF; lepTriggerSF[1]*=(trigSF-trigSFUnc); lepTriggerSF[2]*=(trigSF+trigSFUnc);
		}
	    }

	  Int_t ntops(0);
	  float ptsf(1.0);
	  for(Int_t igen=0; igen<ev.ngtop; igen++)
	    {
	      if(abs(ev.gtop_id[igen])!=6) continue;
	      ntops++;
	      ptsf *= TMath::Exp(0.156-0.00137*ev.gtop_pt[igen]);
	    }
	  if(ptsf>0 && ntops==2)
	    {
	      ptsf=TMath::Sqrt(ptsf);
	      topPtWgt[1]=1./ptsf;
	      topPtWgt[2]=ptsf;
	    }
	  
	  //update pileup weights, if found
	  if(puWgtGr.size())
	    {
	      puWeight[0]=puWgtGr[0]->Eval(ev.putrue);  
	      puWeight[1]=puWgtGr[1]->Eval(ev.putrue); 
	      puWeight[2]=puWgtGr[2]->Eval(ev.putrue);
	    }
	  
	  //update nominal event weight
	  float norm( normH ? normH->GetBinContent(1) : 1.0);
	  wgt=lepTriggerSF[0]*lepSelSF[0]*puWeight[0]*norm;
	  if(ev.ttbar_nw>0) wgt*=ev.ttbar_w[0];
	}

      //nominal selection control histograms
      int nJetsCat=TMath::Min((int)(bJets.size()+lightJets.size()),(int)4);
      int nBtagsCat=TMath::Min((int)(bJets.size()),(int)2);
      std::vector<TString> catsToFill(2,Form("%dj",nJetsCat));
      catsToFill[1]+= Form("%dt",nBtagsCat);
      allPlots["puwgtctr"]->Fill(0.,puWeight[0]!=0 ? wgt/puWeight[0] : 0.);
      allPlots["puwgtctr"]->Fill(1.,wgt);
      std::map<Int_t,Float_t>::iterator rIt=lumiMap.find(ev.run);	  
      for(size_t icat=0; icat<2; icat++)
	{
	  TString tag=catsToFill[icat];
	  if(rIt!=lumiMap.end()) 
	    {
	      Int_t runCtr=std::distance(lumiMap.begin(),rIt);
	      allPlots["ratevsrun_"+tag]->Fill(runCtr,1.e+6/rIt->second);
	    }

	  allPlots["lpt_"+tag]->Fill(isZ ? dilp4.Pt() : lp4.Pt(),wgt);
	  allPlots["leta_"+tag]->Fill(isZ ? dilp4.Eta() : lp4.Eta(),wgt);
	  allPlots["lsip3d_"+tag]->Fill(ev.l_ip3dsig[lepIdx],wgt);
	  allPlots["passSIP3d_"+tag]->Fill(isZPassingSIP3d,wgt);
	  allPlots["lreliso_"+tag]->Fill(ev.l_relIso[lepIdx],wgt);
	  allPlots["jpt_"+tag]->Fill(ev.j_pt[leadingJetIdx],wgt);
	  allPlots["jeta_"+tag]->Fill(fabs(ev.j_eta[leadingJetIdx]),wgt);
	  allPlots["csv_"+tag]->Fill(ev.j_csv[ leadingJetIdx ],wgt);
	  allPlots["nvtx_"+tag]->Fill(ev.nvtx,wgt);
	  allPlots["nvtxdn_"+tag]->Fill(ev.nvtx,wgt*puWeight[2]/puWeight[0]);
	  allPlots["nvtxup_"+tag]->Fill(ev.nvtx,wgt*puWeight[2]/puWeight[0]);
	  allPlots["rho_"+tag]->Fill(ev.rho,wgt);
	  allPlots["metpt_"+tag]->Fill(ev.met_pt[0],wgt);
	  allPlots["metphi_"+tag]->Fill(ev.met_phi[0],wgt);
	  allPlots["mt_"+tag]->Fill(mt,wgt);
	  allPlots["mttbar_"+tag]->Fill(visSystem.M(),wgt);
	  allPlots["ht_"+tag]->Fill(htsum,wgt);
	  allPlots["alpha_"+tag]->Fill(alpha,wgt);
	  if(alpha<0.3) allPlots["RMPF_"+tag]->Fill(RMPF,wgt);
	  if(icat==0) allPlots["nbtags_"+tag]->Fill(bJets.size(),wgt);
	  
	  if(bJets.size())
	    {
	      float mlb=(bJets[0]+(isZ ? dilp4 : lp4)).M();		  
	      float drlb=bJets[0].DeltaR( (isZ ? dilp4 : lp4) );
	      if(bJets.size()>1){
		float mlb2=(bJets[1]+(isZ ? dilp4 : lp4)).M();
		float drlb2=bJets[1].DeltaR( (isZ ? dilp4 : lp4) );
		if(mlb2<mlb)
		      {
			mlb=mlb2;
			drlb=drlb2;
		      }
		  }
	      allPlots["minmlb_"+tag]->Fill(mlb,wgt);		 		 
	      allPlots["drlb_"+tag]->Fill(drlb,wgt);		 		 
	    }	  
	}
      
      //ANALYSIS WITH SYSTEMATICS
      if(!runSysts) continue;
      
      //gen weighting systematics
      if(normH)
	{
	  float mlb(bJets.size() ? (bJets[0]+lp4).M() : 0.);
	  if(bJets.size()>1) mlb=TMath::Min( (float) mlb, (float)(bJets[1]+lp4).M() );
	  
	  for(Int_t igen=0; igen<nGenSysts; igen++)
	    {
	      for(size_t icat=0; icat<2; icat++)
		{
		  float newWgt = wgt*ev.ttbar_w[igen]/ev.ttbar_w[0];
		  
		  //for signal we only consider shapes and acceptance effects as it will be fit
		  if(isTTbar) 
		    newWgt *= normH->GetBinContent(igen+2)/normH->GetBinContent(2);

		  TString tag=catsToFill[icat];	 
		  all2dPlots["metptshapes_"+tag+"_gen"]->Fill(ev.met_pt[0],igen,newWgt);
		  all2dPlots["minmlbshapes_"+tag+"_gen"]->Fill(mlb,igen,newWgt);
		  if(alpha<0.3) all2dPlots["RMPFshapes_"+tag+"_gen"]->Fill(RMPF,igen,newWgt);
		  if(icat==0) all2dPlots["nbtagsshapes_"+tag+"_gen"]->Fill(bJets.size(),igen,newWgt);
		}
	    }
	}

      //experimental systematics
      for(size_t ivar=0; ivar<expSysts.size(); ivar++)
	{
	  TString varName=expSysts[ivar];
	  bool updateJEn(ivar<=jecUncSrcs.size());
	  bool updateBtag(varName.Contains("tagEff"));

	  //update jet selection, if needed
	  for(int isign=0; isign<2; isign++)
	    {
	      //weight systematics
	      float newWgt(wgt);
	      if(varName=="Pileup" && puWeight[0]!=0) newWgt*=(isign==0 ? puWeight[1]/puWeight[0] : puWeight[2]/puWeight[0]);
	      if(
		 (varName=="MuTrigger" && (lid==13 || lid==1300 || lid==2113)) ||
		 (varName=="EleTrigger" && (lid==11 || lid==1100 || lid==2111))
		 )
		newWgt *= (isign==0 ? lepTriggerSF[1]/lepTriggerSF[0] : lepTriggerSF[2]/lepTriggerSF[0]);
	      if(
		 (varName=="MuEfficiency" && (lid==13 ||  lid==1300 || lid==2113)) ||
		 (varName=="EleEfficiency" && (lid==11 || lid==1100 || lid==2111))
		 )
		newWgt *= (isign==0 ? lepSelSF[1]/lepSelSF[0] : lepSelSF[2]/lepSelSF[0]);
	      if(
		 varName=="topPt"
		 )
		newWgt *= topPtWgt[1+isign];

	      //lepton scale systematics
	      TLorentzVector varlp4(lp4),varl2p4(0,0,0,0),vardilp4(0,0,0,0);
	      if(isZ)
		varl2p4.SetPtEtaPhiM(ev.l_pt[tightLeptonsIso[1]],ev.l_eta[tightLeptonsIso[1]],ev.l_phi[tightLeptonsIso[1]],ev.l_mass[tightLeptonsIso[1]]);		
	      if(
		 (varName=="MuScale" && ( lid==13 || lid==1300 || lid==2113)) ||
		 (varName=="EleScale" && (lid==11 || lid==1100 || lid==2111))
		 )
		{
		  varlp4 = (1.0+(isign==0?-1.:1.)*getLeptonEnergyScaleUncertainty(lid,lp4.Pt(),lp4.Eta()) ) *lp4;
		  if(isZ) varl2p4 = (1.0+(isign==0?-1.:1.)*getLeptonEnergyScaleUncertainty(lid,varl2p4.Pt(),varl2p4.Eta()) ) *varl2p4;
		}
	      if(varlp4.Pt()<30 || fabs(varlp4.Eta())>2.1) continue;
	      if(isZ) 
		{
		  if(varl2p4.Pt()<30 || fabs(varl2p4.Eta())>2.1) continue;
		  vardilp4=varlp4+varl2p4;
		  if(fabs(vardilp4.M()-91)>10) continue;
		  if(vardilp4.Pt()<30) continue;
		}

	      //jets
	      std::vector<TLorentzVector> varBJets,varLightJets;
	      TLorentzVector jetDiff(0,0,0,0),jetSum(0,0,0,0);
	      if(!(updateJEn || updateBtag)) { varBJets=bJets; varLightJets=lightJets; }
	      else
		{
		  for (size_t ij=0; ij<resolvedJetIdx.size();ij++)
		    {
		      int k(resolvedJetIdx[ij]);
		      int jflav( abs(ev.j_hadflav[k]) ),jflavForJES( abs(ev.j_flav[k]) );
		      
		      //check kinematics
		      TLorentzVector jp4;
		      jp4.SetPtEtaPhiM(ev.j_pt[k],ev.j_eta[k],ev.j_phi[k],ev.j_mass[k]);
		      //jp4=updateJES(jp4,ev.j_rawsf[k],ev.j_area[k],ev.rho,ev.nvtx,jetCorr);
		      jetDiff -= jp4;
		      
		      //smear jet energy resolution for MC
		      float genJet_pt(0);
		      if(ev.j_g[k]>-1) genJet_pt=ev.g_pt[ ev.j_g[k] ];
		      if(!ev.isData && genJet_pt>0) 
			{
			  std::vector<float> jerSmear=getJetResolutionScales(jp4.Pt(),jp4.Eta(),genJet_pt);
			  if(varName=="JER") jp4 *= jerSmear[isign+1];
			  else               jp4 *= jerSmear[0];
			}
		      
		      //update jet energy correction
		      if(updateJEn && varName!="JER")
			{
			  jecUncs[ivar]->setJetEta(fabs(jp4.Eta()));
			  jecUncs[ivar]->setJetPt(jp4.Pt());
			  if(
			     (jflavForJES==4 && varName=="FlavorPureCharm") || 
			     (jflavForJES==5 && varName=="FlavorPureBottom") || 
			     (jflavForJES==21 && varName=="FlavorPureGluon") || 
			     (jflavForJES<4 && varName=="FlavorPureQuark") || 
			     !varName.Contains("FlavorPure") 
			     )
			    {
			      float jecuncVal=jecUncs[ivar]->getUncertainty(true);
			      jp4 *= (1.0+(isign==0?-1.:1.)*jecuncVal);
			    }
			}
		      
		      jetDiff += jp4;
		      jetSum += jp4;

		      //re-apply kinematics
		      if(jp4.Pt()<30) continue;
		      if(fabs(jp4.Eta()) > 2.4) continue;
		      
		      //b-tag
		      float csv = ev.j_csv[k];	  
		      bool isBTagged(csv>0.800);
		      if(!ev.isData)
			{
			  float jptForBtag(jp4.Pt()>1000. ? 999. : jp4.Pt()), jetaForBtag(fabs(jp4.Eta()));
			  float expEff(1.0), jetBtagSF(1.0);
			  if(jflav==4) 
			    { 
			      expEff        = expBtagEff["c"]->Eval(jptForBtag); 
			      int idx(0);
			      if(varName=="CtagEff") idx=(isign==0 ? 1 : 2);
			      jetBtagSF  = sfbReaders[idx]->eval( BTagEntry::FLAV_C, jetaForBtag, jptForBtag);
			      jetBtagSF *= expEff>0 ? expBtagEffPy8["c"]->Eval(jptForBtag)/expBtagEff["c"]->Eval(jptForBtag) : 0.;
			    }
			  else if(jflav==5)
			    { 
			      expEff=expBtagEff["b"]->Eval(jptForBtag); 
			      int idx(0);
			      if(varName=="BtagEff") idx=(isign==0 ? 1 : 2);
			      jetBtagSF=sfbReaders[idx]->eval( BTagEntry::FLAV_B, jetaForBtag, jptForBtag);
			      jetBtagSF *= expEff>0 ? expBtagEffPy8["b"]->Eval(jptForBtag)/expBtagEff["b"]->Eval(jptForBtag) : 0.;
			    }
			  else
			    {
			      expEff=expBtagEff["udsg"]->Eval(jptForBtag);
			      int idx(0);
			      if(varName=="LtagEff") idx=(isign==0 ? 1 : 2);
			      jetBtagSF=sflReaders[idx]->eval( BTagEntry::FLAV_UDSG, jetaForBtag, jptForBtag);
			      jetBtagSF *= expEff>0 ? expBtagEffPy8["udsg"]->Eval(jptForBtag)/expBtagEff["udsg"]->Eval(jptForBtag) : 0.;
			    }
	      
			  //updated b-tagging decision with the data/MC scale factor
			  myBTagSFUtil.modifyBTagsWithSF(isBTagged,    jetBtagSF,     expEff);
			}

		      //save jet
		      if(isBTagged) varBJets.push_back(jp4);
		      else          varLightJets.push_back(jp4);
		    }
		}
	      
	      //MET variations
	      TLorentzVector varMet(met);
	      varMet += jetDiff;
	      varMet += (varlp4-lp4);
	      if(varName=="UncMET") 
		{
		  TLorentzVector uncMet(met+varlp4+jetSum);
		  uncMet *= (1.0 + (isign==0?-1.:1.)*0.1);
		  varMet = uncMet-varlp4-jetSum;
		}
	      varMet.SetPz(0.); varMet.SetE(varMet.Pt());
	      // float varmt(computeMT(varlp4,varMet) );

	      //ready to fill histograms
	      float mlb(varBJets.size() ? (varBJets[0]+varlp4).M() : 0.);
	      if(varBJets.size()>1) mlb=TMath::Min( (float) mlb, (float)(varBJets[1]+varlp4).M() );
	      
	      //balancing variable
	      float varRMPF(0);
	      TVector2 metT(varMet.Px(),varMet.Py());
	      TVector2 visT(isZ ? vardilp4.Px() : varlp4.Px(), isZ ? vardilp4.Py() : varlp4.Py());
	      varRMPF=1.0+(metT*visT)/visT.Mod2();
	      
	      //update categories
	      int nvarJetsCat=TMath::Min((int)(varBJets.size()+varLightJets.size()),(int)4);
	      int nvarBtagsCat=TMath::Min((int)(varBJets.size()),(int)2);
	      std::vector<TString> varcatsToFill(2,Form("%dj",nvarJetsCat));
	      varcatsToFill[1]+= Form("%dt",nvarBtagsCat);
	      
	      for(size_t icat=0; icat<2; icat++)
                {
		  TString tag=varcatsToFill[icat];
		  all2dPlots["metptshapes_"+tag+"_exp"]->Fill(varMet.Pt(),2*ivar+isign,newWgt);
		  all2dPlots["minmlbshapes_"+tag+"_exp"]->Fill(mlb,2*ivar+isign,newWgt);
		  all2dPlots["RMPFshapes_"+tag+"_exp"]->Fill(varRMPF,2*ivar+isign,newWgt);
		  if(icat!=0) continue;
		  all2dPlots["nbtagsshapes_"+tag+"_exp"]->Fill(varBJets.size(),2*ivar+isign,newWgt);
		}
	    }
	}
    }
  
  //close input file
  f->Close();

  //save histos to file  
  TString selPrefix("");  
  if(flavourSplitting!=NOFLAVOURSPLITTING) selPrefix=Form("%d_",flavourSplitting);
  TString baseName=gSystem->BaseName(outname); 
  TString dirName=gSystem->DirName(outname);
  TFile *fOut=TFile::Open(dirName+"/"+selPrefix+baseName,"RECREATE");
  fOut->cd();
  for (auto& it : allPlots)  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  for (auto& it : all2dPlots)  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  fOut->Close();
}

