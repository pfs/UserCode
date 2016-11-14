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
#include "TopLJets2015/TopAnalysis/interface/LeptonEfficiencyWrapper.h"
#include "TopLJets2015/TopAnalysis/interface/TOP-16-019.h"
#include "TopLJets2015/TopAnalysis/interface/BtagUncertaintyComputer.h"

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"


#include <vector>
#include <iostream>
#include <algorithm>

#include "TMath.h"
#include "TKey.h"

using namespace std;


//
void RunTop16019(TString filename,
		 TString outname,
		 Int_t channelSelection, 
		 Int_t chargeSelection, 
		 FlavourSplitting flavourSplitting,
		 TH1F *normH, 
		 Bool_t runSysts,
		 TString era)
{

  if(filename.Contains("SingleElectron") || filename.Contains("SingleMuon"))
    {
      cout << "Bailing out from analysing " << filename << endl;
      return;
    }

  bool isTTbar( filename.Contains("_TTJets") );

  //LEPTON EFFICIENCIES
  LeptonEfficiencyWrapper lepEffH(filename.Contains("Data13TeV"),era);
  bool hardCodedLES(era=="era2015");


  //prepare output
  TopWidthEvent_t twev;
  TString baseName=gSystem->BaseName(outname); 
  TString dirName=gSystem->DirName(outname);
  TFile *fOut=TFile::Open(dirName+"/"+baseName,"RECREATE");
  fOut->cd();
  TTree *outT=new TTree("twev","twev");
  createTopWidthEventTree(outT,twev);
  outT->SetDirectory(fOut);

  //READ TREE FROM FILE
  MiniEvent_t ev;
  TFile *f = TFile::Open(filename);
  TH1 *puTrue=(TH1 *)f->Get("analysis/putrue");
  puTrue->SetDirectory(0);
  puTrue->Scale(1./puTrue->Integral());
  TTree *t = (TTree*)f->Get("analysis/data");
  attachToMiniEventTree(t,ev,true);
  Int_t nentries(t->GetEntriesFast());
  t->GetEntry(0);
  bool requireEETriggers(false);
  if(ev.isData && filename.Contains("DoubleEG"))       requireEETriggers=true;
  bool requireMMTriggers(false);
  if(ev.isData && filename.Contains("DoubleMuon"))     requireMMTriggers=true;
  bool requireEMTriggers(false);
  if(ev.isData && filename.Contains("MuonEG"))         requireEMTriggers=true;
  //bool vetoEtrigger(false),vetoMtrigger(false);
  //if(ev.isData && filename.Contains("SingleElectron")) vetoMtrigger=true;


  //PILEUP WEIGHTING
  std::vector<TGraph *>puWgtGr;
  if(!ev.isData) puWgtGr=getPileupWeights(era,puTrue);
    
  //B-TAG CALIBRATION
  std::vector<BTagCalibrationReader *> sfbReaders, sflReaders;
  std::map<TString, TGraphAsymmErrors *> expBtagEff, expBtagEffPy8;
  BTagSFUtil myBTagSFUtil;
  if(!ev.isData)
    {
      TString btagUncUrl(era+"/btagSFactors.csv");
      gSystem->ExpandPathName(btagUncUrl);
      BTagCalibration btvcalib("csvv2", btagUncUrl.Data());
      sfbReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "mujets", "central") );
      sflReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "incl", "central") );
      sfbReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "mujets", "up") );
      sflReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "incl", "up") );
      sfbReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "mujets", "down") );
      sflReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "incl", "down") );

      TString btagEffExpUrl(era+"/expTageff.root");
      gSystem->ExpandPathName(btagEffExpUrl);
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

  //jet energy uncertainties
  TString jecUncUrl(era+"/jecUncertaintySources_AK4PFchs.txt");
  gSystem->ExpandPathName(jecUncUrl);
  JetCorrectorParameters *jecParam = new JetCorrectorParameters(jecUncUrl.Data(),"Total");
  JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty( *jecParam );

  //for data only get the lumi per run map                                                                                                                                                                         
  std::map<Int_t,Float_t> lumiMap;
  if(ev.isData) lumiMap=lumiPerRun(era);

  //BOOK HISTOGRAMS
  std::map<TString, TH1 *> allPlots;
  addGenScanCounters(allPlots,f);
  std::map<TString, TH2 *> all2dPlots;
  allPlots["puwgtctr"] = new TH1F("puwgtctr","Weight sums",4,0,4);

  std::vector<TString> lfsVec = { "EE", "EM", "MM" };
  for(size_t ilfs=0; ilfs<lfsVec.size(); ilfs++)   
    { 
      TString tag(lfsVec[ilfs]);
      allPlots["nvtx_"+tag]  = new TH1F("nvtx_"+tag,";Vertex multiplicity;Events",30,0,30);
      allPlots["rho_"+tag]  = new TH1F("rho_"+tag,";#rho;Events",30,0,30);
      allPlots["mll_"+tag]  = new TH1F("mll_"+tag,";Dilepton invariant mass [GeV];Events",20,0,200);
      if(lumiMap.size()) 
	{
	  allPlots["ratevsrun_"+tag] = new TH1F("ratevsrun_"+tag,";Run number; Events/pb",lumiMap.size(),0,lumiMap.size());
	  Int_t runCtr(0);
          for(std::map<Int_t,Float_t>::iterator it=lumiMap.begin(); it!=lumiMap.end(); it++,runCtr++)
            allPlots["ratevsrun_"+tag]->GetXaxis()->SetBinLabel(runCtr+1,Form("%d",it->first));
	}

      for(int i=0; i<2; i++)
	{
	  TString pf(Form("l%d",i));	  
	  for (int j=0; j<2; j++)
            {
              if(j==1) pf+= "b2";
	      allPlots[pf+"pt_"+tag]  = new TH1F(pf+"pt_"+tag,";Lepton p_{t} [GeV];Events",50,0,250);
	      allPlots[pf+"eta_"+tag]  = new TH1F(pf+"eta_"+tag,";Lepton pseudo-rapidity;Events",50,0,2.5);
	    }
	}
      allPlots["njets_"+tag]  = new TH1F("njets_"+tag,";Jet multiplicity;Events",5,0,5);
      allPlots["nbtags_"+tag]  = new TH1F("nbtags_"+tag,";b-tag multiplicity;Events",5,0,5);
      for(int i=0; i<6; i++)
	{
	  TString pf(Form("j%d",i));
	  for (int j=0; j<2; j++)
	    {
	      if(j==1) pf+= "b2";
	      allPlots[pf+"pt_"+tag]  = new TH1F(pf+"pt_"+tag,";Jet transverse momentum [GeV];Events",50,0,250);
	      allPlots[pf+"eta_"+tag]  = new TH1F(pf+"eta_"+tag,";Jet pseudo-rapidity;Events",50,0,4.7);
	    }
	}
      for(int ibj=0; ibj<2; ibj++)
	{
	  for(int k=0; k<2; k++)
	    {
	      TString pf(Form("b%d",ibj));
	      if(k>0) pf += "ch";
	      allPlots[pf+"const_"+tag]     = new TH1F(pf+"const_"+tag,     ";Constituent multiplicity;Events",50,0,50);
	      allPlots[pf+"pullm_"+tag]     = new TH1F(pf+"pullm_"+tag,     ";Pull magnitude;Events",20,0,0.05);
	      allPlots[pf+"pullangle_"+tag]   = new TH1F(pf+"pullangle_"+tag,   ";Pull angle [rad];Events",20,-3.16,3.16);
	    }
	}
    }
  for (auto& it : allPlots)   { it.second->Sumw2(); it.second->SetDirectory(0); }
  for (auto& it : all2dPlots) { it.second->Sumw2(); it.second->SetDirectory(0); }

  //LOOP OVER EVENTS
  for (Int_t iev=0;iev<nentries;iev++)
    {
      t->GetEntry(iev);
      resetTopWidthEvent(twev);
      if(iev%10000==0) printf ("\r [%3.0f/100] done",100.*(float)(iev)/(float)(nentries));

      //select good leptons
      std::vector<TLorentzVector> leptons;
      std::vector<int> selLeptons,selLeptonsId;
      for(int il=0; il<ev.nl; il++)
	{
	  bool passTightKin(ev.l_pt[il]>20 && fabs(ev.l_eta[il])<2.5);
	  bool passTightId(ev.l_id[il]==13 ? (ev.l_pid[il]>>1)&0x1  : (ev.l_pid[il]>>2)&0x1);
	  float relIso(ev.l_relIso[il]);
	  bool passTightIso( ev.l_id[il]==13 ? relIso<0.15 : (ev.l_pid[il]>>1)&0x1 );
	  if(passTightKin && passTightId && passTightIso) 
	    {
	      selLeptons.push_back(il);
	      selLeptonsId.push_back(ev.l_id[il]);
	      TLorentzVector lp4;
	      lp4.SetPtEtaPhiM(ev.l_pt[il],ev.l_eta[il],ev.l_phi[il],ev.l_mass[il]);
	      leptons.push_back(lp4);
	    }
	}
      
      //check if triggers have fired
      bool hasEETrigger(((ev.elTrigger>>1)&0x1)!=0 || ((ev.elTrigger>>4)&0x1)!=0);
      bool hasMMTrigger(((ev.muTrigger>>2)&0x3)!=0);
      bool hasEMTrigger(((ev.elTrigger>>2)&0x3)!=0);
      if(!ev.isData)
	{ 
	  hasEETrigger=true;
	  hasMMTrigger=true;
	  hasEMTrigger=true;
	}
      else
	{
	  if(requireEETriggers) { hasMMTrigger=false; hasEMTrigger=false; }
	  if(requireMMTriggers) { hasEETrigger=false; hasEMTrigger=false; }
	  if(requireEMTriggers) { hasEETrigger=false; hasMMTrigger=false; }
	}

      //decide the channel
      TString chTag("");
      if(selLeptons.size()>=2)
	{
	  if(abs(ev.l_id[ selLeptons[0] ]*ev.l_id[ selLeptons[1] ])==11*13      && hasEMTrigger) chTag="EM";
	  else if(abs(ev.l_id[ selLeptons[0] ]*ev.l_id[ selLeptons[1] ])==13*13 && hasMMTrigger) chTag="MM";
	  else if(abs(ev.l_id[ selLeptons[0] ]*ev.l_id[ selLeptons[1] ])==11*11 && hasEETrigger) chTag="EE";
	}


      //select jets
      std::vector<int> genJetsFlav,genJetsHadFlav, btagStatus;
      std::vector<TLorentzVector> jets,genJets;
      Int_t nbtags=0;
      std::vector<JetPullInfo_t> bJetPulls;
      for (int k=0; k<ev.nj;k++)
	{
	  //check kinematics
	  TLorentzVector jp4;
	  jp4.SetPtEtaPhiM(ev.j_pt[k],ev.j_eta[k],ev.j_phi[k],ev.j_mass[k]);

	  //cross clean with leptons
	  bool overlapsWithLepton(false);
	  for(size_t il=0; il<leptons.size(); il++)
	    {
	      if(jp4.DeltaR(leptons[il])>0.4) continue;
	      overlapsWithLepton=true;
	    }
	  if(overlapsWithLepton) continue;

	  //smear jet energy resolution for MC
	  float genJet_pt(0);
	  if(ev.j_g[k]>-1) genJet_pt=ev.g_pt[ ev.j_g[k] ];
	  if(!ev.isData && genJet_pt>0) 
	    {
	      float jerSmear=getJetResolutionScales(jp4.Pt(),jp4.Eta(),genJet_pt)[0];
	      jp4 *= jerSmear;
	    }

	  //b-tag
	  float csv = ev.j_csv[k];	  
	  bool isBTagged(csv>0.800),isBTaggedUp(isBTagged),isBTaggedDown(isBTagged);
	  if(!ev.isData)
	    {
	      float jptForBtag(jp4.Pt()>1000. ? 999. : jp4.Pt()), jetaForBtag(fabs(jp4.Eta()));
	      float expEff(1.0), jetBtagSF(1.0), jetBtagSFUp(1.0), jetBtagSFDown(1.0);

	      if(abs(ev.j_hadflav[k])==4) 
		{ 	
		  expEff    = expBtagEff["c"]->Eval(jptForBtag); 
		  jetBtagSF = sfbReaders[0]->eval( BTagEntry::FLAV_C, jetaForBtag, jptForBtag);
		  jetBtagSF *= expEff>0 ? expBtagEffPy8["c"]->Eval(jptForBtag)/expBtagEff["c"]->Eval(jptForBtag) : 0.;
		  jetBtagSFUp = sfbReaders[1]->eval( BTagEntry::FLAV_C, jetaForBtag, jptForBtag);
		  jetBtagSFUp *= expEff>0 ? expBtagEffPy8["c"]->Eval(jptForBtag)/expBtagEff["c"]->Eval(jptForBtag) : 0.;
		  jetBtagSFDown = sfbReaders[2]->eval( BTagEntry::FLAV_C, jetaForBtag, jptForBtag);
		  jetBtagSFDown *= expEff>0 ? expBtagEffPy8["c"]->Eval(jptForBtag)/expBtagEff["c"]->Eval(jptForBtag) : 0.;
		}
	      else if(abs(ev.j_hadflav[k])==5) 
		{ 
		  expEff    = expBtagEff["b"]->Eval(jptForBtag); 
		  jetBtagSF = sfbReaders[0]->eval( BTagEntry::FLAV_B, jetaForBtag, jptForBtag);
		  jetBtagSF *= expEff>0 ? expBtagEffPy8["b"]->Eval(jptForBtag)/expBtagEff["b"]->Eval(jptForBtag) : 0.;
		  jetBtagSFUp = sfbReaders[1]->eval( BTagEntry::FLAV_B, jetaForBtag, jptForBtag);
		  jetBtagSFUp *= expEff>0 ? expBtagEffPy8["b"]->Eval(jptForBtag)/expBtagEff["b"]->Eval(jptForBtag) : 0.;
		  jetBtagSFDown = sfbReaders[2]->eval( BTagEntry::FLAV_B, jetaForBtag, jptForBtag);
		  jetBtagSFDown *= expEff>0 ? expBtagEffPy8["b"]->Eval(jptForBtag)/expBtagEff["b"]->Eval(jptForBtag) : 0.;
		}
	      else
		{
		  expEff    = expBtagEff["udsg"]->Eval(jptForBtag);
                  jetBtagSF = sflReaders[0]->eval( BTagEntry::FLAV_UDSG, jetaForBtag, jptForBtag);
		  jetBtagSF *= expEff> 0 ? expBtagEffPy8["udsg"]->Eval(jptForBtag)/expBtagEff["udsg"]->Eval(jptForBtag) : 0.;
                  jetBtagSFUp = sflReaders[1]->eval( BTagEntry::FLAV_UDSG, jetaForBtag, jptForBtag);
		  jetBtagSFUp *= expEff> 0 ? expBtagEffPy8["udsg"]->Eval(jptForBtag)/expBtagEff["udsg"]->Eval(jptForBtag) : 0.;
                  jetBtagSFDown = sflReaders[2]->eval( BTagEntry::FLAV_UDSG, jetaForBtag, jptForBtag);
		  jetBtagSFDown *= expEff> 0 ? expBtagEffPy8["udsg"]->Eval(jptForBtag)/expBtagEff["udsg"]->Eval(jptForBtag) : 0.;
		}
	      
	      //updated b-tagging decision with the data/MC scale factor
	      myBTagSFUtil.modifyBTagsWithSF(isBTagged,      jetBtagSF,      expEff);
	      myBTagSFUtil.modifyBTagsWithSF(isBTaggedUp,    jetBtagSFUp,    expEff);
	      myBTagSFUtil.modifyBTagsWithSF(isBTaggedDown,  jetBtagSFDown,  expEff);
	    }

	  //consider only jets above 30 GeV
	  if(jp4.Pt()<30) continue;

	  //mc truth for this jet
	  Int_t hadFlav=ev.j_hadflav[k];
	  Int_t flav=ev.j_flav[k];
	  TLorentzVector gjp4(0,0,0,0);
	  if(ev.j_g[k]>=0)
	    {
	      int gidx=ev.j_g[k];
	      gjp4.SetPtEtaPhiM( ev.g_pt[gidx], ev.g_eta[gidx], ev.g_phi[gidx], ev.g_m[gidx] );
	    }

	  jets.push_back(jp4);
	  if(fabs(jp4.Eta())>2.5) { isBTagged=false; isBTaggedUp=false; isBTaggedDown=false; }
	  int btagStatusWord(isBTagged | (isBTaggedUp<<1) | (isBTaggedDown<<2));
	  
	  btagStatus.push_back(btagStatusWord);
	  genJets.push_back(gjp4);
	  genJetsFlav.push_back(flav); 
	  genJetsHadFlav.push_back(hadFlav);
	  nbtags += isBTagged;
	  if(isBTagged) bJetPulls.push_back( getPullVector(ev,k) );	  
	}

      //
      //event weight
      //
      float wgt(1.0);
      std::vector<float> puWgts(3,1.0),topPtWgts(2,1.0);
      EffCorrection_t lepSelCorrWgt(1.0,0.0), triggerCorrWgt(1.0,0.0);
      if(!ev.isData)
	{
	  //MC normalization weight
	  float norm( normH ? normH->GetBinContent(1) : 1.0);

	  //top pt
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
              topPtWgts[0]=1./ptsf;
              topPtWgts[1]=ptsf;
            }


	  //account for pu weights and effect on normalization
	  allPlots["puwgtctr"]->Fill(0.,1.0);
	  for(size_t iwgt=0; iwgt<3; iwgt++)
	    {
	      puWgts[iwgt]=puWgtGr[iwgt]->Eval(ev.putrue);  
	      allPlots["puwgtctr"]->Fill(iwgt+1,puWgts[iwgt]);
	    }
	
	  if(chTag!="")
	    {
	      //trigger/id+iso efficiency corrections
	      triggerCorrWgt=lepEffH.getTriggerCorrection(selLeptonsId,leptons);
	      for(size_t il=0; il<2; il++)
		{
		  EffCorrection_t selSF=lepEffH.getOfflineCorrection(selLeptonsId[il],leptons[il].Pt(),leptons[il].Eta());
		  lepSelCorrWgt.second = sqrt( pow(lepSelCorrWgt.first*selSF.second,2)+pow(lepSelCorrWgt.second*selSF.first,2));
		  lepSelCorrWgt.first *= selSF.first;
		}
	    }
	  
	  //update nominal event weight
	  wgt=triggerCorrWgt.first*lepSelCorrWgt.first*puWgts[0]*norm;
	  if(ev.ttbar_nw>0) wgt*=ev.ttbar_w[0];
	}

      //all done... physics
      //preselection
      if(chTag=="") continue;
      if(leptons[0].Pt()<30 && leptons[1].Pt()<30) continue;
      if(fabs(leptons[0].Eta())>2.1 && fabs(leptons[1].Eta())>2.1) continue;

      //nominal selection control histograms
      allPlots["nvtx_"+chTag]->Fill(ev.nvtx,wgt);
      allPlots["rho_"+chTag]->Fill(ev.rho,wgt);
      float mll((leptons[0]+leptons[1]).M());
      allPlots["mll_"+chTag]->Fill(mll,wgt);
      std::map<Int_t,Float_t>::iterator rIt=lumiMap.find(ev.run);
      if(rIt!=lumiMap.end())
	{
	  Int_t runCtr=std::distance(lumiMap.begin(),rIt);
	  allPlots["ratevsrun_"+chTag]->Fill(runCtr,1.e+6/rIt->second);
	}

      if(mll<12) continue;
      allPlots["njets_"+chTag]->Fill(jets.size(),wgt);

      if(jets.size()<2) continue;

      //save leptons
      twev.nl=TMath::Min(2,(int)leptons.size());
      for(int il=0; il<twev.nl; il++)
	{
	  TString pf(Form("l%d",il));
	  allPlots[pf+"pt_"+chTag]->Fill(leptons[il].Pt(),wgt);
	  allPlots[pf+"eta_"+chTag]->Fill(fabs(leptons[il].Eta()),wgt);
	  if(nbtags>1)
	    {
	      allPlots[pf+"b2pt_"+chTag]->Fill(leptons[il].Pt(),wgt);
	      allPlots[pf+"b2eta_"+chTag]->Fill(fabs(leptons[il].Eta()),wgt);
	    }

	  twev.l_pt[il]=leptons[il].Pt();
	  twev.l_eta[il]=leptons[il].Eta();
	  twev.l_phi[il]=leptons[il].Phi();
	  twev.l_m[il]=leptons[il].M();
	  twev.l_id[il]=ev.l_id[ selLeptons[il] ];
	  twev.l_les[il]=getLeptonEnergyScaleUncertainty(twev.l_id[il],twev.l_pt[il],twev.l_eta[il]);
	  if(!hardCodedLES && twev.l_id[il]==13) twev.l_les[il]=TMath::Abs(1-ev.l_scaleUnc[ selLeptons[il] ]);

	  for(Int_t ig=0; ig<ev.ng; ig++)
	    {
	      if(abs(ev.g_id[ig])!=ev.l_id[ selLeptons[il] ]) continue;
	      TLorentzVector glp4;
	      glp4.SetPtEtaPhiM( ev.g_pt[ig], ev.g_eta[ig], ev.g_phi[ig], ev.g_m[ig]);
	      if(glp4.DeltaR( leptons[il] ) > 0.3) continue;
	      twev.gl_id[il]=ev.g_id[ig];
	      twev.gl_pt[il]=ev.g_pt[ig];
	      twev.gl_eta[il]=ev.g_eta[ig];
	      twev.gl_phi[il]=ev.g_phi[ig];
	      twev.gl_m[il]=ev.g_m[ig];
	    }
	}
      
      //save jets
      twev.nj=jets.size();
      for(int ij=0; ij<(int)jets.size(); ij++)
	{
	  TString pf(Form("j%d",ij));
	  if(ij<6)
	    {
	      allPlots[pf+"pt_"+chTag]->Fill(jets[ij].Pt(),wgt);
	      allPlots[pf+"eta_"+chTag]->Fill(fabs(jets[ij].Eta()),wgt);
	      if(nbtags>1)
		{
		  allPlots[pf+"b2pt_"+chTag]->Fill(jets[ij].Pt(),wgt);
		  allPlots[pf+"b2eta_"+chTag]->Fill(fabs(jets[ij].Eta()),wgt);
		}

	    }
	  twev.j_pt[ij]=jets[ij].Pt();
	  twev.j_eta[ij]=jets[ij].Eta();
	  twev.j_phi[ij]=jets[ij].Phi();
	  twev.j_m[ij]=jets[ij].M();
	  twev.j_btag[ij]=btagStatus[ij];
	  twev.gj_flav[ij]=genJetsFlav[ij];
	  twev.gj_hadflav[ij]=genJetsHadFlav[ij];
	  twev.gj_pt[ij]=genJets[ij].Pt();
	  twev.gj_eta[ij]=genJets[ij].Eta();
	  twev.gj_phi[ij]=genJets[ij].Phi();
	  twev.gj_m[ij]=genJets[ij].M();	  
	  twev.j_jer[ij]=1.0;
	  if(twev.gj_pt[ij]>0)
	    {
	      std::vector<float> jerSmear=getJetResolutionScales(twev.j_pt[ij],twev.j_eta[ij],twev.gj_pt[ij]);
	      twev.j_jer[ij]=(jerSmear[0]>0 ? jerSmear[1]/jerSmear[0] : 1.0);
	    }
	  twev.j_jes[ij]=1.0;
	  if(jecUnc)
	    {
	      jecUnc->setJetEta(twev.j_eta[ij]);
	      jecUnc->setJetPt(twev.j_pt[ij]);
	      twev.j_jes[ij]=jecUnc->getUncertainty(true);
	    }
	  
	}

      allPlots["nbtags_"+chTag]->Fill(nbtags,wgt);
      if(nbtags>1)
	{
	  for(Int_t ibj=0; ibj<2; ibj++)
	    {
	      TString pf(Form("b%d",ibj));
	      allPlots[pf+"const_"+chTag]->Fill(bJetPulls[ibj].n,wgt);
	      allPlots[pf+"pullm_"+chTag]->Fill(bJetPulls[ibj].pull.Mod(),wgt);
	      allPlots[pf+"pullangle_"+chTag]->Fill(TMath::ATan2(bJetPulls[ibj].pull.Px(),bJetPulls[ibj].pull.Py()),wgt);
	      allPlots[pf+"chconst_"+chTag]->Fill(bJetPulls[ibj].nch,wgt);
	      allPlots[pf+"chpullm_"+chTag]->Fill(bJetPulls[ibj].chPull.Mod(),wgt);
	      allPlots[pf+"chpullangle_"+chTag]->Fill(TMath::ATan2(bJetPulls[ibj].chPull.Px(),bJetPulls[ibj].chPull.Py()),wgt);
	    }
	}

      if(chTag=="MM") twev.cat=13*13;
      if(chTag=="EM") twev.cat=11*13;
      if(chTag=="EE") twev.cat=11*11;
      twev.nw=9;
      twev.weight[0]=wgt;
      twev.weight[1]=wgt*puWgts[1]/puWgts[0];
      twev.weight[2]=wgt*puWgts[2]/puWgts[0];
      twev.weight[3]=wgt*(triggerCorrWgt.first+triggerCorrWgt.second)/triggerCorrWgt.first;
      twev.weight[4]=wgt*(triggerCorrWgt.first-triggerCorrWgt.second)/triggerCorrWgt.first;
      twev.weight[5]=wgt*(lepSelCorrWgt.first+lepSelCorrWgt.second)/lepSelCorrWgt.first;
      twev.weight[6]=wgt*(lepSelCorrWgt.first-lepSelCorrWgt.second)/lepSelCorrWgt.first;
      twev.weight[7]=wgt*topPtWgts[0];
      twev.weight[8]=wgt*topPtWgts[1];
      if(ev.ttbar_nw>0)
	{
	  twev.nw+=ev.ttbar_nw;
	  for(int iw=1; iw<=ev.ttbar_nw; iw++) twev.weight[8+iw]=wgt*ev.ttbar_w[iw]/ev.ttbar_w[0];
	}
      twev.nt=0;
      twev.met_pt=ev.met_pt[0];
      twev.met_phi=ev.met_phi[0];
      if(ev.ngtop>0)
	{
	  for(int i=0; i<ev.ngtop; i++)
	    {
	      int absid(abs(ev.gtop_id[i]));
	      if(absid!=6 && absid!=1000006 && absid!=1000022) continue;
	      twev.t_pt[twev.nt]=ev.gtop_pt[i];
	      twev.t_eta[twev.nt]=ev.gtop_eta[i];
	      twev.t_phi[twev.nt]=ev.gtop_phi[i];
	      twev.t_m[twev.nt]=ev.gtop_m[i];
	      twev.t_id[twev.nt]=ev.gtop_id[i];	     
	      twev.nt++;
	      if(twev.nt>10) break;
	    }

	  for(Int_t ig=0; ig<ev.ng; ig++)
	    {
	      int absid=abs(ev.g_id[ig]);
	      if(absid!=5 && absid!=13 && absid!=11) continue;
	      twev.t_pt[twev.nt]=ev.g_pt[ig];
	      twev.t_eta[twev.nt]=ev.g_eta[ig];
	      twev.t_phi[twev.nt]=ev.g_phi[ig];
	      twev.t_m[twev.nt]=ev.g_m[ig];
	      twev.t_id[twev.nt]=ev.g_id[ig];	     
	      twev.nt++;
	      if(twev.nt>10) break;
	    }
	}
      outT->Fill();
    }
  
  //close input file
  //f->Close();

  //save histos to file  
  fOut->cd();
  outT->Write();
  for (auto& it : allPlots)  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  for (auto& it : all2dPlots)  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  fOut->Close();
}



//
void createTopWidthEventTree(TTree *t,TopWidthEvent_t &twev)
{
  //event category
  t->Branch("cat", &twev.cat,"cat/I");

  //event weights
  t->Branch("nw",  &twev.nw, "nw/I");
  t->Branch("weight",  twev.weight, "weight[nw]/F");

  //met
  t->Branch("met_pt",  &twev.met_pt, "met_pt/F");
  t->Branch("met_phi",  &twev.met_phi, "met_phi/F");

  //leptons
  t->Branch("nl",  &twev.nl, "nl/I");
  t->Branch("l_pt",  twev.l_pt ,  "l_pt[nl]/F");
  t->Branch("l_eta", twev.l_eta , "l_eta[nl]/F");
  t->Branch("l_phi", twev.l_phi , "l_phi[nl]/F");
  t->Branch("l_m",   twev.l_m ,   "l_m[nl]/F");
  t->Branch("l_les",   twev.l_les ,   "l_les[nl]/F");
  t->Branch("l_id",   twev.l_id ,   "l_id[nl]/I");
  t->Branch("gl_pt",  twev.gl_pt ,  "gl_pt[nl]/F");
  t->Branch("gl_eta", twev.gl_eta , "gl_eta[nl]/F");
  t->Branch("gl_phi", twev.gl_phi , "gl_phi[nl]/F");
  t->Branch("gl_m",   twev.gl_m ,   "gl_m[nl]/F");
  t->Branch("gl_id",  twev.gl_id ,  "gl_id[nl]/I");

  //jets
  t->Branch("nj",  &twev.nj, "nj/I");
  t->Branch("j_pt",  twev.j_pt ,  "j_pt[nj]/F");
  t->Branch("j_eta", twev.j_eta , "j_eta[nj]/F");
  t->Branch("j_phi", twev.j_phi , "j_phi[nj]/F");
  t->Branch("j_m",   twev.j_m ,   "j_m[nj]/F");
  t->Branch("j_jer",   twev.j_jer ,   "j_jer[nj]/F");
  t->Branch("j_jes",   twev.j_jes ,   "j_jes[nj]/F");
  t->Branch("j_btag",   twev.j_btag ,   "j_btag[nj]/I");
  t->Branch("gj_pt",  twev.gj_pt ,  "gj_pt[nj]/F");
  t->Branch("gj_eta", twev.gj_eta , "gj_eta[nj]/F");
  t->Branch("gj_phi", twev.gj_phi , "gj_phi[nj]/F");
  t->Branch("gj_m",   twev.gj_m ,   "gj_m[nj]/F");
  t->Branch("gj_flav",  twev.gj_flav ,  "gj_flav[nj]/I");
  t->Branch("gj_hadflav",  twev.gj_hadflav ,  "gj_hadflav[nj]/I");

  //mc truth
  t->Branch("nt",  &twev.nt, "nt/I");
  t->Branch("t_pt",  twev.t_pt ,  "t_pt[nt]/F");
  t->Branch("t_eta", twev.t_eta , "t_eta[nt]/F");
  t->Branch("t_phi", twev.t_phi , "t_phi[nt]/F");
  t->Branch("t_m",   twev.t_m ,   "t_m[nt]/F");
  t->Branch("t_id",  twev.t_id ,  "t_id[nt]/I");
}

//
void resetTopWidthEvent(TopWidthEvent_t &twev)
{
  twev.cat=0;   twev.nw=0;   twev.nl=0;   twev.nj=0;   twev.nt=0;
  twev.met_pt=0; twev.met_phi=0;
  for(int i=0; i<10; i++) twev.weight[i]=0;
  for(int i=0; i<2; i++) { twev.l_pt[i]=0;   twev.l_eta[i]=0;   twev.l_phi[i]=0;   twev.l_m[i]=0; twev.l_id[i]=0; twev.l_les[i]=0; twev.gl_pt[i]=0;   twev.gl_eta[i]=0;   twev.gl_phi[i]=0;   twev.gl_m[i]=0; twev.gl_id[i]=0; }
  for(int i=0; i<50; i++) { twev.j_pt[i]=0;   twev.j_eta[i]=0;   twev.j_phi[i]=0;   twev.j_m[i]=0; twev.j_btag[i]=0; twev.j_jer[i]=0; twev.j_jes[i]=0; twev.gj_pt[i]=0;   twev.gj_eta[i]=0;   twev.gj_phi[i]=0;   twev.gj_m[i]=0; twev.gj_flav[i]=0; twev.gj_hadflav[i]=0; } 
  for(int i=0; i<10; i++) { twev.t_pt[i]=0;   twev.t_eta[i]=0;   twev.t_phi[i]=0;   twev.t_m[i]=0; twev.t_id[i]=0; }
}

//
void addGenScanCounters(std::map<TString, TH1 *> &plotColl,TFile *fIn)
{
  TH1 *normH=(TH1 *)fIn->Get("analysis/fidcounter0");
  TIter nextkey( fIn->GetDirectory("analysis")->GetListOfKeys() );
  TKey *key;
  while ( (key = (TKey*)nextkey())) {    
    TObject *obj = key->ReadObj();
    TString name(obj->GetName());
    if(!name.Contains("mstop")) continue;
    plotColl[name]=(TH1 *)obj->Clone();
    plotColl[name]->SetDirectory(0);
    if(normH) plotColl[name]->SetBinContent(1,normH->GetBinContent(1)/plotColl[name]->GetBinContent(1));
  }
}
