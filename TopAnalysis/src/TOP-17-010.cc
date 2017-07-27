#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>

#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TopLJets2015/TopAnalysis/interface/TOP-17-010.h"
#include "TopLJets2015/TopAnalysis/interface/LeptonEfficiencyWrapper.h"
#include "TopLJets2015/TopAnalysis/interface/CorrectionTools.h"

#include <vector>
#include <iostream>
#include <algorithm>

#include "TMath.h"
#include "TKey.h"

using namespace std;


//
void RunTop17010(TString filename,
		 TString outname,
		 Int_t channelSelection, 
		 Int_t chargeSelection, 
		 SelectionTool::FlavourSplitting flavourSplitting,
		 TH1F *normH, 
		 Bool_t runSysts,
		 TString era)
{

  //check file type from name
  bool isDataFile(filename.Contains("Data"));
  if(isDataFile) runSysts=false;
  bool isTTbar( filename.Contains("_TTJets") );

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
  TH1 *genPU=(TH1 *)f->Get("analysis/putrue");
  TH1 *triggerList=(TH1 *)f->Get("analysis/triggerList");
  TTree *t = (TTree*)f->Get("analysis/data");
  attachToMiniEventTree(t,ev,true);
  Int_t nentries(t->GetEntriesFast());

  //EVENT SELECTION WRAPPER
  SelectionTool evsel(filename,false,triggerList);

  //CORRECTIONS
  std::vector<RunPeriod_t> runPeriods=getRunPeriods(era);

  //lumi
  TH1F *ratevsrunH=0;
  std::map<Int_t,Float_t> lumiMap;
  if( isDataFile )  
    {
      std::pair<std::map<Int_t,Float_t>, TH1F *> result=parseLumiInfo(era);
      lumiMap   = result.first;
      ratevsrunH = result.second;
    }

  //pileup
  std::map<TString, std::vector<TGraph *> > puWgtGr;
  if( !isDataFile ) puWgtGr=getPileupWeightsMap(era,genPU);

  //b-fragmentation
  float bfragWgt(1.0);

  //lepton efficiencies
  LeptonEfficiencyWrapper lepEffH(isDataFile,era);

  //b-tagging
  BTagSFUtil myBTagSFUtil;
  std::map<TString, std::map<BTagEntry::JetFlavor, BTagCalibrationReader *> > btvsfReaders = getBTVcalibrationReadersMap(era, BTagEntry::OP_MEDIUM);
  std::map<BTagEntry::JetFlavor, TGraphAsymmErrors *> expBtagEff = readExpectedBtagEff(era);

  //JET ENERGY UNCERTAINTIES    
  TString jecUncUrl(era+"/Summer16_23Sep2016V4_MC_UncertaintySources_AK4PFchs.txt");
  gSystem->ExpandPathName(jecUncUrl);
  std::vector<TString> jecUncSrcs;
  std::vector<JetCorrectionUncertainty*> jecUncs;
  if(runSysts)
    {
      jecUncSrcs.push_back("AbsoluteStat");        
      jecUncSrcs.push_back("AbsoluteScale");        
      jecUncSrcs.push_back("AbsoluteMPFBias");        
      jecUncSrcs.push_back("Fragmentation");        
      jecUncSrcs.push_back("SinglePionECAL");  
      jecUncSrcs.push_back("SinglePionHCAL"); 
      jecUncSrcs.push_back("FlavorPureGluon"); 
      jecUncSrcs.push_back("FlavorPureQuark");
      jecUncSrcs.push_back("FlavorPureCharm"); 
      jecUncSrcs.push_back("FlavorPureBottom");
      jecUncSrcs.push_back("TimePtEta");
      jecUncSrcs.push_back("RelativeJEREC1");  
      jecUncSrcs.push_back("RelativeJEREC2"); 
      jecUncSrcs.push_back("RelativeJERHF");
      jecUncSrcs.push_back("RelativePtBB");    
      jecUncSrcs.push_back("RelativePtEC1");  
      jecUncSrcs.push_back("RelativePtEC2");   
      jecUncSrcs.push_back("RelativePtHF");  
      jecUncSrcs.push_back("RelativeBal");  
      jecUncSrcs.push_back("RelativeFSR");
      jecUncSrcs.push_back("RelativeStatFSR");
      jecUncSrcs.push_back("RelativeStatEC"); 
      jecUncSrcs.push_back("RelativeStatHF");
      jecUncSrcs.push_back("PileUpDataMC");    
      jecUncSrcs.push_back("PileUpPtRef");    
      jecUncSrcs.push_back("PileUpPtBB");     
      jecUncSrcs.push_back("PileUpPtEC1");      
      jecUncSrcs.push_back("PileUpPtEC2");      
      jecUncSrcs.push_back("PileUpPtHF");    
    }
  for(size_t i=0; i<jecUncSrcs.size(); i++)
    {
      JetCorrectorParameters *p = new JetCorrectorParameters(jecUncUrl.Data(), jecUncSrcs[i].Data());
      jecUncs.push_back( new JetCorrectionUncertainty(*p) );
    }

  //b-fragmentation weights
  std::map<TString, TGraph*> fragWeights=getBFragmentationWeights(era);
  std::map<TString, std::map<int, double> > semilepBRwgts=getSemilepBRWeights(era);
    
  //BOOK HISTOGRAMS
  std::map<TString, TH1 *> allPlots;
  addGenScanCounters(allPlots,f);
  std::map<TString, TH2 *> all2dPlots;
  allPlots["puwgtctr"] = new TH1F("puwgtctr","Weight sums",4,0,4);

  std::vector<TString> lfsVec = { "EE", "EM", "MM" };
  for(size_t ilfs=0; ilfs<lfsVec.size(); ilfs++)   
    { 
      TString tag(lfsVec[ilfs]);

      if(ratevsrunH) allPlots["ratevsrun_"+tag]=(TH1 *)ratevsrunH->Clone("ratevsrun_"+tag);    

      allPlots["nvtx_"+tag]  = new TH1F("nvtx_"+tag,";Vertex multiplicity;Events",30,0,30);
      allPlots["rho_"+tag]  = new TH1F("rho_"+tag,";#rho;Events",30,0,30);
      allPlots["mll_"+tag]  = new TH1F("mll_"+tag,";Dilepton invariant mass [GeV];Events",20,0,200);

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
    }
  for (auto& it : allPlots)   { it.second->Sumw2(); it.second->SetDirectory(0); }
  for (auto& it : all2dPlots) { it.second->Sumw2(); it.second->SetDirectory(0); }

  //LOOP OVER EVENTS
  for (Int_t iev=0;iev<nentries;iev++)
    {
      TString period("");
      t->GetEntry(iev);
      resetTopWidthEvent(twev);
      if(iev%10000==0) printf ("\r [%3.0f/100] done",100.*(float)(iev)/(float)(nentries));

      //assign a run period and correct the event accordingly
      float puWgt(1.0),puWgtUp(1.0),puWgtDn(1.0),topptsf(1.0);
      if(!ev.isData)
	{
	  period=assignRunPeriod(runPeriods);
	  puWgt   = puWgtGr[period][0]->Eval(ev.g_pu);
	  puWgtUp = puWgtGr[period][1]->Eval(ev.g_pu);
	  puWgtDn = puWgtGr[period][2]->Eval(ev.g_pu);

	  //top pt weighting
	  if(isTTbar)
	    {
	      for(Int_t igen=0; igen<ev.ngtop; igen++)
		{
		  if(abs(ev.gtop_id[igen])!=6) continue;
		  topptsf *= TMath::Exp(0.156-0.00137*ev.gtop_pt[igen]);
		}
	    }
	}

      //EVENT SELECTION
      TString chTag = evsel.flagFinalState(ev);
      if(chTag!="EE" && chTag!="EM" && chTag!="MM") continue;
      std::vector<Particle> &leptons=evsel.getSelLeptons();

      //select jets
      Int_t nbtags=0;
      std::vector<int> genJetsFlav,genJetsHadFlav, btagStatus;
      std::vector<TLorentzVector> jets,genJets;
      for (int k=0; k<ev.nj;k++)
	{
	  //check kinematics
	  TLorentzVector jp4;
	  jp4.SetPtEtaPhiM(ev.j_pt[k],ev.j_eta[k],ev.j_phi[k],ev.j_mass[k]);

	  //cross clean with leptons
	  bool overlapsWithLepton(false);
	  for(size_t il=0; il<leptons.size(); il++)
	    {
	      if(jp4.DeltaR(leptons[il].p4())>0.4) continue;
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
	  bool isBTagged(csv>0.8484),isBTaggedUp(isBTagged),isBTaggedDown(isBTagged);
	  if(!ev.isData)
	    {
	      float jptForBtag(jp4.Pt()>1000. ? 999. : jp4.Pt()), jetaForBtag(fabs(jp4.Eta()));
	      float expEff(1.0), jetBtagSF(1.0), jetBtagSFUp(1.0), jetBtagSFDown(1.0);

	      BTagEntry::JetFlavor hadFlav=BTagEntry::FLAV_UDSG; 
	      if(abs(ev.j_hadflav[k])==4) hadFlav=BTagEntry::FLAV_C; 
	      if(abs(ev.j_hadflav[k])==5) hadFlav=BTagEntry::FLAV_B;  
	     
	      expEff    = expBtagEff[hadFlav]->Eval(jptForBtag); 
	      //float py8corr(expEff>0 ? expBtagEffPy8[hadFlav]->Eval(jptForBtag)/expBtagEff[hadFlav]->Eval(jptForBtag) : 0.);
              
	      jetBtagSF = btvsfReaders[period][hadFlav]->eval_auto_bounds( "central", hadFlav, jetaForBtag, jptForBtag);
	      //jetBtagSF *= py8corr;
	      myBTagSFUtil.modifyBTagsWithSF(isBTagged,      jetBtagSF,      expEff);

	      jetBtagSFUp = btvsfReaders[period][hadFlav]->eval_auto_bounds( "up", hadFlav, jetaForBtag, jptForBtag);
	      //jetBtagSFUp *= py8corr;
	      myBTagSFUtil.modifyBTagsWithSF(isBTaggedUp,    jetBtagSFUp,    expEff);

	      jetBtagSFDown = btvsfReaders[period][hadFlav]->eval_auto_bounds( "down", hadFlav, jetaForBtag, jptForBtag);
	      //jetBtagSFDown *= py8corr;
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
	}

      //
      //event weight
      //
      float wgt(1.0);
      EffCorrection_t eSelCorrWgt(1.0,0.0),mSelCorrWgt(1.0,0.0), triggerCorrWgt(1.0,0.0);
      if(!ev.isData)
	{
	  //MC normalization weight
	  float norm( normH ? normH->GetBinContent(1) : 1.0);

	  //account for pu weights and effect on normalization
	  allPlots["puwgtctr"]->Fill(0.,1.0);
	  allPlots["puwgtctr"]->Fill(1.,puWgt);
	  allPlots["puwgtctr"]->Fill(2.,puWgtUp);
	  allPlots["puwgtctr"]->Fill(3.,puWgtDn);
	
          //trigger/id+iso efficiency corrections
          triggerCorrWgt=lepEffH.getTriggerCorrection(leptons,period);
          for(size_t il=0; il<2; il++)
            {
              EffCorrection_t selSF=lepEffH.getOfflineCorrection(leptons[il].id(),leptons[il].pt(),leptons[il].eta(),period);
              if(abs(leptons[il].id())==11)
                {
                  eSelCorrWgt.second = sqrt( pow(eSelCorrWgt.first*selSF.second,2)+pow(eSelCorrWgt.second*selSF.first,2));
                  eSelCorrWgt.first *= selSF.first;
                }
              else
                {
                  mSelCorrWgt.second = sqrt( pow(mSelCorrWgt.first*selSF.second,2)+pow(mSelCorrWgt.second*selSF.first,2));
                  mSelCorrWgt.first *= selSF.first;
                }
            }
          
          bfragWgt=computeBFragmentationWeight(ev,fragWeights["centralFrag"]);

	  //update nominal event weight
	  wgt=triggerCorrWgt.first*eSelCorrWgt.first*mSelCorrWgt.first*puWgt*bfragWgt*norm;
	  if(ev.g_nw>0) wgt*=ev.g_w[0];
	}
      
      //all done... physics
      //preselection
      if(chTag=="") continue;
      if(leptons[0].pt()<25 && leptons[1].pt()<25) continue;
      if(fabs(leptons[0].eta())>2.5 || fabs(leptons[1].eta())>2.5) continue;
      float mll((leptons[0].p4()+leptons[1].p4()).M());
      if(mll<12) continue;

      //nominal selection control histograms
      allPlots["nvtx_"+chTag]->Fill(ev.nvtx,wgt);
      allPlots["rho_"+chTag]->Fill(ev.rho,wgt);
      allPlots["mll_"+chTag]->Fill(mll,wgt);
      std::map<Int_t,Float_t>::iterator rIt=lumiMap.find(ev.run);
      if(rIt!=lumiMap.end() && ratevsrunH) allPlots["ratevsrun_"+chTag]->Fill(std::distance(lumiMap.begin(),rIt),1./rIt->second);
      allPlots["njets_"+chTag]->Fill(jets.size(),wgt);

      if(jets.size()<2) continue;

      //save leptons
      twev.nl=TMath::Min(2,(int)leptons.size());
      for(int il=0; il<twev.nl; il++)
	{
	  TString pf(Form("l%d",il));
	  allPlots[pf+"pt_"+chTag]->Fill(leptons[il].pt(),wgt);
	  allPlots[pf+"eta_"+chTag]->Fill(fabs(leptons[il].eta()),wgt);
	  if(nbtags>1)
	    {
	      allPlots[pf+"b2pt_"+chTag]->Fill(leptons[il].pt(),wgt);
	      allPlots[pf+"b2eta_"+chTag]->Fill(fabs(leptons[il].eta()),wgt);
	    }

	  twev.l_pt[il]=leptons[il].pt();
	  twev.l_eta[il]=leptons[il].eta();
	  twev.l_phi[il]=leptons[il].phi();
	  twev.l_m[il]=leptons[il].mass();
	  twev.l_id[il]=leptons[il].id();
	  twev.l_les[il]=ev.l_scaleUnc[leptons[il].originalReference()]/leptons[il].p4().E();
	  for(Int_t ig=0; ig<ev.ng; ig++)
	    {
	      if(abs(ev.g_id[ig])!=leptons[il].id()) continue;
	      TLorentzVector glp4;
	      glp4.SetPtEtaPhiM( ev.g_pt[ig], ev.g_eta[ig], ev.g_phi[ig], ev.g_m[ig]);
	      if(glp4.DeltaR( leptons[il].p4() ) > 0.3) continue;
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
	  twev.j_jes.push_back( std::vector<Float_t>() ); //[ij].clear();
	  if(jecUncs.size())
	    {
              for(size_t iju=0; iju<jecUncs.size(); iju++)
                {
                  jecUncs[iju]->setJetEta(twev.j_eta[ij]);
                  jecUncs[iju]->setJetPt(twev.j_pt[ij]);
                  Float_t unc=jecUncs[iju]->getUncertainty(true);
                  if(jecUncSrcs[iju]=="FlavorPureGluon" && abs(twev.gj_flav[ij])!=21) unc=0;
                  if(jecUncSrcs[iju]=="FlavorPureQuark" && (abs(twev.gj_flav[ij])>=4 || abs(twev.gj_flav[ij])==0)) unc=0;
                  if(jecUncSrcs[iju]=="FlavorPureCharm" && abs(twev.gj_flav[ij])!=4) unc=0;
                  if(jecUncSrcs[iju]=="FlavorPureBottom"&& abs(twev.gj_flav[ij])!=5) unc=0;
                  twev.j_jes[ij].push_back( unc );
                }
            }
	}

      allPlots["nbtags_"+chTag]->Fill(nbtags,wgt);

      if(chTag=="MM") twev.cat=13*13;
      if(chTag=="EM") twev.cat=11*13;
      if(chTag=="EE") twev.cat=11*11;
      twev.nw=15;
      twev.weight[0]=wgt;
      twev.weight[1]=wgt*puWgtUp/puWgt;
      twev.weight[2]=wgt*puWgtDn/puWgt;
      twev.weight[3]=wgt*(triggerCorrWgt.first+triggerCorrWgt.second)/triggerCorrWgt.first;
      twev.weight[4]=wgt*(triggerCorrWgt.first-triggerCorrWgt.second)/triggerCorrWgt.first;
      twev.weight[5]=wgt*(eSelCorrWgt.first+eSelCorrWgt.second)/eSelCorrWgt.first;
      twev.weight[6]=wgt*(eSelCorrWgt.first-eSelCorrWgt.second)/eSelCorrWgt.first;
      twev.weight[7]=wgt*(mSelCorrWgt.first+mSelCorrWgt.second)/mSelCorrWgt.first;
      twev.weight[8]=wgt*(mSelCorrWgt.first-mSelCorrWgt.second)/mSelCorrWgt.first;
      twev.weight[9]=wgt*topptsf;
      twev.weight[10]=wgt*computeSemilepBRWeight(ev,semilepBRwgts["semilepbrDown"]);
      twev.weight[11]=wgt*computeSemilepBRWeight(ev,semilepBRwgts["semilepbrUp"]);
      twev.weight[12]=wgt*computeBFragmentationWeight(ev,fragWeights["downFrag"])/bfragWgt;
      twev.weight[13]=wgt*computeBFragmentationWeight(ev,fragWeights["upFrag"])/bfragWgt;
      twev.weight[14]=wgt*computeBFragmentationWeight(ev,fragWeights["PetersonFrag"])/bfragWgt;
      if(ev.g_nw>0)
	{
	  twev.nw+=ev.g_nw;
	  for(int iw=1; iw<=ev.g_nw; iw++) twev.weight[14+iw]=wgt*ev.g_w[iw]/ev.g_w[0];
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
  //for(size_t i=0; i<50; i++) twev.j_jes.push_back( std::vector<Float_t>() );
  t->Branch("j_jes",   &twev.j_jes);
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
  for(int i=0; i<20; i++) twev.weight[i]=0;
  for(int i=0; i<2; i++) { twev.l_pt[i]=0;   twev.l_eta[i]=0;   twev.l_phi[i]=0;   twev.l_m[i]=0; twev.l_id[i]=0; twev.l_les[i]=0; twev.gl_pt[i]=0;   twev.gl_eta[i]=0;   twev.gl_phi[i]=0;   twev.gl_m[i]=0; twev.gl_id[i]=0; }
  for(int i=0; i<50; i++) { twev.j_pt[i]=0;   twev.j_eta[i]=0;   twev.j_phi[i]=0;   twev.j_m[i]=0; twev.j_btag[i]=0; twev.j_jer[i]=0;  twev.gj_pt[i]=0;   twev.gj_eta[i]=0;   twev.gj_phi[i]=0;   twev.gj_m[i]=0; twev.gj_flav[i]=0; twev.gj_hadflav[i]=0; }  twev.j_jes.clear(); //[i].clear(); } 
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
