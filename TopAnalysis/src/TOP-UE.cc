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
#include "TopLJets2015/TopAnalysis/interface/TOP-16-019.h"
#include "TopLJets2015/TopAnalysis/interface/LeptonEfficiencyWrapper.h"
#include "TopLJets2015/TopAnalysis/interface/TOP-UE.h"
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
void RunTopUE(TString filename,
	      TString outname,
	      Int_t channelSelection, 
	      Int_t chargeSelection, 
	      FlavourSplitting flavourSplitting,
	      TH1F *normH, 
	      Bool_t runSysts,
	      TString era)
{

  //check file type from name
  if(filename.Contains("SingleElectron") || filename.Contains("SingleMuon"))
    {
      cout << "Bailing out from analysing " << filename << endl;
      return;
    }
  if(filename.Contains("Data"))
    {
      runSysts=false;
    }
  bool isTTbar( filename.Contains("_TTJets") );

  //LEPTON EFFICIENCIES
  LeptonEfficiencyWrapper lepEffH(filename.Contains("Data13TeV"),era);

  //prepare output
  TopUE_t tue;
  TFile *fOut=TFile::Open(outname,"RECREATE");
  fOut->cd();
  TTree *outT=new TTree("tue","tue");
  createTopUETree(outT,tue);
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

  //BOOK CONTROL HISTOGRAMS
  std::map<TString, TH1 *> allPlots;
  addGenScanCounters(allPlots,f);

  allPlots["puwgtctr"] = new TH1F("puwgtctr","Weight sums",4,0,4);

  std::vector<TString> lfsVec = { "EE", "EM", "MM" };
  for(size_t ilfs=0; ilfs<lfsVec.size(); ilfs++)   
    { 
      TString tag(lfsVec[ilfs]);
      allPlots["nvtx_"+tag]   = new TH1F("nvtx_"+tag,";Vertex multiplicity;Events",40,0,40);
      allPlots["mll_"+tag]    = new TH1F("mll_"+tag,";Dilepton invariant mass [GeV];Events",50,0,400);
      allPlots["ptpos_"+tag]   = new TH1F("ptpos_"+tag,";Lepton transverse momentum [GeV];Events",50,20,200);
      allPlots["ptll_"+tag]    = new TH1F("ptll_"+tag,";Dilepton transverse momentum [GeV];Events",50,0,200);
      allPlots["sumpt_"+tag]   = new TH1F("sumpt_"+tag,";Transverse momentum sum [GeV];Events",50,40,300);
      allPlots["njets_"+tag]  = new TH1F("njets_"+tag,";Jet multiplicity;Events",4,2,6);
      allPlots["nbtags_"+tag] = new TH1F("nbtags_"+tag,";b-tag multiplicity;Events",5,0,5);
    }
  for (auto& it : allPlots)   { it.second->Sumw2(); it.second->SetDirectory(0); }

  //LOOP OVER EVENTS
  for (Int_t iev=0;iev<nentries;iev++)
    {
      t->GetEntry(iev);
      resetTopUE(tue);
      if(iev%10000==0) printf("\r [%3.0f/100] done",100.*(float)(iev)/(float)(nentries));

      //MC TRUTH
      float norm(1.0);
      std::vector<float> puWgts(3,1.0),topPtWgts(2,1.0);
      if(!ev.isData)
	{
	  //MC normalization weight
	  norm = ( normH ? normH->GetBinContent(1) : 1.0);

	  //top truth
	  bool passGenSel(true);
	  TLorentzVector ttGen(0,0,0,0), pseudottGen(0,0,0,0);
	  std::vector<TLorentzVector> pseudobJets,pseudoLeptons;
	  Int_t ntops(0),posLepton(0);
          float ptsf(1.0);
          for(Int_t igen=0; igen<ev.ngtop; igen++)
            {
	      TLorentzVector tp4(0,0,0,0);
	      tp4.SetPtEtaPhiM( ev.gtop_pt[igen], ev.gtop_eta[igen], ev.gtop_phi[igen], ev.gtop_m[igen] );

	      //test PID
	      if(abs(ev.gtop_id[igen])==6000)  
		{
		  pseudottGen += tp4;
		}
	      if(abs(ev.gtop_id[igen])==5000)  
		{ 
		  pseudobJets.push_back(tp4); 
		  if(tp4.Pt()<25 || fabs(tp4.Eta())>2.5) passGenSel=false; 
		}
	      if(abs(ev.gtop_id[igen])==11000 || abs(ev.gtop_id[igen])==13000)  
		{
		  pseudoLeptons.push_back(tp4);
		  if(ev.gtop_id[igen]<0) posLepton=(pseudoLeptons.size()-1);
		  if(tp4.Pt()<20 || fabs(tp4.Eta())>2.5) passGenSel=false; 
		}
              if(abs(ev.gtop_id[igen])==6)
		{
		  ntops++;
		  ttGen += tp4;
		  ptsf *= TMath::Exp(0.156-0.00137*tp4.Pt());
		}
	    }
	  	  
	  //save ttbar and pseudo-ttbar kinematics
	  tue.parton_ptttbar=ttGen.Pt();
	  tue.parton_phittbar=ttGen.Phi();
	  tue.gen_ptttbar=pseudottGen.Pt();
	  tue.gen_phittbar=pseudottGen.Phi();
	  if(pseudoLeptons.size()>1)
	    {
	      tue.gen_mll    = (pseudoLeptons[0]+pseudoLeptons[1]).M();
	      tue.gen_ptpos  = pseudoLeptons[posLepton].Pt();
	      tue.gen_phipos = pseudoLeptons[posLepton].Phi();
	      tue.gen_ptll   = (pseudoLeptons[0]+pseudoLeptons[1]).Pt();
	      tue.gen_phill  = (pseudoLeptons[0]+pseudoLeptons[1]).Phi();
	      tue.gen_sumpt  = pseudoLeptons[0].Pt()+pseudoLeptons[1].Pt();
	      tue.gen_dphill = TMath::Abs(pseudoLeptons[0].DeltaPhi(pseudoLeptons[1]));
	    }
	  

	  //further check if it two jets/two leptons are present and set the flag at generator level
	  if(pseudoLeptons.size()<2 || pseudobJets.size()<2) passGenSel=false;
	  tue.gen_passSel=passGenSel;

	  //count jets at generator level pt>25 |eta|<2.5
	  tue.gen_nb=0;
	  tue.gen_nj=0;
	  for(int igj=0; igj<ev.ng; igj++)
	    {
	      TLorentzVector gjp4(0,0,0,0);
	      gjp4.SetPtEtaPhiM( ev.g_pt[igj], ev.g_eta[igj], ev.g_phi[igj], ev.g_m[igj] );		      
	      if(gjp4.Pt()<25 || fabs(gjp4.Eta())>2.5) continue;
	      bool veto(false);
	      for(size_t il=0; il<pseudoLeptons.size(); il++)
		{
		  if(gjp4.DeltaR(pseudoLeptons[il])<0.1) veto=true;
		}
	      if(veto) continue;
	      if(abs(ev.g_id[igj])==5) tue.gen_nb++;
	      else                     tue.gen_nj++;
	    }

	  //save charged gen particles not associated with the b's or leptons
	  for(int igpf=0; igpf<ev.ngpf; igpf++)
	    {
	      if(ev.gpf_c[igpf]==0) continue;
	      if(fabs(ev.gpf_eta[igpf])>2.5) continue;

	      TLorentzVector gpfp4(0,0,0,0);
	      gpfp4.SetPtEtaPhiM( ev.gpf_pt[igpf], ev.gpf_eta[igpf], ev.gpf_phi[igpf], ev.gpf_m[igpf] );

	      //check if matched to pseudo-top decay products
	      bool matchedToPseudoTop(false);
	      int genJetIdx=ev.gpf_g[igpf];
	      if(genJetIdx>0)
		{
		  int genJetPID=ev.g_id[ genJetIdx ];
		  if( abs(genJetPID)==5 )
		    {
		      TLorentzVector gjp4(0,0,0,0);
		      gjp4.SetPtEtaPhiM( ev.g_pt[genJetIdx], ev.g_eta[genJetIdx], ev.g_phi[genJetIdx], ev.g_m[genJetIdx] );		      
		      for(size_t ib=0; ib<pseudobJets.size(); ib++)
			{
			  if(gjp4.DeltaR(pseudobJets[ib])<0.1) matchedToPseudoTop=true;
			}
		    }
		}
	      for(size_t il=0; il<pseudoLeptons.size(); il++)
		{
		  if(gpfp4.DeltaR(pseudoLeptons[il])<0.1) matchedToPseudoTop=true;
		}
	      if(matchedToPseudoTop) continue;
	      
	      tue.gen_pt[tue.gen_n]=gpfp4.Pt();
	      tue.gen_eta[tue.gen_n]=gpfp4.Eta();
	      tue.gen_phi[tue.gen_n]=gpfp4.Phi();
	      tue.gen_n++;
	    }
	  

	  //other weights
	  //top pt weights
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
	}
      
      //
      //RECO LEVEL analysis
      //

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
	      lp4.SetPtEtaPhiM(ev.l_pt[il],ev.l_eta[il],ev.l_phi[il],0.);
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


      //select the charged PF candidates (veto neutrals and charged within 0.01 of the selected leptons)
      std::vector<std::pair<int,TLorentzVector> > selTracks;
      for(int ipf = 0; ipf < ev.npf; ipf++) {
	if(ev.pf_c[ipf]==0) continue;
	if(ev.pf_pt[ipf]<0.5) continue;
	if(fabs(ev.pf_eta[ipf])>2.5) continue;
	TLorentzVector tkP4(0,0,0,0);
	tkP4.SetPtEtaPhiM(ev.pf_pt[ipf],ev.pf_eta[ipf],ev.pf_phi[ipf],0.);
	if(leptons.size()>0 && tkP4.DeltaR(leptons[0])<0.01) continue;
	if(leptons.size()>1 && tkP4.DeltaR(leptons[1])<0.01) continue;
	selTracks.push_back(std::pair<int,TLorentzVector>(ipf,tkP4) );
      }
      std::vector<int> selTracksBjetAssoc(selTracks.size(),0);

      //jet selection for different JES/JER/btag variations
      size_t maxVariations(runSysts ? 9 : 1);
      for(size_t ivar=0; ivar<maxVariations; ivar++)
	{
	  TLorentzVector met(0,0,0,0);
	  met.SetPtEtaPhiM(ev.met_pt[0],0,ev.met_phi[0],0);

	  TLorentzVector jetDiff(0,0,0,0);
	  std::vector< std::pair<int,TLorentzVector> > bJets, otherJets;
	  for (int k=0; k<ev.nj;k++)
	    {
	      //check kinematics
	      TLorentzVector jp4;
	      jp4.SetPtEtaPhiM(ev.j_pt[k],ev.j_eta[k],ev.j_phi[k],ev.j_mass[k]);
	      jetDiff += jp4;
	      
	      //cross clean with leptons
	      bool overlapsWithLepton(false);
	      for(size_t il=0; il<leptons.size(); il++)
		{
		  if(jp4.DeltaR(leptons[il])>0.4) continue;
		  overlapsWithLepton=true;
		}
	      if(overlapsWithLepton) continue;

	      //vary jet energy scale
	      if(jecUnc && (ivar==1 || ivar==2) ) 
		{
		  jecUnc->setJetEta(jp4.Eta());
		  jecUnc->setJetPt(jp4.Pt());
		  float jecUncVal=jecUnc->getUncertainty(true);
		  if(ivar==1)      jp4 *= (1+jecUncVal);
		  else if(ivar==2) jp4 *= (1-jecUncVal);
		  jetDiff -= jp4;
		}

	      //smear jet energy resolution for MC
	      float genJet_pt(0);
	      if(ev.j_g[k]>-1) genJet_pt=ev.g_pt[ ev.j_g[k] ];
	      if(!ev.isData && genJet_pt>0) 
		{
		  std::vector<float> jerSmear=getJetResolutionScales(jp4.Pt(),jp4.Eta(),genJet_pt);
		  if(ivar==3)      jp4 *= jerSmear[1];
		  else if(ivar==4) jp4 *= jerSmear[2];
		  else             jp4 *= jerSmear[0];
		  jetDiff -= jp4;
		}
	      
	      //b-tag
	      float csv = ev.j_csv[k];	  
	      bool isBTagged(csv>0.800);
	      if(!ev.isData && (ivar==5 || ivar==6 || ivar==7 || ivar==8) )
		{
		  float jptForBtag(jp4.Pt()>1000. ? 999. : jp4.Pt()), jetaForBtag(fabs(jp4.Eta()));
		  float expEff(1.0), jetBtagSF(1.0);
		  
		  if(abs(ev.j_hadflav[k])==4) 
		    { 	
		      int sfIdx(0);
		      if(ivar==5) sfIdx=1;
		      if(ivar==6) sfIdx=2;
		      expEff    = expBtagEff["c"]->Eval(jptForBtag); 
		      jetBtagSF = sfbReaders[sfIdx]->eval( BTagEntry::FLAV_C, jetaForBtag, jptForBtag);
		      jetBtagSF *= expEff>0 ? expBtagEffPy8["c"]->Eval(jptForBtag)/expBtagEff["c"]->Eval(jptForBtag) : 0.;
		    }
		  else if(abs(ev.j_hadflav[k])==5) 
		    { 
		      expEff    = expBtagEff["b"]->Eval(jptForBtag); 
		      int sfIdx(0);
		      if(ivar==5) sfIdx=1;
		      if(ivar==6) sfIdx=2;
		      jetBtagSF = sfbReaders[sfIdx]->eval( BTagEntry::FLAV_B, jetaForBtag, jptForBtag);
		      jetBtagSF *= expEff>0 ? expBtagEffPy8["b"]->Eval(jptForBtag)/expBtagEff["b"]->Eval(jptForBtag) : 0.;
		    }
		  else
		    {
		      expEff    = expBtagEff["udsg"]->Eval(jptForBtag);
		      int sfIdx(0);
		      if(ivar==7) sfIdx=1;
		      if(ivar==8) sfIdx=2;
		      jetBtagSF = sflReaders[sfIdx]->eval( BTagEntry::FLAV_UDSG, jetaForBtag, jptForBtag);
		      jetBtagSF *= expEff> 0 ? expBtagEffPy8["udsg"]->Eval(jptForBtag)/expBtagEff["udsg"]->Eval(jptForBtag) : 0.;
		    }
		  
		  //updated b-tagging decision with the data/MC scale factor
		  myBTagSFUtil.modifyBTagsWithSF(isBTagged,      jetBtagSF,      expEff);
		}
	      
	      
	      std::pair<int, TLorentzVector> selJet(k,jp4);
	      
	      //consider only b-jets above 30 GeV and fiducial in tracker
	      if(isBTagged && jp4.Pt()>30 && abs(jp4.Eta())<2.5)
		{
		  bJets.push_back(selJet);
		  
		  //check if the tracks will be associated to a b-jet
		  if(bJets.size()<=2)
		    {
		      for(size_t ipf=0; ipf<selTracks.size(); ipf++)
			{
			  if(ev.pf_j[ selTracks[ipf].first ] != k) continue;
			  selTracksBjetAssoc[ipf] |= ( 1 << ivar );
			}
		    }
		}
	      else if(jp4.Pt()>30 && abs(jp4.Eta())<2.5)
		{
		  otherJets.push_back(selJet);
		}
	    }
	  
	  //flag if passes alternative selection
	  tue.passSel |= ( (bJets.size()>=2)<<ivar );
	  
	  //count number of jets
	  tue.nb[ivar]=bJets.size();
	  tue.nj[ivar]=otherJets.size()+bJets.size();
	  
	  //reconstruct ttbar flight direction
	  if(bJets.size()>1 && leptons.size()>1)
	    {
	      met += jetDiff;
	      TLorentzVector rec_tt(leptons[0]+leptons[1]+bJets[0].second+bJets[1].second+met);
	      tue.ptttbar[ivar]=rec_tt.Pt();
	      tue.phittbar[ivar]=rec_tt.Phi();	  
	    }
	}

      //finalise dilepton requirements and reset pass selection flags
      //dilepton present, m(ll)>12, at least one with pt>30 and both within |eta|<2.5
      if(chTag=="") tue.passSel=0;
      else
	{
	  //finalize lepton selection
	  float mll=(leptons[0]+leptons[1]).M();
	  if(mll<12) tue.passSel=0;
	  if(leptons[0].Pt()<30 && leptons[1].Pt()<30)  tue.passSel=0;
	  if(fabs(leptons[0].Eta())>2.5 || fabs(leptons[1].Eta())>2.5) tue.passSel=0;
	  
	  //compute variations of the kinematics
	  size_t maxLepVariations(runSysts ? 5 : 1);	  
	  for(size_t ivar=0; ivar<maxLepVariations; ivar++)
	    {
	      std::vector<TLorentzVector> varLeptons;
	      int posLepton(0);
	      for(size_t il=0; il<2; il++)
		{
		  //determine positive lepton
		  if(ev.l_charge[ selLeptons[il] ]>0) posLepton=il;

		  //vary energy if needed
		  float les(0.0);
		  if(abs(selLeptonsId[il])==11 && (ivar==1 || ivar==2)) les=getLeptonEnergyScaleUncertainty(11,leptons[il].Pt(),leptons[il].Eta());
		  if(abs(selLeptonsId[il])==13 && (ivar==3 || ivar==4)) les=TMath::Abs(1-ev.l_scaleUnc[ selLeptons[il] ]);
		  float sign((ivar==1||ivar==3)?+1.0:-1.0);
		  varLeptons.push_back( leptons[il]*(1.0+sign*les) );
		}
	      
	      tue.mll[ivar]    = (varLeptons[0]+varLeptons[1]).M();
	      tue.ptpos[ivar]  = varLeptons[posLepton].Pt();
	      tue.phipos[ivar] = varLeptons[posLepton].Phi();
	      tue.ptll[ivar]   = (varLeptons[0]+varLeptons[1]).Pt();
	      tue.phill[ivar]  = (varLeptons[0]+varLeptons[1]).Phi();
	      tue.sumpt[ivar]  = varLeptons[0].Pt()+varLeptons[1].Pt();
	      tue.dphill[ivar] = TMath::Abs(varLeptons[0].DeltaPhi(varLeptons[1]));
	    }
	}
      
      //at this point we know if the event passes either the GEN level 
      //or one of the RECO level selections
      //if it fails both we don't care about it
      if(tue.passSel==0 && tue.gen_passSel==0) continue;

      //save PF cands
      tue.n=selTracks.size();
      for(size_t ipf=0; ipf<selTracks.size(); ipf++)
	{
	  tue.pt[ipf]         = selTracks[ipf].second.Pt();
	  tue.eta[ipf]        = selTracks[ipf].second.Eta();
	  tue.phi[ipf]        = selTracks[ipf].second.Phi();
	  tue.isInBFlags[ipf] = selTracksBjetAssoc[ipf];
	}

      //event weight
      float wgt(1.0);
      EffCorrection_t lepSelCorrWgt(1.0,0.0), triggerCorrWgt(1.0,0.0);
      if(!ev.isData)
	{
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
      
      
      //nominal selection control histograms
      if( (tue.passSel & 0x1) && chTag!="" )
	{
	  allPlots["nvtx_"+chTag]->Fill(ev.nvtx,wgt);
	  allPlots["mll_"+chTag]->Fill(tue.mll[0],wgt);
	  allPlots["ptpos_"+chTag]->Fill(tue.ptpos[0],wgt);
	  allPlots["ptll_"+chTag]->Fill(tue.ptll[0],wgt);
	  allPlots["sumpt_"+chTag]->Fill(tue.sumpt[0],wgt);
	  allPlots["nbtags_"+chTag]->Fill(tue.nb[0],wgt);
	  if(tue.nb[0]>1) allPlots["njets_"+chTag]->Fill(tue.nj[0],wgt);
	}
      
      //finalize ntuple
      tue.cat=0;
      if(chTag=="MM") tue.cat=13*13;
      if(chTag=="EM") tue.cat=11*13;
      if(chTag=="EE") tue.cat=11*11;
      tue.run=ev.run;
      tue.event=ev.event;
      tue.lumi=ev.lumi;
      tue.nvtx=ev.nvtx;
      if(runSysts)
	{
	  tue.nw=9;
	  tue.weight[0]=wgt;
	  tue.weight[1]=wgt*puWgts[1]/puWgts[0];
	  tue.weight[2]=wgt*puWgts[2]/puWgts[0];
	  tue.weight[3]=wgt*(triggerCorrWgt.first+triggerCorrWgt.second)/triggerCorrWgt.first;
	  tue.weight[4]=wgt*(triggerCorrWgt.first-triggerCorrWgt.second)/triggerCorrWgt.first;
	  tue.weight[5]=wgt*(lepSelCorrWgt.first+lepSelCorrWgt.second)/lepSelCorrWgt.first;
	  tue.weight[6]=wgt*(lepSelCorrWgt.first-lepSelCorrWgt.second)/lepSelCorrWgt.first;
	  tue.weight[7]=wgt*topPtWgts[0];
	  tue.weight[8]=wgt*topPtWgts[1];
	  if(ev.ttbar_nw>0)
	    {
	      tue.nw+=ev.ttbar_nw;
	      for(int iw=1; iw<=ev.ttbar_nw; iw++) tue.weight[8+iw]=wgt*ev.ttbar_w[iw]/ev.ttbar_w[0];
	    }
	}
      else
	{
	  tue.nw=1;
	  tue.weight[0]=1;
	}
   
      //all done, save it
      outT->Fill();
    }
  
  //save histos to file  
  fOut->cd();
  outT->Write();
  for (auto& it : allPlots)  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  fOut->Close();
}



//
void createTopUETree(TTree *t,TopUE_t &tue)
{
  //event category
  t->Branch("run",         &tue.run,         "run/I");
  t->Branch("event",       &tue.event,       "event/I");
  t->Branch("lumi",        &tue.lumi,        "lumi/I");
  t->Branch("cat",         &tue.cat,         "cat/I");
  t->Branch("gen_passSel", &tue.gen_passSel, "gen_passSel/I");
  t->Branch("passSel",     &tue.passSel,     "passSel/I");
  t->Branch("gen_nj",      &tue.gen_nj,      "gen_nj/I");
  t->Branch("gen_nb",      &tue.gen_nb,      "gen_nb/I");
  t->Branch("nj",           tue.nj,          "nj[9]/I");
  t->Branch("nb",           tue.nb,          "nb[9]/I");
  t->Branch("nvtx",        &tue.nvtx,        "nvtx/I");

  //event weights
  t->Branch("nw",     &tue.nw, "nw/I");
  t->Branch("weight",  tue.weight, "weight[nw]/F");

  //ptttbar
  t->Branch("parton_ptttbar",      &tue.parton_ptttbar,   "parton_ptttbar/F");
  t->Branch("parton_phittbar",     &tue.parton_phittbar,  "parton_phittbar/F");
  t->Branch("gen_ptttbar",         &tue.gen_ptttbar,      "gen_pttbar/F");
  t->Branch("gen_phittbar",        &tue.gen_phittbar,     "gen_phittbar/F");
  t->Branch("ptttbar",              tue.ptttbar,          "ptttbar[9]/F");
  t->Branch("phittbar",             tue.phittbar,         "phittbar[9]/F");

  //leptonic quantities
  t->Branch("ptpos",      tue.ptpos ,       "ptpos[5]/F");
  t->Branch("phipos",     tue.phipos ,      "phipos[5]/F");
  t->Branch("ptll",       tue.ptll ,        "ptll[5]/F");
  t->Branch("phill",      tue.phill ,       "phill[5]/F");
  t->Branch("mll",        tue.mll ,         "mll[5]/F");
  t->Branch("sumpt",      tue.sumpt ,       "sumpt[5]/F");
  t->Branch("dphill",     tue.dphill ,      "dphill[5]/F");
  t->Branch("gen_ptpos",  &tue.gen_ptpos ,  "gen_ptpos/F");
  t->Branch("gen_phipos", &tue.gen_phipos , "gen_phipos/F");
  t->Branch("gen_ptll",   &tue.gen_ptll ,   "gen_ptll/F");
  t->Branch("gen_phill",  &tue.gen_phill ,  "gen_phill/F");
  t->Branch("gen_mll",    &tue.gen_mll ,    "gen_mll/F");
  t->Branch("gen_sumpt",  &tue.gen_sumpt ,  "gen_sumpt/F");
  t->Branch("gen_dphill", &tue.gen_dphill , "gen_dphill/F");

  //charged particles
  t->Branch("n",          &tue.n,            "n/I");
  t->Branch("pt",          tue.pt ,          "pt[n]/F");
  t->Branch("eta",         tue.eta ,         "eta[n]/F");
  t->Branch("phi",         tue.phi ,         "phi[n]/F");
  t->Branch("isInBFlags",  tue.isInBFlags ,  "isInBFlags[n]/I");

  //gen charged particles
  t->Branch("gen_n",   &tue.gen_n,     "gen_n/I");
  t->Branch("gen_pt",   tue.gen_pt ,   "gen_pt[gen_n]/F");
  t->Branch("gen_eta",  tue.gen_eta ,  "gen_eta[gen_n]/F");
  t->Branch("gen_phi",  tue.gen_phi ,  "gen_phi[gen_n]/F");
}

//
void resetTopUE(TopUE_t &tue)
{
  //dummy event header
  tue.run=-1;  tue.lumi=0;  tue.event=0;              
  
  //reset selection flags
  tue.cat=0;   
  tue.gen_passSel=0;       
  tue.passSel=0; 

  //reset weights
  tue.nw=0;      
  
  //reset all MC truth
  tue.parton_ptttbar=0;       
  tue.parton_phittbar=0;
  tue.gen_ptttbar=0; 
  tue.gen_phittbar=0;
  tue.gen_ptpos=0;
  tue.gen_phipos=0;
  tue.gen_ptll=0;
  tue.gen_mll=0;
  tue.gen_sumpt=0;
  tue.gen_dphill=0;
  tue.gen_nj=0;             
  tue.gen_nb=0;
  
  //reset particle counters
  tue.n=0; 
  tue.gen_n=0; 
}
