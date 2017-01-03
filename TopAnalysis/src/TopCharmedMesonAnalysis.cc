#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>

#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"
#include "TopLJets2015/TopAnalysis/interface/ObjectTools.h"
#include "TopLJets2015/TopAnalysis/interface/TopCharmedMesonAnalysis.h"
#include "TopLJets2015/TopAnalysis/interface/LeptonEfficiencyWrapper.h"
#include "TopLJets2015/TopAnalysis/interface/BtagUncertaintyComputer.h"

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"


#include <vector>
#include <iostream>

#include "TMath.h"
#include "TopQuarkAnalysis/TopTools/interface/MEzCalculator.h"

using namespace std;

//
void RunTopCharmedMesonAnalysis(TString filename,
				TString outname,
				Int_t channelSelection, 
				Int_t chargeSelection, 
				FlavourSplitting flavourSplitting,
				TH1F *normH, 
				Bool_t runSysts,
				TString era,
				Bool_t debug)
{
  if(debug) cout << "in RunTopCharmedMesonAnalysis" << endl;

  bool isTTbar( filename.Contains("_TTJets") );

  //LEPTON EFFICIENCIES
  LeptonEfficiencyWrapper lepEffH(filename.Contains("Data13TeV"),era);
  

  //READ TREE FROM FILE
  MiniEvent_t ev;
  //TopWidthEvent_t ev;
  TFile *f = TFile::Open(filename);
  //TTree *t = (TTree*)f->Get("twev");
  TH1 *puTrue=(TH1 *)f->Get("analysis/putrue");
  puTrue->SetDirectory(0);
  puTrue->Scale(1./puTrue->Integral());
  TTree *t = (TTree*)f->Get("analysis/data");
  attachToMiniEventTree(t,ev,true);
  Int_t nentries(t->GetEntriesFast());
  t->GetEntry(0);
  if(ev.isData) runSysts=false;
  bool requireEletrigger(false);
  if(ev.isData && filename.Contains("SingleElectron")) requireEletrigger=true;
  bool requireMutrigger(false);
  if(ev.isData && filename.Contains("SingleMuon"))     requireMutrigger=true;
  bool requireEETriggers(false);
  if(ev.isData && filename.Contains("DoubleEG"))       requireEETriggers=true;
  bool requireMMTriggers(false);
  if(ev.isData && filename.Contains("DoubleMuon"))     requireMMTriggers=true;
  bool requireEMTriggers(false);
  if(ev.isData && filename.Contains("MuonEG"))         requireEMTriggers=true;

  cout << "...producing " << outname << " from " << nentries << " events" << (runSysts ? " syst variations will be considered" : "") << endl;

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


  //JET ENERGY SCALE: https://twiki.cern.ch/twiki/bin/view/CMS/JECUncertaintySources#Summer15_uncertainties
  //jet energy uncertainties
  //TString jecUncUrl(era+"/jecUncertaintySources_AK4PFchs.txt");
  //gSystem->ExpandPathName(jecUncUrl);
  //JetCorrectorParameters *jecParam = new JetCorrectorParameters(jecUncUrl.Data(),"Total");
  //JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty( *jecParam );
  
  //LIST OF SYSTEMATICS
  
  //BOOK HISTOGRAMS
  std::map<TString, TH1 *> allPlots;
  //addGenScanCounters(allPlots,f); FIXME
  std::vector<TString> lfsVec = { "_all", "_e", "_ee", "_em", "_mm", "_m" }; 
  allPlots["puwgtctr"] = new TH1F("puwgtctr","Weight sums",2,0,2);

  for(int i = 0; i < (int)lfsVec.size(); i++) {
    TString tag(lfsVec[i]);
    allPlots["nlp"+tag]     = new TH1F("nlp"+tag,";N_{l};Events" ,4,0.,4.);
    allPlots["lp_pt"+tag] = new TH1F("lp_pt"+tag,";Leading Lepton P_{T} [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["l2p_pt"+tag] = new TH1F("l2p_pt"+tag,";Sub-leading Lepton P_{T} [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["dilp_pt"+tag] = new TH1F("dilp_pt"+tag,";Lepton P_{T} [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["dilp_m"+tag] = new TH1F("dilp_m"+tag,";M_{ll} [GeV];Events / 10 GeV", 20, 0,200);
    allPlots["charge"+tag] = new TH1F("charge"+tag,";Charge(l_{1}*l_{2});Events", 5,-2,2);
    allPlots["dR"+tag] = new TH1F("dR"+tag,";dR;Events / 0.05", 20,0.0,1.);

    allPlots["MET"+tag] = new TH1F("MET"+tag,";MET [GeV];Events / 20 GeV", 10,0,200);

    allPlots["nj"+tag]     = new TH1F("nj"+tag,";N_{jets} (P_{T} > 30 GeV);Events" ,8,2.,10.);
    allPlots["nbj"+tag]     = new TH1F("nbj"+tag,";N_{b-jets} (CSV > 0.8);Events" ,4,1.,5.);
    allPlots["j_pt"+tag] = new TH1F("j_pt"+tag,";Leading light Jet P_{T} [GeV];Events / 20 GeV", 15, 0,300);
    allPlots["bj_pt"+tag] = new TH1F("bj_pt"+tag,";Leading b Jet P_{T} [GeV];Events / 20 GeV", 15, 0,300);
    

    allPlots["npf"+tag]     = new TH1F("npf"+tag,";Track multiplicity;Jets" ,20,0.,20.);
    allPlots["pfid"+tag]     = new TH1F("pfid"+tag,";Track id;Jets" ,440,-220.,220.);
    allPlots["massJPsi"+tag]     = new TH1F("massJPsi"+tag,";M(ll) [GeV];Events / 20 MeV" ,100,2,4);
    allPlots["massJPsiK"+tag]     = new TH1F("massJPsiK"+tag,";M(llK) [GeV];Events / 15 MeV" ,100,4.5,6);
     
    /*
      FIXME
      allPlots["massD0"+tag]     = new TH1F("massD0"+tag,";M_{D^{0}};Events / 0.01 GeV" ,30,1.7,2.0);
      allPlots["massD0_lep"+tag]     = new TH1F("massD0_lep"+tag,";M_{K#pi};Events / 0.01" ,30,1.7,2.0);
      allPlots["massD0_mu"+tag]     = new TH1F("massD0_mu"+tag,";M_{K#pi};Events / 0.01 GeV" ,20,1.7,2.0);
      allPlots["massD0_e"+tag]     = new TH1F("massD0_ele"+tag,";M_{K#pi};Events / 0.01 GeV" ,30,1.7,2.0);
      allPlots["massDsmD0loose"+tag]     = new TH1F("massDsmD0loose"+tag,";M_{K#pi};Events / 0.05 GeV" ,6,1.7,2.);
      allPlots["massDsmD0"+tag]     = new TH1F("massDsmD0"+tag,";M_{K#pi};Events / 0.05 GeV" ,6,1.7,2.);
      allPlots["massDs"+tag]     = new TH1F("massDs"+tag,";M_{D^{*}};Events / 0.1 GeV / 0.01" ,30,1.7,2.0);
      allPlots["pi_pt"+tag] = new TH1F("pi_pt"+tag,";#pi^{#pm} P_{T} [GeV];Events / 10 GeV", 15, 0,150);
      
    */
  }


  for (auto& it : allPlots)   { it.second->Sumw2(); it.second->SetDirectory(0); }
  //for (auto& it : all2dPlots) { it.second->Sumw2(); it.second->SetDirectory(0); }
  //LOOP OVER EVENTS
  for (Int_t iev=0;iev<nentries;iev++)
    {
      t->GetEntry(iev);
      if(iev%5000==0) printf ("\r [%3.0f/100] done",100.*(float)(iev)/(float)(nentries));

      //account for pu weights and effect on normalization
      float puWeight(1.0);
      if(!ev.isData) 
	{
	  puWeight=puWgtGr[0]->Eval(ev.putrue);  
	  allPlots["puwgtctr"]->Fill(0.,1.0);
	  allPlots["puwgtctr"]->Fill(1.,puWeight);
	}

      //select 1 good lepton
      //cout << "entering lepton selection" << endl;
      std::vector<int> leptonsIdx,leptonsId;;
      std::vector<TLorentzVector> leptons;
      for(int il=0; il<ev.nl; il++)
	{
          //cout << "in lepton selection" << endl;
	  bool passTightKin(ev.l_pt[il]>30 && fabs(ev.l_eta[il])<2.5);
	  bool passTightId(ev.l_id[il]==13 ? (ev.l_pid[il]>>1)&0x1  : (ev.l_pid[il]>>2)&0x1);	  
	  float relIso(ev.l_relIso[il]);
	  bool passIso( ev.l_id[il]==13 ? relIso<0.15 : (ev.l_pid[il]>>1)&0x1 ); 
	  
	  if(!passTightKin || !passTightId || !passIso) continue;
	  leptonsIdx.push_back(il);
          leptonsId.push_back( ev.l_id[il] );
          TLorentzVector lp4;
          lp4.SetPtEtaPhiM(ev.l_pt[il],ev.l_eta[il],ev.l_phi[il],ev.l_mass[il]);
          leptons.push_back(lp4);
	}
      if(debug) cout << "lepton selection DONE" << endl;

      //check if triggers have fired
      bool hasEETrigger(((ev.elTrigger>>1)&0x1)!=0 || ((ev.elTrigger>>4)&0x1)!=0);
      bool hasMMTrigger(((ev.muTrigger>>2)&0x3)!=0);
      bool hasEMTrigger(((ev.elTrigger>>2)&0x3)!=0);
      bool hasMuTrigger((ev.muTrigger & 0x3)!=0);
      bool hasEleTrigger((ev.elTrigger & 0x1)!=0);
      if(!ev.isData)
	{	 
	  hasMuTrigger=true;
	  hasEleTrigger=true;
	  hasEETrigger=true;
	  hasMMTrigger=true;
	  hasEMTrigger=true;
	}
      else
	{
	  if(requireMutrigger  && !hasMuTrigger) continue;
	  if(requireEletrigger && !hasEleTrigger) continue;
	  if(requireEETriggers && !hasEETrigger) continue;
	  if(requireMMTriggers && !hasMMTrigger) continue;
	  if(requireEMTriggers && !hasEMTrigger) continue;
	}

      //decide the channel
      if(debug) cout << "decide channel" << endl;
      TString chTag("");
      if(leptonsIdx.size()==1 )
	{
	  if(abs(ev.l_id[ leptonsIdx[0] ])==11 && hasEleTrigger) chTag="e";
	  if(abs(ev.l_id[ leptonsIdx[0] ])==13 && hasMuTrigger)  chTag="m";
          if(debug) cout << "found 1 tight lepton ch=" << chTag << endl;
	}
      else if(leptonsIdx.size()>=2)
	{
	  int chId=abs(ev.l_id[ leptonsIdx[0] ])*abs(ev.l_id[ leptonsIdx[1] ]);
	  if(chId==11*11 && hasEETrigger) chTag="ee";
	  else if(chId==13*13 && hasMMTrigger) chTag="mm";
	  else if(chId==11*13 && hasEMTrigger) chTag="em";
          if(debug) cout << "found at least 2 tight leptons ch=" << chTag << endl;
	}
      if(debug) 
	{
	  if(chTag=="") cout << "NO LEPTONS!!" << endl;
	  else cout << "decide channel DONE" << endl;
	}
      if(chTag=="") continue;
      chTag = "_"+chTag;

      //select jets
      std::vector<Jet> bJetsVec,allJetsVec;
      for (int k=0; k<ev.nj;k++)
	{
	  //check kinematics
	  TLorentzVector jp4;
	  jp4.SetPtEtaPhiM(ev.j_pt[k],ev.j_eta[k],ev.j_phi[k],ev.j_mass[k]);

	  //cross clean with respect to leptons 
          bool overlapsWithLepton(false);
          for(size_t il=0; il<leptons.size(); il++) {
            if(jp4.DeltaR(leptons[il])>0.4) continue;
	    overlapsWithLepton=true;
          }
          if(overlapsWithLepton) continue;
          if(debug) cout << "Overlap with lepton DONE" << endl;

	  //smear jet energy resolution for MC
	  float genJet_pt(0);
	  if(ev.j_g[k]>-1) genJet_pt=ev.g_pt[ ev.j_g[k] ];
	  if(!ev.isData && genJet_pt>0) 
	    {
	      float jerSmear=getJetResolutionScales(jp4.Pt(),jp4.Eta(),genJet_pt)[0];
	      jp4 *= jerSmear;
	    }

	  // re-inforce kinematics cuts
	  if(jp4.Pt()<30) continue;
	  if(fabs(jp4.Eta()) > 2.4) continue;
	  
	  //b-tag
	  if(debug) cout << "Starting b-tagging" << endl;
	  float csv = ev.j_csv[k];	  
	  bool isBTagged(csv>0.800);
	  if(!ev.isData)
	    {
	      float jptForBtag(jp4.Pt()>1000. ? 999. : jp4.Pt()), jetaForBtag(fabs(jp4.Eta()));
	      float expEff(1.0), jetBtagSF(1.0);
	      if(abs(ev.j_hadflav[k])==4) 
		{ 
		  expEff    = expBtagEff["c"]->Eval(jptForBtag); 
		  jetBtagSF = sfbReaders[0]->eval( BTagEntry::FLAV_C, jetaForBtag, jptForBtag);
		  jetBtagSF *= expEff>0 ? expBtagEffPy8["c"]->Eval(jptForBtag)/expBtagEff["c"]->Eval(jptForBtag) : 0.;
		}
	      else if(abs(ev.j_hadflav[k])==5) 
		{ 
		  expEff    = expBtagEff["b"]->Eval(jptForBtag); 
		  jetBtagSF = sfbReaders[0]->eval( BTagEntry::FLAV_B, jetaForBtag, jptForBtag);
		  jetBtagSF *= expEff>0 ? expBtagEffPy8["b"]->Eval(jptForBtag)/expBtagEff["b"]->Eval(jptForBtag) : 0.;
		}
	      else
		{
		  expEff    = expBtagEff["udsg"]->Eval(jptForBtag);
                  jetBtagSF = sflReaders[0]->eval( BTagEntry::FLAV_UDSG, jetaForBtag, jptForBtag);
		  jetBtagSF *= expEff> 0 ? expBtagEffPy8["udsg"]->Eval(jptForBtag)/expBtagEff["udsg"]->Eval(jptForBtag) : 0.;
		}
	      
	      //updated b-tagging decision with the data/MC scale factor
	      myBTagSFUtil.modifyBTagsWithSF(isBTagged,    jetBtagSF,     expEff);
	    }
	  if(debug) cout << "b-tagging DONE" << endl;
	  
	  //save jet
          Jet tmpj(jp4, csv, k);
	  for(int ipf = 0; ipf < ev.npf; ipf++) {
	    if(ev.pf_j[ipf] != k) continue;
	    if(ev.pf_c[ipf]==0) continue;
	    TLorentzVector tkP4(0,0,0,0);
	    tkP4.SetPtEtaPhiM(ev.pf_pt[ipf],ev.pf_eta[ipf],ev.pf_phi[ipf],0.);
	    tmpj.addTrack(tkP4,ev.pf_id[ipf]);
	  }
	  tmpj.sortTracksByPt();
	  
          if(isBTagged) bJetsVec.push_back(tmpj);
          allJetsVec.push_back(tmpj);
	}

      //sort by Pt
      sort(bJetsVec.begin(),    bJetsVec.end(),   Jet::sortJetsByPt);
      sort(allJetsVec.begin(),  allJetsVec.end(), Jet::sortJetsByPt);


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
	      triggerCorrWgt=lepEffH.getTriggerCorrection(leptonsId,leptons);
	      for(size_t il=0; il<leptons.size(); il++)
		{
		  if(il>1) break;
		  EffCorrection_t selSF=lepEffH.getOfflineCorrection(leptonsId[il],leptons[il].Pt(),leptons[il].Eta());
		  lepSelCorrWgt.second = sqrt( pow(lepSelCorrWgt.first*selSF.second,2)+pow(lepSelCorrWgt.second*selSF.first,2));
		  lepSelCorrWgt.first *= selSF.first;
		}
	    }
	  
	  //update nominal event weight
	  wgt=triggerCorrWgt.first*triggerCorrWgt.first*lepSelCorrWgt.first*puWgts[0]*norm;
	  if(ev.ttbar_nw>0) wgt*=ev.ttbar_w[0];
	}


      //request at lest 1 b-jet (+ at least 3 jets in total for l+jets)
      //apply Z+low mass DY/quarkonia veto for same flavour channels      
      if(bJetsVec.size()==0) continue;
      if(leptons.size()>=2)
	{
	  float mll=(leptons[0]+leptons[1]).M();
	  float ptll=(leptons[0]+leptons[1]).Pt();
	  if(mll<20) continue;
	  if(chTag!="em" && fabs(mll-91)<15) continue;
	  allPlots["dilp_pt"+chTag]->Fill(ptll,wgt);
	  allPlots["dilp_m"+chTag]->Fill(mll,wgt);
	  allPlots["charge"+chTag]->Fill(ev.l_charge[leptonsIdx[0]]*ev.l_charge[leptonsIdx[1]],wgt);
	}
      else
	{
	  if(allJetsVec.size()<3) continue;
	  allPlots["charge"+chTag]->Fill(ev.l_charge[leptonsIdx[0]],wgt);
	}

      //base selection histograms
      allPlots["MET"+chTag]->Fill(ev.met_pt[0],wgt);
      allPlots["MET_all"]->Fill(ev.met_pt[0],wgt);
      allPlots["j_pt"+chTag]->Fill(allJetsVec[0].getVec().Pt(),wgt);
      allPlots["j_pt_all"]->Fill(allJetsVec[0].getVec().Pt(),wgt);
      allPlots["bj_pt"+chTag]->Fill(bJetsVec[0].getVec().Pt(),wgt);
      allPlots["bj_pt_all"]->Fill(bJetsVec[0].getVec().Pt(),wgt);
      allPlots["nj"+chTag]->Fill(allJetsVec.size(),wgt);
      allPlots["nj_all"]->Fill(allJetsVec.size(),wgt);
      allPlots["nbj"+chTag]->Fill(bJetsVec.size(),wgt);
      allPlots["nbj_all"]->Fill(bJetsVec.size(),wgt);
      allPlots["nlp"+chTag]->Fill(leptons.size(),wgt);
      allPlots["nlp_all"]->Fill(leptons.size(),wgt);
      for(size_t il=0; il<leptons.size(); il++) {       
	if(il>1) break;
	
	TString leptonPF(il == 0 ? "lp" : "l2p");
        allPlots[leptonPF+"_pt"+chTag]->Fill(leptons[il].Pt(),wgt);        
        allPlots[leptonPF+"_pt_all"]->Fill(leptons[il].Pt(),wgt);

        for(size_t ij=0; ij<allJetsVec.size(); ij++) {
          allPlots["dR"+chTag]->Fill(allJetsVec[ij].getVec().DeltaR(leptons[il]),wgt);
          allPlots["dR_all"]->Fill(allJetsVec[ij].getVec().DeltaR(leptons[il]),wgt);
        }
      }

      //charmed resonance analysis : use only jets with CSV>CSVL, up to two per event
      sort(allJetsVec.begin(),    allJetsVec.end(),     Jet::sortJetsByCSV);
      for(size_t ij=0; ij<allJetsVec.size(); ij++)
	{
	  if(ij>1) continue;
	  if(allJetsVec[ij].getCSV()<0.) continue;
	  
	  std::vector<IdTrack> &tracks=allJetsVec[ij].getTracks();
          allPlots["npf"+chTag]->Fill(tracks.size(),wgt);
          allPlots["npf_all"]->Fill(tracks.size(),wgt);
	  for(size_t itk=0; itk<tracks.size(); itk++)
            {
	      allPlots["pfid"+chTag]->Fill(tracks[itk].second,wgt);
	      allPlots["pfid_all"]->Fill(tracks[itk].second,wgt);
	    }
  
	  //J/Psi
	  if(debug) cout << "starting J/Psi->mumu" << endl;
	  float gMassMu(0.1057),gMassK(0.4937);//,gMassPi(0.1396);
	  std::vector<TLorentzVector> pfmuCands,kaonCands;
	  for(size_t itk=0; itk<tracks.size(); itk++)
	    {
	      int abstkid(tracks[itk].second);
	      if(abstkid == 13)
		{
		  TLorentzVector muP4(0,0,0,0);
		  muP4.SetPtEtaPhiM( tracks[itk].first.Pt(), tracks[itk].first.Eta(),
				     tracks[itk].first.Phi(), gMassMu);
		  pfmuCands.push_back( muP4 );
		}
	      if(abstkid == 211)
		{
		  TLorentzVector kP4(0,0,0,0);
		  kP4.SetPtEtaPhiM( tracks[itk].first.Pt(), tracks[itk].first.Eta(),
				    tracks[itk].first.Phi(), gMassK);
		  kaonCands.push_back( kP4 );
		}
	    }
	  if(pfmuCands.size()>1)
	    {
	      float mass12( (pfmuCands[0]+pfmuCands[1]).M() );
	      float mass123( kaonCands.size()>0 ? (pfmuCands[0]+pfmuCands[1]+kaonCands[0]).M() : -1);
              allPlots["massJPsi"+chTag]->Fill(mass12,wgt);
	      allPlots["massJPsi_all"]->Fill(mass12,wgt);
	      if(mass12>3 && mass12<3.2 && mass123>0)
		{
		  allPlots["massJPsiK"+chTag]->Fill(mass123,wgt);
		  allPlots["massJPsiK_all"]->Fill(mass123,wgt);
		}
	    }
	  if(debug) cout << "J/Psi DONE" << endl;
	  
	  continue; //FIXME
	  
        // //D0 and D* 
        // if(debug) cout << "Starting D0 and D*" << endl;
        // nstart = firstTrackIndex(jetindex,&tracks);
        // if((tracks.size() - nstart) < 3) continue;
        // for(int i = nstart; i < nstart+2; i++)
        // //for(int i = 0; i < 3; i++)
        //   for(int j = i+1; j < nstart+2; j++) {
        //     int tk1 = get<0>(tracks.at(i));
        //     int tk2 = get<0>(tracks.at(j));
        //     /*
        //     int tk1 = i;
        //     int tk2 = j;
        //     */
        //     //if(ev.pf_j[tk1] != jetindex) continue;
        //     //if(ev.pf_j[tk2] != jetindex) continue;

        //     allPlots["npf"+chTag+"_meson"]->Fill(ev.npf,wgt);
        //     allPlots["npf"+chTag+"_meson_no_weight"]->Fill(ev.npf,1);
        //     for(unsigned int it = 0; it < jetPt.size(); it++)
        //       bpt = get<2>(jetPt.at(it)) > bpt ? get<2>(jetPt.at(it)) : bpt;
        //     allPlots["bj_pt"+chTag+"_meson"]->Fill(bpt,wgt);
        //     allPlots["bj_pt"+chTag+"_meson_no_weight"]->Fill(bpt,1);
        //     allPlots["pfid"+chTag+"_meson"]->Fill(ev.pf_id[tk1],wgt);
        //     allPlots["nbj"+chTag+"_meson"]->Fill(1,wgt);
        //     allPlots["nbj"+chTag+"_meson_no_weight"]->Fill(1,1);

        //     //opposite sign
        //     if(ev.pf_id[tk1]*ev.pf_id[tk2] != -211*211) continue;

        //     const float gMassK  = 0.4937;
        //     const float gMassPi = 0.1396;
          
        //     p_track1.SetPtEtaPhiM(ev.pf_pt[tk1], ev.pf_eta[tk1], ev.pf_phi[tk1], gMassPi);
        //     p_track2.SetPtEtaPhiM(ev.pf_pt[tk2], ev.pf_eta[tk2], ev.pf_phi[tk2], gMassK);
        //     float mass12 = (p_track1+p_track2).M();
        //     allPlots["dR"+chTag+"_meson"]->Fill(p_track1.DeltaR(p_track2), wgt);
        //     allPlots["dR"+chTag+"_meson_no_weight"]->Fill(p_track1.DeltaR(p_track2), 1);

        //     //if (mass12>1.65 && mass12<2.0)
        //     if (mass12>1.7 && mass12<2.0) {
        //       allPlots["massD0"+chTag]->Fill(mass12,wgt);
        //       allPlots["massD0"+chTag+"_no_weight"]->Fill(mass12,1);
        //     }

        //     //looking for lepton
        //     if(debug) cout << "third lepton" << endl;
        //     //for(int tk3 = 0; tk3 < ev.npf; tk3++)
        //     for(int k = 0; k < (int)tracks.size(); k++) {
        //       int tk3 = get<0>(tracks.at(k));
        //       //if(ev.pf_j[tk3] != jetindex) continue;
        //       if(tk3 == tk1) continue;
        //       if(tk3 == tk2) continue;
        //       if(debug) cout << "third lepton possible" << endl;
            
        //       if(abs(ev.pf_id[tk3]) != 13 && abs(ev.pf_id[tk3]) != 11) continue;
        //       if(debug) cout << "third lepton found" << endl;

        //       if(ev.pf_id[tk2]/abs(ev.pf_id[tk2]) == -ev.pf_id[tk3]/abs(ev.pf_id[tk3])) {
        //         //Kaon and lepton have same charge
        //         //correct mass assumption
        //         if(debug) cout << "correct mass assumption" << endl;
        //         allPlots["masslep"+chTag]->Fill(mass12,wgt);
        //         allPlots["masslep"+chTag+"_no_weight"]->Fill(mass12,1);

        //         if(abs(ev.pf_id[tk3]) == 13)
        //           allPlots["massmu"+chTag]->Fill(mass12,wgt);
        //         if(abs(ev.pf_id[tk3]) == 11)
        //           allPlots["masse"+chTag]->Fill(mass12,wgt);
        //         if(abs(ev.pf_id[tk3]) == 13)
        //           allPlots["massmu"+chTag+"_no_weight"]->Fill(mass12,1);
        //         if(abs(ev.pf_id[tk3]) == 11)
        //           allPlots["masse"+chTag+"_no_weight"]->Fill(mass12,1);

        //       }
        //     }
        //     //looking for pion
        //     if(debug) cout << "D*->pi+D0" << endl;
        //     //for(int tk3 = 0; tk3 < ev.npf; tk3++)
        //     for(int k = 0; k < (int)tracks.size(); k++) {
        //       int tk3 = get<0>(tracks.at(k));
        //       //if(ev.pf_j[tk3] != jetindex) continue;
        //       if(tk3 == tk1) continue;
        //       if(tk3 == tk2) continue;

        //       if(abs(ev.pf_id[tk3]) != 211) continue;
        //       if(debug) cout << "Pion found" << endl;

        //       TLorentzVector p_track3, p_cand;
        //       p_track3.SetPtEtaPhiM(ev.pf_pt[tk3], ev.pf_eta[tk3], ev.pf_phi[tk3], gMassPi);
        //       allPlots["pi_pt"+chTag]->Fill(p_track3.Pt(),wgt);
        //       allPlots["pi_pt"+chTag+"_no_weight"]->Fill(p_track3.Pt(),1);
        //       if( ev.pf_id[tk2]/abs(ev.pf_id[tk2]) == -ev.pf_id[tk3]/abs(ev.pf_id[tk3]) ) {
        //         // Kaon and pion have opposite charges
        //         // I.e. correct mass assumption
        //         if(debug) cout << "correct mass assumption" << endl;
                
        //         p_cand = p_track1+p_track2+p_track3;
        //         allPlots["massDs"+chTag]->Fill(p_cand.M(), wgt);
        //         allPlots["massDs"+chTag+"_no_weight"]->Fill(p_cand.M(), 1);

        //         if(abs(mass12-1.864) < 0.10) { // mass window cut
        //           TLorentzVector p_jet;
        //           p_jet.SetPtEtaPhiM(ev.j_pt[jetindex], ev.j_eta[jetindex], ev.j_phi[jetindex], 0.);

        //           //float hardpt = std::max(ev.pf_pt[tk3], std::max(ev.pf_pt[tk1], ev.pf_pt[tk2]));
        //           //float softpt = std::min(ev.pf_pt[tk3], std::min(ev.pf_pt[tk1], ev.pf_pt[tk2]));
        //           float deltam = p_cand.M() - mass12;

        //           allPlots["massDsmD0loose"+chTag]->Fill(deltam, wgt);
        //           allPlots["massDsmD0loose"+chTag+"_no_weight"]->Fill(deltam, 1);
        //           if(abs(mass12-1.864) < 0.05) { // tighter mass window cut
        //               //FillCharmTree(413,  jetindex, tk1, gMassPi, tk2, gMassK, tk3, gMassPi);
        //               //FillCharmTree(-413, jetindex, deltam, p_cand, p_jet, hardpt, softpt);
        //               allPlots["massDsmD0"+chTag]->Fill(deltam, wgt);
        //               allPlots["massDsmD0"+chTag+"_no_weight"]->Fill(deltam, 1);
        //           }
        //         }
        //       }
        //     }
        //   }
        // if(debug) cout << "D0 and D* DONE" << endl;
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
  if(debug) cout << "writing histograms" << endl;
  for (auto& it : allPlots)  { 
    if(debug) cout << it.second->GetName() << endl;
    if(debug) cout << it.second->GetEntries() << endl;
    it.second->SetDirectory(fOut); it.second->Write(); 
    fOut->cd();
  }
  if(debug) cout << "writing histograms DONE" << endl;
  fOut->Close();
}
