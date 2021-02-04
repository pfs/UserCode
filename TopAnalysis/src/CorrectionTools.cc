#include "TopLJets2015/TopAnalysis/interface/CorrectionTools.h"
#include "TH2F.h"
#include "TSystem.h"
//
std::map<Int_t,Float_t> lumiPerRun(TString era)
{
  return parseLumiInfo(era).first; 
}

//
std::pair<std::map<Int_t,Float_t>, TH1F *> parseLumiInfo(TString era)
{
  std::map<Int_t,Float_t> lumiMap;
  TH1F *countH=0;
  std::pair<std::map<Int_t,Float_t>, TH1F *> toReturn(lumiMap,countH);

  //read out the values from the histogram stored in lumisec.root
  TFile *inF=TFile::Open(Form("%s/lumisec.root",era.Data()),"READ");
  if(inF==0) return toReturn;
  if(inF->IsZombie()) return toReturn;
  TH2F *h=(TH2F *)inF->Get("lumisec_inc");
  int nruns(h->GetNbinsX());
  countH=new TH1F("ratevsrun","ratevsrun;Run;Events/pb",nruns,0,nruns);
  countH->SetDirectory(0);
  for(int xbin=1; xbin<=nruns; xbin++)
    {
      TString run=h->GetXaxis()->GetBinLabel(xbin);
      lumiMap[run.Atoi()]=h->GetBinContent(xbin);
      countH->GetXaxis()->SetBinLabel(xbin,run);
    }
  inF->Close();

  toReturn.first=lumiMap;
  toReturn.second=countH;
  return toReturn;
};

//
std::vector<RunPeriod_t> getRunPeriods(TString era)
{
  //init the conditons of the run
  std::vector<RunPeriod_t> periods;
  if(era.Contains("era2016"))
    {
      periods.push_back( RunPeriod_t("BCDEF", 19323.4) ); 
      periods.push_back( RunPeriod_t("GH",    16551.4) );
    }
  return periods;
}

//
TString assignRunPeriod(const std::vector<RunPeriod_t> & runPeriods, TRandom * rand)
{
  float totalLumi(0.0);
  for (auto periodLumi : runPeriods) totalLumi += periodLumi.second;

  //generate randomly in the total lumi range to pick one of the periods
  float pickLumi( rand != nullptr ? rand -> Uniform(totalLumi) : gRandom -> Uniform(totalLumi) );
  float testLumi(0.0); 
  int iLumi(0);
  for (auto periodLumi : runPeriods) 
    {
      testLumi += periodLumi.second;
      if (pickLumi < testLumi * ( 1 + 1E-9)) 
	{
	  //printf("pickLumi %.9f testLumi %.9f\n", pickLumi, testLumi);
	  break;
	}
      else ++ iLumi;
    }
  //return the period
  return runPeriods[iLumi].first;
}


//
std::vector<TGraph *> getPileupWeights(TString era, TH1 * genPU, TString period)
{
  std::vector<TGraph *> puWgtGr;
  if(genPU == nullptr) return puWgtGr;

  if(genPU -> GetNbinsX() == 1000) genPU -> Rebin(10);
  genPU -> Scale(1.0/genPU -> Integral());

  //readout the pileup weights and take the ratio of data/MC
  TString puWgtUrl(era + "/pileupWgts" + period + ".root");
  gSystem -> ExpandPathName(puWgtUrl);
  TFile * fIn = TFile::Open(puWgtUrl);
  for(size_t ind = 0; ind < 3; ind ++)
    {
      TString grName("pu_nom");
      if (ind == 1) grName = "pu_down";
      if (ind == 2) grName = "pu_up";
      TGraph * puData = (TGraph *) fIn -> Get(grName);
      Float_t totalData = puData -> Integral();
      TH1 * tmp = (TH1 *) genPU -> Clone("tmp");
      for(Int_t xbin = 1; xbin <= tmp -> GetXaxis() -> GetNbins(); xbin ++)
        {
          Float_t yexp = genPU -> GetBinContent(xbin);
          Double_t xobs, yobs;
          puData -> GetPoint(xbin - 1, xobs, yobs);
          tmp -> SetBinContent(xbin, yexp > 0 ? yobs/(totalData*yexp) : 0.0 );
        }
      TGraph * gr = new TGraph(tmp);
      grName.ReplaceAll("pu", "puwgts");
      gr->SetName(period + grName);
      puWgtGr.push_back( gr );
      tmp -> Delete();
    }
  return puWgtGr;
}


//
std::map<TString, std::vector<TGraph *> > getPileupWeightsMap(TString era, TH1 * genPU)
{
  std::map<TString, std::vector<TGraph *> > puWgtGr;
  if (genPU == nullptr) return puWgtGr;

  std::vector<RunPeriod_t> periods = getRunPeriods(era);

  for (auto period : periods) 
    {
      puWgtGr[period.first] = getPileupWeights(era, genPU, period.first);
    }
  return puWgtGr;
}


//apply jet energy resolutions (scaling method)
void smearJetEnergies(MiniEvent_t & ev, std::string option) 
{
  if (ev.isData) return;
  
  for (int kind = 0; kind < ev.nj; kind ++) 
    {
      TLorentzVector jp4;
      jp4.SetPtEtaPhiM(ev.j_pt[kind], ev.j_eta[kind], ev.j_phi[kind], ev.j_mass[kind]);

      //smear jet energy resolution for MC
      float genJet_pt(0.0);
      if (ev.j_g[kind] > -1) genJet_pt = ev.g_pt[ ev.j_g[kind] ];
      if (genJet_pt > 0.0) 
	{
	  smearJetEnergy_old(jp4, genJet_pt, option);
	  ev.j_pt[kind]   = jp4.Pt();
	  ev.j_eta[kind]  = jp4.Eta();
	  ev.j_phi[kind]  = jp4.Phi();
	  ev.j_mass[kind] = jp4.M();
	}
    }
}

//apply jet energy resolutions (hybrid method)
void smearJetEnergies(MiniEvent_t & ev, JME::JetResolution * jer, std::string option) 
{
  if(ev.isData) return;

  TRandom * random = new TRandom3(0); // random seed
  
  for (int kind = 0; kind < ev.nj; kind++) 
    {
      TLorentzVector jp4;
      jp4.SetPtEtaPhiM(ev.j_pt[kind], ev.j_eta[kind], ev.j_phi[kind], ev.j_mass[kind]);

      //smear jet energy resolution for MC
      float genJet_pt(0.0);
      if(ev.j_g[kind] > -1) 
	genJet_pt = ev.g_pt[ ev.j_g[kind] ];
      //scaling method for matched jets
      if(genJet_pt > 0.0) 
	{
	  smearJetEnergy_old(jp4, genJet_pt, option);
	  ev.j_pt[kind]   = jp4.Pt();
	  ev.j_eta[kind]  = jp4.Eta();
	  ev.j_phi[kind]  = jp4.Phi();
	  ev.j_mass[kind] = jp4.M();
	}
      //stochastic smearing for unmatched jets
      else 
	{
	  double jet_resolution = jer -> getResolution({{JME::Binning::JetPt, ev.j_pt[kind]}, {JME::Binning::JetEta, ev.j_eta[kind]}, {JME::Binning::Rho, ev.rho}});
	  smearJetEnergyStochastic(jp4, random, jet_resolution, option);
	  ev.j_pt[kind]   = jp4.Pt();
	  ev.j_eta[kind]  = jp4.Eta();
	  ev.j_phi[kind]  = jp4.Phi();
	  ev.j_mass[kind] = jp4.M();
	}
    }
  
  delete random;
}

//
void smearJetEnergy_old(TLorentzVector & jp4, float genJet_pt, std::string option)
{
  int smearIdx(0);
  if (option == "up") smearIdx   = 1;
  if (option == "down") smearIdx = 2;
  float jerSmear = getJetResolutionScales_old(jp4.Pt(), jp4.Eta(), genJet_pt)[smearIdx];
  jp4 *= jerSmear;
}

//
void smearJetEnergyStochastic(TLorentzVector & jp4, TRandom * random, double resolution, std::string option)
{
  int smearIdx(0);
  if (option == "up")   smearIdx = 1;
  if (option == "down") smearIdx = 2;
  const float jerSmear = getJetResolutionScales_old(jp4.Pt(), jp4.Eta(), 0.0)[smearIdx];
  const float jerFactor = 1.0 + random -> Gaus(0.0, resolution) * sqrt(std::max(pow(jerSmear, 2) - 1.0, 0.0));
  jp4 *= jerFactor;
}

//see working points in https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation80XReReco
void addBTagDecisions(MiniEvent_t & ev, float wp, float wpl) 
{
  for (int kind = 0; kind < ev.nj; kind++) 
    {
      if (ev.j_hadflav[kind] >= 4) ev.j_btag[kind] = (ev.j_csv[kind] > wp);
      else                      ev.j_btag[kind] = (ev.j_csv[kind] > wpl);
    }
}


//details in https://twiki.cern.ch/twiki/bin/view/CMS/BTagCalibration
void updateBTagDecisions_old(MiniEvent_t & ev, 
			 std::map<BTagEntry::JetFlavor, BTagCalibrationReader *> & btvsfReaders,
			 std::map<BTagEntry::JetFlavor, TGraphAsymmErrors*> & expBtagEff, 
			 std::map<BTagEntry::JetFlavor, TGraphAsymmErrors*> & expBtagEffPy8, 
			 BTagSFUtil_old * myBTagSFUtil, 
			 std::string optionbc, 
			 std::string optionlight) 
{
  for (int kind = 0; kind < ev.nj; kind ++) 
    {
      TLorentzVector jp4;
      jp4.SetPtEtaPhiM(ev.j_pt[kind], ev.j_eta[kind], ev.j_phi[kind], ev.j_mass[kind]);

      bool isBTagged(ev.j_btag[kind]);
      if(!ev.isData) 
	{
	  const float jptForBtag(jp4.Pt() > 1000.0 ? 999.0 : jp4.Pt());
	  const float jetaForBtag(fabs(jp4.Eta()));
	  float expEff(1.0);
	  float jetBtagSF(1.0);
      
	  BTagEntry::JetFlavor hadFlav = BTagEntry::FLAV_UDSG;
	  std::string option = optionlight;
	  if (abs(ev.j_hadflav[kind]) == 4) { hadFlav=BTagEntry::FLAV_C; option = optionbc; }
	  if (abs(ev.j_hadflav[kind]) == 5) { hadFlav=BTagEntry::FLAV_B; option = optionbc; }

	  expEff     = expBtagEff[hadFlav] -> Eval(jptForBtag); 
	  jetBtagSF  = btvsfReaders[hadFlav] -> eval_auto_bounds( option, hadFlav, jetaForBtag, jptForBtag);
	  jetBtagSF *= expEff > 0 ? expBtagEffPy8[hadFlav] -> Eval(jptForBtag)/expBtagEff[hadFlav] -> Eval(jptForBtag) : 0.0;
      
	  //updated b-tagging decision with the data/MC scale factor
	  myBTagSFUtil -> modifyBTagsWithSF(isBTagged, jetBtagSF, expEff);
	  ev.j_btag[kind] = isBTagged;
	}
    }
}

//details in https://twiki.cern.ch/twiki/bin/view/CMS/BTagCalibration
std::map<BTagEntry::JetFlavor, BTagCalibrationReader *> getBTVcalibrationReaders(TString era,
										BTagEntry::OperatingPoint btagOP, 
										TString period)
{
  //start the btag calibration
  TString btagUncUrl(era + "/btagSFactors.csv");
  if (era.Contains("era2016")) btagUncUrl = era + "/CSVv2_Moriond17_" + period + ".csv";
  gSystem -> ExpandPathName(btagUncUrl);
  BTagCalibration btvcalib("csvv2", btagUncUrl.Data());

  //start calibration readers for b,c and udsg separately including the up/down variations
  std::map<BTagEntry::JetFlavor, BTagCalibrationReader *> btvCalibReaders;
  btvCalibReaders[BTagEntry::FLAV_B]    =  new BTagCalibrationReader(btagOP, "central", {"up", "down"});
  btvCalibReaders[BTagEntry::FLAV_B]    -> load(btvcalib, BTagEntry::FLAV_B, "mujets");
  btvCalibReaders[BTagEntry::FLAV_C]    =  new BTagCalibrationReader(btagOP, "central", {"up", "down"});
  btvCalibReaders[BTagEntry::FLAV_C]    -> load(btvcalib, BTagEntry::FLAV_C, "mujets");
  btvCalibReaders[BTagEntry::FLAV_UDSG] =  new BTagCalibrationReader(btagOP, "central", {"up", "down"});
  btvCalibReaders[BTagEntry::FLAV_UDSG] -> load(btvcalib, BTagEntry::FLAV_UDSG, "incl");

  //all done
  return btvCalibReaders;
}

//
std::map<TString, std::map<BTagEntry::JetFlavor,BTagCalibrationReader *> > getBTVcalibrationReadersMap(TString era,
												       BTagEntry::OperatingPoint btagOP)
{
  std::map<TString, std::map<BTagEntry::JetFlavor,BTagCalibrationReader *> > btvCalibReadersMap;
  std::vector<RunPeriod_t> periods=getRunPeriods(era);
  for (auto period : periods) {
    btvCalibReadersMap[period.first] = getBTVcalibrationReaders(era, btagOP, period.first);
  }

  return btvCalibReadersMap;
}


//the expections are created with the script scripts/saveExpectedBtagEff.py (cf README)
std::map<BTagEntry::JetFlavor, TGraphAsymmErrors *> readExpectedBtagEff(TString era,TString btagExpPostFix)
{
  //open up the ROOT file with the expected efficiencies
  TString btagEffExpUrl(era+"/expTageff.root");
  btagEffExpUrl.ReplaceAll(".root",btagExpPostFix+".root");
  gSystem->ExpandPathName(btagEffExpUrl);
  TFile *beffIn=TFile::Open(btagEffExpUrl);
  //read efficiency graphs
  std::map<BTagEntry::JetFlavor, TGraphAsymmErrors *> expBtagEff;
  expBtagEff[BTagEntry::FLAV_B]=(TGraphAsymmErrors *)beffIn->Get("b");
  expBtagEff[BTagEntry::FLAV_C]=(TGraphAsymmErrors *)beffIn->Get("c");
  expBtagEff[BTagEntry::FLAV_UDSG]=(TGraphAsymmErrors *)beffIn->Get("udsg");
  beffIn->Close();

  //all done
  return expBtagEff;
}


// See https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECUncertaintySources#Main_uncertainties_2016_80X
void applyJetCorrectionUncertainty(MiniEvent_t &ev, JetCorrectionUncertainty *jecUnc, TString jecVar, TString direction) {
  for (int k = 0; k < ev.nj; k++) {
    if ((jecVar == "FlavorPureGluon"  and not (ev.j_flav[k] == 21 or ev.j_flav[k] == 0)) or
        (jecVar == "FlavorPureQuark"  and not (abs(ev.j_flav[k]) <= 3 and abs(ev.j_flav[k]) != 0)) or
        (jecVar == "FlavorPureCharm"  and not (abs(ev.j_flav[k]) == 4)) or
        (jecVar == "FlavorPureBottom" and not (abs(ev.j_flav[k]) == 5)))
      continue;
    
    TLorentzVector jp4;
    jp4.SetPtEtaPhiM(ev.j_pt[k],ev.j_eta[k],ev.j_phi[k],ev.j_mass[k]);
    applyJetCorrectionUncertainty(jp4,jecUnc,direction);
    ev.j_pt[k]   = jp4.Pt(); 
    ev.j_eta[k]  = jp4.Eta();
    ev.j_phi[k]  = jp4.Phi();
    ev.j_mass[k] = jp4.M();
  }
}

//
void applyJetCorrectionUncertainty(TLorentzVector &jp4,JetCorrectionUncertainty *jecUnc,TString direction)
{
    jecUnc->setJetPt(jp4.Pt());
    jecUnc->setJetEta(jp4.Eta());
    double scale = 1.;
    if (direction == "up")
      scale += jecUnc->getUncertainty(true);
    else if (direction == "down")
      scale -= jecUnc->getUncertainty(false);
    
    jp4 *= scale;
}

//b fragmentation
std::map<TString, TGraph*> getBFragmentationWeights_old(TString era) 
{
  std::map<TString, TGraph*> bfragMap;

  TString bfragWgtUrl(era + "/bfragweights.root");
  gSystem -> ExpandPathName(bfragWgtUrl);
  TFile * fIn = TFile::Open(bfragWgtUrl);
  bfragMap["upFrag"]       = (TGraph *) fIn -> Get("upFrag");
  bfragMap["centralFrag"]  = (TGraph *) fIn -> Get("centralFrag");
  bfragMap["downFrag"]     = (TGraph *) fIn -> Get("downFrag");
  bfragMap["PetersonFrag"] = (TGraph *) fIn -> Get("PetersonFrag");
  return bfragMap;
}

double computeBFragmentationWeight_old(MiniEvent_t &ev, TGraph* wgtGr) 
{
  double weight = 1.0;
  for (int ind = 0; ind < ev.ng; ind ++) 
    {
      if (abs(ev.g_id[ind]) == 5) weight *= wgtGr -> Eval(ev.g_xb[ind]);
    }
  return weight;
}

std::map<TString, std::map<int, double> > getSemilepBRWeights_old(TString era) 
{
  std::map<TString, TGraph*>                bfragMap;
  std::map<TString, std::map<int, double> > brMap;

  TString bfragWgtUrl(era + "/bfragweights.root");
  gSystem -> ExpandPathName(bfragWgtUrl);
  TFile * fIn = TFile::Open(bfragWgtUrl);
  bfragMap["semilepbrUp"] = (TGraph *)   fIn -> Get("semilepbrUp");
  bfragMap["semilepbrDown"] = (TGraph *) fIn -> Get("semilepbrDown");
  
  for (auto const& gr : bfragMap) 
    {
      for (int ind = 0; ind < gr.second -> GetN(); ++ ind) 
	{
	  double x,y;
	  gr.second->GetPoint(ind, x, y);
	  brMap[gr.first][x] = y;
	}
    }
  
  return brMap;
}

double computeSemilepBRWeight_old(MiniEvent_t & ev, std::map<int, double> corr, int pid, bool useabs) 
{
  double weight = 1.0;
  for (int ind = 0; ind < ev.ng; ind ++) 
    {
      if (!ev.g_isSemiLepBhad[ind]) 
	continue;
      if (corr.count(ev.g_bid[ind]) == 0) 
	continue;
      if (!useabs and (pid == 0 or pid == ev.g_bid[ind])) 
	weight *= corr[ev.g_bid[ind]];
      else if (useabs and (pid == 0 or pid == abs(ev.g_bid[ind]))) 
	{
	  weight *= (corr[ev.g_bid[ind]] + corr[-ev.g_bid[ind]])/2.0;
	}
    }
  return weight;
}

void applyTrackingEfficiencySF(MiniEvent_t &ev, double sf) {
  if(ev.isData) return;
  
  TRandom* random = new TRandom3(0); // random seed

  for (int k = 0; k < ev.npf; k++) {
    if (random->Rndm() > sf) {
      //make sure that particle does not pass any cuts
      ev.pf_pt[k]  = 1e-20;
      ev.pf_m[k]   = 1e-20;
      ev.pf_eta[k] = 999.;
      ev.pf_c[k]   = 0;
    }
  }
  
  delete random;
}

