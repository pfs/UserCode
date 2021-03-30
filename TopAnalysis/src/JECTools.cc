#include "TopLJets2015/TopAnalysis/interface/JECTools.h"
#include "TSystem.h"
#include "TFile.h"
#include "TString.h"
#include "TH2F.h"

JECTools::JECTools(TString era)
{
}


//
JECTools::JECTools(TString era, std::string jecVar) : 
  era_(era),
  jetCorr_(nullptr),
  rand_(new TRandom3(0))
{
  TString url(era_ + "/Summer16_25nsV1_MC_EtaResolution_AK4PFchs.txt");
  if (era_.Contains("2017")) url = era_ + "/Fall17_V3_MC_EtaResolution_AK4PFchs.txt";
  gSystem -> ExpandPathName(url);
  jer_ = new JME::JetResolution(url.Data());

  url = era_ + "/Summer16_25nsV1_MC_SF_AK4PFchs.txt";
  if (era_.Contains("2017")) url = era_ + "/Fall17_V3_MC_SF_AK4PFchs.txt";
  gSystem -> ExpandPathName(url);
  jerSF_ = new JME::JetResolutionScaleFactor(url.Data());

  const TString jecUncUrl(era + "/Summer16_23Sep2016V4_MC_UncertaintySources_AK4PFchs.txt");
  gSystem -> ExpandPathName(jecUncUrl);

  JetCorrectorParameters * jecParam                   = new JetCorrectorParameters(jecUncUrl.Data(), jecVar);
  JetCorrectionUncertainty * jecUnc                   = new JetCorrectionUncertainty( *jecParam );
  jecUncs_[jecVar] = jecUnc;
}

//
void JECTools::startFactorizedJetEnergyCorrector(bool isMC) {

  std::cout << "[startFactorizedJetEnergyCorrector][WARNING] don't forget to update with the jet correction files for " << era_ << std::endl;

  //order matters: L1 -> L2 -> L3 (-> Residuals)
  gSystem->ExpandPathName(era_);
  std::vector<std::string> jetCorFiles;
  TString pf("Summer15_25nsV7_");
  pf += (isMC ? "MC" : "DATA");
  std::cout << "Jet energy corrections from " << era_+"/"+pf+"_*_AK4PFchs.txt" << std::endl;
  jetCorFiles.push_back((era_+"/"+pf+"_L1FastJet_AK4PFchs.txt").Data());
  jetCorFiles.push_back((era_+"/"+pf+"_L2Relative_AK4PFchs.txt").Data());
  jetCorFiles.push_back((era_+"/"+pf+"_L3Absolute_AK4PFchs.txt").Data());
  if(!isMC) jetCorFiles.push_back((era_+"/"+pf+"_L2L3Residual_AK4PFchs.txt").Data());
  
  //init the parameters for correction
  std::vector<JetCorrectorParameters> corSteps;
  for(size_t i=0; i<jetCorFiles.size(); i++) corSteps.push_back(JetCorrectorParameters(jetCorFiles[i]));
  
  //return the corrector
  jetCorr_=new FactorizedJetCorrector(corSteps);
}

//apply jet energy resolutions (scaling method)
void JECTools::smearJetEnergies(MiniEvent_t & ev, Variation option) 
{
  if(ev.isData) return;
  
  for (int kind = 0; kind < ev.nj; kind ++) 
    {
      TLorentzVector jp4;
      jp4.SetPtEtaPhiM(ev.j_pt[kind], ev.j_eta[kind], ev.j_phi[kind], ev.j_mass[kind]);

      //smear jet energy resolution for MC
      float genJet_pt(0.0);
      if (ev.j_g[kind] > -1) genJet_pt = ev.g_pt[ ev.j_g[kind] ];
      if (genJet_pt > 0.0) 
	{
	  jp4 = getSmearedJet(jp4, genJet_pt, ev.rho, option);
	  ev.j_pt[kind]   = jp4.Pt();
	  ev.j_eta[kind]  = jp4.Eta();
	  ev.j_phi[kind]  = jp4.Phi();
	  ev.j_mass[kind] = jp4.M();
	}
    }
}

//apply jet energy resolutions (hybrid method)
TLorentzVector JECTools::getSmearedJet(TLorentzVector &jp4, float genJet_pt,float rho,Variation option,float jerVarPartial) {

  TLorentzVector smearedJet(jp4);
  float jerSmear(1.0);
  JME::JetParameters jparam = { {JME::Binning::JetPt, jp4.Pt()}, 
                                {JME::Binning::JetEta, jp4.Eta()}, 
                                {JME::Binning::Rho, rho} };
  float jet_resolution = jer_ -> getResolution(jparam);
  float jer_sf = jerSF_ -> getScaleFactor(jparam, option);
  if (jerVarPartial >= 0 && jerVarPartial <= 1.0) 
    {
      const float jer_nom_sf = jerSF_ -> getScaleFactor(jparam, Variation::NOMINAL);
      float delta_sf = (jer_sf - jer_nom_sf) * jerVarPartial;    
      jer_sf = jer_nom_sf + delta_sf;
    }
  
  //use stochasting smearing for unmatched jets
  if (genJet_pt <= 0.0)
    {
      const float sigma = jet_resolution * std::sqrt(std::max(float(jer_sf * jer_sf - 1.0), float(0.0)));
      jerSmear = std::max(float(1.0 + rand_ -> Gaus(0.0, sigma)), float(0.0));
    }
  else 
    {           
      float dPt = jp4.Pt() - genJet_pt;
      jerSmear = 1.0 + (jer_sf - 1.0) * dPt / jp4.Pt();
    }
  smearedJet *= jerSmear;  
  return smearedJet;
}

//
void JECTools::applyJetCorrectionUncertainty(MiniEvent_t & ev, TString jecVar, Variation direction) 
{
  
  for (int kind = 0; kind < ev.nj; kind ++) 
    {
      if ((jecVar == "FlavorPureGluon"  and not (ev.j_flav[kind] == 21 or ev.j_flav[kind] == 0)) or
	  (jecVar == "FlavorPureQuark"  and not (abs(ev.j_flav[kind]) <= 3 and abs(ev.j_flav[kind]) != 0)) or
	  (jecVar == "FlavorPureCharm"  and not (abs(ev.j_flav[kind]) == 4)) or
	  (jecVar == "FlavorPureBottom" and not (abs(ev.j_flav[kind]) == 5)))
	continue;
    
      TLorentzVector jp4;
      jp4.SetPtEtaPhiM(ev.j_pt[kind], ev.j_eta[kind], ev.j_phi[kind], ev.j_mass[kind]);
      jp4 = getShiftedJet(jp4, jecVar, direction);
      ev.j_pt[kind]   = jp4.Pt(); 
      ev.j_eta[kind]  = jp4.Eta();
      ev.j_phi[kind]  = jp4.Phi();
      ev.j_mass[kind] = jp4.M();
    }
}

//
TLorentzVector JECTools::getShiftedJet(TLorentzVector & jp4, TString jecVar, Variation direction) 
{
  JetCorrectionUncertainty * jecUnc = jecUncs_[jecVar];
  jecUnc -> setJetPt(jp4.Pt());
  jecUnc -> setJetEta(jp4.Eta());
  float scale = 1.0;
  if(direction == Variation::UP) scale += jecUnc -> getUncertainty(true);
  else                           scale -= jecUnc -> getUncertainty(false);
  TLorentzVector shiftedJet(jp4);
  jp4 *= scale;
  return jp4; 
}

//b fragmentation
void JECTools::startBFragmentationWeights() {

  TString bfragWgtUrl(era_+"/bfragweights.root");
  gSystem->ExpandPathName(bfragWgtUrl);
  TFile *fIn=TFile::Open(bfragWgtUrl);
  bfragMap_["upFrag"] = (TGraph *)fIn->Get("upFrag");
  bfragMap_["centralFrag"] = (TGraph *)fIn->Get("centralFrag");
  bfragMap_["downFrag"] = (TGraph *)fIn->Get("downFrag");
  bfragMap_["PetersonFrag"] = (TGraph *)fIn->Get("PetersonFrag");
}

//
float JECTools::computeBFragmentationWeight(MiniEvent_t &ev, TString wgtName) {
  TGraph *wgtGr=bfragMap_[wgtName];
  float weight = 1.;
  for (int i = 0; i < ev.ng; i++) {
    if (abs(ev.g_id[i])==5) weight *= wgtGr->Eval(ev.g_xb[i]);
  }
  return weight;
}

//
void JECTools::startSemilepBRWeights() {
  std::map<TString, TGraph*> bfragMap;

  TString bfragWgtUrl(era_+"/bfragweights.root");
  gSystem->ExpandPathName(bfragWgtUrl);
  TFile *fIn=TFile::Open(bfragWgtUrl);
  bfragMap["semilepbrUp"] = (TGraph *)fIn->Get("semilepbrUp");
  bfragMap["semilepbrDown"] = (TGraph *)fIn->Get("semilepbrDown");
  
  double x,y;
  for (auto const& gr : bfragMap) {
    for (int i = 0; i < gr.second->GetN(); ++i) {      
      gr.second->GetPoint(i,x,y);
      semilepBRwgts_[gr.first][x] = y;
    }
  }
}

//
float JECTools::computeSemilepBRWeight(MiniEvent_t &ev, TString var,int pid, bool useabs) {
  float weight = 1.;
  for (int i = 0; i < ev.ng; i++) {
    if (!ev.g_isSemiLepBhad[i]) continue;
    if (semilepBRwgts_[var].count(ev.g_bid[i]) == 0) continue;
    if (!useabs and (pid == 0 or pid == ev.g_bid[i])) weight *= semilepBRwgts_[var][ev.g_bid[i]];
    else if (useabs and (pid == 0 or pid == abs(ev.g_bid[i]))) {
      weight *= (semilepBRwgts_[var][ev.g_bid[i]]+semilepBRwgts_[var][-ev.g_bid[i]])/2.;
    }
  }
  return weight;
}

JECTools::~JECTools()
{
  /*
  if (jetCorr_)
    {
      delete jetCorr_;
      jetCorr_ = nullptr;
    }
  if (jecParam_)
    {
      delete jecParam_;
      jecParam_ = nullptr;
    }
  if (jer_)
    {
      delete jer_;
      jer_ = nullptr;
    }
  if (jerSF_)
    {
      delete jerSF_;
      jerSF_ = nullptr;
    }
  if (rand_)
    {
      delete rand_;
      rand_ = nullptr;
    }
  for (std::map<TString, JetCorrectionUncertainty *>::iterator it = jecUncs_.begin(); it != jecUncs_.end(); it ++)
    {
      delete it -> second;
    }
  jecUncs_.clear();
  */
}
