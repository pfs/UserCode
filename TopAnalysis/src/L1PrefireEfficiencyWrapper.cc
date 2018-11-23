#include "TopLJets2015/TopAnalysis/interface/L1PrefireEfficiencyWrapper.h"

#include "TFile.h"
#include "TSystem.h"

#include <iostream>


using namespace std;

//
L1PrefireEfficiencyWrapper::L1PrefireEfficiencyWrapper(bool isData,TString era)
{
  if(isData) return;
  init(era);
}

//
void L1PrefireEfficiencyWrapper::init(TString era)
{
  
  era_=2016;
  cout << "[L1PrefireEfficiencyWrapper]" << endl
       << "\tStarting efficiency scale factors for 2016!" << endl
       << "\tThis needs to be fixed once 2017 maps are available" << endl;

  era=era.ReplaceAll("era2017","era2016");
  TString url(era+"/Map_Jet_L1FinOReff_bxm1_looseJet_SingleMuon_Run2016B-H.root");
  gSystem->ExpandPathName(url);
  TFile *fIn=TFile::Open(url);
  effMapsH_["jet_singlemuon"]=(TEfficiency *)fIn->Get("prefireEfficiencyMap")->Clone("jet_singlemuon");
  fIn->Close();
      
  url=era+"/Map_Jet_L1FinOReff_bxm1_looseJet_JetHT_Run2016B-H.root";
  gSystem->ExpandPathName(url);
  fIn=TFile::Open(url);
  effMapsH_["jet_jetht"]=(TEfficiency *)fIn->Get("prefireEfficiencyMap")->Clone("jet_jetht");
  fIn->Close();        
}


//
EffCorrection_t L1PrefireEfficiencyWrapper::getJetBasedCorrection(std::vector<Jet> jets)
{
  EffCorrection_t corr(1.0,0.0);
  if(effMapsH_.size()==0) return corr;

  //iterate up to two highest-pT jets
  for(size_t i=0; i<min(jets.size(),size_t(2)); i++){
    TEfficiency *eff=0;
    if (era_==2016) {
      eff=jets[i].Pt()<300 ? effMapsH_["jet_singlemuon"] : effMapsH_["jet_jetht"];
    }
    if(!eff) continue;
    Int_t ibin=eff->FindFixBin(fabs(jets[i].Eta()),jets[i].Pt());
    Float_t effVal(eff->GetEfficiency(ibin));
    corr.first *= (1-effVal);
    corr.second += pow(0.5*(eff->GetEfficiencyErrorLow(ibin)+eff->GetEfficiencyErrorUp(ibin)),2);
    corr.second += pow(0.2*effVal,2);    
  }
  corr.second=sqrt(corr.second);

  
  return corr;
}

//
L1PrefireEfficiencyWrapper::~L1PrefireEfficiencyWrapper()
{
}
