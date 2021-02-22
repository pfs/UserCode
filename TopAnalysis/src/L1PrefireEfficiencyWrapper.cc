#include "TopLJets2015/TopAnalysis/interface/L1PrefireEfficiencyWrapper.h"

#include "TFile.h"
#include "TKey.h"
#include "TSystem.h"

#include <iostream>


using namespace std;

//
L1PrefireEfficiencyWrapper::L1PrefireEfficiencyWrapper(bool isData,TString era)
{
  if(isData) 
    return;
  init(era);
}

//
void L1PrefireEfficiencyWrapper::init(TString era)
{

  cout << "[L1PrefireEfficiencyWrapper] with efficiencies for " << era << endl;

  TString url(era.Contains("2016") ?
              era + "/L1prefiring_jetpt_2016BtoH.root" :
              era + "/L1prefiring_jetpt_2017BtoF.root");
  gSystem -> ExpandPathName(url);
  TFile * fIn = TFile::Open(url);
  TString key(era.Contains("2016") ?
              "L1prefiring_jetpt_2016BtoH":
              "L1prefiring_jetpt_2017BtoF");
  effMapsH_["jet"] = (TH2 *) fIn -> Get(key);
  effMapsH_["jet"] -> SetDirectory(nullptr);
  fIn -> Close();
      
  url = (era.Contains("2016") ?
         era + "/L1prefiring_photonpt_2016BtoH.root":
         era + "/L1prefiring_photonpt_2017BtoF.root");
  gSystem -> ExpandPathName(url);
  fIn = TFile::Open(url);
  key = (era.Contains("2016") ?
	 "L1prefiring_photonpt_2016BtoH":
	 "L1prefiring_photonpt_2017BtoF");
  effMapsH_["photon"] = (TH2 *) fIn -> Get(key);
  effMapsH_["photon"] -> SetDirectory(nullptr);
  fIn -> Close();        
}


//
EffCorrection_t L1PrefireEfficiencyWrapper::getCorrection(std::vector<Jet> & jets, std::vector<Particle> & photons, bool byMax)
{
  EffCorrection_t corr(1.0, 0.0);

  //iterate up to two highest-pT jets
  Double_t ptmax(effMapsH_["jet"] -> GetYaxis() -> GetXmax());
  for(size_t ind = 0; ind < min(jets.size(), size_t(2)); ind ++)
    {
      if (fabs(jets[ind].Eta()) > 3.0 ) 
	continue;
      if(jets[ind].Pt() < 30.0) 
	continue;    
      const Int_t ibin = effMapsH_["jet"] -> FindBin(fabs(jets[ind].Eta()), min(jets[ind].Pt(), ptmax));
      const Float_t effVal(effMapsH_["jet"] -> GetBinContent(ibin));
      const Float_t effUnc(effMapsH_["jet"] -> GetBinError(ibin));
      if (!byMax) 
	{
	  corr.first *= (1.0 - effVal);
	  corr.second += pow(effUnc, 2) + pow(0.2*effVal, 2);    
	}
      else if(1.0 - effVal < corr.first) 
	{
	  corr.first  = (1.0 - effVal);
	  corr.second = pow(effUnc,2) + pow(0.2*effVal,2);
	}
    }

  //photons
  ptmax = (effMapsH_["photon"] -> GetYaxis() -> GetXmax());
  for(size_t ind = 0; ind < photons.size(); ind ++)
    {
      if (fabs(photons[ind].Eta()) > 3.0 ) 
	continue;
      if (photons[ind].Pt() < 20.0) 
	continue;
      const Int_t ibin = effMapsH_["photon"] -> FindBin(fabs(photons[ind].Eta()), min(photons[ind].Pt(), ptmax));
      const Float_t effVal(effMapsH_["photon"] -> GetBinContent(ibin));
      const Float_t effUnc(effMapsH_["photon"] -> GetBinError(ibin));
      if(!byMax) 
	{
	  corr.first *= (1.0 - effVal);
	  corr.second += pow(effUnc, 2) + pow(0.2 * effVal, 2);    
	}
      else if(1.0 - effVal < corr.first) 
	{
	  corr.first  = (1-effVal);
	  corr.second = pow(effUnc, 2) + pow(0.2 * effVal, 2);
	}
    }

  //finalize error
  corr.second = sqrt(corr.second);  
  return corr;
}

//
EffCorrection_t L1PrefireEfficiencyWrapper::getCorrection(std::vector<Jet> & jets, bool byMax)
{
  EffCorrection_t corr(1.0, 0.0);

  //iterate up to two highest-pT jets
  Double_t ptmax(effMapsH_["jet"] -> GetYaxis() -> GetXmax());
  for(size_t ind = 0; ind < min(jets.size(), size_t(2)); ind ++)
    {
      if(fabs(jets[ind].Eta()) > 3.0 ) 
	continue;
      if(jets[ind].Pt() < 30.0) 
	continue;    
      const Int_t ibin = effMapsH_["jet"] -> FindBin(fabs(jets[ind].Eta()), min(jets[ind].Pt(), ptmax));
      const Float_t effVal(effMapsH_["jet"] -> GetBinContent(ibin));
      const Float_t effUnc(effMapsH_["jet"] -> GetBinError(ibin));
      if (!byMax) 
	{
	  corr.first *= (1.0 - effVal);
	  corr.second += pow(effUnc, 2) + pow(0.2 * effVal, 2);    
	}
      else if(1.0 - effVal < corr.first) 
	{
	  corr.first  = (1.0 - effVal);
	  corr.second = pow(effUnc, 2) + pow(0.2 * effVal, 2);
	}
    }
  
  //finalize error
  corr.second = sqrt(corr.second);  
  return corr;
}

//
L1PrefireEfficiencyWrapper::~L1PrefireEfficiencyWrapper()
{
}
