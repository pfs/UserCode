#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"
#include "TSystem.h"
#include "TFile.h"
#include "TString.h"
#include "TH2F.h"

//
JetPullInfo_t getPullVector( MiniEvent_t &ev, int ijet)
{
  JetPullInfo_t result;
  result.n=0; result.nch=0;
  result.pull=TVector2(0,0);
  result.chPull=TVector2(0,0);
  
  //re-reconstruct the jet direction with the charged tracks
  TLorentzVector jet(0,0,0,0);
  jet.SetPtEtaPhiM(ev.j_pt[ijet], ev.j_eta[ijet], ev.j_phi[ijet], ev.j_mass[ijet]);
  TLorentzVector chargedJet(0,0,0,0);
  TLorentzVector constituent(0,0,0,0);
  std::vector<std::pair<TLorentzVector,bool> > allConstituents;
  unsigned int nCharged = 0;
  for(Int_t idx = 0; idx<ev.npf; idx++)
    {
      if(ev.pf_j[idx]!=ijet) continue;
      constituent.SetPtEtaPhiM( ev.pf_pt[idx], ev.pf_eta[idx], ev.pf_phi[idx], ev.pf_m[idx]);
      bool isCharged(abs(ev.pf_id[idx])==11 ||
		     abs(ev.pf_id[idx])==13 ||
		     abs(ev.pf_id[idx])==211 );
      allConstituents.push_back(std::make_pair(constituent,isCharged) );
      if(isCharged)
	{
	  chargedJet += constituent;
	  ++nCharged;      
	}
    }
  result.n=(Int_t) allConstituents.size();
  result.nch=nCharged;

  //stop here if <2 charged
  if( nCharged < 2 ) return result;

  //compute the pull
  double jetPt        = jet.Pt(),        jetPhi=jet.Phi(),                 jetRapidity=jet.Rapidity();
  double jetPtCharged = chargedJet.Pt(), jetPhiCharged = chargedJet.Phi(), jetRapidityCharged = chargedJet.Rapidity();
  TVector2 r(0,0);
  TVector2 pullAll(0,0);
  TVector2 pullCharged(0,0);
  for(size_t idx = 0; idx<allConstituents.size(); ++idx)
    {
      TLorentzVector &cp4=allConstituents[idx].first;
      bool &isCharged=allConstituents[idx].second;
      double constituentPt       = cp4.Pt();
      double constituentPhi      = cp4.Phi();
      double constituentRapidity = cp4.Rapidity();
      r.Set( constituentRapidity - jetRapidity, TVector2::Phi_mpi_pi( constituentPhi - jetPhi ) );
      pullAll += ( constituentPt / jetPt ) * r.Mod() * r;
      //calculate TVector using only charged tracks
      if( isCharged )
	r.Set( constituentRapidity - jetRapidityCharged, TVector2::Phi_mpi_pi( constituentPhi - jetPhiCharged ) );
      pullCharged += ( constituentPt / jetPtCharged ) * r.Mod() * r;
    }
  
  result.pull=pullAll;
  result.chPull=pullCharged;
  return result;
}


//
Float_t computeMT(TLorentzVector &a, TLorentzVector &b)
{
  return TMath::Sqrt(2*a.Pt()*b.Pt()*(1-TMath::Cos(a.DeltaPhi(b))));
}


// 
float getLeptonEnergyScaleUncertainty(int l_id,float l_pt,float l_eta)
{
  float unc(0.02);
  
  // electron uncertainties for 8 TeV cf. AN-14-145   
  if(abs(l_id)==11 || abs(l_id)==1100 || abs(l_id)==2111)
    {
      float par0(-2.27e-02), par1(-7.01e-02), par2(-3.71e-04);
      if (fabs(l_eta) > 0.8 && fabs(l_eta)<1.5)
        {
          par0 = -2.92e-02;
          par1 = -6.59e-02;
          par2 = -7.22e-04;
        }
      else if(fabs(l_eta)>1.5)
        {
          par0 = -2.27e-02;
          par1 = -7.01e-02;
          par2 = -3.71e-04;
        }
      unc=fabs(par0 * TMath::Exp(par1 * l_pt) + par2);
    }

  return unc;
}

//
FactorizedJetCorrector *getFactorizedJetEnergyCorrector(TString baseDir, bool isMC)
{
  gSystem->ExpandPathName(baseDir);
  
  //order matters: L1 -> L2 -> L3 (-> Residuals)
  std::vector<std::string> jetCorFiles;
  TString pf("Summer15_25nsV7_");
  pf += (isMC ? "MC" : "DATA");
  std::cout << "Jet energy corrections from " << baseDir+"/"+pf+"_*_AK4PFchs.txt" << std::endl;
  jetCorFiles.push_back((baseDir+"/"+pf+"_L1FastJet_AK4PFchs.txt").Data());
  jetCorFiles.push_back((baseDir+"/"+pf+"_L2Relative_AK4PFchs.txt").Data());
  jetCorFiles.push_back((baseDir+"/"+pf+"_L3Absolute_AK4PFchs.txt").Data());
  if(!isMC) jetCorFiles.push_back((baseDir+"/"+pf+"_L2L3Residual_AK4PFchs.txt").Data());
  
  //init the parameters for correction
  std::vector<JetCorrectorParameters> corSteps;
  for(size_t i=0; i<jetCorFiles.size(); i++) corSteps.push_back(JetCorrectorParameters(jetCorFiles[i]));
  
  //return the corrector
  return new FactorizedJetCorrector(corSteps);
}


//Sources :  https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution
std::vector<float> getJetResolutionScales(float pt, float eta, float genjpt)
{
  std::vector<float> res(3,1.0);

  float ptSF(1.0), ptSF_err(0.0);
  if(TMath::Abs(eta)<0.5)       { ptSF=1.109; ptSF_err = 0.008; }
  else if(TMath::Abs(eta)<0.8)  { ptSF=1.138; ptSF_err = 0.013; }
  else if(TMath::Abs(eta)<1.1)  { ptSF=1.114; ptSF_err = 0.013; }
  else if(TMath::Abs(eta)<1.3)  { ptSF=1.123; ptSF_err = 0.024; }
  else if(TMath::Abs(eta)<1.7)  { ptSF=1.084; ptSF_err = 0.011; }  
  else if(TMath::Abs(eta)<1.9)  { ptSF=1.082; ptSF_err = 0.035; }
  else if(TMath::Abs(eta)<2.1)  { ptSF=1.140; ptSF_err = 0.047; }
  else if(TMath::Abs(eta)<2.3)  { ptSF=1.067; ptSF_err = 0.053; }
  else if(TMath::Abs(eta)<2.5)  { ptSF=1.177; ptSF_err = 0.041; }
  else if(TMath::Abs(eta)<2.8)  { ptSF=1.364; ptSF_err = 0.039; }
  else if(TMath::Abs(eta)<3.0)  { ptSF=1.857; ptSF_err = 0.071; }
  else if(TMath::Abs(eta)<3.2)  { ptSF=1.328; ptSF_err = 0.022; }
  else if(TMath::Abs(eta)<5.0)  { ptSF=1.160; ptSF_err = 0.029; }

  res[0] = TMath::Max((Float_t)0.,(Float_t)(genjpt+(ptSF)*(pt-genjpt)))/pt;
  res[1] = TMath::Max((Float_t)0.,(Float_t)(genjpt+(ptSF-ptSF_err)*(pt-genjpt)))/pt;
  res[2] = TMath::Max((Float_t)0.,(Float_t)(genjpt+(ptSF+ptSF_err)*(pt-genjpt)))/pt;
  
  return res;
}

// Histogram tool for automatic creation of 2D uncertainty histograms
HistTool::HistTool(unsigned int nsyst) :
  nsyst_(nsyst)
{}

//
void HistTool::addHist(TString title, TH1* hist) {
  if(hist->InheritsFrom("TH2")) {
    all2dPlots_[title]=(TH2 *)hist;
  }
  else {
    allPlots_[title] = hist;
    if (nsyst_ > 0) {
      all2dPlots_[title+"_syst"] = new TH2F(title+"_syst", hist->GetTitle(), hist->GetNbinsX(), hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax(), nsyst_+1, -0.5, nsyst_+0.5);
      all2dPlots_[title+"_syst"]->SetXTitle(hist->GetXaxis()->GetTitle());
      all2dPlots_[title+"_syst"]->SetYTitle("Variation (0=default)");
    }
  }
}

//
void HistTool::fill(TString title, double value, std::vector<double> weights) {
  if (not allPlots_.count(title)) {
    std::cout << "Histogram " << title << " not registered, not filling." << std::endl;
    return;
  }
  if(allPlots_[title]->InheritsFrom("TH2")) return;
  allPlots_[title]->Fill(value, weights[0]);
  if (nsyst_ > 0) {
    if (weights.size() > nsyst_)
      std::cout << "WARNING: Size of uncertainty weight vector larger than uncertainty histogram size." << std::endl;
    all2dPlots_[title+"_syst"]->Fill(value, 0., weights[0]);
    for (unsigned int i = 1; i < weights.size(); ++i) {
      all2dPlots_[title+"_syst"]->Fill(value, i, weights[0]*weights[i]);
    }
  }
}


