#ifndef _jec_tools_h_
#define _jec_tools_h_

#include <vector>
#include <string>
#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include <CondFormats/JetMETObjects/interface/JetResolutionObject.h>
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "TGraphAsymmErrors.h"
#include "TRandom3.h"

class JECTools
{
 public:

  JECTools(TString era);

  //JES and its uncertainty
  void startFactorizedJetEnergyCorrector(bool isMC);
  FactorizedJetCorrector *getFactorizedJetEnergyCorrector() { return jetCorr_; }
  void applyJetCorrectionUncertainty(MiniEvent_t &ev, TString jecVar, Variation direction);
  TLorentzVector getShiftedJet(TLorentzVector &jp4,TString jecVar,Variation direction);

  //JER smearing and uncertainty
  void smearJetEnergies(MiniEvent_t &ev, Variation option=Variation::NOMINAL);
  TLorentzVector getSmearedJet(TLorentzVector &jp4, float genJet_pt,float rho,Variation option=Variation::NOMINAL);

  //b fragmentation, see https://gitlab.cern.ch/CMS-TOPPAG/BFragmentationAnalyzer
  float computeBFragmentationWeight(MiniEvent_t &ev, TString wgtName);
  float computeSemilepBRWeight(MiniEvent_t &ev, TString var,int pid = 0, bool useabs = true);

 private:
  void startBFragmentationWeights();
  void startSemilepBRWeights();
  TString era_;
  FactorizedJetCorrector *jetCorr_;
  std::map<TString,JetCorrectionUncertainty *> jecUncs_;
  JME::JetResolution *jer_;
  JME::JetResolutionScaleFactor *jerSF_;
  std::map<TString, TGraph*> bfragMap_;
  std::map<TString, std::map<int, float> > semilepBRwgts_;
  TRandom3 *rand_;
};

#endif
