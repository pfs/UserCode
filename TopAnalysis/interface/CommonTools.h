#ifndef _common_tools_h_
#define _common_tools_h_

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"

#include "TVector2.h"

enum FlavourSplitting {NOFLAVOURSPLITTING=0, UDSGSPLITTING=1, CSPLITTING=4, BSPLITTING=5 };

std::map<Int_t,Float_t> lumiPerRun();
Float_t computeMT(TLorentzVector &a, TLorentzVector &b);
FactorizedJetCorrector *getFactorizedJetEnergyCorrector(TString,bool);
std::vector<float> getJetResolutionScales(float pt, float eta, float genjpt);
float getLeptonEnergyScaleUncertainty(int l_id,float l_pt,float l_eta);

struct JetPullInfo_t
{
  Int_t n,nch;
  TVector2 pull,chPull;
};
JetPullInfo_t getPullVector( MiniEvent_t &ev, int ijet);

#endif
