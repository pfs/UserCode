#ifndef _EfficiencyScaleFactorsWrapper_h_
#define _EfficiencyScaleFactorsWrapper_h_

#include "TString.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"

#include <vector>
#include <map>

#include "TopLJets2015/TopAnalysis/interface/ObjectTools.h"

typedef std::pair<float,float> EffCorrection_t;

class EfficiencyScaleFactorsWrapper 
{
 public:
  EfficiencyScaleFactorsWrapper(bool isData,TString era);
  EfficiencyScaleFactorsWrapper(bool isData,TString era,std::map<TString,TString> cfgMap);
  EffCorrection_t getDileptonTriggerCorrection(std::vector<Particle> &leptons);
  EffCorrection_t getTriggerCorrection(std::vector<Particle> leptons={}, 
                                       std::vector<Particle> photons={},
                                       std::vector<Particle> jets={},
                                       TString period = "");
  EffCorrection_t getOfflineCorrection(Particle p,TString period="");
  EffCorrection_t getOfflineCorrection(int pdgId,float pt, float eta,TString period = "");
  ~EfficiencyScaleFactorsWrapper();  
 private:
  void init(TString era);
  int era_;
  std::map<TString,TH2 *> scaleFactorsH_;
  std::map<TString,TGraphAsymmErrors *> scaleFactorsGr_;
  std::map<TString,TString> cfgMap_;
};

#endif
