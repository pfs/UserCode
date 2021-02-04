#ifndef _LeptonEfficiencyWrapper_h_
#define _LeptonEfficiencyWrapper_h_

#include "TString.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"

#include <vector>
#include <map>

#include "TopLJets2015/TopAnalysis/interface/ObjectTools.h"

typedef std::pair<float,float> EffCorrection_t;

class LeptonEfficiencyWrapper 
{
 public:
  LeptonEfficiencyWrapper(bool isData,TString era);
  EffCorrection_t getTriggerCorrection(std::vector<int> &pdgId, std::vector<TLorentzVector> &leptons,TString period = "");
  EffCorrection_t getTriggerCorrection(std::vector<Particle> &leptons,TString period = "");
  EffCorrection_t getOfflineCorrection(int pdgId,float pt,float eta,TString period = "");
  EffCorrection_t getOfflineCorrection(Particle lepton,TString period = "");
  EffCorrection_t getOfflineIsoHFCorrection(int pdgId,float hf);
  EffCorrection_t getTrackingCorrection(int nvtx,TString period = "");
  ~LeptonEfficiencyWrapper();  
 private:
  void init(TString era);
  int era_;
  std::map<TString,TH2 *> lepEffH_;
  std::map<TString,TGraphAsymmErrors *> lepEffGr_;
};

#endif
