#ifndef _L1PrefireEfficiencyWrapper_h_
#define _L1PrefireEfficiencyWrapper_h_

#include "TEfficiency.h"

#include "TopLJets2015/TopAnalysis/interface/EfficiencyScaleFactorsWrapper.h"

class L1PrefireEfficiencyWrapper 
{
 public:
  L1PrefireEfficiencyWrapper(bool isData,TString era);
  EffCorrection_t getJetBasedCorrection(std::vector<Jet> jets={});
  ~L1PrefireEfficiencyWrapper();  
 private:
  void init(TString era);
  int era_;
  std::map<TString,TEfficiency *> effMapsH_;
};

#endif
