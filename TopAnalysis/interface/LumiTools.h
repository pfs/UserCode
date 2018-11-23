#ifndef _lumi_tools_h_
#define _lumi_tools_h_

#include <vector>
#include <string>
#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TGraphAsymmErrors.h"
#include "TRandom3.h"
#include "TH1.h"

class LumiTools
{
 public:
  typedef std::pair<TString,float> RunPeriod_t;
  LumiTools(TString era="era2017",TH1 *genPuH=0);
  TString assignRunPeriod();
  std::map<Int_t,Float_t> lumiPerRun();
  std::vector<Float_t> pileupWeight(Float_t genPu,TString period="");
 private:
  void parseLumiInfo();
  std::map<Int_t,Float_t> lumiPerRun_;
  TH1F *countH_;
  std::vector<RunPeriod_t> runPeriods_;
  void defineRunPeriods();
  //std::map<TString, std::vector<TGraph *> > puWgtGr_;
  std::map<TString, std::vector<TH1 *> > puWgtGr_;
  void parsePileupWeightsMap(TH1 *);
  TString era_;
  TRandom3 rand_;
};

#endif
