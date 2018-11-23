#include "TFile.h"
#include "TH1.h" 
#include "TString.h"
class FakeRateTool{
 public:
  FakeRateTool(TString era, TString fname);
  double getWeight(TString cat, double mjj, double veta);
 private:
  TFile * f;
  double HighMJJFR,  HighMJJFRSyst;

};
