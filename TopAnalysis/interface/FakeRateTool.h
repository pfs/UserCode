#include "TFile.h"
#include "TH1.h" 
#include "TString.h"
#include "TF1.h"
class FakeRateTool{
 public:
  FakeRateTool(TString era, TString fname);
  double getWeight(TString cat, double mjj, double veta);
 private:
  TFile * f;
  double HighMJJFR,  HighMJJFRSyst;
  TH1D * LH, * HHe, * HLe, * He,  * HHb, * HLb, * Hb;
  TF1 * fitL, *fitHe, *fitHb;
};
