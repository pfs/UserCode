#include "TopLJets2015/TopAnalysis/interface/FakeRateTool.h"
FakeRateTool::FakeRateTool(TString era, TString fname){
  f = TFile::Open(era+"/"+fname);
  HighMJJFR = 0.014331;  
  HighMJJFRSyst = 0.006;
}

double FakeRateTool::getWeight(TString cat, double mjj, double veta){
  if (cat.Contains("HighMJJ")) return HighMJJFR;
  TString det = "EB";
  if (! ( fabs(veta) < 1.442) ) det = "EE";
  cat = cat+"_"+det;
  TH1D * h = (TH1D*)f->Get(cat);
  int bin = h->GetXaxis()->FindBin(mjj);
  return h->GetBinContent(bin);  
}
