#include "TopLJets2015/TopAnalysis/interface/FakeRateTool.h"
#include <iostream>
#include "TF1.h"
using namespace std;
FakeRateTool::FakeRateTool(TString era, TString fname){
  f = TFile::Open(era+"/"+fname);
  HighMJJFR = 0.014331;  
  HighMJJFRSyst = 0.006;
  LH = (TH1D*)f->Get("LowVPtHighMJJ_EB"); 
  LH->Fit("pol1");
  fitL = (TF1*)LH->GetListOfFunctions()->FindObject("pol1");
  HHb = (TH1D*)f->Get("HighVPtHighMJJ_EB"); 
  HHe = (TH1D*)f->Get("HighVPtHighMJJ_EE");
  HLb = (TH1D*)f->Get("HighVPtLowMJJ_EB");
  HLe = (TH1D*)f->Get("HighVPtLowMJJ_EE");
  Hb = (TH1D*)HHb->Clone("H_EB");
  Hb->Add(HLb);
  Hb->Fit("pol1");
  fitHb = (TF1*)Hb->GetListOfFunctions()->FindObject("pol1");
  He = (TH1D*)HHe->Clone("H_EE");
  He->Add(HLe);
  He->Fit("pol1");
  fitHe = (TF1*)He->GetListOfFunctions()->FindObject("pol1");

}

double FakeRateTool::getWeight(TString cat, double mjj, double veta){
   /*****************************/
   /* For VPt-binned            */
   /*****************************/
//   if (cat.Contains("LowVPt")) return 0.003294;
//   TString det = "EB";
//   if (! ( fabs(veta) < 1.4442) ) det = "EE";
//   cat ="HighVPt_"+det;
//   TH1D * h = (TH1D*)f->Get(cat);
//   int bin = h->GetXaxis()->FindBin(mjj);
  


  /*****************************/
  /* For MJJ-binned            */
  /*****************************/
  //  if (cat.Contains("HighMJJ")) return HighMJJFR;
  if (cat.Contains("LowVPt")) return fitL->Eval(mjj);
  TString det = "EB";
  if (! ( fabs(veta) < 1.4442) ) det = "EE";
  if (cat.Contains("HighVPt") && det == "EB")  return fitHb->Eval(mjj);
  if (cat.Contains("HighVPt") && det == "EE")  return fitHe->Eval(mjj);
  // int bin = h->GetXaxis()->FindBin(mjj);
  // return h->GetBinContent(bin);  
  return 0.;
}


