#include "TopLJets2015/TopAnalysis/interface/GeneratorTools.h"
#include "TH1F.h"
#include <iostream>

using namespace std;

//
std::vector< WeightSysts_t > getWeightSysts(TFile *f){
  std::vector< WeightSysts_t > systsOfInterest;

  TH1 *h=(TH1 *) f->Get("analysis/generator_initrwgt");
  if(h==0) return systsOfInterest;

  for(Int_t xbin=1; xbin<=h->GetNbinsX(); xbin++){
    TString label=h->GetXaxis()->GetBinLabel(xbin);
    if(label.Length()==0) continue;
    if(label.Contains("mur=0.5 muf=1"))   systsOfInterest.push_back( WeightSysts_t("muRdn",    xbin-1) );
    if(label.Contains("mur=2 muf=1"))     systsOfInterest.push_back( WeightSysts_t("muRup",    xbin-1) );
    if(label.Contains("mur=1 muf=0.5"))   systsOfInterest.push_back( WeightSysts_t("muFdn",    xbin-1) );
    if(label.Contains("mur=1 muf=2"))     systsOfInterest.push_back( WeightSysts_t("muFup",    xbin-1) );
    if(label.Contains("mur=0.5 muf=0.5")) systsOfInterest.push_back( WeightSysts_t("muRmuFdn", xbin-1) );
    if(label.Contains("mur=2 muf=2"))     systsOfInterest.push_back( WeightSysts_t("muRmuFup", xbin-1) );
    if(label.Contains("NNPDF31_nnlo_hessian_pdfas")) {
      Int_t start=label.Index("Member ")+7;
      TString id=label(start,3);
      if(id.EndsWith(" ")) id=label(start,2);
      if(id.EndsWith("o")) id=label(start,1);
      systsOfInterest.push_back( WeightSysts_t("PDF"+id,xbin-1) );
    }
  }

  return systsOfInterest;
}
