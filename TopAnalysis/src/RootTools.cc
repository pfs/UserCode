#include "TGraph.h"
#include "TTree.h"
#include "TString.h"
#include "TROOT.h"
#include "TMath.h"

//
TGraph* ll2dContourPlot(TTree *t, TString x, TString y, Double_t pmin, Double_t pmax)
{

  t->Draw(y+":"+x, "quantileExpected == 1");
  TGraph *bestFit = (TGraph*) gROOT->FindObject("Graph")->Clone("bestfit");

  t->Draw(y+":"+x, Form("%f <= quantileExpected && quantileExpected <= %f && quantileExpected != 1",pmin,pmax));  
  TGraph *gr = (TGraph*) gROOT->FindObject("Graph")->Clone();

  //re-order the points
  Double_t x0 = bestFit->GetX()[0];
  Double_t y0 = bestFit->GetY()[0];
  Double_t *xi = gr->GetX(), *yi = gr->GetY();
  int n = gr->GetN();
  for (int i = 0; i < n; ++i) { xi[i] -= x0; yi[i] -= y0; }
  gr->Sort(&TGraph::CompareArg);
  for (int i = 0; i < n; ++i) { xi[i] += x0; yi[i] += y0; }

  bestFit->Delete();
  return gr;
}
