#include "TopLJets2015/TopAnalysis/interface/HistTool.h"
#include "TSystem.h"
#include "TFile.h"
#include "TString.h"
#include "TH2F.h"
#include <iostream>

// Histogram tool for automatic creation of 2D uncertainty histograms
HistTool::HistTool(unsigned int nsyst) :
  nsyst_(nsyst)
{}

//
void HistTool::addHist(TString title, TH1* hist) {
  if(hist->InheritsFrom("TH2")) {
    all2dPlots_[title]=(TH2 *)hist;
    all2dPlots_[title]->SetDirectory(0);
  }
  else {
    allPlots_[title] = hist;
    allPlots_[title]->SetDirectory(0);
    if (nsyst_ > 0) {
      all2dPlots_[title+"_syst"] = new TH2F(title+"_syst", hist->GetTitle(), hist->GetNbinsX(), hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax(), nsyst_+1, -0.5, nsyst_+0.5);
      all2dPlots_[title+"_syst"]->SetXTitle(hist->GetXaxis()->GetTitle());
      all2dPlots_[title+"_syst"]->SetYTitle("Variation (0=default)");
      all2dPlots_[title+"_syst"]->SetDirectory(0);
    }
  }
}

//
void HistTool::fill(TString title, double value, std::vector<double> weights,std::vector<TString> cats) {
  for(auto &c : cats)
    fill(title,value,weights,c);
}

//
void HistTool::fill(TString title, double value, std::vector<double> weights,TString cat) {
  
  if (not allPlots_.count(title)) {
    //std::cout << "Histogram " << title << " not registered, not filling." << std::endl;
    return;
  }

  if(allPlots_[title]->InheritsFrom("TH2")) return;

  //category specific plot, init if needed
  if(cat!=""){
    TString newTitle=cat+"_"+title;
    if(not allPlots_.count(newTitle)) {
      //std::cout << "Histogram " << title << " for cat=" << cat << " not yet started, adding now." << std::endl;
      allPlots_[newTitle]=(TH1 *)allPlots_[title]->Clone(newTitle);
      allPlots_[newTitle]->SetDirectory(0);
      allPlots_[newTitle]->Reset("ICE");
    }
    title=newTitle;
  }

  allPlots_[title]->Fill(value, weights[0]);
  if (nsyst_ > 0) {
    if (weights.size() > nsyst_)
      std::cout << "WARNING: Size of uncertainty weight vector larger than uncertainty histogram size." << std::endl;
    all2dPlots_[title+"_syst"]->Fill(value, 0., weights[0]);
    for (unsigned int i = 1; i < weights.size(); ++i) {
      all2dPlots_[title+"_syst"]->Fill(value, i, weights[0]*weights[i]);
    }
  }
}


void HistTool::fill2D(TString title, double valueX, double valueY, std::vector<double> weights,std::vector<TString> cats) {
  for(auto &c : cats)
    fill2D(title,valueX,valueY,weights,c);
}


void HistTool::fill2D(TString title, double valueX, double valueY, std::vector<double> weights,TString cat) {
  
  if (not all2dPlots_.count(title)) {
    std::cout << "Histogram " << title << " not registered, not filling." << std::endl;
    return;
  }


  //category specific plot, init if needed
  if(cat!=""){
    TString newTitle=cat+"_"+title;
    if(not all2dPlots_.count(newTitle)) {
      std::cout << "Histogram " << title << " for cat=" << cat << " not yet started, adding now." << std::endl;
      all2dPlots_[newTitle]=(TH2 *)all2dPlots_[title]->Clone(newTitle);
      all2dPlots_[newTitle]->SetDirectory(0);
      all2dPlots_[newTitle]->Reset("ICE");
    }
    title=newTitle;
  }

  all2dPlots_[title]->Fill(valueX,valueY, weights[0]);
}

