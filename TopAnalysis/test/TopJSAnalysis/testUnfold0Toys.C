// Author: Stefan Schmitt
// DESY, 14.10.2008

//  Version 17.6, in parallel to changes in TUnfold
//
//  History:
//    Version 17.5, in parallel to changes in TUnfold
//    Version 17.4, in parallel to changes in TUnfold
//    Version 17.3, in parallel to changes in TUnfold
//    Version 17.2, in parallel to changes in TUnfold
//    Version 17.1, in parallel to changes in TUnfold
//    Version 17.0, updated for using the classes TUnfoldDensity, TUnfoldBinning
//    Version 16.1, parallel to changes in TUnfold
//    Version 16.0, parallel to changes in TUnfold
//    Version 15, with automated L-curve scan
//    Version 14, with changes in TUnfoldSys.cxx
//    Version 13, include test of systematic errors
//    Version 12, catch error when defining the input
//    Version 11,  print chi**2 and number of degrees of freedom
//    Version 10,  with bug-fix in TUnfold.cxx
//    Version 9,  with bug-fix in TUnfold.cxx and TUnfold.h
//    Version 8,  with bug-fix in TUnfold.cxx and TUnfold.h
//    Version 7,  with bug-fix in TUnfold.cxx and TUnfold.h
//    Version 6a, fix problem with dynamic array allocation under windows
//    Version 6, bug-fixes in TUnfold.C
//    Version 5, replace main() by testUnfold1()
//    Version 4, with bug-fix in TUnfold.C
//    Version 3, with bug-fix in TUnfold.C
//    Version 2, with changed ScanLcurve() arguments
//    Version 1, remove L curve analysis, use ScanLcurve() method instead
//    Version 0, L curve analysis included here

/*
  This file is part of TUnfold.

  TUnfold is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  TUnfold is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with TUnfold.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <string> 
#include <sstream>

#define SSTR( x ) static_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

#include <TError.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TFitter.h>
#include <TF1.h>
#include <TStyle.h>
#include <TVector.h>
#include <TGraph.h>
#include <TFile.h>
#include <TLegend.h>

#include "TUnfoldDensity.h"

// #define VERBOSE_LCURVE_SCAN

using namespace std;

///////////////////////////////////////////////////////////////////////
// 
//  Test program for the classes TUnfold and related
//
//  (1) Generate Monte Carlo and Data events
//      The events consist of
//        signal
//        background
//
//      The signal is a resonance. It is generated with a Breit-Wigner,
//      smeared by a Gaussian
//
//  (2) Unfold the data. The result is:
//      The background level
//      The shape of the resonance, corrected for detector effects
//
//      Systematic errors from the MC shape variation are included
//      and propagated to the result
//
//  (3) fit the unfolded distribution, including the correlation matrix
//
//  (4) save six plots to a file testUnfold1.ps
//        1  2  3
//        4  5  6
//      1: 2d-plot of the matrix decsribing the migrations
//      2: generator-level distributions
//             blue: unfolded data, total errors
//             green: unfolded data, statistical errors
//             red: generated data
//             black: fit to green data points
//      3: detector level distributions
//             blue: unfoldede data, folded back through the matrix
//             black: Monte Carlo (with wrong peal position)
//             blue: data
//      4: global correlation coefficients
//      5: chi**2 as a function of log(tau)
//           the star indicates the final choice of tau
//      6: the L curve
//
///////////////////////////////////////////////////////////////////////

TRandom *rnd=0;

TH2 *gHistInvEMatrix;

TVirtualFitter *gFitter=0;

void divideByBinWidth(TH1* hist) {
  for(Int_t i=1; i<=hist->GetNbinsX(); i++) {
     hist->SetBinContent(i,hist->GetBinContent(i)/hist->GetBinWidth(i));
     hist->SetBinError(i,hist->GetBinError(i)/hist->GetBinWidth(i));
  }
}
void normalize(TH1* hist) {
  hist->Scale(1./hist->Integral());
}

int testUnfold0Toys(TString observable = "mult", TString flavor = "all", int nToys = 100)
{
  // switch on histogram errors
  TH1::SetDefaultSumw2();

  // show fit result
  gStyle->SetOptFit(1111);
  
  TString basepath = "/afs/cern.ch/work/m/mseidel/TopAnalysis/CMSSW_8_0_26_patch1/src/TopLJets2015/TopAnalysis/unfolding/fill/";
  
  //============================================
  // generate MC distribution
  //
  TFile file(basepath+"MC13TeV_TTJets.root", "READ");
  
  //TString observable  = "tau21";
          observable += "_charged";
  bool reg = false;
  
  TH1D *histMgenMC = (TH1D*) file.Get(observable+"_"+flavor+"_responsematrix_px");
  //histMgenMC->Scale(1./histMgenMC->Integral());
  TH1D *histMdetMC = (TH1D*) file.Get(observable+"_"+flavor+"_responsematrix_py");
  //histMdetMC->Scale(1./histMdetMC->Integral());
  TH2D *histMdetGenMC = (TH2D*) file.Get(observable+"_"+flavor+"_responsematrix");
  TH1D* histMdetMCsig = histMdetGenMC->ProjectionY("histMdetNonGenMC", 1, -1);
  TH2D* histMdetGenMCSig = (TH2D*) histMdetGenMC->Clone("histMdetGenMCSig");
  
  double lumi = 16551.;
  
  double dataMCSF = lumi*832.;
  histMgenMC->Scale(dataMCSF);
  histMdetMC->Scale(dataMCSF);
  histMdetMCsig->Scale(dataMCSF);
  
  double integral = histMdetGenMC->Integral(0,-1,0,-1);
  //std::cout << "integral: " << integral << std::endl;
  double nExpectedEvents = integral*lumi*832.;
  std::vector<double> wheel = { 0. };
  for (int i = 0; i < histMdetGenMC->GetNcells()+1; ++i) {
    //std::cout << wheel.back() << std::endl;
    wheel.push_back(wheel.back() + histMdetGenMC->GetBinContent(i));
  }
  
  //histMdetGenMC->Scale(1./histMdetGenMC->Integral());
  //histMdetGenMC->Scale(1000000.);
  
  //Set non-reco bins to 0
  //for (int i = 0; i < histMdetGenMCSig->GetNbinsX()+2; ++i) histMdetGenMCSig->SetBinContent(i, 0, 0.);
  //Set non-gen bins to 0
  for (int i = 0; i < histMdetGenMCSig->GetNbinsY()+2; ++i) histMdetGenMCSig->SetBinContent(0, i, 0.);
  //Set non-reco non-gen bin to 0
  //histMdetGenMC->SetBinContent(0, 0, 36500000.);
  
  TRandom3 random(1);
  
  TH2F* hPull = new TH2F("pull", "pull;bin;pull", histMgenMC->GetNbinsX()+2, 0, histMgenMC->GetNbinsX()+2, 50, -5, 5);
  std::vector<std::vector<double> > pseudoresults;
  
  TFile outfile("unfolding/toys/"+observable+"_"+flavor+"_toys.root", "RECREATE");
  
  for (int t = 0; t < nToys; ++t) {
    std::cout<<"toy number "<<t<<std::endl;
    //============================================
    // generate toy distribution
    //
    TH2D* histMdetGenToy = (TH2D*) histMdetGenMC->Clone("histMdetGenToy");
    histMdetGenToy->Scale(dataMCSF);
    histMdetGenToy->Reset();
    
    int drawn = 0;
    while (drawn < nExpectedEvents) {
      double pick = random.Uniform(integral);
      //std::cout << pick << ":";
      auto low = std::lower_bound(wheel.begin(), wheel.end(), pick);
      int position = low - wheel.begin() - 1;
      //std::cout << position << " ";
      //int position = -1;
      //for (unsigned int i = 0; i < wheel.size(); ++i) {
      //  if (pick < wheel[i]) break;
      //  ++position;
      //}
      //std::cout << position << " ";
      histMdetGenToy->SetBinContent(position, histMdetGenToy->GetBinContent(position) + 1);
      ++drawn;
    }
    for (int i = 0; i < histMdetGenToy->GetNcells()+2; ++i) {
      histMdetGenToy->SetBinError(i, sqrt(histMdetGenToy->GetBinContent(i)));
    }
    //histMdetGenToy->FillRandom(histMdetGenMC, nExpectedEvents);    
    
    TH1D *histMdetData = histMdetGenToy->ProjectionY("histMdetData");
    TH1D *histMgenToy = histMdetGenToy->ProjectionX("histMgenToy");

    //double dataMCSF = 36500.; //histMdetData->Integral()/histMdetMC->Integral();
    //dataMCSF = histMdetData->Integral()/histMdetMC->Integral()*100;
    //dataMCSF = 100000.;
    
    //if (not DATA) {
    //  histMdetData->Scale(dataMCSF*832.); // For MC test
    //}
    //histMdetData->Scale(100.);
    //histMdetData->SetBinContent(0, 0.);
    
    //for (int i = 0; i < histMdetData->GetNbinsX()+2; ++i) cout << histMdetData->GetBinContent(i) << " +/- " << histMdetData->GetBinError(i) << endl;
    
    //for (int i = 0; i < histMdetNonGenMC->GetNbinsX()+2; ++i) cout << histMdetNonGenMC->GetBinContent(i)*dataMCSF << endl;
    
    // backgrounds
    std::map<TString, double> backgrounds;
    backgrounds["MC13TeV_TTJets.root"] = 832.; // reco but not generated
    TH1D *histMdetDataBGSubtracted = (TH1D*) histMdetData->Clone();

    TFile bkgfile(basepath+"MC13TeV_TTJets.root", "READ");
    TH2D *histMdetGenMCbkg = (TH2D*) bkgfile.Get(observable+"_"+flavor+"_responsematrix");
    TH1D *histMdetNonGenMCbkg = histMdetGenMCbkg->ProjectionY("histMdetNonGenMCbkg_MC13TeV_TTJets.root", 0, 0);
    TH1D *histMdetNonGenMCall = histMdetGenMCbkg->ProjectionY("histMdetNonGenMCall_MC13TeV_TTJets.root");
    //double SF = histMdetDataBGSubtracted->Integral(1, -1) * histMdetNonGenMCall->Integral(1,-1);
    //std::cout << "SF=" << SF << std::endl;
    //histMdetDataBGSubtracted->Add(histMdetNonGenMCbkg, -SF*832.);
    //*
    for (int i = 0; i < histMdetDataBGSubtracted->GetNbinsX()+2; ++i) {
      double sf = histMdetNonGenMCbkg->GetBinContent(i) / histMdetNonGenMCall->GetBinContent(i);
      histMdetDataBGSubtracted->SetBinContent(i, (1.-sf)*histMdetDataBGSubtracted->GetBinContent(i));
      if (std::isnan(histMdetDataBGSubtracted->GetBinContent(i)))
        histMdetDataBGSubtracted->SetBinContent(i, 0);
      histMdetDataBGSubtracted->SetBinError(i, (1.-sf)*histMdetDataBGSubtracted->GetBinError(i));
      if (std::isnan(histMdetDataBGSubtracted->GetBinError(i)))
        histMdetDataBGSubtracted->SetBinError(i, 0);
      
    }
    //*/
    /*
    for (auto background : backgrounds) {
      TFile bkgfile(basepath+background.first, "READ");
      TH2D *histMdetGenMCbkg = (TH2D*) bkgfile.Get(observable+"_"+flavor+"_responsematrix");
      TH1D *histMdetNonGenMCbkg = histMdetGenMCbkg->ProjectionY("histMdetNonGenMCbkg_"+background.first, 0, 0);
      histMdetDataBGSubtracted->Add(histMdetNonGenMCbkg, -dataMCSF);
    }
    //*/
    //for (int i = 0; i < histMdetDataBGSubtracted->GetNbinsX()+2; ++i) {
    //  std::cout << histMdetDataBGSubtracted->GetBinContent(i) << " +/- " << histMdetDataBGSubtracted->GetBinError(i) << std::endl;
    //}
    

    
    //histMdetData->SetBinContent(0, 1000000.);


    //=========================================================================
    // set up the unfolding
    // define migration matrix
    TUnfold::ERegMode regMode = TUnfold::kRegModeCurvature;
    TUnfold::EConstraint constraint = TUnfold::kEConstraintArea;
    TUnfoldDensity::EDensityMode densityFlags=TUnfoldDensity::kDensityModeBinWidth;
    
    TUnfoldDensity unfold(histMdetGenMCSig,TUnfold::kHistMapOutputHoriz,regMode,constraint,densityFlags);

    // define input and bias scame
    // do not use the bias, because MC peak may be at the wrong place
    // watch out for error codes returned by the SetInput method
    // errors larger or equal 10000 are fatal:
    // the data points specified as input are not sufficient to constrain the
    // unfolding process
    if(unfold.SetInput(histMdetDataBGSubtracted)>=10000) {
      std::cout<<"Unfolding result may be wrong\n";
    }

    //========================================================================
    // the unfolding is done here
    //
    // scan L curve and find best point
    //Int_t nScan=100;
    //// use automatic L-curve scan: start with taumin=taumax=0.0
    //Double_t tauMin=1e-10;
    //Double_t tauMax=1e-2;
    //TSpline *rhoScan;
    //TUnfoldDensity::EScanTauMode tauflag = TUnfoldDensity::kEScanTauRhoAvg;

    //int iBest = unfold.ScanTau(nScan,tauMin,tauMax,&rhoScan,tauflag);
    //
    //// create graphs with one point to visualize best choice of tau
    //double tau[1], rho[1];
    
    //rhoScan->GetKnot(iBest, tau[0], rho[0]);
    //TGraph *bestRho = new TGraph(1,tau,rho);
    //double opt_tau = unfold.GetTau();
    //cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  TAU,iBest --> " << tau[0] << "  " << opt_tau << "  " << rho[0] << endl;

    //==========================================================================
    // print some results
    //
    //std::cout<<"toy number "<<t<<std::endl;
    //std::cout<<"tau="<<unfold.GetTau()<<"\n";
    //std::cout<<"chi**2="<<unfold.GetChi2A()<<"+"<<unfold.GetChi2L()
    //         <<" / "<<unfold.GetNdf()<<"\n";
    //std::cout<<"chi**2(sys)="<<unfold.GetChi2Sys()<<"\n";
    
    // Add systematic uncertainties
    //std::vector<TString> uncertainties;
    //uncertainties.push_back("MC13TeV_TTJets_fsrup.root");
    //uncertainties.push_back("MC13TeV_TTJets_fsrdn.root");
    //uncertainties.push_back("MC13TeV_TTJets_hdampup.root");
    //uncertainties.push_back("MC13TeV_TTJets_hdampdn.root");
    //uncertainties.push_back("MC13TeV_TTJets_isrup.root");
    //uncertainties.push_back("MC13TeV_TTJets_isrdn.root");
    //uncertainties.push_back("MC13TeV_TTJets_ueup.root");
    //uncertainties.push_back("MC13TeV_TTJets_uedn.root");
    //uncertainties.push_back("MC13TeV_TTJets_herwig.root");
    //
    //for (TString uncertainty : uncertainties) {
    //  TFile uncfile(basepath+uncertainty, "READ");
    //  TH2D *histMdetGenMCunc = (TH2D*) uncfile.Get(observable+"_"+flavor+"_responsematrix");
    //  for (int i = 0; i<histMdetGenMCunc->GetNbinsY()+2; ++i) histMdetGenMCunc->SetBinContent(0, i, 0.);
    //  unfold.AddSysError(histMdetGenMCunc, uncertainty, TUnfold::kHistMapOutputHoriz, TUnfoldSys::kSysErrModeMatrix);
    //}


    //==========================================================================
    // retreive results into histograms
    
    unfold.DoUnfold(0.);
    //else unfold.DoUnfold(opt_tau);
    
    outfile.cd();
    TH1 *histMunfold=unfold.GetOutput((std::string("Unfolded_")+std::to_string(t)).c_str());
    histMunfold->Write();
    
    for (int i = 0; i < histMunfold->GetNbinsX()+2; ++i) {
      double pull = (histMunfold->GetBinContent(i) - histMgenToy->GetBinContent(i)) / histMunfold->GetBinError(i); // sqrt(pow(histMunfold->GetBinError(i), 2) + pow(histMgenToy->GetBinError(i), 2));
      //std::cout << histMgenMC->GetBinContent(i) << std::endl;
      hPull->Fill(i, pull);
    }
    
    std::vector<double> pseudoresult;
    double integral = 0.;
    for (int i = 1; i < histMunfold->GetNbinsX()+1; ++i) {
      pseudoresult.push_back(histMunfold->GetBinContent(i));
      integral += histMunfold->GetBinContent(i);
    }
    for (int i = 0; i < histMunfold->GetNbinsX(); ++i) {
      pseudoresult[i] = pseudoresult[i]/integral;
    }
    pseudoresults.push_back(pseudoresult);
  }


  //=====================================================================
  // plot some histograms
  TCanvas output("output","output",500,500);
  output.cd();
  
  hPull->FitSlicesY();
  TH1D *pull_1 = (TH1D*) gDirectory->Get("pull_1");
  pull_1->GetYaxis()->SetRangeUser(-1.,1.);
  pull_1->Fit("pol0");
  
  output.SaveAs("unfolding/toys/"+observable+"_"+flavor+"_pull_mean.pdf");
  output.SaveAs("unfolding/toys/"+observable+"_"+flavor+"_pull_mean.png");
  output.SaveAs("unfolding/toys/"+observable+"_"+flavor+"_pull_mean.root");
  
  TH1D *pull_2 = (TH1D*) gDirectory->Get("pull_2");
  pull_2->GetYaxis()->SetRangeUser(0.,2.);
  pull_2->Fit("pol0");
  
  output.SaveAs("unfolding/toys/"+observable+"_"+flavor+"_pull_width.pdf");
  output.SaveAs("unfolding/toys/"+observable+"_"+flavor+"_pull_width.png");
  output.SaveAs("unfolding/toys/"+observable+"_"+flavor+"_pull_width.root");
  
  for (int i = 0; i < pull_1->GetNbinsX()+2; ++i) {
    pull_1->SetBinError(i, pull_2->GetBinContent(i));
  }
  
  hPull->Draw("colz");
  pull_1->SetLineColor(kRed+1);
  pull_1->SetMarkerColor(kRed+1);
  pull_1->SetLineWidth(4);
  pull_1->Draw("e1,same");
  
  output.SaveAs("unfolding/toys/"+observable+"_"+flavor+"_pull.pdf");
  output.SaveAs("unfolding/toys/"+observable+"_"+flavor+"_pull.png");
  output.SaveAs("unfolding/toys/"+observable+"_"+flavor+"_pull.root");
  
  /*
  vector<vector<double>> pseudoresults_T;
  for (unsigned int i = 0; i < pseudoresults[0].size(); i++) {
    std::vector<double> column;
    for (unsigned int j = 0; j < pseudoresults.size(); j++) {
      column.push_back(pseudoresults[j][i]);
    }
    pseudoresults_T.push_back(column);
  }
  vector<vector<double>> cov = compute_covariance_matrix(pseudoresults_T);
  */

  return 0;
}
