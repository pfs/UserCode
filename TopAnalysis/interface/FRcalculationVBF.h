#ifndef __FRcalculationVBF__
#define __FRcalculationVBF__

#include "TH1.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH2.h"
#include "RooDataHist.h"
#include "RooRealVar.h"
#include "RooHistPdf.h"
#include "RooExtendPdf.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include <vector>
#include <iostream>
#include <sstream>



using namespace RooFit;
using namespace std;
void fixZeroBins(TH1D * h){
  for(int i = 0; i < h->GetXaxis()->GetNbins(); i++){
    if(h->GetBinContent(i+1) == 0) h->SetBinContent(i+1, 0.1);
  }
}

TH1D* mergeBin(TH1D * h){
  double x[]={0,500,1000,4000};
  TString name = h->GetName()+TString("_Merged");
  TH1D * ret = new TH1D(name, name, 3, x);
  ret->SetBinContent(1,h->GetBinContent(1));
  ret->SetBinContent(2,h->GetBinContent(2));
  ret->SetBinError(1,h->GetBinError(1));
  ret->SetBinError(2,h->GetBinError(2));
  ret->SetBinContent(3,h->GetBinContent(3)+h->GetBinContent(4));
  ret->SetBinError(3,sqrt(h->GetBinError(3)*h->GetBinError(3)+h->GetBinError(4)*h->GetBinError(4)));  
  return ret;
}

TH1D*  rawFRMaker(TH1D& Num1, TH1D& Den1, TString fname = "Data13TeV_SinglePhoton_2017.root", TString catName = "LowMJJA", TString det = "EB"){
  TFile * f = TFile::Open(fname);
  TH2D * Num = (TH2D*)f->Get(catName+"_tightMjj"+det);
  TH2D * Den = (TH2D*)f->Get(catName+"_looseMjj"+det);

  Num1 = *((TH1D*)Num->ProjectionY("rawNum"));
  Den1 = *((TH1D*)Den->ProjectionY("Den"));
  TH1D * ratio = (TH1D*)Num1.Clone("rawFR");
  ratio->Divide(&Den1);
  return ratio;
}

void makeFile(){
  TFile * f1 = TFile::Open("LowMJJA_FR_EE.root");
  TFile * f2 = TFile::Open("LowMJJA_FR_EB.root");
  TFile * f3 = TFile::Open("HighMJJA_FR_EB.root");
  TH1D * h1 = (TH1D*)f1->Get("fracRatio_Merged");
  h1->SetName("LowMJJ_EE");
  TH1D * h2 = (TH1D*)f2->Get("fracRatio_Merged");
  h2->SetName("LowMJJ_EB");
  TH1D * h3 = (TH1D*)f3->Get("fracRatio");
  h3->SetName("HigMJJ_EB");
  TFile * f = new TFile("fakeRatios.root","recreate");
  f->cd();
  h1->Write();
  h2->Write();
  h3->Write();
  f->Close();
}

void promptEstimator(TString dataname = "Data13TeV_SinglePhoton_2017.root", TString catName = "LowMJJA", TString det = "EB", TString pname = "MC13TeV_GJets.root", TString qcd = "Data13TeV_JetHTQCD_2017.root"){
  double pHigh = 0.00996;
  if (det == "EE")
    pHigh = 0.0271;
  TFile * fD = TFile::Open(dataname);
  TH2D * data = (TH2D*)fD->Get(catName+"_relaxedTightMjj"+det);
  std::vector<TH1D*> dataSihih;
  std::vector<double> trueNData;
  for(int i = 0; i<4; i++){
    std::stringstream name;
    name <<"data_MJJbin"<<i+1;
    dataSihih.push_back((TH1D*)data->ProjectionX(name.str().c_str(),i+1,i+1));
    trueNData.push_back(dataSihih[dataSihih.size()-1]->Integral(1,pHigh));
  }
  TH1D * tmpMJJ = (TH1D*)data->ProjectionY("DataMJJ");
  RooRealVar sihih("sihih","sihih",0,0.1);
  TFile * fP = TFile::Open(pname);
  TH2D * prompt2 = (TH2D*)fP->Get(catName+"_tightMjj"+det);
  TH1D * prompt = (TH1D*)prompt2->ProjectionX("Prompt");  
  prompt->Scale(1/prompt->Integral());
  double pFrac = prompt->Integral(1,prompt->GetXaxis()->FindBin(pHigh)); 
  prompt->Rebin(4);
  prompt->Scale(1/prompt->Integral());
  RooDataHist dP("dP","dP",sihih,Import(*prompt));
  RooHistPdf promptPdf("promptPdf","promptPdf",sihih,dP,0) ;
  TFile * fQ = TFile::Open(qcd);
  TH2D * QCD2 = (TH2D*)fQ->Get(catName+"_allMjj"+det);
  //TH2D * QCD2 = (TH2D*)fD->Get(catName+"_tmpQCDMjj"+det);
  TH1D * QCD = (TH1D*)QCD2->ProjectionX("QCD");
  QCD->Scale(1/QCD->Integral());
  double fFrac = QCD->Integral(1,QCD->GetXaxis()->FindBin(pHigh)); 

  QCD->Rebin(4);
  QCD->Scale(1/QCD->Integral());
  RooDataHist dQ("dQ","dQ",sihih,Import(*QCD));
  RooHistPdf QCDPdf("QCDPdf","QCDPdf",sihih,dQ,0) ;
  std::vector<TCanvas*> c;
  TH1D* nPrompt        = (TH1D*)tmpMJJ->Clone("nPrompt");
  TH1D* nFake          = (TH1D*)tmpMJJ->Clone("nFake");
  TH1D* fakeFraction   = (TH1D*)tmpMJJ->Clone("fakeFraction");
  for(unsigned int i =0; i<dataSihih.size();i++){
    cout <<"Number of entries: "<<dataSihih[i]->GetEntries()<<endl;
    dataSihih[i]->Rebin(4);
    if(dataSihih[i]->GetEntries() == 0) continue;
    fixZeroBins(dataSihih[i]);
    std::stringstream name;
    name <<"dh_MJJbin"<<i+1;
    RooDataHist dh(name.str().c_str(),name.str().c_str(),sihih,Import(*dataSihih[i]));
    name.str("");
    name <<"frac_MJJbin"<<i+1;
    RooRealVar frac(name.str().c_str(),name.str().c_str(),0.8,0,1) ;
    name.str("");
    name <<"norm_MJJbin"<<i+1;
    RooRealVar norm(name.str().c_str(),name.str().c_str(),0.9*dataSihih[i]->GetEntries(),0,dataSihih[i]->GetEntries()) ;

    name.str("");
    name <<"model_MJJbin"<<i+1;
    RooAddPdf model(name.str().c_str(),name.str().c_str(),RooArgList(promptPdf,QCDPdf),RooArgList(frac)) ;
    RooExtendPdf pext("pext", "pext", model, norm);

    pext.chi2FitTo(dh);



    nPrompt->SetBinContent(i+1,pFrac*dataSihih[i]->GetEntries()*frac.getVal());
    nFake->SetBinContent(i+1,fFrac*(1-frac.getVal())*dataSihih[i]->GetEntries());

    fakeFraction->SetBinContent(i+1,fFrac*(1-frac.getVal()));
    fakeFraction->SetBinError(i+1,fFrac*frac.getError());

    cout<<pFrac<<" "<<fFrac<<endl;
    RooPlot* mesframe = sihih.frame() ;
    dh.plotOn(mesframe);
    model.plotOn(mesframe);
    model.plotOn(mesframe,Components(promptPdf),LineStyle(kDashed),LineColor(kGreen)) ;
    model.plotOn(mesframe,Components(QCDPdf),LineStyle(kDashed),LineColor(kRed)) ;
    name.str("");
    name <<"canvas_MJJbin"<<i+1;
    TCanvas * C = new TCanvas(name.str().c_str(),name.str().c_str(),0,0,600,600);
	C->SetLogy();
    C->cd();
    mesframe->Draw();
    c.push_back(C);
  }
  TH1D rawnum, *num, den, *rawRatio, *Ratio;
  rawRatio = rawFRMaker(rawnum,den,dataname,catName,det);
  TH1D* negPrompt = (TH1D*)nPrompt->Clone("negPrompt");
  num = (TH1D*)rawnum.Clone("Num");
  negPrompt->Scale(-1);
  num->Add(negPrompt);
  Ratio = (TH1D*)num->Clone("Ratio");
  Ratio->Divide(&den);
  TH1D* fracRatio = (TH1D*)rawnum.Clone("fracRatio");
  fracRatio->Multiply(fakeFraction);
  TH1D* mfracRatio = mergeBin(fracRatio);
  TH1D* mden =  mergeBin(&den);
  fracRatio->Divide(&den);
  mfracRatio->Divide(mden);


  TFile * out = new TFile(catName+"_FR_"+det+".root","recreate");
  out->cd();
  for(unsigned int i =0; i<c.size();i++){
    c[i]->Write();
  }

  rawnum.Write();
  num->Write();
  den.Write();
  rawRatio->Write();
  Ratio->Write();
  nPrompt->Write();
  nFake->Write();
  tmpMJJ->Write();
  fakeFraction->Write();
  fracRatio->Write();
  mfracRatio->Write();
  QCD->Write();
  prompt->Write();
  out->Close();
  

}

#endif


