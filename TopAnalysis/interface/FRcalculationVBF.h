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
#include <string>

using namespace RooFit;
using namespace std;

TString corrCat(TString dataset, TString cat){
  if(dataset.Contains("JetHT") && cat.Contains("HighVPt")){
    cat = "HighVPtHighMJJA";
  }
  return cat;
}

void fixZeroBins(TH1D * h){
  cout<<h->GetName()<<endl;
  for(int i = 0; i < h->GetXaxis()->GetNbins(); i++){
    cout<<h->GetBinContent(i+1) <<endl;
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

TH1D*  rawFRMaker(TH1D& Num1, TH1D& Den1, TString fname = "Data13TeV_SinglePhoton_2017.root", TString catName = "LowMJJA", TString det = "EB", TString binvar="Mjj"){
  TFile * f = TFile::Open(fname);
  TString tmpCat = corrCat(fname,catName);
  TH2D * Num = (TH2D*)f->Get(tmpCat+"_tight"+binvar+det);
  TH2D * Den = (TH2D*)f->Get(tmpCat+"_loose"+binvar+det);
  if(tmpCat == catName && catName == "HighVPtA"){
    Num = (TH2D*)f->Get("HighVPtHighMJJA_tight"+binvar+det);
    Den = (TH2D*)f->Get("HighVPtHighMJJA_loose"+binvar+det);
    TH2D * tmpNum = (TH2D*)f->Get("HighVPtLowMJJA_tight"+binvar+det);
    TH2D * tmpDen = (TH2D*)f->Get("HighVPtLowMJJA_loose"+binvar+det);
    Num->Add(tmpNum);
    Den->Add(tmpDen);
  }
  
  Num1 = *((TH1D*)Num->ProjectionY("rawNum"));
  //temporary fix for tight definition vs relaxed tight
  //if(catName.Contains("LowVPt")) Num1.Scale(5042.000000/4536.000000);
  Den1 = *((TH1D*)Den->ProjectionY("Den"));
  TH1D * ratio = (TH1D*)Num1.Clone("rawFR");
  ratio->Divide(&Den1);
  return ratio;
  
}

void makeFile(TString binvar){
  const int nCat = 4;
  TString cats[nCat] = {"LowVPtHighMJJ","HighVPtHighMJJ","HighVPtLowMJJ","HighVPt"};
  TString det[2] = {"EB","EE"};
  std::vector<TH1D*> hists;
  for(int iCat = 0; iCat < nCat; iCat++){
    for(int iDet = 0; iDet < 2; iDet++){
      TFile * f = TFile::Open(cats[iCat]+"A_FR_"+det[iDet]+".root");
      if(f){
	TString hname(binvar == "Mjj"? "fracRatio_Merged":"fracRatio");
	TH1D * h = (TH1D*) f->Get(hname);
	h->SetName(cats[iCat]+"_"+det[iDet]);
	hists.push_back(h);
      }
    }
  }

  TFile * f = new TFile("fakeRatios.root","recreate");
  f->cd();
  for(unsigned int iHist = 0; iHist < hists.size(); iHist++)
    hists[iHist]->Write();
  f->Close();
}

void promptEstimator(TString binvar= "Mjj", TString dataname = "Data13TeV_SinglePhoton_2017.root", TString catName = "LowMJJA", TString det = "EB", TString pname = "MC13TeV_GJets.root", TString qcd = "Data13TeV_JetHTQCD_2017.root", float oer = 1.){
  double pHigh = 0.00996;
  if (det == "EE")
    pHigh = 0.0271;
  TFile * fD = TFile::Open(dataname);
  TString tmpCat = corrCat(dataname,catName);
  TH2D * data = 0;
  TString gDataHist("relaxedTight"+binvar),qcdDataHist("tmpQCD"+binvar),gMCHist("tight"+binvar);
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // tmpQCD <--> relaxedTight
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if(tmpCat == catName && catName == "HighVPtA"){
    data = (TH2D*)fD->Get("HighVPtHighMJJA_"+gDataHist+det);   
    TH2D * tmpdata = (TH2D*)fD->Get("HighVPtLowMJJA_"+gDataHist+det);      
    data->Add(tmpdata);
  } else data = (TH2D*)fD->Get(tmpCat+"_"+gDataHist+det);
  if(!data) {   
    cout << "Wrong category name "<< catName<<endl;
    return;
  }
  std::vector<TH1D*> dataSihih;
  std::vector<double> trueNData;
  for(int i = 0; i<data->GetYaxis()->GetNbins(); i++){
    std::stringstream name;
    name <<"data_"+binvar+"bin"<<i+1;
    dataSihih.push_back((TH1D*)data->ProjectionX(name.str().c_str(),i+1,i+1));
    trueNData.push_back(dataSihih[dataSihih.size()-1]->Integral(1,pHigh));
  }
  TH1D * tmpMJJ = (TH1D*)data->ProjectionY("DataMJJ");
  RooRealVar sihih("sihih","sihih",0,0.1);
  TFile * fP = TFile::Open(pname);
  TH2D * prompt2 = (TH2D*)fP->Get(catName+"_"+gMCHist+det);
  if(catName == "HighVPtA"){
    prompt2 = (TH2D*)fP->Get("HighVPtHighMJJA_"+gMCHist+det);
    TH2D * tmpPrompt2 = (TH2D*)fP->Get("HighVPtLowMJJA_"+gMCHist+det);
    prompt2->Add(tmpPrompt2);
  }
  TH1D * prompt = (TH1D*)prompt2->ProjectionX("Prompt");  
  prompt->Scale(1/prompt->Integral());
  double pFrac = prompt->Integral(1,prompt->GetXaxis()->FindBin(pHigh)); 
  prompt->Rebin(4);
  prompt->Scale(1/prompt->Integral());
  RooDataHist dP("dP","dP",sihih,Import(*prompt));
  RooHistPdf promptPdf("promptPdf","promptPdf",sihih,dP,0) ;
  //  if(qcd.Contains("2016") && det == "EE") qcd = dataname;
  TFile * fQ = TFile::Open(qcd);
  tmpCat = corrCat(qcd,catName);
  TH2D * QCD2 = (TH2D*)fQ->Get(tmpCat+"_"+qcdDataHist+det);

  //  TH2D * QCD2 = (TH2D*)fD->Get(tmpCat+"_tmpQCDMjj"+det);

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
    name <<"dh_"+binvar+"bin"<<i+1;
    RooDataHist dh(name.str().c_str(),name.str().c_str(),sihih,Import(*dataSihih[i]));
    name.str("");
    name <<"frac_"+binvar+"bin"<<i+1;
    RooRealVar frac(name.str().c_str(),name.str().c_str(),0.8,0,1.) ;
    name.str("");
    name <<"norm_"+binvar+"bin"<<i+1;
    RooRealVar norm(name.str().c_str(),name.str().c_str(),0.9*dataSihih[i]->GetEntries(),0,dataSihih[i]->GetEntries()) ;

    name.str("");
    name <<"model_"+binvar+"bin"<<i+1;
    RooAddPdf model(name.str().c_str(),name.str().c_str(),RooArgList(promptPdf,QCDPdf),RooArgList(frac)) ;
    RooExtendPdf pext("pext", "pext", model, norm);

    pext.chi2FitTo(dh);



    nPrompt->SetBinContent(i+1,pFrac*dataSihih[i]->GetEntries()*frac.getVal());
    nFake->SetBinContent(i+1,oer*fFrac*(1-frac.getVal())*dataSihih[i]->GetEntries());

    fakeFraction->SetBinContent(i+1,oer*fFrac*(1-frac.getVal()));
    //    if(catName.Contains("LowVPt")) fakeFraction->SetBinContent(i+1,oer*(1-(pFrac*frac.getVal())));
    fakeFraction->SetBinError(i+1,oer*fFrac*frac.getError());
    //    if(catName.Contains("LowVPt")) fakeFraction->SetBinError(i+1,oer*pFrac*frac.getError());

    cout<<pFrac<<" "<<fFrac<<" "<<frac.getVal()<<endl;
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
  rawRatio = rawFRMaker(rawnum,den,dataname,catName,det,binvar);
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


