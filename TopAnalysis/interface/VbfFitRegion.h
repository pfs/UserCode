#ifndef VBFFITREGION_H
#define VBFFITREGION_H

#include "TH1.h"
#include "TDirectory.h"
#include "TString.h"
#include "TFile.h"

#include "HiggsAnalysis/CombinedLimit/interface/RooParametricHist.h"
#include "RooAddition.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooDataHist.h"

#include <sstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace RooFit;
const float nMinInBin = 0.1;
////////////////////////////////////
// Acontainer for histograms, etc //
// of SR or CR                    //
////////////////////////////////////

class VbfFitRegion{
 public:
 VbfFitRegion(TString channel, TString v, TString Hist, int bin, bool SR):chan(channel),boson(v),hist(Hist),nBin(bin),isSR(SR){
    this->setHistograms();
    if(isSR)
      for(int i = 0; i< nBin; i++) {
	binsSigCR.push_back(NULL);
	binsBkgCR.push_back(NULL);
      }
    else
      for(int i = 0; i< nBin; i++) {
	binsSigSR.push_back(NULL);
	binsBkgSR.push_back(NULL);
      }
  };
  ~VbfFitRegion(){};

   VbfFitRegion* clone() const { 
     return new VbfFitRegion(*this); 
   } 
 
   void setBkg(){ 
     TString bkgs[]={"Top+VV","DY","QCD","#gamma+jets"}; 
     TDirectory * dir = (TDirectory*)((TFile::Open("plotter_"+chan+".root"))->Get(chan+boson+"_"+hist)); 
     dir->ls();
     hBkg = (TH1F*)dir->Get(chan+boson+"_"+hist+"_"+bkgs[0]); 
     hBkg->SetNameTitle("Background_"+chan+boson,"Background in "+chan+boson); 
     hBkg->SetLineColor(kBlack);
     hBkg->SetMarkerColor(kBlack);

     hBkgCorr = (TH1F*)hBkg->Clone("BackgroundNLO_"+chan+boson); 
     hBkgCorr->SetTitle("Background (NLO) in "+chan+boson); 

     int nRebin =(int)((double)hBkg->GetXaxis()->GetNbins()/(double)nBin);      
     for(int i = 1; i < 4; i++){ 
       TH1F * tmp = (TH1F*)dir->Get(chan+boson+"_"+hist+"_"+bkgs[i]); 
       if (tmp != NULL) {
	 tmp->SetLineColor(kBlack);
	 tmp->SetMarkerColor(kBlack);
	 hBkg->Add(tmp);
	 if(i == 3){
	   TH1F * tmpCorr = correctBackground(tmp);
	   hBkgCorr->Add(tmpCorr);
	 } else hBkgCorr->Add(tmp);
       }
     } 
     hBkg->Rebin(nRebin);
     hBkgCorr->Rebin(nRebin);
     for(int i = 0; i<nRebin; i++){
       if(hBkg->GetBinContent(i) == 0)
	 hBkg->SetBinContent(i,nMinInBin);
       if(hBkgCorr->GetBinContent(i) == 0)
	 hBkgCorr->SetBinContent(i,nMinInBin);
     }
   } 

   void setSig(){ 
     TString signal = "EWK #gammaJJ"; 
     if(boson == "MM") 
       signal = "EWK ZJJ"; 
     TDirectory * dir = (TDirectory*)((TFile::Open("plotter_"+chan+".root"))->Get(chan+boson+"_"+hist)); 
     hSig = (TH1F*)dir->Get(chan+boson+"_"+hist+"_"+signal); 
     hSig->SetNameTitle("Signal_"+chan+boson,"Signal in "+chan+boson); 
     int nRebin =(int)((double)hSig->GetXaxis()->GetNbins()/(double)nBin); 
     hSig->Rebin(nRebin);     
     for(int i = 0; i<nRebin; i++){
       if(hSig->GetBinContent(i) == 0)
	 hSig->SetBinContent(i,nMinInBin);
     }
   } 

   void setData(TString boson = "A"){ 
     TDirectory * dir = (TDirectory*)((TFile::Open("plotter_"+chan+".root"))->Get(chan+boson+"_"+hist)); 
     hData = (TH1F*)dir->Get(chan+boson+"_"+hist); 
     hData->SetNameTitle("Data_"+chan+boson,"Data in "+chan+boson); 
     int nRebin =(int)((double)hData->GetXaxis()->GetNbins()/(double)nBin); 
     hData->Rebin(nRebin);     
   } 
   void setHistograms(){ 

     this->setBkg(); 
     this->setSig(); 
     this->setData(); 

     TFile* fBkg = new TFile("Background_"+chan+boson+".root","recreate"); 
     fBkg->cd(); 
     hBkg->Write();     
     hBkgCorr->Write();
     fBkg->Close(); 

     TFile* fSig = new TFile("Signal_"+chan+boson+".root","recreate"); 
     fSig->cd(); 
     hSig->Write(); 
     fSig->Close(); 

    TFile* fData = new TFile("Data_"+chan+boson+".root","recreate");
    fData->cd();
    hData->Write();
    fData->Close();
  }
  TH1F * correctBackground(TH1F * hin){
    TFile * f = TFile::Open("tf_plotter.root");
    TH1F * tf = (TH1F*) f->Get(chan+"MM_"+hist+"/"+chan+"_"+hist+"_nlo2lo"); 
    TH1F * hout = (TH1F*)hin->Clone(hin->GetName()+TString("_NLOcorr"));
    for(int i = 0; i< hin->GetXaxis()->GetNbins(); i++){
      float binVal = hin->GetBinCenter(i+1);
      int iBin     = tf->GetXaxis()->FindBin(binVal);
      float cf     = tf->GetBinContent(iBin);
      hout->SetBinContent(i+1,cf* hin->GetBinContent(i+1));
      hout->SetBinError(i+1,cf* hin->GetBinError(i+1));
    }
    return hout;
  }
  RooDataHist * getDataDH(RooRealVar * var){
    return  new RooDataHist("Data"+chan+boson,"Data"+chan+boson,*var,Import(*this->hData));
  }

  RooDataHist * getBkgDH(RooRealVar * var){
    return new RooDataHist("Background"+chan+boson,"Background"+chan+boson,*var,Import(*this->hBkg));
  }
  RooDataHist * getBkgDHNLO(RooRealVar * var){
    return new RooDataHist("BackgroundNLO"+chan+boson,"BackgroundNLO"+chan+boson,*var,Import(*this->hBkgCorr));
  }
  RooDataHist * getSigDH(RooRealVar * var){
    return new RooDataHist("Signal"+chan+boson,"Signal"+chan+boson,*var,Import(*this->hSig));
  }

  RooRealVar * getSigHistNorm(){
    stringstream s;
    s << hSig->Integral();
    return new RooRealVar("Signal"+chan+boson+"_norm", "",atof(s.str().c_str()));
  }
  RooRealVar * getBkgHistNorm(){
    stringstream s;
    s << hBkg->Integral();
    return new RooRealVar("Background"+chan+boson+"_norm", "",atof(s.str().c_str()));
  }
  RooRealVar * getBkgNLOHistNorm(){
    stringstream s;
    s << hBkgCorr->Integral();
    return new RooRealVar("BackgroundNLO"+chan+boson+"_norm", "",atof(s.str().c_str()));
  }
  TH1F * getModelHist(int id = 0){ //0: signal, 1:background, 2:backgroung at NLO
    if (id == 0) return hSig;
    if (id == 1) return hBkg;
    return hBkgCorr;
  }

  std::vector<RooRealVar*> getModelBinsCR(int id = 0){
    if(id == 0) return binsSigCR;
    return binsBkgCR;
  }

  std::vector<RooFormulaVar*> getModelBinsSR(int id = 0){
    if(id == 0) return binsSigSR;
    return binsBkgSR;
  }

  void addModelBinsCR(RooRealVar * var, int id = 0){
    if(id == 0)       binsSigCR.push_back(var);
    else if(id == 1)  binsBkgCR.push_back(var);
  }

  void addModelBinsSR(RooFormulaVar * var, int id = 0){
    if(id == 0)       binsSigSR.push_back(var);
    else if(id == 1)  binsBkgSR.push_back(var);
  }

  void setPH(RooParametricHist * h, int id = 0){
    if(id == 0) sigPH = h;
    if(id == 1) bkgPH = h;
  }

  void setPH_norm(RooAddition * a, int id = 0){
    if(id == 0) sigPH_norm = a;
    if(id == 1) bkgPH_norm = a;
  }

  TString chan, boson, hist;
  int nBin;
  bool isSR;
  TH1F * hSig, * hBkg, * hBkgCorr, * hData;
  RooParametricHist * sigPH, * bkgPH;
  RooAddition * sigPH_norm, * bkgPH_norm;
  std::vector<RooRealVar*>    binsSigCR;
  std::vector<RooFormulaVar*> binsSigSR;
  std::vector<RooRealVar*>    binsBkgCR;
  std::vector<RooFormulaVar*> binsBkgSR;
};

////////////////////////////////////////////
// Simple class to make transfer function //
// in a combine-friendly format           //
////////////////////////////////////////////

class TF{
 public:
 TF(VbfFitRegion * SR, VbfFitRegion * CR, std::vector<std::pair<TString,double> > sNui = {}, std::vector<std::pair<TString,double> > bNui = {}):sr(SR),cr(CR){
    if((SR->chan != CR->chan) || (SR->boson == CR->boson)){
      cout<< "Transfer function is made between the ragions from the same channel and different boson types;"<<endl;
      cout<< "Cannot transfer from "<<CR->chan<<CR->boson<<" to "<<SR->chan<<SR->boson<<"!!!"<<endl;
      return;
    }
    name = "TF_"+CR->chan;
    for(unsigned int i = 0; i<sNui.size(); i++) {
      RooRealVar * tmp = new RooRealVar(name+"_"+sNui[i].first+"_Sig",sNui[i].first+" nuisance parameter of signal",0);
      sigNuisances.push_back(make_pair(tmp,1+fabs(sNui[i].second)));
    }
    for(unsigned int i = 0; i<bNui.size(); i++) {
      RooRealVar * tmp = new RooRealVar(name+"_"+bNui[i].first+"_Bkg",bNui[i].first+" nuisance parameter of background",0);
      bkgNuisances.push_back(make_pair(tmp,1+fabs(bNui[i].second)));
    }
    nBin = sr->nBin;
    if(sr->nBin != cr->nBin) {
      cout<<"Provide histograms with the same binning in SR and CR!"<<endl;
      return;
    }
  };
  ~TF(){};

  // Based on current MC, we can set the TF as constant
  void creatTFHists(bool shapeOnly = false){   
    bkgTF = (TH1F*)sr->hBkg->Clone("BackgroundTF_"+sr->chan);
    bkgTF->SetTitle("BackgroundTF");
    bkgTF->Sumw2();

    TH1F * crBkg = (TH1F*)cr->hBkg->Clone("NormalCRBkg_"+sr->chan);
    if(crBkg->Integral() != 0 )
      crBkg->Scale(1./crBkg->Integral());

    if(shapeOnly){
      if(bkgTF->Integral() != 0)
	bkgTF->Scale(1./bkgTF->Integral());
      bkgTF->Divide(crBkg);
    } else {
      bkgTF->Divide(cr->hBkg);
    }

    bkgETF = (TH1F*)bkgTF->Clone(TString("Err_")+bkgTF->GetName());
    for(int i = 0; i < nBin+1; i++){
      cout << "Background Bin "<<i+1<<": "<<bkgTF->GetBinContent(i)<<" +/- "<<bkgTF->GetBinError(i)<<endl;
      double err = bkgTF->GetBinError(i);
      if(bkgTF->GetBinContent(i) < 0) err = -999;
      bkgETF->SetBinContent(i,err);
      if(bkgTF->GetBinContent(i) == 0)
	err = 0;
      else
	err = fabs(err/bkgTF->GetBinContent(i));
      bkgETF->SetBinError(i,err);
    }
    

    sigTF = (TH1F*)sr->hSig->Clone("SignalTF_"+sr->chan);
    sigTF->SetTitle("SignalTF");
    sigTF->Sumw2();
    TH1F * crSig = (TH1F*)cr->hSig->Clone("NormalCRSig_"+sr->chan);
    if(crSig->Integral() != 0 )
      crSig->Scale(1./crSig->Integral());

    if(shapeOnly){
      if(sigTF->Integral() != 0)
	sigTF->Scale(1./sigTF->Integral());
      sigTF->Divide(crSig);
    } else {
      sigTF->Divide(cr->hSig);
    }

    sigETF = (TH1F*)sigTF->Clone(TString("Err_")+sigTF->GetName());
    for(int i = 0; i < nBin+1; i++){
      cout << "Signal Bin "<<i+1<<": "<<sigTF->GetBinContent(i)<<" +/- "<<sigTF->GetBinError(i)<<endl;
      double err = sigTF->GetBinError(i);
      if(sigTF->GetBinContent(i) < 0) err = -999;;
      sigETF->SetBinContent(i,err);
      if(sigTF->GetBinContent(i) == 0)
	err = 0;
      else
	err = fabs(err/sigTF->GetBinContent(i));
      sigETF->SetBinError(i,err);
    }
    TFile* fTF = new TFile("TransferFactors_"+sr->chan+".root","recreate");
    fTF->cd();
    sigTF->Write();
    sigETF->Write();
    bkgTF->Write();
    bkgETF->Write();
    fTF->Close();    
  }

  void setTFs(){
    this->bkgTForm = this->computeTF(1);
    this->sigTForm = this->computeTF(0);
  }

  
  RooFormulaVar * getTFormula(int id = 0){
    if (id == 0) return sigTForm;
    return bkgTForm;
  }
 private:
  RooFormulaVar * computeTF(int id = 0){ //0: signal 1: background
    TH1F * Hist = NULL;
    TString tfName = "";
    TString prefix = "tmp_";
    if (id == 0) {
      Hist = (TH1F*)sigTF->Clone(prefix+sigTF->GetName());
      tfName = "Signal";
    } else if (id == 1){
      Hist = (TH1F*)bkgTF->Clone(prefix+bkgTF->GetName());
      tfName = "Background";
    } else {
      cout << "Bad TF ID!!"<<endl;
      return 0;
    }
    Hist->Fit("pol0");
    TF1 * f = Hist->GetFunction("pol0");
    RooRealVar * statErr = new RooRealVar("Stat"+name+tfName,"Stat"+name+tfName,0);
    stringstream formule;
    formule.str("");
    formule << f->GetParameter(0) << "*TMath::Power("<<1+(f->GetParError(0)/ f->GetParameter(0))<<",@0)";
    RooArgList * args = new RooArgList();
    args->add(*statErr);
    if (id == 0){
      for(unsigned int iNui = 0; iNui < sigNuisances.size(); iNui++){
	formule << "*TMath::Power("<<sigNuisances[iNui].second<<",@"<<iNui+1<<")";
	args->add(*sigNuisances[iNui].first);
      }
    } else {
      for(unsigned int iNui = 0; iNui < bkgNuisances.size(); iNui++){
	formule << "*TMath::Power("<<bkgNuisances[iNui].second<<",@"<<iNui+1<<")";
	args->add(*bkgNuisances[iNui].first);
      }
    }
    return new RooFormulaVar(name+tfName,"Trasnfer factor for "+name+" in "+tfName,formule.str().c_str(),*args);
  }

  TString name;
  VbfFitRegion * sr, * cr;
  TH1F * bkgTF, * sigTF, * bkgETF, *sigETF; 
  RooFormulaVar * bkgTForm, * sigTForm;
  std::vector<std::pair<RooRealVar*,double> > bkgNuisances, sigNuisances;
  int nBin;
};

#endif
