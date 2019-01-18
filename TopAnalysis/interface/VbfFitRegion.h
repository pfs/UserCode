
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

typedef std::map<TString, std::pair<double, std::pair<double, double> > > MeanErr;
typedef std::map<TString, std::pair <double, double> > YieldsErr;
const int nBackgrounds = 5;
const TString bkgs[nBackgrounds]={"DY","#gamma+jets","QCD","W,Top,VV","#gamma#gamma"}; 
using namespace std;
using namespace RooFit;
const float nMinInBin = 0.1;

/////////////////////////////////
// A container for systematics //
/////////////////////////////////
class systematics{
 public:
 systematics(TString name_):name(name_){};
  ~systematics(){};
  
  YieldsErr getYieldSystematics(double nominal, TH2F* h){ 
    YieldsErr ret;
    if(name.Contains("exp")){
      cout << "Nominal: "<<nominal<<endl;
      for(int iSyst = 0; iSyst < h->GetYaxis()->GetNbins(); iSyst+=2){
	TString label = h->GetYaxis()->GetBinLabel(iSyst+1);
	TString syst = label(0,label.Sizeof()-3);
	double yield = (h->ProjectionX(h->GetName()+label,iSyst+1,iSyst+1))->Integral();
	cout << label<<": "<<yield<<endl;
	yield = 1+((yield - nominal)/nominal);
	std::pair<double,double> var;
	if(label.EndsWith("up")) var.first = yield;
	else if(label.EndsWith("dn"))  var.second = yield;
	label = h->GetYaxis()->GetBinLabel(iSyst+2);
	yield = (h->ProjectionX(h->GetName()+label,iSyst+2,iSyst+2))->Integral();
	yield = 1+((yield - nominal)/nominal);
	if(label.EndsWith("up")) var.first = yield;
	else if(label.EndsWith("dn"))  var.second = yield;

	if ((var.second < 1. && var.first < 1.) || (var.second > 1. && var.first > 1.)) {
	  yield = std::max(var.first,var.second);
	  if (yield < 1) yield = 2 - yield;
	  var.first = yield;
	  var.second = yield;
	}
	if(!label.Contains("JEC") && !label.Contains("JER") && !label.Contains("prefire"))
	  ret[syst] = var;

	if(syst.Contains("trig")) ret[syst] = make_pair(1.03,1.03);
      }
    } else {
      cout << "Expect experimental uncertainties for pure yield variation" <<endl;
      throw std::exception();
    }
    return ret;
  }  

  
  void setYieldSystematicsB(double nominal){ 

    yieldSystBkg = this->getYieldSystematics(nominal,this->bkg);
  }  
  void setYieldSystematicsS(double nominal){ 
    yieldSystSig =  this->getYieldSystematics(nominal,this->sig);
  }  

  std::vector<TH1F*> getShapeSystematics(TH2F * h){
    std::vector<TH1F*> ret;
    TString hName = h->GetName();
    TString process( hName.Contains("Background_")? "Bkg":"Signal");
    if(process == "Signal") rateSystSig.clear();
    else rateSystBkg.clear();
    std::pair<int,int> hasSyst(0,0);
    int iPDF = 0;
    stringstream pdfName;
    for(int iSyst = 0; iSyst < h->GetYaxis()->GetNbins(); iSyst++){
      TString label = h->GetYaxis()->GetBinLabel(iSyst+1);
      pdfName.str("");
      pdfName << "PDF" << iPDF; 
      TString syst = ((label.EndsWith("dn") || label.EndsWith("up")) ? TString(label(0,label.Sizeof()-3)) : TString(pdfName.str().c_str()));
      TString postfix(label.EndsWith("dn")? "Down" : "Up");
      if(syst.Contains("PDF")) {
	postfix = "Up" ;
	iPDF++;
      }
      if(syst.Contains("PDF0")) continue;
      TH1D * H = h->ProjectionX(process+"_"+syst+postfix,iSyst+1,iSyst+1);
      TH1F * tmpH = new TH1F();
      H->Copy(*tmpH);
      if(name.Contains("exp") && !label.Contains("JEC") && !label.Contains("JER") && !label.Contains("prefire")){
	if(process == "Signal") rateSystSig.push_back(tmpH);
	if(process == "Bkg")    rateSystBkg.push_back(tmpH);
	continue;
      } else ret.push_back(tmpH);
      if(syst.Contains("PDF")){
	postfix = "Down";
	ret.push_back((TH1F*)tmpH->Clone(process+"_"+syst+postfix));
      }
      shapeSystMap[syst] = hasSyst;
    }
    return ret;
  }
  void setShapeSystematicsB(){
    if(bkg != NULL)
      shapeSystBkg = this->getShapeSystematics(bkg);
  }
  void setShapeSystematicsS(){
    cout <<"Signal >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" <<endl;
    if(sig != NULL)
      shapeSystSig = this->getShapeSystematics(sig);
  }
  
  void setSystematicsMap(){
    for (auto& b : shapeSystMap) {
      for (unsigned int i = 0; i < shapeSystBkg.size(); i++){
	if (TString(shapeSystBkg[i]->GetName()).Contains(b.first)) {
	  shapeSystMap[b.first].second = 1;
	}
      }
      for (unsigned int i = 0; i < shapeSystSig.size(); i++){
	if (TString(shapeSystSig[i]->GetName()).Contains(b.first)){
	  shapeSystMap[b.first].first = 1;
	}
      }
    }
  }
  void printMap(){
    cout << "Printing ........"<<endl;
    for (auto& b : shapeSystMap) {
	  cout << b.first<<" "<< shapeSystMap[b.first].first << " " <<shapeSystMap[b.first].second<<endl;      
    }
  }
  void setBkg(TString chan, TString boson, TString hist, int nBin, TString year){     
     TDirectory * dir = (TDirectory*)((TFile::Open("plotter_"+chan+year+".root"))->Get(chan+boson+"_"+hist+"_"+name)); 
     dir->ls();
     int initBkg = ((boson == "A" )? 1 : 0);
   
     bkg = (TH2F*)dir->Get(chan+boson+"_"+hist+"_"+name+"_"+bkgs[initBkg]); 
     // Temporary fix for 2016 ============
     if(bkg == NULL) return;
     bkg->SetNameTitle("Background_"+chan+boson+"_syst_"+name,"Background "+name+" syst. in "+chan+boson); 
     bkg->SetLineColor(kBlack);
     bkg->SetMarkerColor(kBlack);

     int nRebin =(int)((double)bkg->GetXaxis()->GetNbins()/(double)nBin);      
     for(int i = 1; i < nBackgrounds; i++){ 
       if(boson == "A") break;
       if (bkgs[i] == "QCD") continue;
       TH2F * tmp = (TH2F*)dir->Get(chan+boson+"_"+hist+"_"+name+"_"+bkgs[i]); 
       if (tmp != NULL) {
	 tmp->SetLineColor(kBlack);
	 tmp->SetMarkerColor(kBlack);
	 bkg->Add(tmp);
       }
     } 
     bkg->RebinX(nRebin);
     for(int i = 0; i<nRebin; i++){
       for(int j = 0;  j < bkg->GetYaxis()->GetNbins();j++){
	 if(bkg->GetBinContent(i,j) == 0)
	   bkg->SetBinContent(i,j,nMinInBin);
       }
     }
  }
  void setSig(TString chan, TString boson, TString hist, int nBin, TString year){ 
    TString signal = "EWK #gammajj"; 
    if(boson == "MM") 
      signal = "EWK Zjj"; 
    TDirectory * dir = (TDirectory*)((TFile::Open("plotter_"+chan+year+".root"))->Get(chan+boson+"_"+hist+"_"+name)); 
    sig = (TH2F*)dir->Get(chan+boson+"_"+hist+"_"+name+"_"+signal); 
    sig->SetNameTitle("Signal_"+chan+boson+"_syst_"+name,"Signal "+name+" syst. in "+chan+boson); 
    int nRebin =(int)((double)sig->GetXaxis()->GetNbins()/(double)nBin); 
    sig->RebinX(nRebin);     
    for(int i = 0; i<nRebin; i++){
      for(int  j = 0;  j < sig->GetYaxis()->GetNbins(); j++){
	if(sig->GetBinContent(i,j) == 0)
	  sig->SetBinContent(i,j,nMinInBin);
      }
    }
  }
  TString name;
  TH2F * sig, * bkg;
  std::vector<TH1F*> shapeSystSig, shapeSystBkg;
  std::vector<TH1F*> rateSystSig, rateSystBkg;
  YieldsErr yieldSystSig, yieldSystBkg;
  std::map<TString, std::pair<int, int> > shapeSystMap;
};



/////////////////////////////////////
// A container for histograms, etc //
// of SR or CR                     //
/////////////////////////////////////

class VbfFitRegion{
 public:
 VbfFitRegion(TString channel, TString v, TString Hist, TString year_, int bin, bool SR):chan(channel),boson(v),hist(Hist),year(year_),nBin(bin),isSR(SR){
    Exp = new systematics("exp");
    Theo = new systematics("th");
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
   
   void setSystematics(){
     Exp->setBkg(chan, boson, hist, nBin, year);
     Theo->setBkg(chan, boson, hist, nBin, year);
     Exp->setSig(chan, boson, hist, nBin, year);
     Theo->setSig(chan, boson, hist, nBin, year);
     Exp->setShapeSystematicsB();
     Exp->setShapeSystematicsS();
     Exp->setSystematicsMap();
     //     Exp->printMap();
     Theo->setShapeSystematicsB();
     Theo->setShapeSystematicsS();
     Theo->setSystematicsMap();
     //     Theo->printMap();
     //Yields
     double normSig = hSig->Integral();
     double normBkg = hBkg->Integral();
     Exp->setYieldSystematicsB(normBkg);
     Exp->setYieldSystematicsS(normSig);
   }

   void setBkg(){ 
     TDirectory * dir = (TDirectory*)((TFile::Open("plotter_"+chan+year+".root"))->Get(chan+boson+"_"+hist)); 
     dir->ls();
     hBkg = (TH1F*)dir->Get(chan+boson+"_"+hist+"_"+bkgs[0]); 
     hBkg->SetNameTitle("Background_"+chan+boson,"Background in "+chan+boson); 
     hBkg->SetLineColor(kBlack);
     hBkg->SetMarkerColor(kBlack);

     hBkgCorr = (TH1F*)hBkg->Clone("BackgroundNLO_"+chan+boson); 
     hBkgCorr->SetTitle("Background (NLO) in "+chan+boson); 

     int nRebin =(int)((double)hBkg->GetXaxis()->GetNbins()/(double)nBin);      
     for(int i = 1; i < nBackgrounds; i++){ 
       if (bkgs[i] == "QCD") continue;
       TH1F * tmp = (TH1F*)dir->Get(chan+boson+"_"+hist+"_"+bkgs[i]); 
       if (tmp != NULL) {
	 tmp->SetLineColor(kBlack);
	 tmp->SetMarkerColor(kBlack);
	 hBkg->Add(tmp);
	 /* if(i == 3){ */
	 /*   TH1F * tmpCorr = correctBackground(tmp); */
	 /*   hBkgCorr->Add(tmpCorr); */
	 /* } else  */
	 hBkgCorr->Add(tmp);
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
     TString signal = "EWK #gammajj"; 
     if(boson == "MM") 
       signal = "EWK Zjj"; 
     TDirectory * dir = (TDirectory*)((TFile::Open("plotter_"+chan+year+".root"))->Get(chan+boson+"_"+hist)); 
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
     TDirectory * dir = (TDirectory*)((TFile::Open("plotter_"+chan+year+".root"))->Get(chan+boson+"_"+hist)); 
     hData = (TH1F*)dir->Get(chan+boson+"_"+hist); 
     hData->SetNameTitle("Data_"+chan+boson,"Data in "+chan+boson); 
     int nRebin =(int)((double)hData->GetXaxis()->GetNbins()/(double)nBin); 
     hData->Rebin(nRebin);     
   } 
   void setHistograms(){ 

     this->setBkg(); 
     this->setSig(); 
     this->setData(); 
     this->setSystematics();
     TFile* fBkg = new TFile("Background_"+chan+boson+year+".root","recreate"); 
     fBkg->cd(); 
     hBkg->Write();     
     hBkgCorr->Write();
     for(unsigned int i = 0; i < (Theo->shapeSystBkg).size(); i++) Theo->shapeSystBkg[i]->Write();
     for(unsigned int i = 0; i < (Exp->shapeSystBkg).size(); i++) Exp->shapeSystBkg[i]->Write();
     fBkg->Close(); 

     TFile* fSig = new TFile("Signal_"+chan+boson+year+".root","recreate"); 
     fSig->cd(); 
     hSig->Write(); 
     for(unsigned int i = 0; i < (Theo->shapeSystSig).size(); i++) Theo->shapeSystSig[i]->Write();
     for(unsigned int i = 0; i < (Exp->shapeSystSig).size(); i++) Exp->shapeSystSig[i]->Write();
     fSig->Close(); 

    TFile* fData = new TFile("Data_"+chan+boson+year+".root","recreate");
    fData->cd();
    hData->Write();
    fData->Close();
  }

  TH1F * correctBackground(TH1F * hin){
    TFile * f = TFile::Open("tf_plotter_"+year+".root");
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

  TString chan, boson, hist, year;
  int nBin;
  bool isSR;
  TH1F * hSig, * hBkg, * hBkgCorr, * hData;
  RooParametricHist * sigPH, * bkgPH;
  RooAddition * sigPH_norm, * bkgPH_norm;
  std::vector<RooRealVar*>    binsSigCR;
  std::vector<RooFormulaVar*> binsSigSR;
  std::vector<RooRealVar*>    binsBkgCR;
  std::vector<RooFormulaVar*> binsBkgSR;
  systematics * Exp, * Theo;
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
      double err = sigTF->GetBinError(i);
      if(sigTF->GetBinContent(i) < 0) err = -999;;
      sigETF->SetBinContent(i,err);
      if(sigTF->GetBinContent(i) == 0)
	err = 0;
      else
	err = fabs(err/sigTF->GetBinContent(i));
      sigETF->SetBinError(i,err);
    }
    TFile* fTF = new TFile("TransferFactors_"+sr->chan+"_"+sr->year+".root","recreate");
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
