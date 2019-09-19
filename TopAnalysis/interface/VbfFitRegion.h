
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
#include "VbfSystematics.h"
#include <sstream>
#include <iostream>
#include <string>
#include <vector>


/////////////////////////////////////
// A container for histograms, etc //
// of SR or CR                     //
/////////////////////////////////////

class VbfFitRegion{
 public:
 VbfFitRegion(TString channel, TString v, TString Hist, TString year_, int bin, bool SR, bool onlyShape=false, bool NLODefult_ = false):chan(channel),boson(v),hist(Hist),year(year_),nBin(bin),isSR(SR),NLODefault(NLODefult_){
    TFile * f_ = TFile::Open("tf_plotter_"+year+".root");
    tf_ = (TH1F*) f_->Get(chan+"MM_"+hist+"/"+chan+"_"+hist+"_Z_data2lo");//nlo2lo");  
    
    /* int nonZero(0); */
    /* for(int b = tf_->GetXaxis()->GetNbins()-1; b > 0; b--){ */
    /*   if(fabs(tf_->GetBinContent(b+1)) > 0.0001) { */
    /* 	nonZero = b+1; */
    /* 	break; */
    /*   } */
    /* } */
    /* bool blind(nonZero > int(0.5*tf_->GetXaxis()->GetNbins()) && nonZero < int(tf_->GetXaxis()->GetNbins()-2)); */
    /* if(blind){ */
    /*   cout <<"last value: "<<tf_->GetXaxis()->GetBinLowEdge(nonZero+1)<<endl; */
    /*   fit_ = new TF1("pol1","pol1",-1.0,tf_->GetXaxis()->GetBinLowEdge(nonZero+1)); */
    /*   tf_->Fit(fit_,"R"); */
    /*   fit_->SetRange(-1.,1.); */
    /*   cout<<fit_->Eval(0.9)<<endl; */
    /* } else tf_->Fit("pol1"); */
    tf_->Fit("pol1");
    fit_ = (TF1*)tf_->GetListOfFunctions()->FindObject("pol1");
    Exp = new systematics("exp", onlyShape, year);
    Theo = new systematics("th", onlyShape, year);
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
     double normSig = hSig->Integral();
     double normBkg = hBkg->Integral();
     Exp->setBkg(chan, boson, hist, nBin, tf_, fit_);
     Theo->setBkg(chan, boson, hist, nBin, tf_, fit_);
     Exp->setSig(chan, boson, hist, nBin);
     Theo->setSig(chan, boson, hist, nBin);
     Exp->setShapeSystematicsB(normBkg,hBkg,0);
     Exp->setShapeSystematicsB(normBkg,hBkgCorr,1);
     Exp->setShapeSystematicsB(normBkg,hBkgCorrLin[0],2);
     Exp->setShapeSystematicsS(normSig, hSig);
     Exp->setSystematicsMap();
     //     Exp->printMap();
     Theo->setShapeSystematicsB(normBkg, hBkg, 0);
     Theo->setShapeSystematicsB(normBkg, hBkgCorr, 1);
     Theo->setShapeSystematicsB(normBkg, hBkgCorrLin[0], 2);
     Theo->setShapeSystematicsS(normSig, hSig);
     Theo->setSystematicsMap();
     

     std::pair<TH1F*,TH1F*> updownI = convert1SDto2SD(hBkg,hBkgCorr);
     updownI.first->SetName("Bkg_NLODown");
     Theo->shapeSystBkg.push_back(updownI.first);
     updownI.second->SetName("Bkg_NLOUp");
     Theo->shapeSystBkg.push_back(updownI.second);
     Theo->shapeSystMap["NLO"]=make_pair(0,1);

     TH1F* tmp = (TH1F*)hBkgCorrLin[1]->Clone("BkgNLO_NLOLin"+chan+year+"Up");
     Theo->shapeSystBkgNLO.push_back(tmp);
     tmp = (TH1F*)hBkgCorrLin[2]->Clone("BkgNLO_NLOLin"+chan+year+"Down");
     Theo->shapeSystBkgNLO.push_back(tmp);
     Theo->shapeSystMap["NLOLin"+chan+year]=make_pair(0,1);
     std::pair<TH1F*,TH1F*> updown = convert1SDto2SD(hBkgCorr,hBkg);
     updown.first->SetName("BkgNLOBinned_NLOBinnedUp");
     Theo->shapeSystBkgNLOBinned.push_back(updown.first);
     updown.second->SetName("BkgNLOBinned_NLOBinnedDown");
     Theo->shapeSystBkgNLOBinned.push_back(updown.second);
     Theo->shapeSystMap["NLOBinned"]=make_pair(0,1);

     //     Theo->printMap();
     //Yields
     Exp->setYieldSystematicsB(normBkg);
     Exp->setYieldSystematicsS(normSig);
     Theo->setYieldSystematicsB(normBkg);
     Theo->setYieldSystematicsS(normSig);
     Theo->setPureAcceptanceScaleS();
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
     hBkgCorrLin.clear();
     hBkgCorrLin.push_back((TH1F*)hBkg->Clone("BackgroundNLOLin_"+chan+boson));
     hBkgCorrLin.push_back((TH1F*)hBkg->Clone("BackgroundNLOLin_"+chan+boson+"Up"));
     hBkgCorrLin.push_back((TH1F*)hBkg->Clone("BackgroundNLOLin_"+chan+boson+"Down"));

     int nRebin =(int)((double)hBkg->GetXaxis()->GetNbins()/(double)nBin);      
    
     for(int i = 1; i < nBackgrounds; i++){ 
       if (bkgs[i] == "QCD") continue;
       if (TString(bkgs[i]).Contains("Fake")) continue;
       TH1F * tmp = (TH1F*)dir->Get(chan+boson+"_"+hist+"_"+bkgs[i]); 
       if (tmp != NULL) {
	 tmp->SetLineColor(kBlack);
	 tmp->SetMarkerColor(kBlack);
	 hBkg->Add(tmp);
	 if(bkgs[i] == "#gamma+jets"){
	   std::vector<TH1F*> tmpCorrVec = correctBackground(tmp,tf_,fit_);
	   hBkgCorr->Add(tmpCorrVec[0]);
	   hBkgCorrLin[0]->Add(tmpCorrVec[1]);
	   hBkgCorrLin[1]->Add(tmpCorrVec[2]);
	   hBkgCorrLin[2]->Add(tmpCorrVec[3]);
	 } else {
	   hBkgCorr->Add(tmp);
	   hBkgCorrLin[0]->Add(tmp);
	   hBkgCorrLin[1]->Add(tmp);
	   hBkgCorrLin[2]->Add(tmp);
	 }
       }
     } 
     

     /* std::vector<TH1F*> tmpCorrVec = correctBackground(hBkg,tf_,fit_); */
     /* hBkgCorrLin.push_back((TH1F*)tmpCorrVec[1]->Clone("BackgroundNLOLin_"+chan+boson)); */
     /* hBkgCorrLin.push_back((TH1F*)tmpCorrVec[2]->Clone("BackgroundNLOLin_"+chan+boson+"Up")); */
     /* hBkgCorrLin.push_back((TH1F*)tmpCorrVec[3]->Clone("BackgroundNLOLin_"+chan+boson+"Down")); */
     /* hBkgCorr = (TH1F*)tmpCorrVec[0]->Clone("BackgroundNLO_"+chan+boson); */
     /* hBkgCorr->SetTitle("Background (NLO) in "+chan+boson); */
     hBkg->Rebin(nRebin);
     hBkgCorr->Rebin(nRebin);
     for(int i = 0; i<nRebin; i++){
       if(hBkg->GetBinContent(i) == 0)
	 hBkg->SetBinContent(i,nMinInBin);
       if(hBkgCorr->GetBinContent(i) == 0)
	 hBkgCorr->SetBinContent(i,nMinInBin);
     }
     hBkgCorr->Scale(hBkg->Integral()/hBkgCorr->Integral());
     for(unsigned int iH = 0; iH<hBkgCorrLin.size(); iH++){
       hBkgCorrLin[iH]->Rebin(nRebin);
       hBkgCorrLin[iH]->Scale(hBkg->Integral()/hBkgCorrLin[iH]->Integral());
     }

     hFake = (TH1F*)dir->Get(chan+boson+"_"+hist+"_"+bkgs[5]);
     if (hFake != NULL){
       hFake->SetNameTitle("Fake_"+chan+boson,"Fake in "+chan+boson);
       hFake->SetLineColor(kBlack);
       hFake->SetMarkerColor(kBlack);
       hFake->Rebin(nRebin);
       for(int i = 0; i<nRebin; i++){
	 if(hFake->GetBinContent(i) == 0)
	   hFake->SetBinContent(i,0.1*nMinInBin);
       }
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
     //     std::vector<TH1F*> tmpCorrVec;
     //for(unsigned int iH = 0; iH<hBkgCorrLin.size(); iH++) hBkgCorrLin[iH]->Write();
     hBkgCorrLin[0]->Write();
     if(hFake != NULL)
       hFake->Write();
     for(unsigned int i = 0; i < (Theo->shapeSystBkg).size(); i++) {
       Theo->shapeSystBkg[i]->Write();
     }
     for(unsigned int i = 0; i < (Theo->shapeSystBkgNLO).size(); i++) Theo->shapeSystBkgNLO[i]->Write();
     for(unsigned int i = 0; i < (Theo->shapeSystBkgNLOBinned).size(); i++) Theo->shapeSystBkgNLOBinned[i]->Write();
     for(unsigned int i = 0; i < (Exp->shapeSystBkg).size(); i++) {
       Exp->shapeSystBkg[i]->Write();
     }
     for(unsigned int i = 0; i < (Exp->shapeSystBkgNLO).size(); i++) Exp->shapeSystBkgNLO[i]->Write();
     for(unsigned int i = 0; i < (Exp->shapeSystBkgNLOBinned).size(); i++) Exp->shapeSystBkgNLOBinned[i]->Write();
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
  bool isSR, NLODefault;
  TH1F * hSig, * hBkg, * hBkgCorr, * hData, * hFake, * tf_;
  std::vector<TH1F*> hBkgCorrLin;
  RooParametricHist * sigPH, * bkgPH;
  RooAddition * sigPH_norm, * bkgPH_norm;
  std::vector<RooRealVar*>    binsSigCR;
  std::vector<RooFormulaVar*> binsSigSR;
  std::vector<RooRealVar*>    binsBkgCR;
  std::vector<RooFormulaVar*> binsBkgSR;
  systematics * Exp, * Theo, * ExpNlo, * TheoNlo;
  TF1 * fit_;
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
