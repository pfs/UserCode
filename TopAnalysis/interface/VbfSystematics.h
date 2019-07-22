
#ifndef VBFSYST_H
#define VBFSYST_H

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
const int nBackgrounds = 6;
const TString bkgs[nBackgrounds]={"DY","#gamma+jets","QCD","W,Top,VV","#gamma#gamma","Fake #gamma"}; 
using namespace std;
using namespace RooFit;
const float nMinInBin = 0.1;


std::pair<TH1F*,TH1F*> convert1SDto2SD(TH1F* nominal, TH1F* systematic){
  std::string hname(systematic->GetName());
  std::string p1 = hname.substr(0,hname.find("Down"));
  if(p1 == hname) p1 = hname.substr(0,hname.find("Up"));
  TString sysName = p1.c_str();//(hname.EndsWith("Up")? hname(0,hname.First("Up")) : (hname.EndsWith("Down")? hname(0,hname.First("Down")) : hname));
  TH1F * up = (TH1F*)nominal->Clone(sysName+"Up");
  TH1F * down = (TH1F*)nominal->Clone(sysName+"Down");
  for(int i = 0; i < nominal->GetXaxis()->GetNbins(); i++){
    double diff = 0.5*fabs(nominal->GetBinContent(i+1) - systematic->GetBinContent(i+1));
    up->SetBinContent(i+1, nominal->GetBinContent(i+1) + diff);
    down->SetBinContent(i+1, nominal->GetBinContent(i+1) - diff);    
  }
  up->Scale(nominal->Integral()/up->Integral());
  down->Scale(nominal->Integral()/down->Integral());
  return std::make_pair(up,down);
}

std::vector<TH1F*> correctBkgLinearFitNLO(TF1 * f, TH1F * hin){
  std::vector<TH1F*> ret;
  double slope = f->GetParameter(1);
  double slopeEr = f->GetParError(1);
  TH1F * h = (TH1F*)hin->Clone(hin->GetName()+TString("_NLOcorrLine"));
  h->Multiply(f);
  //h->Scale(hin->Integral()/h->Integral());
  ret.push_back(h);
  TF1 * fUp = (TF1*)f->Clone("UpFunc");
  fUp->SetParameter(1,slope+slopeEr);
  TH1F * hu = (TH1F*)hin->Clone(hin->GetName()+TString("_NLOcorrLineUp"));    
  hu->Multiply(fUp);
  delete fUp;
  //hu->Scale(hin->Integral()/hu->Integral());
  ret.push_back(hu);
  TF1 * fDown = (TF1*)f->Clone("DownFunc");
  fDown->SetParameter(1,slope-slopeEr);
  TH1F * hd = (TH1F*)hin->Clone(hin->GetName()+TString("_NLOcorrLineDown"));    
  hd->Multiply(fDown);
  //hd->Scale(hin->Integral()/hd->Integral());
  delete fDown;
  ret.push_back(hd);
  return ret;
}


std::vector<TH1F*> correctBackground(TH1F* hin, TH1F* tf, TF1 * fit, bool doLin = true){
  //cout<<fit->GetParameter(0)<<endl;
  std::vector<TH1F*> ret;
  TH1F * hout = (TH1F*)hin->Clone(hin->GetName()+TString("_NLOcorr"));
  for(int i = 0; i< hin->GetXaxis()->GetNbins(); i++){
    float binVal = hin->GetBinCenter(i+1);
    int iBin     = tf->GetXaxis()->FindBin(binVal);
    float cf     = tf->GetBinContent(iBin);
    hout->SetBinContent(i+1,cf* hin->GetBinContent(i+1));
    hout->SetBinError(i+1,cf* hin->GetBinError(i+1));
  } 
  ret.push_back(hout);
  if(doLin){
    std::vector<TH1F*> lin = correctBkgLinearFitNLO(fit,hin);
    ret.insert(ret.end(), lin.begin(), lin.end());
  }
  return ret;
}

std::vector<TH2F*> correctBkg2D(TH2F * hin, TH1F * tf, TF1 * fit){
  TH2F * ret = (TH2F*)hin->Clone(TString("tmp2D")+hin->GetName());
  TH2F * ret2 = (TH2F*)hin->Clone(TString("tmp2DLin")+hin->GetName());
  stringstream s("");
  for(int y = 0; y < hin->GetYaxis()->GetNbins(); y++){
    s.str("");
    s << "tmp1D" << hin->GetName() << y+1;
    TH1F * tmp = (TH1F*)hin->ProjectionX(s.str().c_str(),y+1,y+1);
    std::vector<TH1F*> rw = correctBackground(tmp, tf, fit);
    for(int x = 0; x < rw[0]->GetXaxis()->GetNbins(); x++){
      ret->SetBinContent(x+1,y+1,rw[0]->GetBinContent(x+1));
      ret2->SetBinContent(x+1,y+1,rw[1]->GetBinContent(x+1));
    }
  }
  std::vector<TH2F*> out;
  out.push_back(ret);
  out.push_back(ret2);
  return out;
}


/////////////////////////////////
// A container for systematics //
/////////////////////////////////
class systematics{
 public:
 systematics(TString name_, bool onlyshape = false, TString year_="2017"):name(name_),shapeOnly(onlyshape), year(year_){
    NLO = false;
    this->setJECCorr(); 
    TString hashtag = "";
    if(year == "2017") hashtag = "3129835";
    if(year == "2016") hashtag = "0c522df";
    TFile * genwgt = TFile::Open("$CMSSW_BASE/src/TopLJets2015/TopAnalysis/data/era"+year+"/genweights_"+hashtag+".root");
    this->normH = (TH1F*)genwgt->Get("MC13TeV_"+year+"_EWKAJJ");
  };
  ~systematics(){};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Correlation map for JEC                                                                                 //
// https://docs.google.com/spreadsheets/d/1JZfk78_9SD225bcUuTWVo4i02vwI5FfeVKH-dwzUdhM/edit#gid=1345121349 //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
  int getSystBin(TString syst){
    TString binLabel = "";
    if(syst.Contains("muRmuF")){
      if(syst.Contains("up"))
	binLabel = "mur2muf2";
      else binLabel = "mur0.5muf0.5";
    } else if(syst.Contains("muR")){
      if(syst.Contains("up"))
	binLabel = "mur2muf1";
      else binLabel = "mur0.5muf1";
    } else if(syst.Contains("muF")){
      if(syst.Contains("up"))
	binLabel = "mur1muf2";
      else binLabel = "mur1muf0.5";
    }
    int ret = -1;
    for(int i = 0; i<12;i++){
      if(TString(normH->GetXaxis()->GetBinLabel(i+1)).Contains(binLabel)){
	ret = i+1;
	break;
      }
    }
    return ret;
  }
  std::map<TString,float> JECcorrMap;
  void setJECCorr(){
    JECcorrMap["AbsoluteMPFBias"] = 1 ;
    JECcorrMap["AbsoluteScale"] = 1;
    JECcorrMap["AbsoluteStat"] = 0;
    JECcorrMap["FlavorQCD"] = 1;
    JECcorrMap["Fragmentation"] = 1;
    JECcorrMap["PileUpDataMC"] = 1 ;
    JECcorrMap["PileUpPtBB"] = 1;
    JECcorrMap["PileUpPtEC1"] = 1;
    JECcorrMap["PileUpPtEC2"] = 1;
    JECcorrMap["PileUpPtHF"] = 1;
    JECcorrMap["PileUpPtRef"] = 1;
    JECcorrMap["RelativeFSR"] = 1;
    JECcorrMap["RelativeJEREC1"] = 0;
    JECcorrMap["RelativeJEREC2"] = 0;
    JECcorrMap["RelativeJERHF"] = 1;
    JECcorrMap["RelativePtBB"] = 1;
    JECcorrMap["RelativePtEC1"] = 0;
    JECcorrMap["RelativePtEC2"] = 0;
    JECcorrMap["RelativePtHF"] = 1;
    JECcorrMap["RelativeBal"] = 1;
    JECcorrMap["RelativeSample"] = 0;
    JECcorrMap["RelativeStatEC"] = 0 ;
    JECcorrMap["RelativeStatFSR"] = 0;
    JECcorrMap["RelativeStatHF"] = 0;
    JECcorrMap["SinglePionECAL"] = 1;
    JECcorrMap["SinglePionHCAL"] = 1;
    JECcorrMap["TimePtEta"] = 0;
    ////////////////////////////////////
    // Not in the table. Self-guessed //
    ////////////////////////////////////
    JECcorrMap["FlavorPureBottom"] = 1;
    JECcorrMap["FlavorPureCharm"] = 1;
    JECcorrMap["FlavorPureGluon"] = 1;
    JECcorrMap["FlavorPureQuark"] = 1;
  }
  inline void setNLO(){NLO = true;}
  YieldsErr getPureAcceptanceScale(double nominal, TH2F* h){ 
    YieldsErr ret;
    int iPDF = 0;
    stringstream pdfName;
    nominal = nominal*normH->GetBinContent(1);

    for(int iSyst = 0; iSyst < h->GetYaxis()->GetNbins(); iSyst++){
      TString label = h->GetYaxis()->GetBinLabel(iSyst+1);
      if(!(label.Contains("mu") || label.Contains("Mu") || label.Contains("MU"))) continue;
      pdfName.str("");
      pdfName << "PDF" << iPDF; 
      int nBinsyst = getSystBin(label);
      TString syst = ((label.EndsWith("dn") || label.EndsWith("up")) ? TString(label(0,label.Sizeof()-3)) : TString(pdfName.str().c_str()));
      double yield = (h->ProjectionX(h->GetName()+label,iSyst+1,iSyst+1))->Integral()*normH->GetBinContent(nBinsyst);   
      yield = 1.+((yield - nominal)/nominal);
      std::pair<double,double> var;
      if(label.EndsWith("up")) var.first = yield;
      else if(label.EndsWith("dn"))  var.second = yield;
      for(int jSyst = 0; jSyst < h->GetYaxis()->GetNbins(); jSyst++){
	if(jSyst == iSyst) continue;
	TString label_(h->GetYaxis()->GetBinLabel(jSyst+1));
	TString syst_ = ((label_.EndsWith("dn") || label_.EndsWith("up")) ? TString(label_(0,label_.Sizeof()-3)) : TString(pdfName.str().c_str()));
	if(syst_ != syst) continue;
	nBinsyst = getSystBin(label_);
	yield = (h->ProjectionX(h->GetName()+label_,jSyst+1,jSyst+1))->Integral()*normH->GetBinContent(nBinsyst);
	yield = 1.+((yield - nominal)/nominal);
	if(label.EndsWith("up")) var.second = yield;
	else if(label.EndsWith("dn"))  var.first = yield;
      }
      
      //Temporary
      if(syst.Contains("PDF")) continue;

      ret["rate"+syst] = var;
    }    
    return ret;
  }  
  
  YieldsErr getYieldSystematics(double nominal, TH2F* h){ 
    YieldsErr ret;
    //    if(name.Contains("exp")){
      cout << "Nominal: "<<nominal<<endl;
      int iPDF = 0;
      stringstream pdfName;
      
      for(int iSyst = 0; iSyst < h->GetYaxis()->GetNbins(); iSyst++){
	TString label = h->GetYaxis()->GetBinLabel(iSyst+1);
	if((label.Contains("mu") || label.Contains("Mu") || label.Contains("MU"))) continue;
	pdfName.str("");
	pdfName << "PDF" << iPDF; 
	TString syst = ((label.EndsWith("dn") || label.EndsWith("up")) ? TString(label(0,label.Sizeof()-3)) : TString(pdfName.str().c_str()));
	double yield = (h->ProjectionX(h->GetName()+label,iSyst+1,iSyst+1))->Integral();
	yield = 1.+((yield - nominal)/nominal);
	std::pair<double,double> var;
	if(label.EndsWith("up")) var.first = yield;
	else if(label.EndsWith("dn"))  var.second = yield;
	for(int jSyst = 0; jSyst < h->GetYaxis()->GetNbins(); jSyst++){
	  if(jSyst == iSyst) continue;
	  if(!TString(h->GetYaxis()->GetBinLabel(jSyst+1)).Contains(syst)) continue;
	  yield = (h->ProjectionX(h->GetName()+label,jSyst+1,jSyst+1))->Integral();
	  yield = 1.+((yield - nominal)/nominal);
	  if(label.EndsWith("up")) var.second = yield;
	  else if(label.EndsWith("dn"))  var.first = yield;
	}	
	if ((var.second < 1. && var.first < 1.) || (var.second > 1. && var.first > 1.)) {
	  yield = std::max(var.first,var.second);
	  float uncert = fabs(1 - yield);
	  var.first = fabs(1.0-uncert);
	  var.second = fabs(1.0+uncert);
	}
	if (name.Contains("exp")) {
	  if(label.Contains("JEC")){
	    if (JECcorrMap[syst(0,syst.Sizeof()-4)] == 0 ) syst = syst+year;
	  } else syst = syst+year;
	}
	
	//Temporary
	if(syst.Contains("PDF")) continue;
        if(syst.Contains("trig")) {
	  var.first = 0.97;
	  var.second = 1.03;
	}


	///////////////////////
	// Checks on numbers //
	///////////////////////
	//-- very small systs 
	bool rmSyst = false;
	rmSyst = (fabs(float(var.first - 1.)) < 0.001 && fabs(float(var.second - 1.)) < 0.001);
	rmSyst = (rmSyst || ((float(var.first) == float(1)) && (float(var.second) == float(0))));
	rmSyst = (rmSyst || ((float(var.first) == float(0)) && (float(var.second) == float(1))));
	rmSyst = (rmSyst || ((float(var.first) == float(var.second)) && (float(var.second) == float(0))));
	if(rmSyst){
	  var.first = 1.;
	  var.second = 1.;
	}
	// absolut 0 
        if(float(var.first) == float(0) || float(var.second) == float(0)){
	  float uncert = (float(var.first) != float(0) ? fabs(1-var.first) : fabs(1-var.second));
	  if(uncert < 0.005) {
	    var.first = 1.;
	    var.second = 1.;
	  }
	  var.first = fabs(1-uncert);
	  var.second = fabs(1+uncert);
	}
	
	if(shapeOnly){
	  ret["rate"+syst] = var;
	} else if(name.Contains("exp") && !label.Contains("JEC") && !label.Contains("JER") && !label.Contains("prefire") && !label.Contains("aes"))
	  ret[syst] = var;	
      }
    return ret;
  }  

  
  void setYieldSystematicsB(double nominal){ 
    yieldSystBkg = this->getYieldSystematics(nominal,this->bkg);
    for (auto& s : yieldSystBkg) {
      if (s.first == "ratemuF"){
	yieldSystBkg["ratemuFQCD"] = make_pair(yieldSystBkg["ratemuF"].first,yieldSystBkg["ratemuF"].second);
	yieldSystBkg["ratemuFEWK"] = make_pair(1.0,1.0);
      } else if (s.first == "ratemuR"){
	yieldSystBkg["ratemuRQCD"] = make_pair(yieldSystBkg["ratemuR"].first,yieldSystBkg["ratemuR"].second);
	yieldSystBkg["ratemuREWK"] = make_pair(1.0,1.0);
      }
    }
  }  
  void setYieldSystematicsS(double nominal){ 
    yieldSystSig =  this->getYieldSystematics(nominal,this->sig);
    for (auto& s : yieldSystSig) {
      if (s.first == "ratemuF"){
	yieldSystSig["ratemuFEWK"] = make_pair(yieldSystSig["ratemuF"].first,yieldSystSig["ratemuF"].second);
	yieldSystSig["ratemuFQCD"] = make_pair(1.0,1.0);
      } else if (s.first == "ratemuR"){
	yieldSystSig["ratemuREWK"] = make_pair(yieldSystSig["ratemuR"].first,yieldSystSig["ratemuR"].second);
	yieldSystSig["ratemuRQCD"] = make_pair(1.0,1.0);
      }
    }
  }  
  void setPureAcceptanceScaleS(){ 
    accScaleSig = this->getPureAcceptanceScale(sigAcc1D->Integral(),this->sigAcc);
    for (auto& s : accScaleSig) {
      if (s.first == "ratemuF"){
	accScaleSig["ratemuFEWK"] = make_pair(accScaleSig["ratemuF"].first,accScaleSig["ratemuF"].second);
	accScaleSig["ratemuFQCD"] = make_pair(1.0,1.0);
      } else if (s.first == "ratemuR"){
	accScaleSig["ratemuREWK"] = make_pair(accScaleSig["ratemuR"].first,accScaleSig["ratemuR"].second);
	accScaleSig["ratemuRQCD"] = make_pair(1.0,1.0);
      } else if (s.first == "ratemuRmuF"){
	accScaleSig["ratemuRmuFEWK"] = make_pair(accScaleSig["ratemuRmuF"].first,accScaleSig["ratemuRmuF"].second);
	accScaleSig["ratemuRmuFQCD"] = make_pair(1.0,1.0);
      }
    }
  }  

  std::vector<TH1F*> getShapeSystematics(TH2F * h, double nominal, TH1F* hNom = NULL){
    std::vector<TH1F*> ret;
    TString hName = h->GetName();
    TString process( hName.Contains("Background_")? "Bkg": 
		     (hName.Contains("BackgroundNLO_")? "BkgNLOBinned":
		      (hName.Contains("BackgroundNLOLin_")? "BkgNLO": "Signal")));



    if(process == "Signal") rateSystSig.clear();
    else rateSystBkg.clear();
    std::pair<int,int> hasSyst(0,0);
    int iPDF = 0;
    stringstream pdfName;
    for(int iSyst = 0; iSyst < h->GetYaxis()->GetNbins(); iSyst++){
      TString label = h->GetYaxis()->GetBinLabel(iSyst+1);
      pdfName.str("");
      pdfName << "PDF" << iPDF; 
      TString syst = ((label.EndsWith("dn") || label.EndsWith("up") || label.EndsWith("qg")) ? TString(label(0,label.Sizeof()-3)) : TString(pdfName.str().c_str()));
      if (label.EndsWith("qg")) syst = label;
      if (name.Contains("exp")) {
	if(label.Contains("JEC")){
	  if (JECcorrMap[syst(0,syst.Sizeof()-4)] == 0) syst = syst+year;
	} else if(!label.EndsWith("qg")) syst = syst+year;
      }
      TString postfix((label.EndsWith("dn") || label.EndsWith("qg"))? "Down" : "Up");
      
      ////////////////////////////////////////////////////////////
      // Decorrelating the Scale Unc. between QCD and EWK gjets //
      ////////////////////////////////////////////////////////////
      TString EwkQcdScale = "";
      if (syst.Contains("mu")) EwkQcdScale = ((process == "Signal") ? "EWK": "QCD");
      /////////////////////////////////////
      // Replicating PDF for Up and Down //
      /////////////////////////////////////
      if(syst.Contains("PDF")) {
	postfix = "Down" ;
	iPDF++;
      }
      if(syst.Contains("PDF0")) continue;
      TH1D * H = h->ProjectionX(process+"_"+syst+postfix,iSyst+1,iSyst+1);
      TH1F * tmpH = new TH1F();
      H->Copy(*tmpH);
      if(shapeOnly) {
	//tmpH->Sumw2();
	tmpH->Scale(nominal/tmpH->Integral());
	
	if(!label.EndsWith("qg") || hNom == NULL) ret.push_back(tmpH);
	else {
	  std::pair<TH1F*,TH1F*> UD = convert1SDto2SD(hNom, tmpH);
	  ret.push_back(UD.first);
	  ret.push_back(UD.second);
	}
	if (syst.Contains("mu")) ret.push_back((TH1F*)tmpH->Clone(process+"_"+syst+EwkQcdScale+postfix));
	//cout<<ret[ret.size()-1]->GetName()<<endl;
      } else {
	if(name.Contains("exp") && !label.Contains("JEC") && !label.Contains("JER") && !label.Contains("prefire") && !label.Contains("aes")){
	  if(process == "Signal") rateSystSig.push_back(tmpH);
	  if(process == "Bkg")    rateSystBkg.push_back(tmpH);
	  continue;
	} else {
	  ret.push_back(tmpH);
	  if (syst.Contains("mu")) ret.push_back((TH1F*)tmpH->Clone(process+"_"+syst+EwkQcdScale+postfix));
	}
      }
      //      if(syst.Contains("PDF") || syst.EndsWith("qg")){
      if(syst.Contains("PDF")){
	postfix = "Up";
	ret.push_back((TH1F*)tmpH->Clone(process+"_"+syst+postfix));
	//cout<<ret[ret.size()-1]->GetName()<<endl;
      }
      if(process == "Signal" || process == "Bkg"){ 
	shapeSystMap[syst] = hasSyst;
	if (syst.Contains("mu")) shapeSystMap[syst+EwkQcdScale] = hasSyst;
      }
    }
    return ret;
  }
  void setShapeSystematicsB(double nominal, TH1F* hNom = NULL, int NLO = 0){
    if(bkg != NULL && NLO==0)
      shapeSystBkg = this->getShapeSystematics(bkg, nominal, hNom);
    if(bkgCorr != NULL && NLO==1)
      shapeSystBkgNLOBinned = this->getShapeSystematics(bkgCorr, nominal, hNom);
    if(bkgCorrLin != NULL && NLO==2)
      shapeSystBkgNLO = this->getShapeSystematics(bkgCorrLin, nominal, hNom);
  }

  void setShapeSystematicsS(double nominal, TH1F* hNom = NULL){
    if(sig != NULL)
      shapeSystSig = this->getShapeSystematics(sig, nominal, hNom);
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
  void setBkg(TString chan, TString boson, TString hist, int nBin, TH1F * tf_, TF1 * fit_){     
     TDirectory * dir = (TDirectory*)((TFile::Open("plotter_"+chan+year+".root"))->Get(chan+boson+"_"+hist+"_"+name)); 
     int initBkg = 0; //((boson == "A" )? 1 : 0);
   
     bkg = (TH2F*)dir->Get(chan+boson+"_"+hist+"_"+name+"_"+bkgs[initBkg]); 
     // Temporary fix for 2016 ============
     if(bkg == NULL) return;
     bkg->SetNameTitle("Background_"+chan+boson+"_syst_"+name,"Background "+name+" syst. in "+chan+boson); 
     bkg->SetLineColor(kBlack);
     bkg->SetMarkerColor(kBlack);

     bkgCorr = (TH2F*)bkg->Clone("BackgroundNLO_"+chan+boson+"_syst_"+name);
     bkgCorr->SetTitle("Background (NLO) syst in "+chan+boson);
     bkgCorrLin = (TH2F*)bkg->Clone("BackgroundNLOLin_"+chan+boson+"_syst_"+name);

     int nRebin =(int)((double)bkg->GetXaxis()->GetNbins()/(double)nBin);      
     for(int i = 1; i < nBackgrounds; i++){ 
       //if(boson == "A") break;
       if (bkgs[i] == "QCD" ) continue;
       if (TString(bkgs[i]).Contains("Fake")) continue;
       TH2F * tmp = (TH2F*)dir->Get(chan+boson+"_"+hist+"_"+name+"_"+bkgs[i]); 
       if (tmp != NULL) {
	 tmp->SetLineColor(kBlack);
	 tmp->SetMarkerColor(kBlack);
	 bkg->Add(tmp);
	 if(bkgs[i] == "#gamma+jets"){
	   std::vector<TH2F*> tmpCorrVec = correctBkg2D(tmp,tf_,fit_);
           bkgCorr->Add(tmpCorrVec[0]);
           bkgCorrLin->Add(tmpCorrVec[1]);
         } else {
           bkgCorr->Add(tmp);
           bkgCorrLin->Add(tmp);
         }
       }
     } 
     bkg->RebinX(nRebin);
     bkgCorr->RebinX(nRebin);
     bkgCorrLin->RebinX(nRebin);
     for(int i = 0; i<nRebin; i++){
       for(int j = 0;  j < bkg->GetYaxis()->GetNbins();j++){
	 if(bkg->GetBinContent(i,j) <= 0)
	   bkg->SetBinContent(i,j,nMinInBin);
	 if(bkgCorr->GetBinContent(i,j) <= 0)
	   bkgCorr->SetBinContent(i,j,nMinInBin);
	 if(bkgCorrLin->GetBinContent(i,j) <= 0)
	   bkgCorrLin->SetBinContent(i,j,nMinInBin);
       }
     }
  }
  void setSig(TString chan, TString boson, TString hist, int nBin){ 
    TString signal = "EWK #gammajj"; 
    if(boson == "MM") 
      signal = "EWK Zjj"; 
    TDirectory * dir = (TDirectory*)((TFile::Open("plotter_"+chan+year+".root"))->Get(chan+boson+"_"+hist+"_"+name)); 
    sig = (TH2F*)dir->Get(chan+boson+"_"+hist+"_"+name+"_"+signal); 
    sig->SetNameTitle("Signal_"+chan+boson+"_syst_"+name,"Signal "+name+" syst. in "+chan+boson); 
    if(name == "th"){
      TDirectory * dir2 = (TDirectory*)((TFile::Open("plotterAcc_"+chan+year+".root"))->Get(chan+boson+"_"+hist+"Acc_"+name)); 
      dir2->ls();
      sigAcc = (TH2F*)dir2->Get(chan+boson+"_"+hist+"Acc_"+name+"_"+signal); 
      sigAcc->SetNameTitle("SignalAcc_"+chan+boson+"_syst_"+name,"Signal "+name+" syst. in "+chan+boson); 
      TDirectory * dir3 = (TDirectory*)((TFile::Open("plotterAcc_"+chan+year+".root"))->Get(chan+boson+"_"+hist+"Acc")); 
      sigAcc1D = (TH1F*)dir3->Get(chan+boson+"_"+hist+"Acc_"+signal); 
      sigAcc1D->SetNameTitle("SignalAcc_"+chan+boson+"_syst","Signal Acc. syst. in "+chan+boson); 
    }
    int nRebin =(int)((double)sig->GetXaxis()->GetNbins()/(double)nBin); 
    sig->RebinX(nRebin);
    if(name == "th"){
      sigAcc->RebinX(nRebin);
      sigAcc1D->Rebin(nRebin);
    }
    for(int i = 0; i<nRebin; i++){
      if(name == "th"){
	if(sigAcc1D->GetBinContent(i) <= 0){
	  sigAcc1D->SetBinContent(i,nMinInBin);
	}
      }
      for(int  j = 0;  j < sig->GetYaxis()->GetNbins(); j++){
	if(sig->GetBinContent(i,j) <= 0)
	  sig->SetBinContent(i,j,nMinInBin);
      }
      if(name == "th"){
	for(int  j = 0;  j < sigAcc->GetYaxis()->GetNbins(); j++){
	  if(sigAcc->GetBinContent(i,j) <= 0)
	    sigAcc->SetBinContent(i,j,nMinInBin);
	}
      }
    }
  }
  TString name; 
  bool shapeOnly, NLO;
  TString year;
  TH2F * sig, * bkg, * bkgCorr, * bkgCorrLin;
  TH2F * sigAcc;
  TH1F * normH, * sigAcc1D;
  std::vector<TH1F*> shapeSystSig, shapeSystBkg, shapeSystBkgNLOBinned, shapeSystBkgNLO;
  std::vector<TH1F*> rateSystSig, rateSystBkg;
  YieldsErr accScaleSig, yieldSystSig, yieldSystBkg;
  std::map<TString, std::pair<int, int> > shapeSystMap;
};


#endif
