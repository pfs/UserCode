
//////////////////////////////////////
// Currently the simplest case:     //
//    - Only signal region          //
//    - Hence no transfer factor    //
//////////////////////////////////////

#ifndef WSPROVIDER_H
#define WSPROVIDER_H
#include <iostream>
#include <sstream>
#include <cstdlib>
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooHistPdf.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooBernstein.h"
#include "RooSimultaneous.h"
#include "RooFitResult.h"
#include "RooCategory.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include <vector>
#include <sstream>
#include <iostream>
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooFitResult.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TString.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TColor.h"
#include "TColorWheel.h"
#include "TColorGradient.h"
#include "RooVoigtian.h"
#include "RooCBShape.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooRandom.h"
#include <algorithm>
#include "RooWorkspace.h"
#include <fstream>
#include "RooCategory.h"
#include "RooWorkspace.h"
#include "RooRandom.h"
#include <iostream>
#include <sstream>
#include <vector>
#include "TFile.h"
#include "TArrow.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TString.h"
#include "TStyle.h"
#include "TChain.h"
#include "TH2.h"
#include "TH1.h"
#include "TF1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"
#include "TMath.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TVirtualFitter.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooPolynomial.h"
#include "RooExponential.h"
#include "RooArgList.h"
#include "RooGenericPdf.h"
#include "RooExtendPdf.h"
#include "RooFitResult.h"
#include "RooChebychev.h"
#include "RooDataHist.h"
#include "RooVoigtian.h"
#include "RooPlot.h"
#include "RooArgusBG.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooParametricHist.h"
#include <stdio.h>
#include <math.h>
#include <map>
#include "VbfFitRegion.h"
using namespace std;
using namespace RooFit;



//////////////////////////////////////////////////
// Simple class to make RooParametricHist       //
// All variables can be accessed via this class //
//////////////////////////////////////////////////

class WorkspaceProvider{
 public:
 WorkspaceProvider(TString histName, VbfFitRegion * sr, VbfFitRegion * cr,std::vector<std::pair<TString,double> > stfu = {}, std::vector<std::pair<TString,double> >btfu={}): hist(histName),chan(sr->chan){
    ws      = new RooWorkspace("ws"+chan, "ws"+chan);
    var     = new RooRealVar("var","var",-2,3); // to be tuned later
    ws->import(*var);
    for(unsigned int i = 0; i < stfu.size(); i++)
      sTFUnc.push_back(stfu[i]);
    for(unsigned int i = 0; i < btfu.size(); i++)
      bTFUnc.push_back(btfu[i]);
    SR = sr->clone();
    CR = cr->clone();
    cout <<"Set TF ----------"<<endl;
    this->setTF();
    cout<<"Create PH signal -"<<endl;
    this->createParamHists(0);
    cout<<"Create PH bkg ----"<<endl;
    this->createParamHists(1);

    wsNoTF      = new RooWorkspace("wsNoTF"+chan, "wsNoTF"+chan);
    wsNoTF->import(*var);
  }

  ~WorkspaceProvider(){}

  void setTF(){
    wsTF = new TF(SR, CR, sTFUnc, bTFUnc);
    wsTF->creatTFHists();
    wsTF->setTFs();
  }

  void createParamHists(int id = 0){//0: signal 1:background 
    RooArgList * argsMM = new RooArgList(); //CR
    RooArgList * argsA = new RooArgList();  //SR
    stringstream s;
    for (int i = 1; i < CR->getModelHist(id)->GetXaxis()->GetNbins()+1; i++){
      s.str("");
      s<<CR->getModelHist(id)->GetName()<<"_bin"<<i;
      double cont = CR->getModelHist(id)->GetBinContent(i);
      double err  = 3*cont;
      if(cont < 0){
	cont = 3;
	err  = 10;
      }
      CR->addModelBinsCR(new RooRealVar(s.str().c_str(),s.str().c_str(), cont, 0, err), id);
      argsMM->add(*CR->getModelBinsCR(id)[i-1]);

      s.str("");
      s<<SR->getModelHist(id)->GetName()<<"_bin"<<i;      
      SR->addModelBinsSR(new RooFormulaVar(s.str().c_str(),s.str().c_str(), "@0*@1", RooArgList(*wsTF->getTFormula(id),*CR->getModelBinsCR(id)[i-1])), id);
      argsA->add(*SR->getModelBinsSR(id)[i-1]);     
    }

    TString tmpName = "Signal";
    if(id == 1) tmpName = "Background";
    SR->setPH(new RooParametricHist(tmpName+"_PH_"+SR->chan+SR->boson, tmpName +" PDF in "+SR->boson+" region",*var,*argsA,*SR->hData),id);
    SR->setPH_norm(new RooAddition (tmpName+"_PH_"+SR->chan+SR->boson+"_norm","Total Number of "+tmpName+" events in "+SR->boson+" region",*argsA),id);
    CR->setPH(new RooParametricHist(tmpName+"_PH_"+CR->chan+CR->boson, tmpName +" PDF in "+CR->boson+" region",*var,*argsMM,*CR->hData),id);
    CR->setPH_norm(new RooAddition (tmpName+"_PH_"+CR->chan+CR->boson+"_norm","Total Number of "+tmpName+" events in "+CR->boson+" region",*argsMM),id);
  }

  void import(bool doSignalPH = false){

    // Signal region
    ws->import(*SR->getDataDH(var));
    ws->import(*SR->getBkgDH(var));
    ws->import(*SR->getSigDH(var));
    ws->import(*SR->getSigHistNorm());
    ws->import(*SR->getBkgHistNorm());
    for(unsigned int iTh = 0; iTh < SR->Theo->shapeSystBkg.size(); iTh++){
      TString thName = SR->Theo->shapeSystBkg[iTh]->GetName();
      RooDataHist * thD = new RooDataHist(thName+"DH",thName+"DH",*var,Import(*SR->Theo->shapeSystBkg[iTh]));
      RooHistPdf * thPdf = new RooHistPdf(thName,thName,*var,*thD);
      ws->import(*thD);
      ws->import(*thPdf);
    }
    for(unsigned int iExp = 0; iExp < SR->Exp->shapeSystBkg.size(); iExp++){
      TString exName = SR->Exp->shapeSystBkg[iExp]->GetName();
      RooDataHist * exD = new RooDataHist(exName+"DH",exName+"DH",*var,Import(*SR->Exp->shapeSystBkg[iExp]));
      RooHistPdf * exPdf = new RooHistPdf(exName,exName,*var,*exD);
      ws->import(*exD);
      ws->import(*exPdf);
    }
    for(unsigned int iTh = 0; iTh < SR->Theo->shapeSystSig.size(); iTh++){
      TString thName = SR->Theo->shapeSystSig[iTh]->GetName();
      RooDataHist * thD = new RooDataHist(thName,thName,*var,Import(*SR->Theo->shapeSystSig[iTh]));
      RooHistPdf * thPdf = new RooHistPdf(thName,thName,*var,*thD);
      ws->import(*thD);
      ws->import(*thPdf);
    }
    for(unsigned int iExp = 0; iExp < SR->Exp->shapeSystSig.size(); iExp++){
      TString exName = SR->Exp->shapeSystSig[iExp]->GetName();
      RooDataHist * exD = new RooDataHist(exName,exName,*var,Import(*SR->Exp->shapeSystSig[iExp]));
      RooHistPdf * exPdf = new RooHistPdf(exName,exName,*var,*exD);
      ws->import(*exD);
      ws->import(*exPdf);
    }
    if(doSignalPH){
      ws->import(*SR->sigPH);
      ws->import(*SR->sigPH_norm);
    }
    ws->import(*SR->bkgPH);
    ws->import(*SR->bkgPH_norm);


    cout<<"Start of the control region"<<endl;
    // Control region
    ws->import(*CR->getDataDH(var));
    ws->import(*CR->getBkgDH(var));
    ws->import(*CR->getSigDH(var));
    ws->import(*CR->getSigHistNorm());
    ws->import(*CR->getBkgHistNorm());
    /* for(unsigned int iTh = 0; iTh < CR->Theo->shapeSystBkg.size(); iTh++){ */
    /*   TString thName = CR->Theo->shapeSystBkg[iTh]->GetName(); */
    /*   RooDataHist * thD = new RooDataHist(thName,thName,*var,Import(*CR->Theo->shapeSystBkg[iTh])); */
    /*   ws->import(*thD); */
    /* } */
    /* for(unsigned int iExp = 0; iExp < CR->Exp->shapeSystBkg.size(); iExp++){ */
    /*   TString exName = CR->Exp->shapeSystBkg[iExp]->GetName(); */
    /*   RooDataHist * exD = new RooDataHist(exName,exName,*var,Import(*CR->Exp->shapeSystBkg[iExp])); */
    /*   ws->import(*exD); */
    /* } */
    if(doSignalPH){
      ws->import(*CR->sigPH);
      ws->import(*CR->sigPH_norm);
    }
    ws->import(*CR->bkgPH);
    ws->import(*CR->bkgPH_norm);

    
    ws->import(*wsTF->getTFormula(0));
    ws->import(*wsTF->getTFormula(1));

    TString opt = "SigPH";
    if(!doSignalPH) opt = "NoSigPH";
    TFile * fOut = new TFile("Channel_"+CR->chan+opt+CR->year+".root","recreate");
    fOut->cd();
    ws->Write();
    fOut->Save();
    fOut->Close();


    wsNoTF->import(*SR->getDataDH(var));
    wsNoTF->import(*SR->getBkgDHNLO(var));
    wsNoTF->import(*SR->getSigDH(var));
    RooHistPdf * bkgPdf = new RooHistPdf(SR->getBkgDHNLO(var)->GetName()+TString("Pdf"),"",*var,*SR->getBkgDHNLO(var));
    RooHistPdf * sigPdf = new RooHistPdf(SR->getSigDH(var)->GetName()+TString("Pdf"),"",*var,*SR->getSigDH(var));
    wsNoTF->import(*bkgPdf);
    wsNoTF->import(*sigPdf);
    wsNoTF->import(*SR->getSigHistNorm());
    wsNoTF->import(*SR->getBkgNLOHistNorm());
    for(unsigned int iTh = 0; iTh < SR->Theo->shapeSystBkg.size(); iTh++){
      TString thName = SR->Theo->shapeSystBkg[iTh]->GetName();
      RooDataHist * thD = new RooDataHist(thName+"DH",thName+"DH",*var,Import(*SR->Theo->shapeSystBkg[iTh]));
      RooHistPdf * thPdf = new RooHistPdf(thName,thName,*var,*thD);
      wsNoTF->import(*thD);
      wsNoTF->import(*thPdf);
    }
    for(unsigned int iExp = 0; iExp < SR->Exp->shapeSystBkg.size(); iExp++){
      TString exName = SR->Exp->shapeSystBkg[iExp]->GetName();
      RooDataHist * exD = new RooDataHist(exName+"DH",exName+"DH",*var,Import(*SR->Exp->shapeSystBkg[iExp]));
      RooHistPdf * exPdf = new RooHistPdf(exName,exName,*var,*exD);
      wsNoTF->import(*exD);
      wsNoTF->import(*exPdf);
    }
    for(unsigned int iTh = 0; iTh < SR->Theo->shapeSystSig.size(); iTh++){
      TString thName = SR->Theo->shapeSystSig[iTh]->GetName();
      if(thName.Contains("PDF"))
	cout <<thName <<" "<<SR->boson <<endl;
      RooDataHist * thD = new RooDataHist(thName,thName,*var,Import(*SR->Theo->shapeSystSig[iTh]));
      RooHistPdf * thPdf = new RooHistPdf(thName,thName,*var,*thD);
      wsNoTF->import(*thD);
      wsNoTF->import(*thPdf);
    }
    for(unsigned int iExp = 0; iExp < SR->Exp->shapeSystSig.size(); iExp++){
      TString exName = SR->Exp->shapeSystSig[iExp]->GetName();
      RooDataHist * exD = new RooDataHist(exName,exName,*var,Import(*SR->Exp->shapeSystSig[iExp]));
      RooHistPdf * exPdf = new RooHistPdf(exName,exName,*var,*exD);
      wsNoTF->import(*exD);
      wsNoTF->import(*exPdf);
    }
    TFile * fOut2 = new TFile("Channel_"+CR->chan+"_BkgNLO_"+CR->year+".root","recreate");
    fOut2->cd();
    wsNoTF->Write();
    fOut2->Save();
    fOut2->Close();
  }

  void makeCardNLO(YieldsErr YieldErrors, TString boson){
    TString outname = chan+"_"+boson+"_NLO"+SR->year+".txt";
    TString binName = "signal";
    ofstream myfile;
    myfile.setf(ios_base::fixed);
    myfile.precision(4);
    myfile.open(outname);
    myfile << "Datacard for Signal Region with Gamma+Jets corrected to NLO"<<endl;
    myfile << "imax *  number of categories" << endl;
    myfile << "jmax *  number of samples minus 1" << endl;
    myfile << "kmax *  number of nuisance parameters (sources of systematical uncertainties)" << endl;

    myfile << "\n------------" << endl;
    myfile << "shapes\tSignal\t" <<binName <<"\tSignal_"    <<chan<<boson<<SR->year<<".root\t"    <<SR->hSig->GetName()    <<  "\t$PROCESS_$SYSTEMATIC"<< endl; 
    myfile << "shapes\tBkg\t"    <<binName <<"\tBackground_"<<chan<<boson<<SR->year<<".root\t"    <<SR->hBkgCorr->GetName()<<  "\t$PROCESS_$SYSTEMATIC"<< endl; 
    myfile << "shapes\tdata_obs\t"<<binName<<"\tData_"<<chan<<boson<<SR->year<<".root\t"          <<SR->hData->GetName()   <<  "\t$PROCESS_$SYSTEMATIC"<< endl;
    myfile << "------------" << endl;
    myfile << "bin\t"<<binName << endl;
    myfile << "observation\t-1.0" << endl;
    myfile << "------------" << endl;

    myfile << "bin\t"<<binName<<"\t"<<binName<< endl;
    myfile << "process\tSignal\tBkg" << endl;
    myfile << "process\t0\t1" << endl;
    myfile << "rate\t"<<SR->hSig->Integral()<<"\t"<<SR->hBkgCorr->Integral()<< endl;

    myfile << "------------" << endl;
    myfile << "BkgTheoryRate\tlnN\t-\t1.3"<<endl;
    for (auto& x : YieldErrors) {
      if(x.second.first == x.second.second )
	myfile << x.first << "\tlnN\t" << x.second.first << "\t"<<x.second.first << endl;
      else
	myfile << x.first << "\tlnN\t" << x.second.first << "/" << x.second.second << "\t-" << endl;
    }
    stringstream line;
    for (auto& s : SR->Exp->yieldSystSig) {
      if(boson == "A" && s.first.Contains("mes")) continue;
      if(boson == "MM" && s.first.Contains("aes")) continue;
      line.str("");
      line << s.first << "\tlnN\t";
      if(s.second.first == s.second.second ) {
	if( s.second.second != 1)
	  line << s.second.first << "\t";
	else
	  line <<"-\t";
      } else line << s.second.first << "/" << s.second.second <<"\t" ;
      for (auto& b : SR->Exp->yieldSystBkg) {
	if(s.first != b.first) continue;
	if(b.second.first == b.second.second ) {
	  if( b.second.second != 1)
	    line << b.second.first << "\t";
	  else
	    line <<"-\t";
	} else line << b.second.first << "/" << b.second.second;
      }
      myfile << line.str() << endl;
    }
    for (auto& s : SR->Exp->shapeSystMap) {
      myfile << s.first << "\tshape";
      TString tmp(s.second.first == 1 ? "\t1" : "\t-");
      myfile << tmp;
      tmp = (s.second.second == 1 ? "\t1" : "\t-");
      myfile << tmp << endl;
    }
    for (auto& s : SR->Theo->shapeSystMap) {
      myfile << s.first << "\tshape";
      TString tmp(s.second.first == 1 ? "\t1" : "\t-");
      myfile << tmp;
      tmp = (s.second.second == 1 ? "\t1" : "\t-");
      myfile << tmp << endl;
    }
    myfile.close();
  }

  void makeCard(YieldsErr YieldErrors, TString boson, bool doSignalPH = false, double sigEff=1, double bkgEff=1){
    // No shape uncertainty yet!
    TString opt = "SigPH";
    if(!doSignalPH) opt= "NoSigPH";

    TString outname = chan+"_"+boson+"_"+opt+SR->year+".txt";
    ofstream myfile;
    myfile.setf(ios_base::fixed);
    myfile.precision(4);
    myfile.open(outname);
    TString binName = "signal";
    if(boson == "MM"){
      myfile << "Control Region Datacard -- control category" <<endl;
      binName = "Control";
    } else if(boson == "A")
      myfile << "Signal Region Datacard -- signal category"<<endl;
    myfile << "imax *  number of categories" << endl;
    myfile << "jmax *  number of samples minus 1" << endl;
    myfile << "kmax *  number of nuisance parameters (sources of systematical uncertainties)" << endl;

    myfile << "\n------------" << endl;
    if(boson == TString("A")){
      if(doSignalPH)
	myfile << "shapes\tSignal"<<boson<<"\t"  <<binName<<"\tChannel_"<<chan<<opt<<SR->year<<".root ws"<<chan<<":"<<SR->sigPH->GetName()   << endl; 
      else
	myfile << "shapes\tSignal"<<boson<<"\t"  <<binName<<"\tChannel_"<<chan<<opt<<SR->year<<".root ws"<<chan<<":"<<SR->getSigDH(var)->GetName() << endl; 
      myfile << "shapes\tBkg"<<boson<<"\t"     <<binName<<"\tChannel_"<<chan<<opt<<SR->year<<".root ws"<<chan<<":"<<SR->bkgPH->GetName()   << endl;
    } else if (boson == TString("MM")){
      if(doSignalPH)
	myfile << "shapes\tSignal"<<boson<<"\t"  <<binName<<"\tChannel_"<<chan<<opt<<CR->year<<".root ws"<<chan<<":"<<CR->sigPH->GetName()   << endl; 
      else
	myfile << "shapes\tSignal"<<boson<<"\t"  <<binName<<"\tChannel_"<<chan<<opt<<CR->year<<".root ws"<<chan<<":"<<CR->getSigDH(var)->GetName()   << endl; 
      myfile << "shapes\tBkg"<<boson<<"\t"     <<binName<<"\tChannel_"<<chan<<opt<<CR->year<<".root ws"<<chan<<":"<<CR->bkgPH->GetName()   << endl;
    }
    myfile << "shapes\tdata_obs\t"<<binName<<"\tChannel_"<<chan<<opt<<SR->year<<".root ws"<<chan<<":"<<"Data"<<chan<<boson << endl;
    myfile << "------------" << endl;
    myfile << "bin\t"<<binName << endl;
    myfile << "observation\t-1.0" << endl;
    myfile << "------------" << endl;

    myfile << "bin\t"<<binName<<"\t"<<binName<< endl;
    myfile << "process\tSignal"<<boson<<"\tBkg"<<boson << endl;
    myfile << "process\t0\t1" << endl;
    if(boson == TString("A")){
      myfile << "rate\t"<<SR->hSig->Integral()<<"\t"<<SR->hBkg->Integral()<< endl;
    } else if (boson == TString("MM")){
      myfile << "rate\t"<<CR->hSig->Integral()<<"\t"<<CR->hBkg->Integral()<< endl;
    }
    myfile << "------------" << endl;
    
    for (auto& x : YieldErrors) {
      if(x.second.first == x.second.second )
	myfile << x.first << "\tlnN\t" << x.second.first << "\t"<<x.second.first << endl;
      else
	myfile << x.first << "\tlnN\t" << x.second.first << "/" << x.second.second << "\t-" << endl;
    }
    if(boson == "MM"){
      myfile << "------------" << endl;
      myfile << "# free floating parameters, we do not need to declare them, but its a good idea to "<<endl;
      for(int i = 0; i < SR->nBin; i++){
	if(doSignalPH)
	  myfile << SR->getModelBinsSR(0)[i]->GetName()<<"\tflatParam "<<endl;
	myfile << SR->getModelBinsSR(1)[i]->GetName()<<"\tflatParam "<<endl;
      }
    } else {
      myfile << "StatTF_"<<chan<<"Background\tparam\t0\t1"<<endl;
      if(doSignalPH)
	myfile << "StatTF_"<<chan<<"Signal\tparam\t0\t1"<<endl;
    }
    myfile.close();
  }

  void plotSystSig(){
    plotSysts(SR);
  }
  void plotter(TH1F* nom, std::vector<TH1F*> hists, TString process, TDirectory * sigD, TString dName){
    for(unsigned int iSyst = 0; iSyst < hists.size(); iSyst++){
      sigD->cd();
      TString hName = hists[iSyst]->GetName();
      if(hName.Contains("Down")) continue;
      TString systName = hName(process.Length()+1,hName.Length()-(process.Length()+3));
      for(unsigned int iSyst2 = 0; iSyst2 < hists.size(); iSyst2++){
	if(TString(hists[iSyst2]->GetName()).Contains("Up")) continue;
	if(!TString(hists[iSyst2]->GetName()).Contains(systName)) continue;
	TCanvas * c = new TCanvas(process+"_"+systName, systName, 448,99,500,500);
	c->SetHighLightColor(2);
	c->Range(0,0,1,1);
	c->SetFillColor(0);
	c->SetBorderMode(0);
	c->SetBorderSize(2);
	c->SetTickx(1);
	c->SetTicky(1);
	c->SetLeftMargin(0.16);
	c->SetRightMargin(0.03);
	c->SetTopMargin(0.05);
	c->SetBottomMargin(0.13);
	c->SetFrameFillStyle(0);
	c->SetFrameBorderMode(0);
	TPad *pad1 = new TPad(TString("pad_")+c->GetName(), "",0,0.20,1,1);
	pad1->Draw();
	pad1->cd();
	pad1->Range(-140.4878,-1.348856,737.561,4.924722);
	pad1->SetFillColor(0);
	pad1->SetBorderMode(0);
	pad1->SetBorderSize(2);
	pad1->SetTickx(1);
	pad1->SetTicky(1);
	pad1->SetLeftMargin(0.16);
	pad1->SetRightMargin(0.03);
	pad1->SetTopMargin(0.05);
	pad1->SetBottomMargin(0.13);
	pad1->SetFrameFillStyle(0);
	pad1->SetFrameBorderMode(0);
	pad1->SetFrameFillStyle(0);
	pad1->SetFrameBorderMode(0);
	nom->SetLineColor(kBlack);
	nom->SetFillColor(0);
	nom->SetStats(0);
	nom->SetMarkerColor(kBlack);
	nom->SetMarkerStyle(20);
	nom->SetTitle("Nominal");
	nom->Draw();
	hists[iSyst]->SetLineColor(kRed);
	hists[iSyst]->SetFillColor(0);
	hists[iSyst]->SetMarkerColor(kRed);
	hists[iSyst]->SetMarkerStyle(20);
	hists[iSyst]->SetStats(0);
	hists[iSyst]->SetTitle(systName+" up");
	hists[iSyst]->Draw("sames");
	hists[iSyst2]->SetLineColor(kBlue);
	hists[iSyst2]->SetFillColor(0);
	hists[iSyst2]->SetStats(0);
	hists[iSyst2]->SetMarkerColor(kBlue);
	hists[iSyst2]->SetMarkerStyle(20);
	hists[iSyst2]->SetTitle(systName + " down");
	hists[iSyst2]->Draw("sames");
	pad1->BuildLegend();
	c->cd();
	TPad *   pad2 = new TPad(TString("pad2_")+c->GetName(), "",0,0,1,0.2);
	pad2->Draw();
	pad2->cd();
	pad2->Range(-195.1219,-61.64157,1024.39,129.679);
	pad2->SetFillColor(0);
	pad2->SetBorderMode(0);
	pad2->SetBorderSize(2);
	pad2->SetGridy();
	pad2->SetTicky(1);
	pad2->SetLeftMargin(0.16);
	pad2->SetRightMargin(0.03);
	pad2->SetTopMargin(0.05);
	pad2->SetBottomMargin(0.13);
	pad2->SetFrameFillStyle(0);
	pad2->SetFrameBorderMode(0);
	pad2->SetFrameFillStyle(0);
	pad2->SetFrameBorderMode(0);
	TH1F * tmpU = (TH1F*)hists[iSyst]->Clone();
	tmpU->Divide(nom);
	tmpU->GetYaxis()->SetRangeUser(0.8,1.2);
	tmpU->SetTitle("");
	tmpU->GetYaxis()->SetNdivisions(504,kFALSE);
	tmpU->GetYaxis()->SetLabelSize(0.12);
	tmpU->Draw();
	TH1F * tmpD = (TH1F*)hists[iSyst2]->Clone();
	tmpD->Divide(nom);	
	tmpD->Draw("sames");
	c->cd();
	sigD->cd();
	c->Write();
	c->SaveAs("plots/"+dName+"_"+process+"_"+systName+".pdf");
	c->SaveAs("plots/"+dName+"_"+process+"_"+systName+".png");
      }
    }
  }
  void plotSysts(VbfFitRegion * r){
    TString dName = r->chan+r->boson+r->year;
    const int dir_err = system("mkdir -p plots/"+r->chan+r->boson+r->year);
    TFile * f = new TFile("plots/"+r->chan+r->boson+r->year+".root","recreate");    

    TDirectory * sigD = f->mkdir("Signal"); 
    plotter(r->hSig, r->Exp->shapeSystSig,"Signal",sigD, dName);
    plotter(r->hSig, r->Exp->rateSystSig,"Signal",sigD, dName);
    plotter(r->hSig, r->Theo->shapeSystSig,"Signal",sigD, dName);
    f->cd();

    TDirectory * bkgD = f->mkdir("Background"); 
    plotter(r->hBkg, r->Exp->shapeSystBkg,"Bkg",bkgD, dName);
    plotter(r->hBkg, r->Exp->rateSystBkg,"Bkg",bkgD, dName);
    plotter(r->hBkg, r->Theo->shapeSystBkg,"Bkg",bkgD, dName);
    f->cd();

    f->Close();    
  }
 private:
  TString hist,chan;
  VbfFitRegion * SR, *CR;
  TF * wsTF;
  RooRealVar * var;
  RooWorkspace * ws, * wsNoTF;
  std::vector<std::pair<TString,double> > bTFUnc, sTFUnc;

};

YieldsErr splitter(const std::string &s, char delim1 , char delim2){
  YieldsErr ret;
  int nSyst = count(s.begin(),s.end(),delim1) + 1;
  string mystr = s;
  for(int i = 0; i < nSyst; i++){
    int pos1 = mystr.find(delim1);
    string syst = mystr.substr(0,pos1-1);
    string name = syst.substr(0,syst.find(delim2));
    syst = syst.substr(syst.find(delim2)+1);
    double down = (double)atof(syst.substr(0,syst.find(delim2)).c_str());
    syst = syst.substr(syst.find(delim2)+1);
    double up = (double)atof(syst.c_str());
    ret[name]=make_pair(down,up);
  }
  cout << ret.size()<<endl;
  return ret;
}
#endif
