#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TStyle.h"

#include<iostream>
#include <vector>

struct Acceptance_t
{
  float acceptance,stat,qcdScale,pdf,alphaS;
};

enum TTbarSample {AMCATNLOFXFX,AMCATNLOFXFX_HERWIG,M169v5,M172v5,M175v5,PSUP,PSDOWN};
Acceptance_t acceptanceEst(TTbarSample sample,Int_t njets,TString baseDir="/store/cmst3/user/psilva/LJets2015/64217e8");

Acceptance_t acceptanceEst(TTbarSample sample,Int_t njets,TString baseDir)
{
  float brCorr(1.0);
  std::vector<TString> files;
  if(sample==AMCATNLOFXFX)
    {
      brCorr=0.9883;
      for(Int_t i=0; i<15; i++) files.push_back(Form("root://eoscms//eos/cms/%s/MC13TeV_TTJets_amcatnloFXFX/MergedMiniEvents_%d.root",baseDir.Data(),i));
    }
  if(sample==AMCATNLOFXFX_HERWIG)
    {
      brCorr=0.9883;
      for(Int_t i=0; i<8; i++) files.push_back(Form("root://eoscms//eos/cms/%s/MC13TeV_TTJets_amcatnloFXFX_herwig/MergedMiniEvents_%d.root",baseDir.Data(),i));
    }
  if(sample==M169v5)
    for(Int_t i=0; i<5; i++) files.push_back(Form("root://eoscms//eos/cms/%s/MC13TeV_TTJets_m169v5/MergedMiniEvents_%d.root",baseDir.Data(),i));
  if(sample==M172v5)
    for(Int_t i=0; i<5; i++) files.push_back(Form("root://eoscms//eos/cms/%s/MC13TeV_TTJets/MergedMiniEvents_%d.root",baseDir.Data(),i));
  if(sample==PSUP)
    for(Int_t i=0; i<5; i++) files.push_back(Form("root://eoscms//eos/cms/%s/MC13TeV_TTJets_scaleup/MergedMiniEvents_%d.root",baseDir.Data(),i));
  if(sample==PSDOWN)
    for(Int_t i=0; i<5; i++) files.push_back(Form("root://eoscms//eos/cms/%s/MC13TeV_TTJets_scaledown/MergedMiniEvents_%d.root",baseDir.Data(),i));
  if(sample==M175v5)
    for(Int_t i=0; i<5; i++) files.push_back(Form("root://eoscms//eos/cms/%s/MC13TeV_TTJets_m175v5/MergedMiniEvents_%d.root",baseDir.Data(),i));
  
  //add results from all files
  TH1 *hcounter=0,*hnormcounter=0,*hsystlist=0;
  for(size_t i=0; i<files.size(); i++)
    {
      TFile *inF=TFile::Open(files[i]);
      if(inF==0) continue;
      if(hcounter==0) 
	{
	  hcounter = (TH1*)inF->Get(Form("analysis/fidcounter%d",njets))->Clone("hcounter");
	  hcounter->SetDirectory(0);
	  hnormcounter = (TH1*)inF->Get("analysis/fidcounter0")->Clone("hnormcounter");
	  hnormcounter->SetDirectory(0);
	  hsystlist = (TH1*)inF->Get("analysis/generator_initrwgt")->Clone("hsystlist");
	  hsystlist->SetDirectory(0);
	}
      else
	{
	  hcounter->Add((TH1*)inF->Get(Form("analysis/fidcounter%d",njets))->Clone("hcounter") );
	  hnormcounter->Add((TH1F*)inF->Get("analysis/fidcounter0"));
	  hsystlist->Add( (TH1*)inF->Get("analysis/generator_initrwgt")->Clone("hsystlist"));
	}
    }

  //report on final result
  float acceptance=hcounter->GetBinContent(2)/hnormcounter->GetBinContent(1)*brCorr;
  float acceptanceUnc=hcounter->GetBinError(2)/hnormcounter->GetBinContent(1)*brCorr;

  std::vector<float>qcdVars(6,0);
  TGraph *pdfVariations=new TGraph;
  std::vector<float> alphaSvars;
  for(int xbin=3; xbin<=hcounter->GetNbinsX(); xbin++)
    {
      float acceptanceVar=hcounter->GetBinContent(xbin)/hnormcounter->GetBinContent(xbin);
      if(acceptanceVar==0) continue;
      float deltaAcceptance=acceptanceVar-acceptance;
      TString varTitle=hsystlist->GetXaxis()->GetBinLabel(xbin);
      if(varTitle.Contains("PDF set = 260"))
	{
	  pdfVariations->SetPoint(pdfVariations->GetN(),deltaAcceptance,deltaAcceptance);
	}
      if(varTitle.Contains("PDF set = 265000") || varTitle.Contains("PDF set = 266000"))
	{
	  alphaSvars.push_back(deltaAcceptance);
	}
      else if (varTitle.Contains("hdamp=mt=172.5"))
	{
	  Float_t muR(1.0), muF(1.0);
	  if(varTitle.Contains("muR=2")) muR=2.0;
	  if(varTitle.Contains("muR=0.5")) muR=0.5;
	  if(varTitle.Contains("muF=2")) muF=2;
	  if(varTitle.Contains("muF=0.5")) muF=0.5;	  
	  if(muF==1.0 && muR==0.5) qcdVars[0]=deltaAcceptance;
	  if(muF==1.0 && muR==2.0) qcdVars[1]=deltaAcceptance;
	  if(muF==0.5 && muR==1.0) qcdVars[2]=deltaAcceptance;
	  if(muF==2.0 && muR==1.0) qcdVars[3]=deltaAcceptance;
	  if(muF==0.5 && muR==0.5) qcdVars[4]=deltaAcceptance;
	  if(muF==2.0 && muR==2.0) qcdVars[5]=deltaAcceptance;
	}
    }
  
  Float_t scaleUnc=TMath::Sqrt(
			       pow(TMath::Max(fabs(qcdVars[0]),fabs(qcdVars[1])),2)+
			       pow(TMath::Max(fabs(qcdVars[2]),fabs(qcdVars[3])),2)+
			       pow(TMath::Max(fabs(qcdVars[4]),fabs(qcdVars[5])),2)
			       );
  Float_t pdfUnc=pdfVariations->GetRMS();
  Float_t alphaSUnc=alphaSvars.size()==2 ? TMath::Max(fabs(alphaSvars[0]),fabs(alphaSvars[1])) : 0;
  //Float_t totalUnc(TMath::Sqrt(pow(acceptanceUnc,2)+pow(scaleUnc,2)+pow(pdfUnc,2)+pow(alphaSUnc,2)));
  Float_t totalUnc(acceptanceUnc+scaleUnc+pdfUnc+alphaSUnc);
  
  cout << "Inclusive jet multiplicity >=" << njets << endl
       << "\t Acc       =   " << acceptance << endl
       << "\t Stat      +/- " << acceptanceUnc << endl
       << "\t QCD scale +/- " << scaleUnc << endl
       << "\t PDF       +/- " << pdfUnc << endl
       << "\t alphaS    +/- " << alphaSUnc << endl
       << "\t --------------------------" << endl
       << "\t Total     +/- " << totalUnc << endl;

  Acceptance_t toReturn;
  toReturn.acceptance=acceptance;
  toReturn.stat=acceptanceUnc;
  toReturn.qcdScale=scaleUnc;
  toReturn.pdf=pdfUnc;
  toReturn.alphaS=alphaSUnc;
  return toReturn;
}

void scanAcceptanceEstimateFor(TTbarSample sample,TString baseDir="/store/cmst3/user/psilva/LJets2015/64217e8")
{
  TGraphErrors *grStat=new TGraphErrors();
  TGraphErrors *grScales=new TGraphErrors();
  TGraphErrors *grTotal=new TGraphErrors();
  for(Int_t ij=1; ij<5; ij++)
    {
      Acceptance_t result=acceptanceEst(sample,ij,baseDir);
      Int_t np=grStat->GetN();
      grStat->SetPoint(np,ij,result.acceptance);
      grStat->SetPointError(np,0.5,result.stat);
      grScales->SetPoint(np,ij,result.acceptance);
      grScales->SetPointError(np,0.5,result.qcdScale);
      //Float_t totalUnc(TMath::Sqrt(pow(result.stat,2)+pow(result.qcdScale,2)+pow(result.pdf,2)+pow(result.alphaS,2)));
      Float_t totalUnc(result.stat+result.qcdScale+result.pdf+result.alphaS);
      grTotal->SetPoint(np,ij,result.acceptance);
      grTotal->SetPointError(np,0.5,totalUnc);
    }

  grStat->SetMarkerStyle(20);
  grStat->SetFillStyle(0);
  grStat->SetTitle("Stat");
  grStat->SetName("Stat");
  grScales->SetMarkerStyle(1);
  grScales->SetFillStyle(1001);
  grScales->SetFillColor(kGray+1);
  grScales->SetTitle("QCD scale");
  grScales->SetName("scales");
  grTotal->SetMarkerStyle(1);
  grTotal->SetFillStyle(1001);
  grTotal->SetFillColor(kGray);
  grTotal->SetTitle("Total");
  grTotal->SetName("Total");

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TCanvas *c=new TCanvas("c","c",500,500);
  c->SetTopMargin(0.05);
  c->SetRightMargin(0.05);
  c->SetLeftMargin(0.12);
  c->SetBottomMargin(0.1);
  grTotal->Draw("a2");
  grScales->Draw("2");
  grStat->Draw("p");
  grTotal->GetXaxis()->SetTitle("Inclusive jet multiplicity");
  grTotal->GetYaxis()->SetTitle("Acceptance");
  grTotal->GetXaxis()->SetNdivisions(5);
  TLegend *leg=c->BuildLegend();
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->SetHeader("#bf{CMS} #it{simulation}");
  c->Modified();
  c->Update();
}

