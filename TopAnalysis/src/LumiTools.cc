#include "TopLJets2015/TopAnalysis/interface/LumiTools.h"
#include "TFile.h"
#include "TH2F.h"
#include "TF1.h"
#include "TSystem.h"

#include <iostream>

using namespace std;

//
LumiTools::LumiTools(TString era,TH1 *genPU):
  countH_(0),
  era_(era)
{
  parseLumiInfo();
  defineRunPeriods();
  parsePileupWeightsMap(genPU);
}

//
std::map<Int_t,Float_t> LumiTools::lumiPerRun()
{
  return lumiPerRun_;
}

//
void LumiTools::parseLumiInfo()
{
  //read out the values from the histogram stored in lumisec.root
  TFile *inF=TFile::Open(Form("%s/lumisec.root",era_.Data()),"READ");
  if(inF==0) return;
  if(inF->IsZombie()) return;
  TH2F *h=(TH2F *)inF->Get("lumisec_inc");
  int nruns(h->GetNbinsX());
  countH_=new TH1F("ratevsrun","ratevsrun;Run;Events/pb",nruns,0,nruns);
  countH_->SetDirectory(0);
  for(int xbin=1; xbin<=nruns; xbin++)
    {
      TString run=h->GetXaxis()->GetBinLabel(xbin);
      lumiPerRun_[run.Atoi()]=h->GetBinContent(xbin);
      countH_->GetXaxis()->SetBinLabel(xbin,run);
    }
  inF->Close();
};

//
void LumiTools::defineRunPeriods()
{
  runPeriods_.clear();
  if(era_.Contains("era2017") || era_.Contains("era2016"))
    {
      runPeriods_.push_back(std::pair<TString,float> ("",1.0));
    }
}

//
TString LumiTools::assignRunPeriod()
{
  if(runPeriods_.size()==0) return "";
  if(runPeriods_.size()==1) return runPeriods_[0].first;

  float totalLumi(0.);
  for (auto periodLumi : runPeriods_) totalLumi += periodLumi.second;

  //generate randomly in the total lumi range to pick one of the periods
  float pickLumi(rand_.Uniform(totalLumi));
  float testLumi(0); 
  int iLumi(0);
  for (auto periodLumi : runPeriods_) {
    testLumi += periodLumi.second;
    if (pickLumi < testLumi) break;
    else ++iLumi;
  }

  //return the period
  return runPeriods_[iLumi].first;
}

//
void LumiTools::parsePileupWeightsMap(TH1 *genPU)
{  
  if(genPU==0) return;
  float totalExp(genPU->Integral());
  if(totalExp<=0) return;
  genPU->Scale(1./totalExp);
  for (auto period : runPeriods_) {

    puWgtGr_[period.first]=std::vector<TH1 *>(3,0);
    //puWgtGr_[period.first]=std::vector<TGraph *>(3,0);

    //readout the pileup weights and take the ratio of data/MC
    TString puWgtUrl(era_+"/pileupWgts"+period.first+".root");
    gSystem->ExpandPathName(puWgtUrl);
    TFile *fIn=TFile::Open(puWgtUrl);
    for(size_t i=0; i<3; i++)
      {
        TString grName("pu_nom");
        if(i==1) grName="pu_down";
        if(i==2) grName="pu_up";
        TGraph *puData=(TGraph *)fIn->Get(grName);
        Float_t totalData=puData->Integral();
        TH1 *tmp=(TH1 *)genPU->Clone("tmp");
        for(Int_t xbin=1; xbin<=tmp->GetXaxis()->GetNbins(); xbin++)
          {
            Float_t yexp=genPU->GetBinContent(xbin);
            Double_t xobs,yobs;
            puData->GetPoint(xbin-1,xobs,yobs);
            float wgt(yexp>0 ? float(yobs)/(totalData*yexp) : 0.);
            if(isinf(wgt) || isnan(wgt)) wgt=0.;
            tmp->SetBinContent(xbin,wgt);
          }
        
        grName.ReplaceAll("pu","puwgts");
        tmp->SetName(period.first+grName);
        tmp->SetDirectory(0);
        puWgtGr_[period.first][i]=tmp;

        //TGraph *gr=new TGraph(tmp);
        //grName.ReplaceAll("pu","puwgts");
        //gr->SetName(period.first+grName);
        //puWgtGr_[period.first][i]=gr;
        //tmp->Delete();
      }
  }
}

//
std::vector<Float_t> LumiTools::pileupWeight(Float_t genPu,TString period)
{
  using namespace std;
  std::vector<Float_t> toReturn(3,1.0);
  if(puWgtGr_.find(period)==puWgtGr_.end()) return toReturn;
  for(size_t i=0; i<3; i++)
    {      
      TH1 *h=puWgtGr_[period][i];
      float minPu( h->GetXaxis()->GetXmin() ), maxPu( h->GetXaxis()->GetXmax()-0.01 );
      float puForEval=TMath::Max(TMath::Min(genPu,maxPu),minPu);
      toReturn[i]=h->GetBinContent(puForEval);
      //toReturn[i]=puWgtGr_[period][i]->Eval(genPu);
    }

  return toReturn;
}
