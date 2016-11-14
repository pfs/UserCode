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

enum TTbarSample {PY8,HWpp,PSUP,PSDN};
void acceptanceEst(TTbarSample sample,TString baseDir="~/work/LJets-5TeV/analysis_mu",bool addElectron=false);

void acceptanceEst(TTbarSample sample,TString baseDir,bool addElectron)
{
  float brCorr(1.0);
  TString file("");  
  if(sample==PY8)  file=baseDir+"/MCTTNominal_v2.root";
  if(sample==HWpp) file=baseDir+"/MCTTHerwig_v1.root";
  if(sample==PSUP) file=baseDir+"/MCTTScaleUp_v2.root";
  if(sample==PSDN) file=baseDir+"/MCTTScaleDown_v2.root";
					
  //read results from file
  TFile *inF=TFile::Open(file);
  if(inF==0) return;
  TH1 *hcounter = (TH1*)inF->Get("fidcounter")->Clone("hcounter");
  hcounter->SetDirectory(0);
  TH1 *wgtCounter=(TH1*)inF->Get("wgtcounter")->Clone("hnormcounter");
  wgtCounter->SetDirectory(0);
  inF->Close();

  if(addElectron)
    {
      file=file.ReplaceAll("analysis_mu","analysis_e");
      inF=TFile::Open(file);
      if(inF){
	hcounter->Add((TH1*)inF->Get("fidcounter")->Clone("hcounter") );
	//wgtCounter->Add( (TH1*)inF->Get("wgtcounter")->Clone("hnormcounter") );
	inF->Close();
      }
    }
  hcounter->Divide(wgtCounter);

  
  //report on final result
  float acceptance=hcounter->GetBinContent(1);
  float acceptanceUnc=hcounter->GetBinError(1);
  cout << "\t Acc       =   " << acceptance << endl
       << "\t Stat      +/- " << acceptanceUnc << endl;
    
  //systematic variations
  float totalUnc(acceptanceUnc);
  if(sample==PY8)
    {

      std::vector<float>qcdVars(6,0);
      TGraph *pdfVariations=new TGraph;
      std::vector<float> alphaSvars;
      for(int xbin=2; xbin<=111; xbin++)
	{
	  float acceptanceVar=hcounter->GetBinContent(xbin);
	  if(acceptanceVar==0) continue;
	  float deltaAcceptance=acceptanceVar-acceptance;
	  if(xbin>=10 && xbin<=109)
	    pdfVariations->SetPoint(pdfVariations->GetN(),deltaAcceptance,deltaAcceptance);
	  if(xbin==110 || xbin==111)
	    alphaSvars.push_back(deltaAcceptance);
	  if(xbin==2) qcdVars[3]=deltaAcceptance;
	  if(xbin==3) qcdVars[2]=deltaAcceptance;
	  if(xbin==4) qcdVars[1]=deltaAcceptance;
	  if(xbin==7) qcdVars[0]=deltaAcceptance;
	  if(xbin==5) qcdVars[5]=deltaAcceptance;
	  if(xbin==9) qcdVars[4]=deltaAcceptance;
	}
      
      Float_t scaleUnc=TMath::Sqrt(
				   pow(TMath::Max(fabs(qcdVars[0]),fabs(qcdVars[1])),2)+
				   pow(TMath::Max(fabs(qcdVars[2]),fabs(qcdVars[3])),2)+
				   pow(TMath::Max(fabs(qcdVars[4]),fabs(qcdVars[5])),2)
				   );
      Float_t pdfUnc=pdfVariations->GetRMS();
      Float_t alphaSUnc=alphaSvars.size()==2 ? TMath::Max(fabs(alphaSvars[0]),fabs(alphaSvars[1])) : 0;
      totalUnc+=scaleUnc+pdfUnc+alphaSUnc;

      cout << "\t QCD scale +/- " << scaleUnc << endl
	   << "\t PDF       +/- " << pdfUnc << endl
	   << "\t alphaS    +/- " << alphaSUnc << endl;
    }
  cout << "\t --------------------------" << endl
       << "\t Total     +/- " << totalUnc << endl;
}

