#include "TopLJets2015/TopAnalysis/interface/LeptonEfficiencyWrapper.h"

#include "TFile.h"
#include "TSystem.h"

#include <iostream>


using namespace std;

//
LeptonEfficiencyWrapper::LeptonEfficiencyWrapper(bool isData,TString era)
{
  if(isData) return;
  init(era);
}

//
void LeptonEfficiencyWrapper::init(TString era)
{
  //2015 dataset
  if(era.Contains("era2015"))
    {
      era_=2015;
      TString lepEffUrl(era+"/muonEfficiencies.root");
      gSystem->ExpandPathName(lepEffUrl);
      TFile *fIn=TFile::Open(lepEffUrl);
      lepEffH_["m_sel"]=(TH2 *)fIn->Get("m_sel");
      lepEffH_["m_sel"]->SetDirectory(0);
      lepEffH_["m_singleleptrig"]=(TH2 *)fIn->Get("m_trig");
      lepEffH_["m_singleleptrig"]->SetDirectory(0);
      fIn->Close();
      
      lepEffUrl=era+"/electronEfficiencies.root";
      gSystem->ExpandPathName(lepEffUrl);
      fIn=TFile::Open(lepEffUrl);
      lepEffH_["e_sel"]=(TH2 *)fIn->Get("EGamma_SF2D");
      lepEffH_["e_sel"]->SetDirectory(0);
      fIn->Close();
    }

  //2016 dataset
  else
    {
      era_=2016;
      TString lepEffUrl(era+"/SingleMuonTrigger_Z_RunBCD_prompt80X_7p65.root");
      gSystem->ExpandPathName(lepEffUrl);
      TFile *fIn=TFile::Open(lepEffUrl);
      lepEffH_["m_singleleptrig"]=(TH2 *)fIn->Get("IsoMu22_OR_IsoTkMu22_PtEtaBins_Run273158_to_274093/efficienciesDATA/abseta_pt_DATA")->Clone();
      lepEffH_["m_singleleptrig"]->SetDirectory(0);
      fIn->Close();

      lepEffUrl=era+"/muonTrackingEffs.root";
      gSystem->ExpandPathName(lepEffUrl);
      fIn=TFile::Open(lepEffUrl);
      lepEffGr_["m_tk"]=(TGraphAsymmErrors *)fIn->Get("ratio_eta");
      fIn->Close();

      lepEffUrl=era+"/MuonID_Z_RunBCD_prompt80X_7p65.root";
      gSystem->ExpandPathName(lepEffUrl);
      fIn=TFile::Open(lepEffUrl);      
      lepEffH_["m_sel"]=(TH2 *)fIn->Get("MC_NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio")->Clone();
      lepEffH_["m_sel"]->SetDirectory(0);
      fIn->Close();

      lepEffUrl=era+"/MuonIso_Z_RunBCD_prompt80X_7p65.root";
      gSystem->ExpandPathName(lepEffUrl);
      fIn=TFile::Open(lepEffUrl);      
      TH2 *isoH=(TH2 *)fIn->Get("MC_NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1/abseta_pt_ratio");
      for(Int_t xbin=1; xbin<=(lepEffH_["m_sel"])->GetNbinsX(); xbin++)
	for(Int_t ybin=1; ybin<=(lepEffH_["m_sel"])->GetNbinsY(); ybin++)
	  {
	    float sfid(lepEffH_["m_sel"]->GetBinContent(xbin,ybin)), sfiso(isoH->GetBinContent(xbin,ybin));
	    float sfidUnc(lepEffH_["m_sel"]->GetBinError(xbin,ybin)), sfisoUnc(isoH->GetBinError(xbin,ybin));
	    float sf(sfid*sfiso), sfUnc(sqrt(pow(sfid*sfisoUnc,2)+pow(sfidUnc*sfiso,2)));
	    lepEffH_["m_sel"]->SetBinContent(xbin,ybin,sf);
	    lepEffH_["m_sel"]->SetBinError(xbin,ybin,sfUnc);
	  }
      fIn->Close();
      
      lepEffUrl=era+"/egammaEff_tight_SF2D.root";
      gSystem->ExpandPathName(lepEffUrl);
      fIn=TFile::Open(lepEffUrl);
      lepEffH_["e_sel"]=(TH2 *)fIn->Get("EGamma_SF2D")->Clone();
      lepEffH_["e_sel"]->SetDirectory(0);     
      fIn->Close();
    }
}

//
EffCorrection_t LeptonEfficiencyWrapper::getTriggerCorrection(std::vector<int> pdgId,std::vector<TLorentzVector> leptons)
{
  EffCorrection_t corr(1.0,0.0);
  if(pdgId.size()<1 || leptons.size()<1) return corr;
  if(era_==2015)
    {
      if(leptons.size()>=2)
	{
	  int cat=abs(pdgId[0]*pdgId[1]);
	  if(cat==13*13)      { corr.first=0.894; corr.second=0.894; }
	  else if(cat==11*13) { corr.first=0.931; corr.second=0.931; }
	  else if(cat==11*11) { corr.first=0.930; corr.second=0.930; }
	}
      else
	{
	  TString hname(abs(pdgId[0])==11 ? "e" : "m");
	  hname += "_singleleptrig";

	  TH2 *h=lepEffH_[hname];
	  float minEtaForEff( h->GetXaxis()->GetXmin() ), maxEtaForEff( h->GetXaxis()->GetXmax()-0.01 );
	  float etaForEff=TMath::Max(TMath::Min(float(fabs(leptons[0].Eta())),maxEtaForEff),minEtaForEff);
	  Int_t etaBinForEff=h->GetXaxis()->FindBin(etaForEff);
	  
	  float minPtForEff( h->GetYaxis()->GetXmin() ), maxPtForEff( h->GetYaxis()->GetXmax()-0.01 );
	  float ptForEff=TMath::Max(TMath::Min(float(leptons[0].Pt()),maxPtForEff),minPtForEff);
	  Int_t ptBinForEff=h->GetYaxis()->FindBin(ptForEff);
	  
	  corr.first=h->GetBinContent(etaBinForEff,ptBinForEff);
	  corr.second=h->GetBinContent(etaBinForEff,ptBinForEff);
	}
    }
  else
    {
      //cf .https://indico.cern.ch/event/532751/contributions/2170250/subcontributions/196785/attachments/1287403/1915605/Talk_TSG_0806.pdf
      if(leptons.size()>=2)
	{
	  float leadPt(TMath::Max(leptons[0].Pt(),leptons[1].Pt())), trailPt(TMath::Min(leptons[0].Pt(),leptons[1].Pt()));
	  int cat=abs(pdgId[0]*pdgId[1]);
	  if(cat==13*13)      
	    { 
	      if(leadPt<40)
		{
		  if(trailPt<40)      { corr.first=0.908; corr.second=0.039; }
		}
	      else if(leadPt<70)
		{
		  if(trailPt<40)      { corr.first=0.897; corr.second=0.033; }
		  else if(trailPt<70) { corr.first=0.859; corr.second=0.054; }
		}
	      else 
		{
		  if(trailPt<40)      { corr.first=0.851; corr.second=0.057; }
		  else if(trailPt<70) { corr.first=0.870; corr.second=0.050; }
		  else                { corr.first=0.769; corr.second=0.174; }
		}
	    }
	  else if(cat==11*13)
	    { 
	      if(abs(pdgId[0])==11) { leadPt=leptons[0].Pt(); trailPt=leptons[1].Pt(); }
	      else                  { leadPt=leptons[1].Pt(); trailPt=leptons[0].Pt(); }
	      if(leadPt<40)
                {
                  if(trailPt<40)      { corr.first=0.868; corr.second=0.032; }
                  else if(trailPt<70) { corr.first=0.874; corr.second=0.036; }
                  else                { corr.first=0.933; corr.second=0.060; }
		}
	      else if(leadPt<70)
		{
                  if(trailPt<40)      { corr.first=0.863; corr.second=0.032; }
                  else if(trailPt<70) { corr.first=0.902; corr.second=0.034; }
                  else                { corr.first=0.768; corr.second=0.070; }
                }
	      else
		{
                  if(trailPt<40)      { corr.first=0.901; corr.second=0.042; }
                  else if(trailPt<70) { corr.first=0.838; corr.second=0.052; }
                  else                { corr.first=0.950; corr.second=0.05; }
		}
	    }
	  else if(cat==11*11) 
	    {
	      if(leadPt<40)
                {
                  if(trailPt<40)      { corr.first=0.853; corr.second=0.056; }
		}
	      else if(leadPt<70)
		{
                  if(trailPt<40)      { corr.first=0.954; corr.second=0.026; }
                  else if(trailPt<70) { corr.first=0.944; corr.second=0.042; }
                }
	      else
		{
                  if(trailPt<40)      { corr.first=0.918; corr.second=0.046; }
                  else if(trailPt<70) { corr.first=0.984; corr.second=0.037; }
                  else                { corr.first=0.941; corr.second=0.122; }
		}
	    }

	  //add a 2% uncertainty
	  corr.second=sqrt(0.02*0.02+corr.second);
	}
      else 
	{
	  TString hname(abs(pdgId[0])==11 ? "e" : "m");
	  hname += "_singleleptrig";

	  if( lepEffH_.find(hname)!=lepEffH_.end() )
	    {
	      TH1 *h=lepEffH_[hname];
	      float minEtaForEff( h->GetXaxis()->GetXmin() ), maxEtaForEff( h->GetXaxis()->GetXmax()-0.01 );
	      float etaForEff=TMath::Max(TMath::Min(float(fabs(leptons[0].Eta())),maxEtaForEff),minEtaForEff);
	      Int_t etaBinForEff=h->GetXaxis()->FindBin(etaForEff);

	      float minPtForEff( h->GetYaxis()->GetXmin() ), maxPtForEff( h->GetYaxis()->GetXmax()-0.01 );
	      float ptForEff=TMath::Max(TMath::Min(float(leptons[0].Pt()),maxPtForEff),minPtForEff);
	      Int_t ptBinForEff=h->GetYaxis()->FindBin(ptForEff);

	      corr.first=h->GetBinContent(etaBinForEff,ptBinForEff);
	      corr.second=h->GetBinError(etaBinForEff,etaBinForEff);
	    }
	}
    }

  return corr;
}

//
EffCorrection_t LeptonEfficiencyWrapper::getOfflineCorrection(int pdgId,float pt,float eta)
{
  EffCorrection_t corr(1.0,0.0);

  //update correction from histo, if found
  TString idstr(abs(pdgId)==11 ? "e" : "m");
  TString hname(idstr);
  hname+="_sel";
  if( lepEffH_.find(hname)!=lepEffH_.end() )
    {
      TH2 *h=lepEffH_[hname];
      float minEtaForEff( h->GetXaxis()->GetXmin() ), maxEtaForEff( h->GetXaxis()->GetXmax()-0.01 );
      float etaForEff=TMath::Max(TMath::Min(float(fabs(eta)),maxEtaForEff),minEtaForEff);
      Int_t etaBinForEff=h->GetXaxis()->FindBin(etaForEff);

      float minPtForEff( h->GetYaxis()->GetXmin() ), maxPtForEff( h->GetYaxis()->GetXmax()-0.01 );
      float ptForEff=TMath::Max(TMath::Min(pt,maxPtForEff),minPtForEff);
      Int_t ptBinForEff=h->GetYaxis()->FindBin(ptForEff);
      
      corr.first=h->GetBinContent(etaBinForEff,ptBinForEff);
      corr.second=h->GetBinError(etaBinForEff,ptBinForEff);
      
      //tracking efficiency (if available)
      hname=idstr+"_tk";
      if(lepEffGr_.find(hname)!=lepEffGr_.end())
	{
	  Double_t x(0.),xdiff(9999.),y(0.);
	  float tkEffSF(1.0),tkEffSFUnc(0);
	  for(Int_t ip=0; ip<lepEffGr_["m_tk"]->GetN(); ip++)
	    {
	      lepEffGr_["m_tk"]->GetPoint(ip,x,y);
	      float ixdiff(TMath::Abs(eta-x));
	      if(ixdiff>xdiff) continue;
	      xdiff=ixdiff;
	      tkEffSF=y;
	      tkEffSFUnc=lepEffGr_["m_tk"]->GetErrorY(ip);
	    }
	  corr.second = sqrt(pow(tkEffSFUnc*corr.first,2)+pow(tkEffSF*corr.second,2));
	  corr.first  = corr.first*tkEffSF;
	}
    }

  return corr;
}

LeptonEfficiencyWrapper::~LeptonEfficiencyWrapper()
{
  //for(auto& it : lepEffH_) it.second->Delete();
}
