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
  // 2016 pPb/Pbp dataset
   else if(era.Contains("era2016pPb") || era.Contains("era2016Pbp"))
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

      lepEffUrl=era+"/muonIsoHFEffs.root";
      gSystem->ExpandPathName(lepEffUrl);
      fIn=TFile::Open(lepEffUrl);
      lepEffGr_["m_isoHF"]=(TGraphAsymmErrors *)fIn->Get("Graph");
      fIn->Close();

      lepEffUrl=era+"/electronRECOEfficiencies.root";
      gSystem->ExpandPathName(lepEffUrl);
      fIn=TFile::Open(lepEffUrl);
      lepEffH_["e_reco"]=(TH2 *)fIn->Get("EGamma_SF2D");
      lepEffH_["e_reco"]->SetDirectory(0);
      fIn->Close();

      lepEffUrl=era+"/electronIDEfficiencies.root";
      gSystem->ExpandPathName(lepEffUrl);
      fIn=TFile::Open(lepEffUrl);
      lepEffH_["e_sel"]=(TH2 *)fIn->Get("EGamma_SF2D");
      lepEffH_["e_sel"]->SetDirectory(0);
      fIn->Close();

      lepEffUrl=era+"/electronHLTEfficiencies.root";
      gSystem->ExpandPathName(lepEffUrl);
      fIn=TFile::Open(lepEffUrl);
      lepEffH_["e_singleleptrig"]=(TH2 *)fIn->Get("EGamma_SF2D")->Clone();
      lepEffH_["e_singleleptrig"]->SetDirectory(0);     
      fIn->Close();
      
      lepEffUrl=era+"/electronIsoHFEffs.root";
      gSystem->ExpandPathName(lepEffUrl);
      fIn=TFile::Open(lepEffUrl);
      lepEffGr_["e_isoHF"]=(TGraphAsymmErrors *)fIn->Get("Graph");
      fIn->Close();

    }  
  //2016 dataset
  else if(era.Contains("era2016"))
    {
      era_=2016;
      std::vector<TString> periods = { "BCDEF", "GH" };
      bool onlyGH = false;
      if(era.Contains("era2016GH")) {
        periods = { "GH" };
        onlyGH = true;
        era.Remove(era.Length()-2, 2);
      }
      for (auto period : periods) {
        //MUONS
        TString lepEffUrl(era+"/SingleMuonTriggerEfficienciesAndSF_"+period+".root");
        gSystem->ExpandPathName(lepEffUrl);
        TFile *fIn=TFile::Open(lepEffUrl);
        lepEffH_["m_singleleptrig"+period]=(TH2 *)fIn->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio")->Clone();
        lepEffH_["m_singleleptrig"+period]->SetDirectory(0);
        fIn->Close();

        lepEffUrl=era+"/MuonTracking_EfficienciesAndSF_"+period+".root";
        gSystem->ExpandPathName(lepEffUrl);
        fIn=TFile::Open(lepEffUrl);
        //global muon
        lepEffGr_["m_tk"+period]=(TGraphAsymmErrors *)fIn->Get("ratio_eff_eta3_dr030e030_corr");
        lepEffGr_["m_tk_aeta"+period]=(TGraphAsymmErrors *)fIn->Get("ratio_eff_aeta_dr030e030_corr");
        lepEffGr_["m_tk_nvtx"+period]=(TGraphAsymmErrors *)fIn->Get("ratio_eff_vtx_dr030e030_corr");
        //tracker-only (used for charged tracks correction and ucnertainty)
        lepEffGr_["m_tk0"+period]=(TGraphAsymmErrors *)fIn->Get("ratio_eff_eta3_tk0_dr030e030_corr");
        lepEffGr_["m_tk0_aeta"+period]=(TGraphAsymmErrors *)fIn->Get("ratio_eff_aeta_tk0_dr030e030_corr");
        lepEffGr_["m_tk0_nvtx"+period]=(TGraphAsymmErrors *)fIn->Get("ratio_eff_vtx_tk0_dr030e030_corr");
        fIn->Close();

        lepEffUrl=era+"/MuonIdEfficienciesAndSF_"+period+".root";
        gSystem->ExpandPathName(lepEffUrl);
        fIn=TFile::Open(lepEffUrl);      
        lepEffH_["m_sel"+period]=(TH2 *)fIn->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio")->Clone();
        lepEffH_["m_sel"+period]->SetDirectory(0);
        fIn->Close();

        lepEffUrl=era+"/MuonIsoEfficienciesAndSF_"+period+".root";
        gSystem->ExpandPathName(lepEffUrl);
        fIn=TFile::Open(lepEffUrl);      
        TH2 *isoH=(TH2 *)fIn->Get("TightISO_TightID_pt_eta/abseta_pt_ratio");
        for(Int_t xbin=1; xbin<=(lepEffH_["m_sel"+period])->GetNbinsX(); xbin++)
          for(Int_t ybin=1; ybin<=(lepEffH_["m_sel"+period])->GetNbinsY(); ybin++)
            {
              float sfid(lepEffH_["m_sel"+period]->GetBinContent(xbin,ybin)), sfiso(isoH->GetBinContent(xbin,ybin));
              float sfidUnc(lepEffH_["m_sel"+period]->GetBinError(xbin,ybin)), sfisoUnc(isoH->GetBinError(xbin,ybin));
              float sf(sfid*sfiso), sfUnc(sqrt(pow(sfid*sfisoUnc,2)+pow(sfidUnc*sfiso,2)));
              lepEffH_["m_sel"+period]->SetBinContent(xbin,ybin,sf);
              lepEffH_["m_sel"+period]->SetBinError(xbin,ybin,sfUnc);
            }
        fIn->Close();
        
        //ELECTRONS
        lepEffUrl=era+"/SingleElectron_TriggerSF_Run2016"+period+"_v3_prel.root";
        gSystem->ExpandPathName(lepEffUrl);
        fIn=TFile::Open(lepEffUrl);
        lepEffH_["e_singleleptrig"+period]=(TH2 *)fIn->Get("SF")->Clone();
        lepEffH_["e_singleleptrig"+period]->SetDirectory(0);
        fIn->Close();
        
        lepEffUrl=era+"/ElectronReco_egammaEffi.txt_EGM2D.root";
        if(onlyGH)
          lepEffUrl=era+"/ElectronReco_egammaEffi.txt_EGM2D_GH.root";
        gSystem->ExpandPathName(lepEffUrl);
        fIn=TFile::Open(lepEffUrl);
        lepEffH_["e_rec"+period]=(TH2 *)fIn->Get("EGamma_SF2D")->Clone();
        lepEffH_["e_rec"+period]->SetDirectory(0);     
        fIn->Close();
        
        lepEffUrl=era+"/ElectronIdTight_egammaEffi.txt_EGM2D.root";
        //TODO
        //if(onlyGH)
        //  lepEffUrl=era+"/ElectronIdTight_egammaEffi.txt_EGM2D_GH.root";
        gSystem->ExpandPathName(lepEffUrl);
        fIn=TFile::Open(lepEffUrl);      
	lepEffH_["e_sel"+period]=(TH2 *)fIn->Get("EGamma_SF2D")->Clone();
        lepEffH_["e_sel"+period]->SetDirectory(0);
        fIn->Close();
      }
    }
}

//
EffCorrection_t LeptonEfficiencyWrapper::getTriggerCorrection(std::vector<Particle> &leptons,TString period)
{
  EffCorrection_t corr(1.0,0.0);
  if(leptons.size()<1) return corr;
  if(era_==2015)
    {
      if(leptons.size()>=2)
        {
          int cat=abs(leptons[0].id()*leptons[1].id());
          if(cat==13*13)      { corr.first=0.894; corr.second=0.894; }
          else if(cat==11*13) { corr.first=0.931; corr.second=0.931; }
          else if(cat==11*11) { corr.first=0.930; corr.second=0.930; }
        }
      else
        {
          TString hname(abs(leptons[0].id())==11 ? "e" : "m");
          hname += "_singleleptrig"+period;

          TH2 *h=lepEffH_[hname];
          float minEtaForEff( h->GetXaxis()->GetXmin() ), maxEtaForEff( h->GetXaxis()->GetXmax()-0.01 );
          float etaForEff=TMath::Max(TMath::Min(float(fabs(leptons[0].eta())),maxEtaForEff),minEtaForEff);
          Int_t etaBinForEff=h->GetXaxis()->FindBin(etaForEff);
          
          float minPtForEff( h->GetYaxis()->GetXmin() ), maxPtForEff( h->GetYaxis()->GetXmax()-0.01 );
          float ptForEff=TMath::Max(TMath::Min(float(leptons[0].pt()),maxPtForEff),minPtForEff);
          Int_t ptBinForEff=h->GetYaxis()->FindBin(ptForEff);
          
          corr.first=h->GetBinContent(etaBinForEff,ptBinForEff);
          corr.second=h->GetBinContent(etaBinForEff,ptBinForEff);
        }
    }
  else if(era_==2016)
    {
      //cf AN 2016/392 v2 Figs 3(ee) 6(mm) 23 (em)
      if(leptons.size()>=2)
        {
          float leadPt(TMath::Max(leptons[0].pt(),leptons[1].pt())), trailPt(TMath::Min(leptons[0].pt(),leptons[1].pt()));
          int cat=abs(leptons[0].id()*leptons[1].id());
          if(cat==13*13)      
            { 
              if(leadPt<40)
          {
            if(trailPt<40)      { corr.first=0.995; corr.second=0.006; }
          }
              else if(leadPt<70)
          {
            if(trailPt<40)      { corr.first=1.001; corr.second=0.004; }
            else if(trailPt<70) { corr.first=0.998; corr.second=0.003; }
          }
              else 
          {
            if(trailPt<40)      { corr.first=0.983; corr.second=0.006; }
            else if(trailPt<70) { corr.first=0.997; corr.second=0.004; }
            else                { corr.first=0.999; corr.second=0.007; }
          }
            }
          else if(cat==11*13)
            { 
              if(abs(leptons[0].id())==11) { leadPt=leptons[0].pt(); trailPt=leptons[1].pt(); }
              else                         { leadPt=leptons[1].pt(); trailPt=leptons[0].pt(); }
              if(leadPt<40)
                      {
                        if(trailPt<40)      { corr.first=0.980; corr.second=0.009; }
                        else if(trailPt<70) { corr.first=0.979; corr.second=0.010; }
                        else                { corr.first=0.960; corr.second=0.017; }
          }
              else if(leadPt<70)
          {
                        if(trailPt<40)      { corr.first=0.988; corr.second=0.008; }
                        else if(trailPt<70) { corr.first=0.985; corr.second=0.009; }
                        else                { corr.first=0.953; corr.second=0.015; }
                      }
              else
          {
                        if(trailPt<40)      { corr.first=0.991; corr.second=0.017; }
                        else if(trailPt<70) { corr.first=0.988; corr.second=0.011; }
                        else                { corr.first=1.006; corr.second=0.017; }
          }
            }
          else if(cat==11*11) 
            {
              if(leadPt<40)
                      {
                        if(trailPt<40)      { corr.first=0.980; corr.second=0.010; }
          }
              else if(leadPt<70)
          {
                        if(trailPt<40)      { corr.first=0.988; corr.second=0.004; }
                        else if(trailPt<70) { corr.first=1.002; corr.second=0.004; }
                      }
              else
          {
                        if(trailPt<40)      { corr.first=0.989; corr.second=0.005; }
                        else if(trailPt<70) { corr.first=0.990; corr.second=0.005; }
                        else                { corr.first=0.989; corr.second=0.010; }
          }
            }

          //add a 2% uncertainty
          corr.second=sqrt(0.02*0.02+corr.second);
        }
      else 
        {
          TString hname(abs(leptons[0].id())==11 ? "e" : "m");
          hname += "_singleleptrig"+period;

          if( abs(leptons[0].id())==13 and lepEffH_.find(hname)!=lepEffH_.end() )
            {
              TH1 *h=lepEffH_[hname];
              float minEtaForEff( h->GetXaxis()->GetXmin() ), maxEtaForEff( h->GetXaxis()->GetXmax()-0.01 );
              float etaForEff=TMath::Max(TMath::Min(float(fabs(leptons[0].eta())),maxEtaForEff),minEtaForEff);
              Int_t etaBinForEff=h->GetXaxis()->FindBin(etaForEff);

              float minPtForEff( h->GetYaxis()->GetXmin() ), maxPtForEff( h->GetYaxis()->GetXmax()-0.01 );
              float ptForEff=TMath::Max(TMath::Min(float(leptons[0].pt()),maxPtForEff),minPtForEff);
              Int_t ptBinForEff=h->GetYaxis()->FindBin(ptForEff);

              corr.first=h->GetBinContent(etaBinForEff,ptBinForEff);
              corr.second=h->GetBinError(etaBinForEff,ptBinForEff);
            }
          //electron histogram uses eta, not abs(eta)
          else if( abs(leptons[0].id())==11 and lepEffH_.find(hname)!=lepEffH_.end() )
            {
              TH1 *h=lepEffH_[hname];
              float minEtaForEff( h->GetXaxis()->GetXmin() ), maxEtaForEff( h->GetXaxis()->GetXmax()-0.01 );
              float etaForEff=TMath::Max(TMath::Min(float(leptons[0].eta()),maxEtaForEff),minEtaForEff);
              Int_t etaBinForEff=h->GetXaxis()->FindBin(etaForEff);

              float minPtForEff( h->GetYaxis()->GetXmin() ), maxPtForEff( h->GetYaxis()->GetXmax()-0.01 );
              float ptForEff=TMath::Max(TMath::Min(float(leptons[0].pt()),maxPtForEff),minPtForEff);
              Int_t ptBinForEff=h->GetYaxis()->FindBin(ptForEff);

              corr.first=h->GetBinContent(etaBinForEff,ptBinForEff);
              corr.second=h->GetBinError(etaBinForEff,ptBinForEff);
            }
        }
    }

  return corr;
}

//
EffCorrection_t LeptonEfficiencyWrapper::getOfflineCorrection(int pdgId,float pt,float eta,TString period)
{
  EffCorrection_t corr(1.0,0.0);

  //update correction from histo, if found
  TString idstr(abs(pdgId)==11 ? "e" : "m");
  TString hname(idstr);
  hname+="_sel"+period;
  if( lepEffH_.find(hname)!=lepEffH_.end() )
    {
      TH2 *h=lepEffH_[hname];
      float minEtaForEff( h->GetXaxis()->GetXmin() ), maxEtaForEff( h->GetXaxis()->GetXmax()-0.01 );
      float etaForEff;
      if (minEtaForEff >= 0.) //axis is abseta
        etaForEff=TMath::Max(TMath::Min(float(fabs(eta)),maxEtaForEff),minEtaForEff);
      else //axis is signed eta
        etaForEff=TMath::Max(TMath::Min(float(eta),maxEtaForEff),minEtaForEff);
      Int_t etaBinForEff=h->GetXaxis()->FindBin(etaForEff);

      float minPtForEff( h->GetYaxis()->GetXmin() ), maxPtForEff( h->GetYaxis()->GetXmax()-0.01 );
      float ptForEff=TMath::Max(TMath::Min(pt,maxPtForEff),minPtForEff);
      Int_t ptBinForEff=h->GetYaxis()->FindBin(ptForEff);
      
      corr.first=h->GetBinContent(etaBinForEff,ptBinForEff);
      corr.second=h->GetBinError(etaBinForEff,ptBinForEff);
      
      //reco efficiency (if available)
      hname=idstr+"_rec"+period;
      if(lepEffH_.find(hname)!=lepEffH_.end())
        {
	  float minEtaForEff( lepEffH_[hname]->GetXaxis()->GetXmin() ), maxEtaForEff( lepEffH_[hname]->GetXaxis()->GetXmax()-0.01 );
	  float etaForEff( minEtaForEff >= 0. ?
			   TMath::Max(TMath::Min(float(fabs(eta)),maxEtaForEff),minEtaForEff) :
			   TMath::Max(TMath::Min(float(eta),maxEtaForEff),minEtaForEff) );
	  Int_t etaBinForEff=lepEffH_[hname]->GetXaxis()->FindBin(etaForEff);

	  float minPtForEff( lepEffH_[hname]->GetYaxis()->GetXmin() ), maxPtForEff( lepEffH_[hname]->GetYaxis()->GetXmax()-0.01 );
	  float ptForEff=TMath::Max(TMath::Min(pt,maxPtForEff),minPtForEff);
	  Int_t ptBinForEff=lepEffH_[hname]->GetYaxis()->FindBin(ptForEff);
	  
	  float recEff=lepEffH_[hname]->GetBinContent(etaBinForEff,ptBinForEff);
	  float recEffUnc=lepEffH_[hname]->GetBinError(etaBinForEff,ptBinForEff);

	  corr.second=sqrt(pow(recEff*corr.second,2)+pow(recEffUnc*corr.first,2));
	  corr.first*=recEff;
	}

      //tracking efficiency (if available)
      //TODO: add nvtx-dependency? divide by average sf(aeta) then
      hname=idstr+"_tk_aeta"+period;
      if(lepEffGr_.find(hname)!=lepEffGr_.end())
        {
          Double_t x(0.),xdiff(9999.),y(0.);
          float tkEffSF(1.0),tkEffSFUnc(0);
          for(Int_t ip=0; ip<lepEffGr_[hname]->GetN(); ip++)
            {
              lepEffGr_[hname]->GetPoint(ip,x,y);
              float ixdiff(TMath::Abs(fabs(eta)-x));
              if(ixdiff>xdiff) continue;
              xdiff=ixdiff;
              tkEffSF=y;
              tkEffSFUnc=lepEffGr_[hname]->GetErrorY(ip);
            }
          corr.second = sqrt(pow(tkEffSFUnc*corr.first,2)+pow(tkEffSF*corr.second,2));
          corr.first  = corr.first*tkEffSF;
        }
      
      //reco efficiency (if available)
      hname=idstr+"_reco";
      if(lepEffH_.find(hname)!=lepEffH_.end() )
	{
	  TH2 *h=lepEffH_[hname];
	  float minEtaForEff( h->GetXaxis()->GetXmin() ), maxEtaForEff( h->GetXaxis()->GetXmax()-0.01 );
	  float etaForEff=TMath::Max(TMath::Min(float(fabs(eta)),maxEtaForEff),minEtaForEff);
	  Int_t etaBinForEff=h->GetXaxis()->FindBin(etaForEff);
	  
	  float minPtForEff( h->GetYaxis()->GetXmin() ), maxPtForEff( h->GetYaxis()->GetXmax()-0.01 );
	  float ptForEff=TMath::Max(TMath::Min(pt,maxPtForEff),minPtForEff);
	  Int_t ptBinForEff=h->GetYaxis()->FindBin(ptForEff);
	  
	  corr.second = sqrt(pow(h->GetBinError(etaBinForEff,ptBinForEff)*corr.first,2)+pow(h->GetBinError(etaBinForEff,ptBinForEff)*corr.second,2));
	  corr.first  = corr.first*h->GetBinContent(etaBinForEff,ptBinForEff);
	  
	}
    }
  
  return corr;
}

//
EffCorrection_t LeptonEfficiencyWrapper::getOfflineCorrection(Particle lepton,TString period)
{
  return getOfflineCorrection(lepton.id(),lepton.pt(),lepton.eta(), period);
}

//
EffCorrection_t LeptonEfficiencyWrapper::getTriggerCorrection(std::vector<int> &pdgId, std::vector<TLorentzVector> &leptons,TString period)
{
  std::vector<Particle> lepParts;
  for(size_t i=0; i<pdgId.size(); i++)
    lepParts.push_back( Particle(leptons[i],0,pdgId[i],0,0,1) );
  return getTriggerCorrection(lepParts, period);
}

EffCorrection_t LeptonEfficiencyWrapper::getOfflineIsoHFCorrection(int pdgId,float hf)
{
  EffCorrection_t corr(1.0,0.0);

  //update correction from graph, if found
  TString idstr(abs(pdgId)==11 ? "e" : "m");
  TString hname(idstr);
  hname+="_isoHF";
  if(lepEffGr_.find(hname)!=lepEffGr_.end())
    {
      corr.second = 0;
      corr.first  = lepEffGr_[hname]->Eval(hf);
    }

    
  return corr;
}

//
EffCorrection_t LeptonEfficiencyWrapper::getTrackingCorrection(int nvtx, TString period)
{
  EffCorrection_t corr(1.0,0.0);
  
  //nvtx-dependent tracking efficiency correction
  TString hname = "m_tk0_nvtx"+period;
  if(lepEffGr_.find(hname)!=lepEffGr_.end()) {
    corr.first  = lepEffGr_[hname]->Eval(TMath::Min(nvtx,40));
    corr.second = 1.-corr.first;
  }
  
  return corr;
}

LeptonEfficiencyWrapper::~LeptonEfficiencyWrapper()
{
  //for(auto& it : lepEffH_) it.second->Delete();
}
