#include "TopLJets2015/TopAnalysis/interface/EfficiencyScaleFactorsWrapper.h"
#include "TopLJets2015/TopAnalysis/interface/JSONWrapper.h"

#include "TFile.h"
#include "TSystem.h"

#include <iostream>


using namespace std;

//
EfficiencyScaleFactorsWrapper::EfficiencyScaleFactorsWrapper(bool isData,TString era)
{
  if(isData) return;
  init(era);
}

EfficiencyScaleFactorsWrapper::EfficiencyScaleFactorsWrapper(bool isData,TString era,std::map<TString,TString> cfgMap):
  cfgMap_(cfgMap)
{
  if(isData) return;
  init(era);
}

//
void EfficiencyScaleFactorsWrapper::init(TString era)
{
  if(era.Contains("era2017")) era_=2017;
  if(era.Contains("era2016")) era_=2016;

  cout << "[EfficiencyScaleFactorsWrapper]" << endl
       << "\tStarting efficiency scale factors for " << era << endl
       << "\tWarnings: no trigger SFs for any object" << endl
       << "\t          uncertainties returned are of statistical nature only" << endl
       << "\tDon't forget to fix these and update these items!" << endl;

  //PHOTONS
  TString gid(era_==2016? "MVAWP80" : "MVAwp80");
  if(cfgMap_.find("g_id")!=cfgMap_.end()) gid=cfgMap_["g_id"]; 
  TString url(era+"/2017_Photons"+gid+".root");
  if(era_==2016) url=era+"/2016LegacyReReco_Photon"+gid+".root";
  gSystem->ExpandPathName(url);
  TFile *fIn=TFile::Open(url);
  scaleFactorsH_["g_id"]=(TH2 *)fIn->Get("EGamma_SF2D")->Clone();
  scaleFactorsH_["g_id"]->SetDirectory(0);
  fIn->Close();
  
  url=era+"/egammaEffi.txt_EGM2D.root";
  gSystem->ExpandPathName(url);
  fIn=TFile::Open(url);
  scaleFactorsH_["g_rec"]=(TH2 *)fIn->Get("EGamma_SF2D")->Clone();
  scaleFactorsH_["g_rec"]->SetDirectory(0);     
  fIn->Close();
  
  //MUONS
  if(era_==2016){
    url=era+"/MuonTracking_EfficienciesAndSF_BCDEF.root";
    gSystem->ExpandPathName(url);
    fIn=TFile::Open(url);
    scaleFactorsGr_["m_tk"]=(TGraphAsymmErrors *)fIn->Get("ratio_eff_aeta_dr030e030_corr");
    fIn->Close();

    url=era+"/MuonTracking_EfficienciesAndSF_GH.root";
    gSystem->ExpandPathName(url);
    fIn=TFile::Open(url);
    scaleFactorsGr_["m_tkGH"]=(TGraphAsymmErrors *)fIn->Get("ratio_eff_aeta_dr030e030_corr");
    fIn->Close();


  } else {
    std::cout << "No tracking effciency for 2017 muons!!!" << endl;
  }

  TString mid("TightID");
  if(cfgMap_.find("m_id")!=cfgMap_.end()) mid=cfgMap_["m_id"];
  if(era_==2016) {
    url=era+"/RunBCDEF_SF_MuID.root";
    gSystem->ExpandPathName(url); 
    fIn=TFile::Open(url);
    scaleFactorsH_["m_id"]=(TH2F *)fIn->Get("NUM_"+mid+"_DEN_genTracks_eta_pt")->Clone();
    scaleFactorsH_["m_id"]->SetDirectory(0);
    fIn->Close();

    url=era+"/RunGH_SF_MuID.root";
    gSystem->ExpandPathName(url);
    fIn=TFile::Open(url);
    scaleFactorsH_["m_idGH"]=(TH2F *)fIn->Get("NUM_"+mid+"_DEN_genTracks_eta_pt")->Clone();
    scaleFactorsH_["m_idGH"]->SetDirectory(0);
    fIn->Close();
  } else {
    url=era+"/RunBCDEF_SF_MuID.root"; 
    fIn=TFile::Open(url);
    scaleFactorsH_["m_id"]=(TH2F *)fIn->Get("NUM_"+mid+"_DEN_genTracks_pt_abseta")->Clone();
    scaleFactorsH_["m_id"]->SetDirectory(0);
    fIn->Close();
  }

  TString miso("TightRelIso"),mid4iso(mid);
  if(cfgMap_.find("m_iso")!=cfgMap_.end())    miso=cfgMap_["m_iso"];
  if(cfgMap_.find("m_id4iso")!=cfgMap_.end()) mid4iso=cfgMap_["m_id4iso"];
  if(era_==2016) { 
    url=era+"/RunBCDEF_SF_MuISO.root";
    gSystem->ExpandPathName(url);
    fIn=TFile::Open(url);
    scaleFactorsH_["m_iso"]=(TH2F *)fIn->Get("NUM_"+miso+"_DEN_"+mid4iso+"_eta_pt")->Clone();
    scaleFactorsH_["m_iso"]->SetDirectory(0);
    fIn->Close();
    url=era+"/RunGH_SF_MuISO.root";
    gSystem->ExpandPathName(url);
    fIn=TFile::Open(url);
    scaleFactorsH_["m_isoGH"]=(TH2F *)fIn->Get("NUM_"+miso+"_DEN_"+mid4iso+"_eta_pt")->Clone();
    scaleFactorsH_["m_isoGH"]->SetDirectory(0);
    fIn->Close();
  }else {
    url=era+"/RunBCDEF_SF_MuISO.root";
    gSystem->ExpandPathName(url);
    fIn=TFile::Open(url);
    scaleFactorsH_["m_iso"]=(TH2F *)fIn->Get("NUM_"+miso+"_DEN_"+mid4iso+"_pt_abseta")->Clone();
    scaleFactorsH_["m_iso"]->SetDirectory(0);
    fIn->Close();
  }

  //ELECTRONS
  url=era+"/egammaEffi.txt_EGM2D.root";
  gSystem->ExpandPathName(url);
  fIn=TFile::Open(url);
  scaleFactorsH_["e_rec"]=(TH2 *)fIn->Get("EGamma_SF2D")->Clone();
  scaleFactorsH_["e_rec"]->SetDirectory(0);     
  fIn->Close();
  
  TString eid("MVA80");
  if(cfgMap_.find("e_id")!=cfgMap_.end()) eid=cfgMap_["e_id"]; 
  url=era+"/2017_Electron"+eid+".root";
  if(era_==2016) url=era+"/2016LegacyReReco_Electron"+eid+".root";
  gSystem->ExpandPathName(url);
  fIn=TFile::Open(url);      
  scaleFactorsH_["e_id"]=(TH2 *)fIn->Get("EGamma_SF2D")->Clone();
  scaleFactorsH_["e_id"]->SetDirectory(0);
  fIn->Close();
  
}

//
EffCorrection_t EfficiencyScaleFactorsWrapper::getDileptonTriggerCorrection(std::vector<Particle> &leptons){

  int dilCode=abs(leptons[0].id()*leptons[1].id());
  EffCorrection_t corr(1.0,0.0);
  if(era_==2016) {
    //values from AN 2016/392 (v3) 
    //as scale factors are approximately constant as function of number of jets take:
    //   - central value from events with 2 jets
    //   - uncertainty as 2 jets unc. + max. difference wrt to other jet mult
    //the offline ids to use must be tight ones
    float leta(fabs(leptons[0].Eta()));
    if(dilCode==11*13) {
      if(leta<0.09)      { corr.first=0.995; corr.second=TMath::Sqrt(0.002*0.002+0.005*0.005); }
      else if (leta<1.2) { corr.first=0.996; corr.second=TMath::Sqrt(0.003*0.003+0.011*0.011); }
      else if (leta<2.1) { corr.first=0.994; corr.second=TMath::Sqrt(0.003*0.003+0.019*0.019); }
      else               { corr.first=0.992; corr.second=TMath::Sqrt(0.011*0.011+0.016*0.016); } //0 jets has too large error
    }
    if(dilCode==11*11) {
      if(leta<0.09)      { corr.first=0.991; corr.second=TMath::Sqrt(0.002*0.002+0.027*0.027); }
      else if (leta<1.2) { corr.first=0.994; corr.second=TMath::Sqrt(0.004*0.004+0.018*0.018); }
      else if (leta<2.1) { corr.first=0.992; corr.second=TMath::Sqrt(0.004*0.004+0.025*0.025); }
      else               { corr.first=0.982; corr.second=TMath::Sqrt(0.014*0.014+0.029*0.029); }
    }
    if(dilCode==13*13) {
      if(leta<0.09)      { corr.first=0.994; corr.second=TMath::Sqrt(0.002*0.003+0.015*0.015); }
      else if (leta<1.2) { corr.first=0.994; corr.second=TMath::Sqrt(0.003*0.003+0.021*0.021); }
      else if (leta<2.1) { corr.first=0.989; corr.second=TMath::Sqrt(0.002*0.002+0.015*0.015); }
      else               { corr.first=0.977; corr.second=TMath::Sqrt(0.009*0.009+0.014*0.014); } //0 jets has too large error
    }
    
  }

  return corr;
}



//
EffCorrection_t EfficiencyScaleFactorsWrapper::getTriggerCorrection(std::vector<Particle> leptons, 
                                                                    std::vector<Particle> photons,
                                                                    std::vector<Particle> jets,
                                                                    TString period)
{
  EffCorrection_t corr(1.0,0.0);
  if(era_==2017)
    {
      if(leptons.size()==1)
        {
          TString hname(abs(leptons[0].id())==11 ? "e_trig" : "m_trig");          
          if( abs(leptons[0].id())==13 and scaleFactorsH_.find(hname)!=scaleFactorsH_.end() )
            {
              TH1 *h=scaleFactorsH_[hname];
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
          else if( abs(leptons[0].id())==11 and scaleFactorsH_.find(hname)!=scaleFactorsH_.end() )
            {
              TH1 *h=scaleFactorsH_[hname];
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
EffCorrection_t EfficiencyScaleFactorsWrapper::getOfflineCorrection(Particle p,TString period)
{
  return getOfflineCorrection(p.id(),p.pt(),p.eta(),period);
}

//
EffCorrection_t EfficiencyScaleFactorsWrapper::getOfflineCorrection(int pdgId,float pt,float eta,TString period)
{
  EffCorrection_t corr(1.0,0.01);

  TString corrSteps[]={"rec","tk","id","iso"};
  for(size_t icor=0; icor<sizeof(corrSteps)/sizeof(TString); icor++) {

    TString idstr("e");
    if(abs(pdgId)==13) idstr="m";
    if(abs(pdgId)==22) idstr="g";
    TString hname(idstr+"_"+corrSteps[icor]);
    if(abs(pdgId)==13) hname += period;
    
    Float_t iSF(1.0), iSFUnc(0.0);
    if(scaleFactorsGr_.find(hname)!=scaleFactorsGr_.end() ) {
      Double_t x(0.),xdiff(9999.),y(0.);
      for(Int_t ip=0; ip<scaleFactorsGr_[hname]->GetN(); ip++)
        {
          scaleFactorsGr_[hname]->GetPoint(ip,x,y);
          float ixdiff(TMath::Abs(fabs(eta)-x));
          if(ixdiff>xdiff) continue;
          xdiff=ixdiff;
          iSF=y;
          iSFUnc=scaleFactorsGr_[hname]->GetErrorY(ip);
        }
    }
    else if(scaleFactorsH_.find(hname)!=scaleFactorsH_.end() ) {

      TH2 *h=scaleFactorsH_[hname];
      Double_t xval(eta), yval(pt);
      if(idstr=="m") {
        if(era_==2017) {
          xval=pt;
          yval=fabs(eta);
        }
      }
      xval=TMath::Max( TMath::Min(xval,h->GetXaxis()->GetXmax()-0.001), h->GetXaxis()->GetXmin() );
      Int_t xbin=h->GetXaxis()->FindBin(xval);
      yval=TMath::Max( TMath::Min(yval,h->GetYaxis()->GetXmax()-0.001), h->GetYaxis()->GetXmin() );
      Int_t ybin=h->GetYaxis()->FindBin(yval);      

      iSF=h->GetBinContent(xbin,ybin);
      iSFUnc=h->GetBinError(xbin,ybin);
    }
   
    corr.second = sqrt(pow(iSFUnc*corr.first,2)+pow(iSF*corr.second,2));
    corr.first  = corr.first*iSF;
  }
     
  //
  return corr;
}

//
EfficiencyScaleFactorsWrapper::~EfficiencyScaleFactorsWrapper()
{
}
