#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"
#include "TSystem.h"


//
Float_t computeMT(TLorentzVector &a, TLorentzVector &b)
{
  return TMath::Sqrt(2*a.Pt()*b.Pt()*(1-TMath::Cos(a.DeltaPhi(b))));
}

//
std::map<Int_t,Float_t> lumiPerRun()
{
  std::map<Int_t,Float_t> lumiMap;
  lumiMap[254231]=   32626.916   ;
  lumiMap[254232]=   108539.723  ;
  lumiMap[254790]=  11333305.274 ;
  lumiMap[254852]=   898963.031  ;
  lumiMap[254879]=  1798475.919  ;
  lumiMap[254906]=  1565208.707  ;
  lumiMap[254907]=  1070993.080  ;
  lumiMap[254914]=   923324.411  ;
  lumiMap[256630]=  1019427.537  ;
  lumiMap[256673]=    5821.004   ;
  lumiMap[256674]=   97107.612   ;
  lumiMap[256675]=  7631339.155  ;
  lumiMap[256676]=  9586678.621  ;
  lumiMap[256677]=  16208124.083 ;
  lumiMap[256801]=  9289181.921  ;
  lumiMap[256842]=   17564.969   ;
  lumiMap[256843]=  39192996.677 ;
  lumiMap[256866]=   60179.433   ;
  lumiMap[256867]=  4778327.656  ;
  lumiMap[256868]=  23626060.836 ;
  lumiMap[256869]=  1613257.519  ;
  lumiMap[256926]=  1585513.104  ;
  lumiMap[256941]=  9153369.805  ;
  lumiMap[257461]=  3273371.101  ;
  lumiMap[257531]=  8952857.360  ;
  lumiMap[257599]=  5277913.939  ;
  lumiMap[257613]=  80288701.786 ;
  lumiMap[257614]=   898910.938  ;
  lumiMap[257645]=  66251074.200 ;
  lumiMap[257682]=  14059859.130 ;
  lumiMap[257722]=   874139.924  ;
  lumiMap[257723]=  6416461.542  ;
  lumiMap[257735]=   576143.428  ;
  lumiMap[257751]=  28892223.256 ;
  lumiMap[257804]=   225829.957  ;
  lumiMap[257805]=  18191777.239 ;
  lumiMap[257816]=  25831347.642 ;
  lumiMap[257819]=  16070065.308 ;
  lumiMap[257968]=  17947956.702 ;
  lumiMap[257969]=  41437510.749 ;
  lumiMap[258129]=  6161039.580  ;
  lumiMap[258136]=  3833715.336  ;
  lumiMap[258157]=  4130426.007  ;
  lumiMap[258158]= 112150208.043 ;
  lumiMap[258159]=  27041879.753 ;
  lumiMap[258177]= 112357734.179 ;
  lumiMap[258211]=  6899616.879  ;
  lumiMap[258213]=  12447784.863 ;
  lumiMap[258214]=  16299123.425 ;
  lumiMap[258215]=   443760.789  ;
  lumiMap[258287]=  14271300.581 ;
  lumiMap[258403]=  16554699.075 ;
  lumiMap[258425]=  10948640.280 ;
  lumiMap[258426]=   808721.923  ;
  lumiMap[258427]=  8497851.929  ;
  lumiMap[258428]=  12440664.974 ;
  lumiMap[258432]=   298695.064  ;
  lumiMap[258434]=  32645147.197 ;
  lumiMap[258440]=  47654602.747 ;
  lumiMap[258444]=  2208821.299  ;
  lumiMap[258445]=  17379231.195 ;
  lumiMap[258446]=  7906567.040  ;
  lumiMap[258448]=  37636207.590 ;
  lumiMap[258655]=   412374.500  ;
  lumiMap[258656]=  27561949.634 ;
  lumiMap[258694]=  16613108.138 ;
  lumiMap[258702]=  31593447.906 ;
  lumiMap[258703]=  33749411.575 ;
  lumiMap[258705]=  8215733.522  ;
  lumiMap[258706]=  56015291.210 ;
  lumiMap[258712]=  36912048.837 ;
  lumiMap[258713]=  10868729.417 ;
  lumiMap[258714]=  4462940.479  ;
  lumiMap[258741]=  4899047.520  ;
  lumiMap[258742]=  65372682.457 ;
  lumiMap[258745]=  22816248.664 ;
  lumiMap[258749]=  48011842.080 ;
  lumiMap[258750]=  15311166.469 ;
  lumiMap[259626]=  11503939.036 ;
  lumiMap[259637]=  15843833.799 ;
  lumiMap[259681]=  2006428.466  ;
  lumiMap[259683]=  7733152.101  ;
  lumiMap[259685]=  55748876.683 ;
  lumiMap[259686]=  27125232.494 ;
  lumiMap[259721]=  12400448.429 ;
  lumiMap[259809]=  14370193.633 ;
  lumiMap[259810]=  9903086.201  ;
  lumiMap[259811]=  7470396.336  ;
  lumiMap[259813]=   746162.774  ;
  lumiMap[259817]=   362610.422  ;
  lumiMap[259818]=  13130237.492 ;
  lumiMap[259820]=  12560062.290 ;
  lumiMap[259821]=  16180451.962 ;
  lumiMap[259822]=  32721804.046 ;
  lumiMap[259861]=  6561060.297  ;
  lumiMap[259862]=  45860217.938 ;
  lumiMap[259884]=  6731111.093  ;
  lumiMap[259890]=  9701207.990  ;
  lumiMap[259891]=  9603195.320  ;
  lumiMap[260373]=  10920147.469 ;
  lumiMap[260424]=  66688251.029 ;
  lumiMap[260425]=  23599504.405 ;
  lumiMap[260426]=  43930543.476 ;
  lumiMap[260427]=  15969446.707 ;
  lumiMap[260431]=  35126694.498 ;
  lumiMap[260532]=  69073584.559 ;
  lumiMap[260533]=  1195476.609  ;
  lumiMap[260534]=  32043973.431 ;
  lumiMap[260536]=  14466413.325 ;
  lumiMap[260538]=  22368836.359 ;
  lumiMap[260541]=  1829959.151  ;
  lumiMap[260575]=  1721667.572  ;
  lumiMap[260576]=  16664531.028 ;
  lumiMap[260577]=  8251536.906  ;
  lumiMap[260593]=  35893405.704 ;
  lumiMap[260627]= 178937353.997 ;

  return lumiMap;
};

// 
float getLeptonEnergyScaleUncertainty(int l_id,float l_pt,float l_eta)
{
  float unc(0.02);
  
  // electron uncertainties for 8 TeV cf. AN-14-145   
  if(abs(l_id)==11 || abs(l_id)==1100 || abs(l_id)==2111)
    {
      float par0(-2.27e-02), par1(-7.01e-02), par2(-3.71e-04);
      if (fabs(l_eta) > 0.8 && fabs(l_eta)<1.5)
        {
          par0 = -2.92e-02;
          par1 = -6.59e-02;
          par2 = -7.22e-04;
        }
      else if(fabs(l_eta)>1.5)
        {
          par0 = -2.27e-02;
          par1 = -7.01e-02;
          par2 = -3.71e-04;
        }
      unc=fabs(par0 * TMath::Exp(par1 * l_pt) + par2);
    }

  return unc;
}


//
FactorizedJetCorrector *getFactorizedJetEnergyCorrector(TString baseDir, bool isMC)
{
  gSystem->ExpandPathName(baseDir);
  
  //order matters: L1 -> L2 -> L3 (-> Residuals)
  std::vector<std::string> jetCorFiles;
  TString pf("Summer15_25nsV7_");
  pf += (isMC ? "MC" : "DATA");
  std::cout << "Jet energy corrections from " << baseDir+"/"+pf+"_*_AK4PFchs.txt" << std::endl;
  jetCorFiles.push_back((baseDir+"/"+pf+"_L1FastJet_AK4PFchs.txt").Data());
  jetCorFiles.push_back((baseDir+"/"+pf+"_L2Relative_AK4PFchs.txt").Data());
  jetCorFiles.push_back((baseDir+"/"+pf+"_L3Absolute_AK4PFchs.txt").Data());
  if(!isMC) jetCorFiles.push_back((baseDir+"/"+pf+"_L2L3Residual_AK4PFchs.txt").Data());
  
  //init the parameters for correction
  std::vector<JetCorrectorParameters> corSteps;
  for(size_t i=0; i<jetCorFiles.size(); i++) corSteps.push_back(JetCorrectorParameters(jetCorFiles[i]));
  
  //return the corrector
  return new FactorizedJetCorrector(corSteps);
}


//Sources :  https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution
std::vector<float> getJetResolutionScales(float pt, float eta, float genjpt)
{
  std::vector<float> res(3,1.0);

  float ptSF(1.0), ptSF_err(0.0);
  if(TMath::Abs(eta)<0.8)       { ptSF=1.061; ptSF_err = 0.023; }
  else if(TMath::Abs(eta)<1.3)  { ptSF=1.088; ptSF_err = 0.029; }
  else if(TMath::Abs(eta)<1.9)  { ptSF=1.106; ptSF_err = 0.030; }
  else if(TMath::Abs(eta)<2.5)  { ptSF=1.126; ptSF_err = 0.094; }
  else if(TMath::Abs(eta)<3.0)  { ptSF=1.343; ptSF_err = 0.123; }
  else if(TMath::Abs(eta)<3.2)  { ptSF=1.303; ptSF_err = 0.111; }
  else if(TMath::Abs(eta)<5.0)  { ptSF=1.320; ptSF_err = 0.286; }

  res[0] = TMath::Max((Float_t)0.,(Float_t)(genjpt+(ptSF)*(pt-genjpt)))/pt;
  res[1] = TMath::Max((Float_t)0.,(Float_t)(genjpt+(ptSF-ptSF_err)*(pt-genjpt)))/pt;
  res[2] = TMath::Max((Float_t)0.,(Float_t)(genjpt+(ptSF+ptSF_err)*(pt-genjpt)))/pt;
  
  return res;
}
