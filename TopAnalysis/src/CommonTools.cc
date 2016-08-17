#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"
#include "TSystem.h"


//
JetPullInfo_t getPullVector( MiniEvent_t &ev, int ijet)
{
  JetPullInfo_t result;
  result.n=0; result.nch=0;
  result.pull=TVector2(0,0);
  result.chPull=TVector2(0,0);
  
  //re-reconstruct the jet direction with the charged tracks
  TLorentzVector jet(0,0,0,0);
  jet.SetPtEtaPhiM(ev.j_pt[ijet], ev.j_eta[ijet], ev.j_phi[ijet], ev.j_mass[ijet]);
  TLorentzVector chargedJet(0,0,0,0);
  TLorentzVector constituent(0,0,0,0);
  std::vector<std::pair<TLorentzVector,bool> > allConstituents;
  unsigned int nCharged = 0;
  for(Int_t idx = 0; idx<ev.npf; idx++)
    {
      if(ev.pf_j[idx]!=ijet) continue;
      constituent.SetPtEtaPhiM( ev.pf_pt[idx], ev.pf_eta[idx], ev.pf_phi[idx], ev.pf_m[idx]);
      bool isCharged(abs(ev.pf_id[idx])==11 ||
		     abs(ev.pf_id[idx])==13 ||
		     abs(ev.pf_id[idx])==211 );
      allConstituents.push_back(std::make_pair(constituent,isCharged) );
      if(isCharged)
	{
	  chargedJet += constituent;
	  ++nCharged;      
	}
    }
  result.n=(Int_t) allConstituents.size();
  result.nch=nCharged;

  //stop here if <2 charged
  if( nCharged < 2 ) return result;

  //compute the pull
  double jetPt        = jet.Pt(),        jetPhi=jet.Phi(),                 jetRapidity=jet.Rapidity();
  double jetPtCharged = chargedJet.Pt(), jetPhiCharged = chargedJet.Phi(), jetRapidityCharged = chargedJet.Rapidity();
  TVector2 r(0,0);
  TVector2 pullAll(0,0);
  TVector2 pullCharged(0,0);
  for(size_t idx = 0; idx<allConstituents.size(); ++idx)
    {
      TLorentzVector &cp4=allConstituents[idx].first;
      bool &isCharged=allConstituents[idx].second;
      double constituentPt       = cp4.Pt();
      double constituentPhi      = cp4.Phi();
      double constituentRapidity = cp4.Rapidity();
      r.Set( constituentRapidity - jetRapidity, TVector2::Phi_mpi_pi( constituentPhi - jetPhi ) );
      pullAll += ( constituentPt / jetPt ) * r.Mod() * r;
      //calculate TVector using only charged tracks
      if( isCharged )
	r.Set( constituentRapidity - jetRapidityCharged, TVector2::Phi_mpi_pi( constituentPhi - jetPhiCharged ) );
      pullCharged += ( constituentPt / jetPtCharged ) * r.Mod() * r;
    }
  
  result.pull=pullAll;
  result.chPull=pullCharged;
  return result;
}


//
Float_t computeMT(TLorentzVector &a, TLorentzVector &b)
{
  return TMath::Sqrt(2*a.Pt()*b.Pt()*(1-TMath::Cos(a.DeltaPhi(b))));
}

//
std::map<Int_t,Float_t> lumiPerRun(TString era)
{
  std::map<Int_t,Float_t> lumiMap;
  if(era=="era2015")
    {
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
    }
  else
    {
      lumiMap[ 273158 ] = 51635046.420  ;
      lumiMap[ 273302 ] = 16819239.527  ;
      lumiMap[ 273402 ] = 6792207.559   ;
      lumiMap[ 273403 ] = 1878485.596   ;
      lumiMap[ 273404 ] = 630136.094    ;
      lumiMap[ 273405 ] = 824379.492    ;
      lumiMap[ 273406 ] = 3887742.993   ;
      lumiMap[ 273408 ] = 185921.869    ;
      lumiMap[ 273409 ] = 9725934.152   ;
      lumiMap[ 273410 ] = 2817437.075   ;
      lumiMap[ 273411 ] = 759189.735    ;
      lumiMap[ 273425 ] = 27512155.549  ;
      lumiMap[ 273446 ] = 2202659.483   ;
      lumiMap[ 273447 ] = 26659723.187  ;
      lumiMap[ 273448 ] = 23333828.306  ;
      lumiMap[ 273449 ] = 11742205.666  ;
      lumiMap[ 273450 ] = 33190912.201  ;
      lumiMap[ 273492 ] = 12517186.413  ;
      lumiMap[ 273493 ] = 10691190.039  ;
      lumiMap[ 273494 ] = 8704028.577   ;
      lumiMap[ 273502 ] = 48296078.647  ;
      lumiMap[ 273503 ] = 26175295.146  ;
      lumiMap[ 273554 ] = 27628360.593  ;
      lumiMap[ 273555 ] = 12597573.412  ;
      lumiMap[ 273725 ] = 159616078.734 ;
      lumiMap[ 273728 ] = 5010364.712   ;
      lumiMap[ 273730 ] = 89451996.572  ;
      lumiMap[ 274094 ] = 10059427.765  ;
      lumiMap[ 274146 ] = 5271396.948   ;
      lumiMap[ 274159 ] = 3668736.593   ;
      lumiMap[ 274160 ] = 18560369.363  ;
      lumiMap[ 274161 ] = 43482410.731  ;
      lumiMap[ 274172 ] = 6380273.172   ;
      lumiMap[ 274198 ] = 12694886.888  ;
      lumiMap[ 274199 ] = 69534024.918  ;
      lumiMap[ 274200 ] = 66425877.607  ;
      lumiMap[ 274240 ] = 10328978.619  ;
      lumiMap[ 274241 ] = 137315787.867 ;
      lumiMap[ 274244 ] = 59427904.348  ;
      lumiMap[ 274250 ] = 97146811.386  ;
      lumiMap[ 274251 ] = 67084377.583  ;
      lumiMap[ 274283 ] = 358761.932    ;
      lumiMap[ 274284 ] = 30760422.968  ;
      lumiMap[ 274286 ] = 4337739.639   ;
      lumiMap[ 274314 ] = 9756511.369   ;
      lumiMap[ 274315 ] = 68207093.787  ;
      lumiMap[ 274316 ] = 134065238.812 ;
      lumiMap[ 274317 ] = 310073.517    ;
      lumiMap[ 274319 ] = 28103257.785  ;
      lumiMap[ 274335 ] = 146491732.457 ;
      lumiMap[ 274336 ] = 1790263.505   ;
      lumiMap[ 274337 ] = 1831090.279   ;
      lumiMap[ 274338 ] = 91493706.851  ;
      lumiMap[ 274339 ] = 10377782.305  ;
      lumiMap[ 274344 ] = 71995570.263  ;
      lumiMap[ 274345 ] = 18266582.622  ;
      lumiMap[ 274382 ] = 7622787.273   ;
      lumiMap[ 274387 ] = 59660679.277  ;
      lumiMap[ 274388 ] = 261038800.403 ;
      lumiMap[ 274420 ] = 28799015.579  ;
      lumiMap[ 274421 ] = 55267231.822  ;
      lumiMap[ 274422 ] = 284957505.173 ;
      lumiMap[ 274440 ] = 65179944.181  ;
      lumiMap[ 274441 ] = 65643264.358  ;
      lumiMap[ 274442 ] = 101738861.928 ;
      lumiMap[ 274954 ] = 990679.863    ;
      lumiMap[ 274955 ] = 4903780.935   ;
      lumiMap[ 274968 ] = 198787386.103 ;
      lumiMap[ 274969 ] = 130136824.076 ;
      lumiMap[ 274970 ] = 4426188.047   ;
      lumiMap[ 274971 ] = 100582997.672 ;
      lumiMap[ 274998 ] = 115887069.133 ;
      lumiMap[ 274999 ] = 167346833.683 ;
      lumiMap[ 275000 ] = 16276448.639  ;
      lumiMap[ 275001 ] = 203827216.709 ;
      lumiMap[ 275059 ] = 6237815.233   ;
      lumiMap[ 275066 ] = 18717462.852  ;
      lumiMap[ 275067 ] = 72728810.479  ;
      lumiMap[ 275068 ] = 145890403.924 ;
      lumiMap[ 275073 ] = 69945928.163  ;
      lumiMap[ 275074 ] = 78014694.422  ;
      lumiMap[ 275124 ] = 57282239.944  ;
      lumiMap[ 275125 ] = 153206070.148 ;
      lumiMap[ 275282 ] = 15120897.318  ;
      lumiMap[ 275283 ] = 22765754.482  ;
      lumiMap[ 275284 ] = 12635751.870  ;
      lumiMap[ 275290 ] = 2724125.194   ;
      lumiMap[ 275291 ] = 52618068.750  ;
      lumiMap[ 275292 ] = 15968243.944  ;
      lumiMap[ 275293 ] = 25627096.583  ;
      lumiMap[ 275309 ] = 84180084.990  ;
      lumiMap[ 275310 ] = 228559917.513 ;
      lumiMap[ 275311 ] = 110613086.238 ;
      lumiMap[ 275319 ] = 24989475.511  ;
      lumiMap[ 275337 ] = 70068728.197  ;
      lumiMap[ 275338 ] = 77062388.059  ;
      lumiMap[ 275344 ] = 51398067.798  ;
      lumiMap[ 275345 ] = 60786304.289  ;
      lumiMap[ 275370 ] = 52765296.013  ;
      lumiMap[ 275371 ] = 96379521.539  ;
      lumiMap[ 275375 ] = 210811736.697 ;
      lumiMap[ 275376 ] = 332911666.521 ;
      lumiMap[ 275657 ] = 20104258.744  ;
      lumiMap[ 275658 ] = 65783229.906  ;
      lumiMap[ 275659 ] = 3242784.743   ;
      lumiMap[ 275761 ] = 1706814.919   ;
      lumiMap[ 275767 ] = 722032.655    ;
      lumiMap[ 275772 ] = 10594354.027  ;
      lumiMap[ 275773 ] = 1209542.133   ;
      lumiMap[ 275774 ] = 56385797.450  ;
      lumiMap[ 275776 ] = 23844819.394  ;
      lumiMap[ 275777 ] = 49161764.584  ;
      lumiMap[ 275778 ] = 46971370.652  ;
      lumiMap[ 275782 ] = 102352560.172 ;
      lumiMap[ 275832 ] = 74243062.618  ;
      lumiMap[ 275833 ] = 47422478.958  ;
      lumiMap[ 275834 ] = 53227355.782  ;
      lumiMap[ 275835 ] = 2121590.895   ;
      lumiMap[ 275836 ] = 189644121.048 ;
      lumiMap[ 275837 ] = 85911906.190  ;
      lumiMap[ 275847 ] = 205316609.057 ;
      lumiMap[ 275886 ] = 7744154.853   ;
      lumiMap[ 275890 ] = 241038016.589 ;
      lumiMap[ 275911 ] = 74878108.405  ;
      lumiMap[ 275912 ] = 53245958.710  ;
      lumiMap[ 275913 ] = 78576586.835  ;
      lumiMap[ 275918 ] = 51039088.804  ;
      lumiMap[ 275920 ] = 62972829.338  ;
      lumiMap[ 275921 ] = 1032540.440   ;
      lumiMap[ 275923 ] = 14915918.861  ;
      lumiMap[ 275931 ] = 10722080.834  ;
      lumiMap[ 275963 ] = 19953501.924  ;
      lumiMap[ 276092 ] = 15083053.906  ;
      lumiMap[ 276097 ] = 101625242.318 ;
      lumiMap[ 276242 ] = 275473847.283 ;
      lumiMap[ 276243 ] = 79349747.546  ;
      lumiMap[ 276244 ] = 136344065.622 ;
      lumiMap[ 276282 ] = 210161541.213 ;
      lumiMap[ 276283 ] = 171192116.106 ;
      lumiMap[ 276315 ] = 34592805.682  ;
      lumiMap[ 276317 ] = 26892256.997  ;
      lumiMap[ 276318 ] = 100279131.747 ;
      lumiMap[ 276355 ] = 6337750.950   ;
      lumiMap[ 276361 ] = 153033538.256 ;
      lumiMap[ 276363 ] = 212058213.970 ;
      lumiMap[ 276384 ] = 125193608.313 ;
      lumiMap[ 276437 ] = 364937772.598 ;
      lumiMap[ 276454 ] = 65824328.594  ;
      lumiMap[ 276458 ] = 34208739.803  ;
      lumiMap[ 276495 ] = 38883684.433  ;
      lumiMap[ 276501 ] = 411670895.729 ;
      lumiMap[ 276502 ] = 87951123.579  ;
      lumiMap[ 276525 ] = 413058119.501 ;
      lumiMap[ 276527 ] = 23093998.094  ;
      lumiMap[ 276528 ] = 39543412.387  ;
      lumiMap[ 276542 ] = 151715866.299 ;
      lumiMap[ 276543 ] = 151987535.534 ;
      lumiMap[ 276544 ] = 22599776.026  ;
      lumiMap[ 276545 ] = 27386707.571  ;
      lumiMap[ 276581 ] = 75636838.564  ;
      lumiMap[ 276582 ] = 153540104.533 ;
      lumiMap[ 276583 ] = 8046310.347   ;
      lumiMap[ 276584 ] = 177293.086    ;
      lumiMap[ 276585 ] = 35841452.444  ;
      lumiMap[ 276586 ] = 99526396.503  ;
      lumiMap[ 276587 ] = 110778975.673 ;
      lumiMap[ 276653 ] = 88382468.945  ;
      lumiMap[ 276655 ] = 165414641.285 ;
      lumiMap[ 276659 ] = 30782867.099  ;
      lumiMap[ 276775 ] = 212273692.403 ;
      lumiMap[ 276776 ] = 242982640.483 ;
      lumiMap[ 276794 ] = 89856377.590  ;
      lumiMap[ 276807 ] = 31086349.677  ;
      lumiMap[ 276808 ] = 161488158.437 ;
      lumiMap[ 276810 ] = 45897433.991  ;
      lumiMap[ 276811 ] = 310487543.429 ;
    }
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
std::vector<TGraph *> getPileupWeights(TString era,TH1 *puTrue)
{
  std::vector<TGraph *>puWgtGr;
  TString puWgtUrl(era+"/pileupWgts.root");
  gSystem->ExpandPathName(puWgtUrl);
  TFile *fIn=TFile::Open(puWgtUrl);
  for(size_t i=0; i<3; i++)
    {
      TString grName("pu_nom");
      if(i==1) grName="pu_down";
      if(i==2) grName="pu_up";
      TGraph *puData=(TGraph *)fIn->Get(grName);
      Float_t totalData=puData->Integral();
      TH1 *tmp=(TH1 *)puTrue->Clone("tmp");
      for(Int_t xbin=1; xbin<=tmp->GetXaxis()->GetNbins(); xbin++)
	{
	  Float_t yexp=puTrue->GetBinContent(xbin);
	  Double_t xobs,yobs;
	  puData->GetPoint(xbin-1,xobs,yobs);
	  tmp->SetBinContent(xbin, yexp>0 ? yobs/(totalData*yexp) : 0. );
	}
      TGraph *gr=new TGraph(tmp);
      grName.ReplaceAll("pu","puwgts");
      gr->SetName(grName);
      puWgtGr.push_back( gr );
      tmp->Delete();
    }
  return puWgtGr;
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