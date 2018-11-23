#include "TopLJets2015/TopAnalysis/interface/FRcalculationVBF.h"
int main( int argc, char** argv )
{
  TString fGammaData("Data13TeV_SinglePhoton_2017.root"), fJetData("Data13TeV_JetHT_2017.root"), fJetQCD("Data13TeV_JetHTQCD_2017.root"), fGammaMC("MC13TeV_GJets.root");
  for (int i=1; i<argc; i++) {
    TString input(argv[i]);
    if ( input=="--fGdata"){
      i++;
      fGammaData = TString(argv[i]);
      continue;
    } else if ( input=="--fJdata"){
      i++;
      fJetData = TString(argv[i]);
      continue;
    } else if( input=="--fJQCD"){
      i++;
      fJetQCD = TString(argv[i]);
      continue;
    } else if( input=="--fGMC"){
      i++;
      fGammaMC = TString(argv[i]);
      continue;
    } 
  }
  //HighMJJ: bias from VBF trigger --> use jet data
  promptEstimator(fJetData,"HighMJJA", "EB", fGammaMC, fJetQCD);

  //LowMJJ: barrel
  promptEstimator(fGammaData,"LowMJJA", "EB", fGammaMC, fJetQCD);
  //LowMJJ: endcap
  promptEstimator(fGammaData,"LowMJJA", "EE", fGammaMC, fJetQCD);

  makeFile();
  return 0;
}



