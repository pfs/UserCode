#include "TopLJets2015/TopAnalysis/interface/FRcalculationVBF.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/TMVARegGui.h"
 
using namespace std;
using namespace TMVA;
int main( int argc, char** argv )
{
  TString fGammaData("Data13TeV_SinglePhoton_2016.root"), fJetData("Data13TeV_JetHT_2016.root"), fJetQCD("Data13TeV_JetHTQCD_2016.root"), fGammaMC("MC13TeV_GJets.root");
  std::string categories;
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
    } else if( input == "--cats"){
      i++;
      categories = std::string(argv[i]);
    }
  }
  std::vector<TString> cats = TMVA::gTools().SplitString( categories, ':' );
  for(unsigned int icat = 0; icat < cats.size(); icat++){
    cout << cats[icat]<<endl;
    //HighMJJ: bias from VBF trigger --> use jet data
    if(cats[icat].Contains("LowVPtHighMJJ")){ 
      promptEstimator(fJetData,cats[icat], "EB", fGammaMC, fJetQCD);
      continue;
    }
    promptEstimator(fGammaData,cats[icat], "EB", fGammaMC, fJetQCD);
    promptEstimator(fGammaData,cats[icat], "EE", fGammaMC, fJetQCD);
  }


  makeFile();
  return 0;
}



