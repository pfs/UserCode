#include "TopLJets2015/TopAnalysis/interface/WorkspaceProvider.h"
#include "TopLJets2015/TopAnalysis/interface/VbfFitRegion.h"
using namespace std;

int main( int argc, char** argv )
{
  TString channel, histname, boson;
  YieldsErr YE;
  double sigEff =1;
  double bkgEff =1;
  int nBin = 4;
  bool doSignalPH = false;
  for (int i=1; i<argc; i++) {
    TString input(argv[i]);
    if ( input=="--Chan"){
      i++;
      channel = TString(argv[i]);
      continue;
    } else if ( input=="--V"){
      i++;
      boson = TString(argv[i]);
      continue;
    } else if( input=="--Hist"){
      i++;
      histname = TString(argv[i]);
      continue;
    } else if( input=="--YieldErr"){
      i++;
      YE = splitter(string(argv[i]),':',',');     
      continue;
    } else if( input=="--sigEff"){
      i++;
      sigEff = (double)atof(TString(argv[i]));
      continue;
    } else if( input=="--bkgEff"){
      i++;
      bkgEff = (double)atof(TString(argv[i]));
      continue;
    } else if( input=="--nBin"){
      i++;
      nBin = (int)atof(TString(argv[i]));
      continue;
    } else if (input=="--doSignalPH"){
      doSignalPH = true;
      continue;
    }
  }


  VbfFitRegion * SR = new VbfFitRegion(channel, TString("A"), histname, nBin, true);
  VbfFitRegion * CR = new VbfFitRegion(channel, TString("MM"), histname, nBin, false);
  
  WorkspaceProvider wsp(histname,SR, CR);
  wsp.import(doSignalPH);
  wsp.makeCard(YE, TString("A"), doSignalPH, sigEff, bkgEff);
  wsp.makeCard(YE, TString("MM"), doSignalPH, sigEff, bkgEff);
  wsp.makeCardNLO(YE, TString("A"));
  return 0;
}
