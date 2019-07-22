#include "TopLJets2015/TopAnalysis/interface/WorkspaceProvider.h"
#include "TopLJets2015/TopAnalysis/interface/VbfFitRegion.h"
using namespace std;

int main( int argc, char** argv )
{
  TString channel, histname, boson, year;
  YieldsErr YE;
  double sigEff =1;
  double bkgEff =1;
  int nBin = 4;
  bool shapeOnly = false;
  bool doSignalPH = false;
  bool NLO(false);
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
    } else if (input=="--year"){
      i++;
      year = TString(argv[i]);
      continue;
    } else if (input=="--shapeOnly"){
      shapeOnly = true;
      continue;
    } else if (input=="--nloDefault"){
      NLO = true;
      continue;
    }
  }
  cout <<" ------------ "<<NLO<<endl;
  if(false) cout << doSignalPH << "\t"<<sigEff << "\t" <<bkgEff<<endl;
  VbfFitRegion * SR = new VbfFitRegion(channel, TString("A"), histname, year, nBin, true, shapeOnly,NLO);
  VbfFitRegion * CR = new VbfFitRegion(channel, TString("MM"), histname, year,nBin, false, shapeOnly,NLO);
  
  WorkspaceProvider wsp(histname,SR, CR);
  // wsp.import(doSignalPH);
  //wsp.makeCard(YE, TString("A"), doSignalPH, sigEff, bkgEff);
  // wsp.makeCard(YE, TString("MM"), doSignalPH, sigEff, bkgEff);
  wsp.makeCardNLO(YE, TString("A"));
  wsp.makeCardNLO(YE, TString("A"), "NLOLin");
  wsp.makeCardNLO(YE, TString("A"), "NLOBinned");
  wsp.plotSystSig();  
  return 0;
}
