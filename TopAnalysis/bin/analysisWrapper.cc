#include <iostream>

#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"
#include "TopLJets2015/TopAnalysis/interface/TOP-16-006.h"
#include "TopLJets2015/TopAnalysis/interface/TOP-16-019.h"
#include "TopLJets2015/TopAnalysis/interface/TOP-16-023.h"
#include "TopLJets2015/TopAnalysis/interface/TOP-UE.h"
#include "TopLJets2015/TopAnalysis/interface/TopCharmedMesonAnalysis.h"
#include "TopLJets2015/TopAnalysis/interface/TOPJetShape.h"

#include "TH1F.h"
#include "TFile.h"

using namespace std;

//
void printHelp()
{
  cout << "analysisWrapper options are:" << endl
       << "\t --in - input file" << endl
       << "\t --out - output file" << endl
       << "\t --channel - channel to analyze" << endl
       << "\t --charge  - charge selection to apply" << endl
       << "\t --flav    - flavour selection to apply" << endl
       << "\t --runSysts - activate running systematics" << endl
       << "\t --era      - era directory to use for corrections, uncertainties" << endl
       << "\t --normTag  - normalization tag" << endl
       << "\t --method   - method to run" << endl;
}

//
int main(int argc, char* argv[])
{
  //get input arguments
  TString in(""),out(""),era(""),normTag(""),method("");
  bool runSysts(false),debug(false);
  int channel(0),charge(0),flav(0);
  for(int i=1;i<argc;i++){
    string arg(argv[i]);
    if(arg.find("--help") !=string::npos)                     { printHelp(); return -1;} 
    else if(arg.find("--runSysts")!=string::npos )            { runSysts=true;  }
    else if(arg.find("--channel")!=string::npos && i+1<argc)  { sscanf(argv[i+1],"%d",&channel); i++;}
    else if(arg.find("--charge")!=string::npos && i+1<argc)   { sscanf(argv[i+1],"%d",&charge); i++;}
    else if(arg.find("--flav")!=string::npos && i+1<argc)     { sscanf(argv[i+1],"%d",&flav); i++;}
    else if(arg.find("--in")!=string::npos && i+1<argc)       { in=argv[i+1]; i++;}
    else if(arg.find("--out")!=string::npos && i+1<argc)      { out=argv[i+1]; i++;}
    else if(arg.find("--debug")!=string::npos)                { debug=true; }
    else if(arg.find("--normTag")!=string::npos && i+1<argc)  { normTag=argv[i+1]; i++;}
    else if(arg.find("--era")!=string::npos && i+1<argc)      { era=argv[i+1]; i++;}
    else if(arg.find("--method")!=string::npos && i+1<argc)   { method=argv[i+1]; i++;}
  }

  //open normalization file
  TH1F *normH=0;
  TFile *normF=TFile::Open(era+"/genweights.root");
  if(normF)
    {
      normH=(TH1F *)normF->Get(normTag);
      if(normH) {
        normH->SetDirectory(0);
        normH->SetTitle(normTag);
      }
      normF->Close();
    }
  if(normH==0)
    {
      cout << "Check normalization file genweights.root in era=" << era 
	   << " and tag (" << normTag << ")" << endl
	   << "Will run without any" << endl;
      printHelp();
      //return -1;
    }
  
  //check input/output
  if(in=="" || out=="")
    {
      cout << "Check input/output=" << in << "/" << out << endl;
      printHelp();
      return -1;
    }

  //check method to run
  if(method=="TOP-16-006::RunTop16006")                                  RunTop16006(in,out,channel,charge,FlavourSplitting(flav),normH,runSysts,era);
  else if(method=="TOP-16-019::RunTop16019")                             RunTop16019(in,out,channel,charge,FlavourSplitting(flav),normH,runSysts,era);
  else if(method=="TOP-16-023::RunTop16023")                             RunTop16023(in,out,channel,charge,FlavourSplitting(flav),normH,runSysts,era);
  else if(method=="TopCharmedMesonAnalysis::RunTopCharmedMesonAnalysis") RunTopCharmedMesonAnalysis(in,out,channel,charge,FlavourSplitting(flav),normH,runSysts,era,debug);
  else if(method=="TOP-UE::RunTopUE")                                    RunTopUE(in,out,channel,charge,FlavourSplitting(flav),normH,runSysts,era);
  else if(method=="TOPJetShape::RunTopJetShape")                         RunTopJetShape(in,out,channel,charge,FlavourSplitting(flav),normH,runSysts,era);
  else
    {
      cout << "Check method=" << method <<endl;
      printHelp();
      return -1;
    }

  //all done
  return 0;
}  





