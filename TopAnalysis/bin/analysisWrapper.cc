#include <iostream>

#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"
#include "TopLJets2015/TopAnalysis/interface/ExclusiveTop.h"
#include "TopLJets2015/TopAnalysis/interface/VBFVectorBoson.h"

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
       << "\t --systVar  - specify single systematic variation" << endl
       << "\t --era      - era directory to use for corrections, uncertainties" << endl
       << "\t --normTag  - normalization tag" << endl
       << "\t --method   - method to run" << endl;
}

//
int main(int argc, char* argv[])
{
  //get input arguments
  TString in(""),out(""),era(""),normTag(""),method("");
  std::string systVar("");
  bool runSysts(false);
  bool debug(false);
  int channel(0),charge(0),flav(0);
  for(int i=1;i<argc;i++){
    string arg(argv[i]);
    if(arg.find("--help") !=string::npos)                     { printHelp(); return -1;} 
    else if(arg.find("--runSysts")!=string::npos )            { runSysts=true;  }
    else if(arg.find("--systVar")!=string::npos && i+1<argc)  { systVar=argv[i+1]; i++;}
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

  if(debug) cout << "Debug mode is active, runSysts=" << runSysts << endl;

  //open normalization file
  TH1F *normH=0;
  TH1F *puH=0;
  TFile *normF=TFile::Open(era+"/genweights.root");
  if(normF)
    {
      normH=(TH1F *)normF->Get(normTag);
      if(normH) 
        normH->SetDirectory(0);
      puH=(TH1F *)normF->Get(normTag+"_pu");
      if(puH)
        puH->SetDirectory(0);    
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
  if(method=="ExclusiveTop::RunExclusiveTop")          RunExclusiveTop(in,out,channel,charge,normH,puH,era,debug);
  else if(method=="VBFVectorBoson::RunVBFVectorBoson") RunVBFVectorBoson(in,out,channel,charge,normH,puH,era,debug);
  else {
      cout << "Check method=" << method <<endl;
      printHelp();
      return -1;
    }

  //all done
  return 0;
}  





