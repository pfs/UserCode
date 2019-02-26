#include <iostream>

#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"
#include "TopLJets2015/TopAnalysis/interface/ExclusiveTop.h"
#include "TopLJets2015/TopAnalysis/interface/TOP-17-010.h"
#include "TopLJets2015/TopAnalysis/interface/VBFVectorBoson.h"
#include "TopLJets2015/TopAnalysis/interface/PhotonTrigEff.h"

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
       << "\t --flag    - job flag to apply" << endl
       << "\t --runSysts - activate running systematics" << endl
       << "\t --systVar  - specify single systematic variation" << endl
       << "\t --era      - era directory to use for corrections, uncertainties" << endl
       << "\t --normTag  - normalization tag" << endl
       << "\t --method   - method to run" << endl
       << "\t --CR       - make CR for fake rate, based on pu jet id" << endl
       << "\t --QCDTemp  - can be true only if is CR. Used to make QCD templates for fake photons" << endl
       << "\t --SRfake   - makes the region on which the fake ratio should be applied" << endl
       << "\t --mvatree  - store selected events in a tree for mva" << endl;
}

//
int main(int argc, char* argv[])
{
  //get input arguments
  TString in(""),out(""),era(""),normTag(""),method(""),genWeights("genweights.root");
  std::string systVar("");
  bool runSysts(false);
  bool debug(false), skimtree(false), CR(false), QCDTemp(false), SRfake(false);
  int channel(0),charge(0),flag(-1);
  float xsec(1.0);
  for(int i=1;i<argc;i++){
    string arg(argv[i]);
    if(arg.find("--help") !=string::npos)                       { printHelp(); return -1;} 
    else if(arg.find("--runSysts")!=string::npos )              { runSysts=true;  }
    else if(arg.find("--systVar")!=string::npos && i+1<argc)    { systVar=argv[i+1]; i++;}
    else if(arg.find("--channel")!=string::npos && i+1<argc)    { sscanf(argv[i+1],"%d",&channel); i++;}
    else if(arg.find("--charge")!=string::npos && i+1<argc)     { sscanf(argv[i+1],"%d",&charge); i++;}
    else if(arg.find("--flag")!=string::npos && i+1<argc)       { sscanf(argv[i+1],"%d",&flag); i++;}
    else if(arg.find("--in")!=string::npos && i+1<argc)         { in=argv[i+1]; i++;}
    else if(arg.find("--out")!=string::npos && i+1<argc)        { out=argv[i+1]; i++;}
    else if(arg.find("--debug")!=string::npos)                  { debug=true; }
    else if(arg.find("--CR")!=string::npos)                     { CR=true; }
    else if(arg.find("--QCDTemp")!=string::npos)                { QCDTemp=true; }
    else if(arg.find("--SRfake")!=string::npos)                 { SRfake=true; }
    else if(arg.find("--mvatree")!=string::npos)                { skimtree=true; }
    else if(arg.find("--normTag")!=string::npos && i+1<argc)    { normTag=argv[i+1]; i++;}
    else if(arg.find("--era")!=string::npos && i+1<argc)        { era=argv[i+1]; i++;}
    else if(arg.find("--method")!=string::npos && i+1<argc)     { method=argv[i+1]; i++;}
    else if(arg.find("--genWeights")!=string::npos && i+1<argc) { genWeights=argv[i+1]; i++;}
    else if(arg.find("--xsec")!=string::npos && i+1<argc)       { sscanf(argv[i+1],"%f",&xsec); i++;}
  }

  if(debug) cout << "Debug mode is active, runSysts=" << runSysts << endl;

  //open normalization file
  TH1F *normH=0;
  TH1F *puH=0;
  TFile *normF=TFile::Open(Form("%s/%s",era.Data(),genWeights.Data()));
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
      cout << "Check normalization file " << genWeights << " in era=" << era 
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
  if(method=="ExclusiveTop::RunExclusiveTop") {
    RunExclusiveTop(in,out,channel,charge,normH,puH,era,debug);
  }
  else if(method=="PhotonTrigEff::RunPhotonTrigEff") {
    RunPhotonTrigEff(in,out,normH,puH,era,debug);
  }
  else if(method=="VBFVectorBoson::RunVBFVectorBoson") {
    VBFVectorBoson myVBF(in,out,normH,puH,era,xsec,debug,CR,QCDTemp,SRfake,skimtree,true);
    myVBF.runAnalysis();
  }
  else if(method=="TOP17010::TOP17010") {
    TOP17010 myTOP17010(in,out,normH,puH,era,flag,debug);
    myTOP17010.runAnalysis();
  }
  else {
    cout << "Check method=" << method <<endl;
    printHelp();
    return -1;
  }
  
  //all done
  return 0;
}  





