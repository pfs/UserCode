#include "TopLJets2015/TopAnalysis/interface/TMVAClassification.h"

using namespace std;
using namespace TMVA;

int main( int argc, char** argv )
{
   // Select methods (don't look at this code - not of interest)
  BDTOptimizer* bdt_options;
  bdt_options = NULL;
  TString methodList,ext,category;
  TString DEFAULT_INFNAME  = "";
  for (int i=1; i<argc; i++) {
    TString regMethod(argv[i]);
    cout<<regMethod<<endl;
    if( regMethod=="--ext"){
	i++;
	ext = TString(argv[i]);
	continue;
    }
    if( regMethod=="--vbf"){
      i++;
      bdt_options = new BDTOptimizer("vbf" , "BDT_VBF" , string(argv[i] ) );
      for( auto i : *bdt_options)
	i.PrintAll(cout);
      methodList += "VBF";
      continue;
    }
    if( regMethod=="--cuts"){
      //i++;
      methodList += "Cuts";
      continue;
    }
    if( regMethod=="--cutsD"){
      //i++;
      methodList += "CutsD";
      continue;
    }
    if( regMethod=="--fisher"){
      //i++;
      methodList += "Fisher";
      continue;
    }
    if( regMethod=="--boostedfisher"){
      //i++;
      methodList += "BoostedFisher";
      continue;
    }
    if( regMethod=="--indir"){
      i++;
      DEFAULT_INFNAME = std::string (*(argv + i)).c_str();
      continue;
    }
    if( regMethod=="--cat"){
      i++;
      category = std::string (*(argv + i)).c_str();
      if(!category.Contains(":")){
	std::cout<<"Provide \"P:Q\" with P in {MM, A} and Q in {VBF, HighPt, HighPtVBF, V1J} " << std::endl;
	return -1;
      }
      continue;
    }
  }
  
  return TMVAClassification(methodList , ext , bdt_options, DEFAULT_INFNAME, category);
}
