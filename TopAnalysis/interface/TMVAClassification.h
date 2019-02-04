#ifndef __TMVAClassification_h_
#define __TMVAClassification_h_

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/MethodBDT.h"
#include "TMVA/Config.h"
#include "TMVA/ResultsClassification.h"
#include "TMVA/ROCCurve.h"
// Authored by Hamed Bakhshiansohi 08/12/2016
// I edited tmva/tmvagui/src/mvas.cxx file
// and changed "x" -> "d" option for kolmogorov test so that the overtraining works
// you need to compile the root again

using namespace std;
using namespace TMVA;

class BDTOptions{
public:
  int ntrees ;
  double minnodesize ;
  int maxdepth ;
  double adaboostbeta ;
  int ncuts ;
  TString name;
  
  BDTOptions(  TString name_ = "" , int ntrees_ = 50 , double minnodesize_ = 2. , int maxdepth_ = 50 , double adaboostbeta_ = 0.5 , int ncuts_ = 10  );
  void Set( string options_ );
  TString to_string( TString name );

  TMVA::MethodBDT* BookMethod( Factory* factory , DataLoader* dataloader, TString ext );
  

  double ROCIntegral;
  double Signal_TrainTest_Chi2 , Signal_TrainTest_Kolmo , Bkg_TrainTest_Kolmo , Bkg_TrainTest_Chi2 ; 
  void CalcEfficiency(TMVA::IMethod * method , TDirectory* dir , TString dsname , TString ext);
  
  double ReadROC(TMVA::IMethod * method, TString ext  );
  double ReadOverTrainingParam(TDirectory* dir , TString methodTitle , TString dsname , bool signal , bool kolmo );

  void PrintAll(std::ostream& out  , bool header = false){
    if( header )
      out << "name" << ","
	   << "ntrees" << ","
	   << "minnodesize" << ","
	   << "maxdepth" << ","
	   << "adaboostbeta" << ","
	   << "ncuts" << ","
	   << "ROCIntegral" << ","
	   << "Signal_TrainTest_Kolmo" << ","
	   << "Signal_TrainTest_Chi2" << ","
	   << "Bkg_TrainTest_Kolmo" << ","
	   << "Bkg_TrainTest_Chi2" << endl;
    else
      out << name << ","
	   << ntrees << ","
	   << minnodesize << ","
	   << maxdepth << ","
	   << adaboostbeta << ","
	   << ncuts << ","
	   << ROCIntegral << ","
	   << Signal_TrainTest_Kolmo << ","
	   << Signal_TrainTest_Chi2 << ","
	   << Bkg_TrainTest_Kolmo << ","
	   << Bkg_TrainTest_Chi2 << endl;
  }
};


class BDTOptimizer : public std::vector<BDTOptions> {
public:
  TString dsName;
  TString MVATitle;
  
  BDTOptimizer( TString dsname , TString MVATitle_ , string options_ );
  void FillRange( vector<string> range_info , vector<int>& vals );
  void FillRange( vector<string> range_info , vector<double>& vals );

  void EvaluateAll( Factory* factory , TDirectory* dir , TString ext){
    if( this->size() == 0 )
      return;
    
    ofstream myfile;
    myfile.open ( MVATitle + ".csv" );
    this->at(0).PrintAll(myfile , true);
    for(auto option : *this){
      option.CalcEfficiency( factory->GetMethod( dsName, option.name+ext ) , dir , dsName, ext ) ;
      option.PrintAll(myfile);
    }
  }
};

class DataLoaderWrapper: public TMVA::DataLoader{
public:
 DataLoaderWrapper(TString name= "defaultName"):DataLoader(name){};
  ~DataLoaderWrapper(){};

  void readInputs(TString fname){
    inputVars.clear();
    ifstream f(fname);
    string line;
    int iLine = 0;
    if (f.is_open()) {
      while (getline(f, line)) {
	iLine++;
	TString linestr(line.c_str());
	int nameIndex = linestr.Index(",");
	TString variable = linestr(0,nameIndex);
	TString sub      = linestr(nameIndex+1, linestr.Length() - nameIndex);
	int qId = -1;
	for(int i = 0; i < sub.Length(); i++){if(! TString(sub(i,1)).IsWhitespace()) {qId = i; break;}}
	TString varName  = sub(qId, sub.Length() - qId);
	inputVars.push_back(make_pair(variable,varName));
      }
      f.close();
    }
  }

  void setVariables(){
    for(unsigned int i = 0; i < inputVars.size(); i++){
      TMVA::DataLoader::AddVariable( inputVars[i].first, inputVars[i].second, "", 'F');
    }
  }
  void setCutOptVars(){
     TMVA::DataLoader::AddVariable( "mjj", "mjj", "", 'F' ) ;
     TMVA::DataLoader::AddVariable( "j_pt[1]", "SubLeadJetPt", "", 'F' ) ;
  };
  void setBestVars(bool isHighMJJ = false, bool excludeCuts = false){
    if (!excludeCuts){
      //this->setCutOptVars();//Not needed for High/Low MJJ categories
      ///////////////////////////////////////////////
      // Best variables w/o cuts on mjj and jet pt //
      ///////////////////////////////////////////////
      // TMVA::DataLoader::AddVariable( "detajj", "detajj", "", 'F' ) ;    
      // TMVA::DataLoader::AddVariable( "ht", "ht", "", 'F' ) ;
      // TMVA::DataLoader::AddVariable("j_c2_02[0]",        "jet_c2_021",       "", 'F' ) ;
      // TMVA::DataLoader::AddVariable("j_gaptd[0]",        "jet_gaptd1",       "", 'F' ) ;

      //Variables for the May25th studies
      /* TMVA::DataLoader::AddVariable( "ht", "ht", "", 'F' ) ;  */
      /* TMVA::DataLoader::AddVariable("j_gawidth[0]",      "jet_gawidth1",     "", 'F' ) ; */
      /* TMVA::DataLoader::AddVariable( "forwardeta", "forwardeta", "", 'F' ) ; */
      /* TMVA::DataLoader::AddVariable("j_c1_05[0]",        "jet_c1_051",       "", 'F' ) ; */
      /* TMVA::DataLoader::AddVariable( "balance", "balance", "", 'F' ) ; */

      //Variables for the HighMJJ category
      TMVA::DataLoader::AddVariable( "ht", "ht", "", 'F' ) ;
      TMVA::DataLoader::AddVariable( "forwardeta", "forwardeta", "", 'F' ) ;
      TMVA::DataLoader::AddVariable( "balance", "balance", "", 'F' ) ;
      TMVA::DataLoader::AddVariable( "dphivj0", "dphivj0", "", 'F' ) ;
      TMVA::DataLoader::AddVariable( "mjj", "mjj", "", 'F' ) ;
      TMVA::DataLoader::AddVariable( "j_qg[0]", "leadjet_qg", "", 'F' ) ;
      //      TMVA::DataLoader::AddVariable("j_c2_02[1]",        "jet_c2_021",       "", 'F' ) ;   
      //      TMVA::DataLoader::AddVariable( "aplanarity", "aplanarity", "", 'F' ) ;
      //      TMVA::DataLoader::AddVariable( "D", "D", "", 'F' ) ;
      //      TMVA::DataLoader::AddVariable( "j_qg[1]", "subleadjet_qg", "", 'F' ) ;

      
      /* if(isHighMJJ) */
      /* 	TMVA::DataLoader::AddVariable( "dphivj1", "dphivj1", "", 'F' ); */
      /* else */
      /* 	TMVA::DataLoader::AddVariable("j_c2_00[0]",        "jet_c2_000",       "", 'F' ) ; */


    } else {}
  };
  void setCentralJetVars(){
     TMVA::DataLoader::AddVariable( "centjy", "centjy", "", 'F' ) ;
     TMVA::DataLoader::AddVariable( "ncentjj", "ncentjj", "", 'F' ) ;
     TMVA::DataLoader::AddVariable( "dphivj2", "dphivj2", "", 'F' ) ;
     TMVA::DataLoader::AddVariable( "dphivj3", "dphivj3", "", 'F' ) ;
  };
  void setLowVPtHighMJJVariables(){
    TMVA::DataLoader::AddVariable( "mjj",              "mjj", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "dphijj",           "dphijj", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "ystar",            "ystar", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "dphibjj",          "dphibjj", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "balance",          "balance", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "subleadj_gawidth", "subleadj_gawidth", "", 'F' ) ;
    /* TMVA::DataLoader::AddVariable( "j_c2_00[0]",       "jet_c2_001",       "", 'F' ) ; */
    /* TMVA::DataLoader::AddVariable( "j_c2_00[1]",       "jet_c2_002",       "", 'F' ) ; */
    TMVA::DataLoader::AddVariable( "j_qg[0]",          "leadjet_qg", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "dphivj0",          "dphivj0", "", 'F' ) ;
    //    TMVA::DataLoader::AddVariable( "dphivj2",          "dphivj2", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "mht",              "mht", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "ht",               "ht", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "aplanarity",       "aplanarity", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "sphericity",       "sphericity", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "circularity",      "circularity", "", 'F' ) ;

  }

  void setHighVPtHighMJJVariables(){
    TMVA::DataLoader::AddVariable( "mjj",              "mjj", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "jjpt",             "jjpt", "", 'F' ) ;
    //    TMVA::DataLoader::AddVariable( "dphijj",           "dphijj", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "ystar",            "ystar", "", 'F' ) ;
    //    TMVA::DataLoader::AddVariable( "dphibjj",          "dphibjj", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "balance",          "balance", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "subleadj_gawidth", "subleadj_gawidth", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "j_qg[0]",          "leadjet_qg", "", 'F' ) ;
    //    TMVA::DataLoader::AddVariable( "j_qg[1]",          "subleadjet_qg", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "dphivj0",          "dphivj0", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "dphivj1",          "dphivj1", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "ht",               "ht", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "sphericity",       "sphericity", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "D",                "D", "", 'F' ) ; 
    TMVA::DataLoader::AddVariable( "circularity",      "circularity", "", 'F' ) ;
  }
  void setHighVPtLowMJJVariables(){

    TMVA::DataLoader::AddVariable( "leadj_pt",         "leadj_pt", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "mjj",              "mjj", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "detajj",           "detajj", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "dphijj",           "dphijj", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "ystar",            "ystar", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "dphibjj",          "dphibjj", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "balance",          "balance", "", 'F' ) ;
    /* TMVA::DataLoader::AddVariable( "j_c2_00[0]",       "jet_c2_001",       "", 'F' ) ; */
    /* TMVA::DataLoader::AddVariable( "j_c2_00[1]",       "jet_c2_002",       "", 'F' ) ; */
    TMVA::DataLoader::AddVariable( "j_qg[1]",          "subleadjet_qg", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "dphivj0",          "dphivj0", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "dphivj1",          "dphivj1", "", 'F' ) ;
    /* TMVA::DataLoader::AddVariable( "dphivj2",          "dphivj2", "", 'F' ) ; */
    /* TMVA::DataLoader::AddVariable( "dphivj3",          "dphivj3", "", 'F' ) ; */
    TMVA::DataLoader::AddVariable( "mht",              "mht", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "sphericity",       "sphericity", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "D",                "D", "", 'F' ) ;     
  }

  void setHighVPtVariables(){
    TMVA::DataLoader::AddVariable( "mjj",              "mjj", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "jjpt",             "jjpt", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "detajj",           "detajj", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "dphijj",           "dphijj", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "ystar",            "ystar", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "relbpt",           "relbpt", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "dphibjj",          "dphibjj", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "balance",          "balance", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "j_c2_00[0]",       "jet_c2_001",       "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "j_c2_00[1]",       "jet_c2_002",       "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "j_qg[0]",          "leadjet_qg", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "j_qg[1]",          "subleadjet_qg", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "dphivj1",          "dphivj1", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "dphivj2",          "dphivj2", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "ht",               "ht", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "isotropy",         "isotropy", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "D",                "D", "", 'F' ) ; 
  }



  void setAllVariables(bool excludeCuts = false){

    TMVA::DataLoader::AddVariable( "forwardeta",       "forwardeta", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "leadj_pt",         "leadj_pt", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "subleadj_pt",      "subleadj_pt", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "mjj",              "mjj", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "jjpt",             "jjpt", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "detajj",           "detajj", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "dphijj",           "dphijj", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "ystar",            "ystar", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "relbpt",           "relbpt", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "dphibjj",          "dphibjj", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "balance",          "balance", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "leadj_gawidth",    "leadj_gawidth", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "subleadj_gawidth", "subleadj_gawidth", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "j_qg[0]",          "leadjet_qg", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "j_qg[1]",          "subleadjet_qg", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "dphivj0",          "dphivj0", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "dphivj1",          "dphivj1", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "mht",              "mht", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "ht",               "ht", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "isotropy",         "isotropy", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "aplanarity",       "aplanarity", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "sphericity",       "sphericity", "", 'F' ) ;
    TMVA::DataLoader::AddVariable( "C",                "C", "", 'F' ) ; 
    TMVA::DataLoader::AddVariable( "D",                "D", "", 'F' ) ; 
    TMVA::DataLoader::AddVariable( "circularity",      "circularity", "", 'F' ) ;
  };
  std::vector<std::pair<TString, TString> > inputVars; 
};

int TMVAClassification( TString myMethodList , TString extention, BDTOptimizer* bdt_options, TString sig, TString bkg, TString category, TString card );
std::vector<std::string> split(const std::string &s, char delim);


//needed to compile fully in cmssw
namespace TMVA{
  namespace TMVAGlob{
    void NormalizeHists( TH1* sig, TH1* bkg ) ;
  }
}

#endif
