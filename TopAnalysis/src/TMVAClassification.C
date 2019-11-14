#include "TopLJets2015/TopAnalysis/interface/TMVAClassification.h"

using namespace std;
using namespace TMVA;


void TMVAGlob::NormalizeHists( TH1* sig, TH1* bkg ) 
{
   if (sig->GetSumw2N() == 0) sig->Sumw2();
   if (bkg && bkg->GetSumw2N() == 0) bkg->Sumw2();
      
   if(sig->GetSumOfWeights()!=0) {
      Float_t dx = (sig->GetXaxis()->GetXmax() - sig->GetXaxis()->GetXmin())/sig->GetNbinsX();
      sig->Scale( 1.0/sig->GetSumOfWeights()/dx );
   }
   if (bkg != 0 && bkg->GetSumOfWeights()!=0) {
     Float_t dx = (bkg->GetXaxis()->GetXmax() - bkg->GetXaxis()->GetXmin())/bkg->GetNbinsX();
     bkg->Scale( 1.0/bkg->GetSumOfWeights()/dx );
   }
}

int TMVAClassification( TString myMethodList , TString extention, BDTOptimizer* bdt_options, TString sigName, TString bkgName, TString category, TString card )
{

  TMVA::MethodBDT* method;
  
   TMVA::Tools::Instance();
   std::map<std::string,int> Use;
   Use["VBF"] = 0;
   Use["Cuts"] = 0;
   Use["CutsD"] = 0;
   Use["Fisher"] = 0;
   Use["BoostedFisher"] = 0;
   TString dlName = "";   

   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++)
	it->second = 0;

      std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
	std::string regMethod(mlist[i]);

	if (Use.find(regMethod) == Use.end()) {
	  std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
	  for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
	  std::cout << std::endl;
	  return 1;
	}
	Use[regMethod] = 1;
	dlName = regMethod+"DL";
	cout<< regMethod<<"\tis set"<<endl;
      }
   }


   TString outfileName( "TMVA_"+ extention  +".root" );

   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
   
   TMVA::Factory *factory = new TMVA::Factory( extention , outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I:AnalysisType=Classification" );
   
   TFile *inputS = TFile::Open(sigName ); // To be fixed
   TTree *signalTree     = (TTree*)inputS->Get("data"); // To be fixed
   
   // std::cout << "Debug:=signal " << data << "  "  << std::endl;
   Double_t signalWeight     = 1.0;
   TString default_w_str = "evtWeight";
   
   // TMVA::DataLoader *dataloader = NULL;
   DataLoaderWrapper * dataloader = NULL;
   TFile *inputB = TFile::Open( bkgName ); //To be fixed
   TTree *background  = (TTree*)inputB->Get("data");

   if(bdt_options){
     cout<<"bdt_options->dsName = "<<bdt_options->dsName<<endl;
     dlName = bdt_options->dsName;
   }
   dataloader = new DataLoaderWrapper(dlName);
   
   std::vector<TString> cats = TMVA::gTools().SplitString( string(category), ':' );
   TString jet = " j_pt[0] > 50 && j_pt[1] > 50 ";
   TString lowMJJ = " mjj < 1000 ";
   TString highMJJ = " mjj > 1000 ";
   // TString looseMJJ = " mjj > 200 ";
   // std::cout << "Debug:=signal event rate " << netries << "  "  << std::endl;
   TString lowVpt = " gamma_pt[0] < 200 ";

   TString cut = "";
   bool isLowMJJ (false), isLowVpt(false), isHighMJJ(false), isHighVpt(false);
   bool isHighPt(false), isVBF(false); // , isLooseVBF(false) , isLooseVPt(false);
   for (unsigned int iCat = 0; iCat < cats.size(); iCat++){
     if (cats[iCat] == "LowVPt")  {cout<<"It is LowVPt!"<<endl; isLowVpt = true;}
     if (cats[iCat] == "HighVPt") {cout<<"It is HighVPt!"<<endl; isHighVpt = true; }
     if (cats[iCat] == "HighMJJ") {cout<<"It is HighMJJ!"<<endl; isHighMJJ = true; }
     if (cats[iCat] == "LowMJJ")  {cout<<"It is LowMJJ!"<<endl; isLowMJJ = true; }
     if (cats[iCat] == "HighPt")  {cout<<"It is HighPt!"<<endl; isHighPt = true; }
     if (cats[iCat] == "VBF")     {cout<<"It is VBF!"<<endl; isVBF = true; }
     //  if (cats[iCat] == "LooseVPt"){cout<<"It is LooseVPt!"<<endl; isLooseVPt = true; }
     //  if (cats[iCat] == "LooseVBF"){cout<<"It is LooseVBF!"<<endl; isLooseVBF = true; }

     cut+= "category." + cats[iCat] + " == 1 ";
     if(iCat < cats.size() - 1)
       cut+= " && ";
     
   }

   bool isLowVptHighMJJ  = (isLowVpt && isHighMJJ);
   //  bool isLowVptLowMJJ  = (isLowVpt && isLowMJJ);
   bool isHighVptHighMJJ = (isHighVpt && isHighMJJ);
   bool isHighVptLowMJJ  = (isHighVpt && isLowMJJ);
  

   isHighVpt = (isHighVptHighMJJ || isHighVptLowMJJ);
   if (isVBF)    isLowVpt = true;
   //  if (isLooseVBF)    isLowVpt = true;
   if (isHighPt) isHighVpt       = true;
   if (isHighVpt && extention.Contains("HighV") && extention.Contains("HighM")) isHighVptHighMJJ = true;
   if (isHighVpt && extention.Contains("HighV") && extention.Contains("LowM")) isHighVptLowMJJ = true;
   if (isLowVpt  && extention.Contains("HighM")) isLowVptHighMJJ = true;
   // if (isLowVpt  && extention.Contains("LowM")) isLowVptLowMJJ = true;
   
   
   if (isVBF){
     cut = cut + " && " + jet + " && " + lowVpt + " && mjj > 200 ";
     if(isLowVptHighMJJ) cut = cut + " && " + highMJJ;
   }
   if (isHighPt){
     cut = cut + " && mjj > 200 && " + jet;
     if (isHighVptHighMJJ) cut = cut + " && " + highMJJ;
     if (isHighVptLowMJJ)  cut = cut + " && " + lowMJJ;
   }
   //Loose region
   /*    if (isLooseVBF){
	cut = cut + " && mjj > 200 && " + jet;
     if (isHighVptHighMJJ) cut = cut + " && " + highMJJ;
     if (isHighVptLowMJJ)  cut = cut + " && " + lowMJJ;
     if (isLowVptHighMJJ) cut = cut + " && " + highMJJ;
     if (isLowVptLowMJJ)  cut = cut + " && " + lowMJJ;
     }*/

      /*   if (isLooseVPt){
     cut = cut + " && mjj > 200 && " + " && mjj<500 &&" + jet;
     if (isHighVptHighMJJ) cut = cut + " && " + highMJJ;
     if (isHighVptLowMJJ)  cut = cut + " && " + lowMJJ;
     if (isLowVptHighMJJ) cut = cut + " && " + highMJJ;
     if (isLowVptLowMJJ)  cut = cut + " && " + lowMJJ;
     }*/
  

   TString cutB = cut;
   int idx = cutB.Index("category.A");
   cutB.Replace(idx,10, "category.MM");
   TCut dyCut  = TCut(cutB);
   TCut mycuts = TCut(cut);
   TCut mycutb = TCut(cut) || dyCut;
   cout << "Cut is : " <<cut<<endl;
   dataloader->readInputs(card);
   /*
   if (Use["VBF"]) {
     if (isLowVptHighMJJ)
       dataloader->setVariables();
     else if (isHighVpt && !isHighVptHighMJJ  && !isHighVptLowMJJ)
       dataloader->setAllVariables(); 
     else if(isHighVptHighMJJ)
       dataloader->setVariables();
       else if(isHighVptLowMJJ)*/

   if (Use["VBF"]) {     

     dataloader->setVariables();
   }
   else if (Use["Cuts"] || Use["CutsD"]) dataloader->setCutOptVars();
   else if (Use["Fisher"] || Use["BoostedFisher"]) dataloader->setBestVars(isLowVptHighMJJ, false);
  
   dataloader->AddSignalTree( signalTree,     signalWeight );
   dataloader->AddBackgroundTree( background, 1 );

   dataloader->SetBackgroundWeightExpression( default_w_str  );
   dataloader->SetSignalWeightExpression( default_w_str );

 
   dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,
					   // "nTrain_Signal=1000:nTrain_Background=5000:SplitMode=Random:NormMode=NumEvents:!V");// To use the info in the Tree
					   // "nTrain_Signal=1000:"+nTrainBkg+":SplitMode=Random:NormMode=NumEvents:!V");// To use the info in the Tree
					   "nTrain_Signal=0:nTest_Signal=0:nTrain_Background=0:nTest_Background=0:SplitMode=Alternate:NormMode=None:!V");// To use the info in the Tree
					   // "SplitMode=Alternate:NormMode=None:!V");// To use the info in the Tree

     
   if( Use["VBF"] ){

     std::cout << "Loading VBF trees" << endl;
     cout<<bdt_options->size()<<endl;
     for(auto bdt_option : *bdt_options){
       method = bdt_option.BookMethod( factory , dataloader, extention );
       cout << "method "<<method->GetName() << endl;
     }
   } else if (Use["Cuts"] ){
     std::cout << "Loading Cuts trees" << endl; 
     factory->BookMethod( dataloader, TMVA::Types::kCuts, "Cuts","!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );
   } else if (Use["CutsD"])
     factory->BookMethod( dataloader, TMVA::Types::kCuts, "CutsD","!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate" );   
   else if (Use["Fisher"]) 
     factory->BookMethod( dataloader, TMVA::Types::kFisher, extention, "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );
   else if (Use["BoostedFisher"]) 
     factory->BookMethod( dataloader, TMVA::Types::kFisher, "BoostedFisher", "H:!V:Boost_Num=20:Boost_Transform=log:Boost_Type=AdaBoost:Boost_AdaBoostBeta=0.2:!Boost_DetailedMonitoring" );
     

     

   
   // --------------------------------------------------------------------------------------------------
   //  Now you can optimize the setting (configuration) of the MVAs using the set of training events
   // STILL EXPERIMENTAL and only implemented for BDT's !
   //
   // factory->OptimizeAllMethods("SigEffAt001","Scan");
   // factory->OptimizeAllMethods("ROCIntegral","FitGA");
   //
   // --------------------------------------------------------------------------------------------------

   // Now you can tell the factory to train, test, and evaluate the MVAs
   //
   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();

   // --------------------------------------------------------------

   // Save the output
   if(bdt_options)
     bdt_options->EvaluateAll( factory , outputFile, extention );
   
   
   outputFile->Close();

   delete factory;
   return 0;
}

std::vector<std::string> split(const std::string &s, char delim) {
  std::vector<std::string> elems;
  
  std::stringstream ss;
  ss.str(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }

  return elems;
}

BDTOptions::BDTOptions(  TString name_ , int ntrees_  , double minnodesize_ , int maxdepth_  , double adaboostbeta_  , int ncuts_  ){
  name = name_;
  ntrees = ntrees_;
  ncuts = ncuts_;
  adaboostbeta = adaboostbeta_;
  maxdepth = maxdepth_;
  minnodesize = minnodesize_;
}

void BDTOptions::Set( string options_ ){
  vector<string> options =  split( options_ , ':' );
  for(auto opt_ : options){
    vector<string> opt =  split( opt_ , '=' ) ;
    if(opt.size() != 2 ){
      cout << opt_ << " is a wrong word" << endl;
      continue;
    }
    if( opt[0] == "nt" )
      ntrees = stoi( opt[1] );
    else if( opt[0] == "mns" )
      minnodesize = stod( opt[1] );
    else if( opt[0] == "md" )
      maxdepth = stoi( opt[1] );
    else if( opt[0] == "abb" )
      adaboostbeta = stod( opt[1] );
    else if( opt[0] == "nc" )
      ncuts = stoi( opt[1] );
    else
      cout << opt[0] << " is a wrong word " << endl;
  }
}
  
TString BDTOptions::to_string( TString name ){
  return TString( name + "_nt" + std::to_string(ntrees)
		  + "_mns" + TString::Format("%.1f",minnodesize)
		  + "_md" + std::to_string(maxdepth)
		  + "_abb" + TString::Format("%.1f" , adaboostbeta)
		  + "_nc" + std::to_string(ncuts) );
}


BDTOptimizer::BDTOptimizer( TString dsname , TString MVATitle_ , string options_ ){
  dsName = dsname ;
  MVATitle = MVATitle_;


  vector<int> ntrees ;
  vector<double> minnodesize ;
  vector<int> maxdepth ;
  vector<double> adaboostbeta ;
  vector<int> ncuts ;

  cout << options_ <<endl;
  
  vector<string> options =  split( options_ , ':' );
  for(auto opt_ : options){
    cout << opt_ <<endl;
    vector<string> opt =  split( opt_ , '=' ) ;
    if(opt.size() != 2 ){
      cout << opt_ << " is a wrong word" << endl;
      continue;
    }
    cout << opt[0] << opt[1] << endl;
    vector<string> range_info = {opt[1]};
    if(opt[1].find( ',' ) != string::npos)
      range_info = split( opt[1] , ',' );
      
    if( opt[0] == "nt" )
      FillRange( range_info , ntrees);
    else if( opt[0] == "mns" )
      FillRange( range_info , minnodesize);
    else if( opt[0] == "md" )
      FillRange( range_info , maxdepth);
    else if( opt[0] == "abb" )
      FillRange( range_info , adaboostbeta);
    else if( opt[0] == "nc" )
      FillRange( range_info , ncuts);
    else
      cout << opt[0] << " is a wrong word " << endl;

    cout << ntrees.size() << " " << minnodesize.size() << " " << maxdepth.size() << " " << adaboostbeta.size() << " " << ncuts.size() << endl;
  }

  for(auto nt : ntrees)
    for(auto mns : minnodesize)
      for(auto md : maxdepth)
	for(auto adabb : adaboostbeta)
	  for(auto nc : ncuts)
	    this->emplace_back( MVATitle + to_string( this->size() ) ,  nt , mns , md , adabb , nc );
	    
  for(auto i : *this)
    i.PrintAll(cout);
	    
}

void BDTOptimizer::FillRange( vector<string> range_info , vector<int>& vals ){
  if(range_info.size() == 1)
    vals.push_back( stoi(range_info[0]) );
  else if(range_info.size() != 3){
    vals.push_back(-1);
    cout << "wrong range info : " << range_info[0] << endl;
  }
  else{
    int from = stoi(range_info[0]);
    int to = stoi(range_info[1]);
    int step = stoi(range_info[2]);

    while( from <= to ){
      vals.push_back( from );
      from += step;
    }
  }
}
void BDTOptimizer::FillRange( vector<string> range_info , vector<double>& vals ){
  if(range_info.size() == 1)
    vals.push_back( stod(range_info[0]) );
  else if(range_info.size() != 3){
    vals.push_back(-1.);
    cout << "wrong range info : " << range_info[0] << endl;
  }
  else{
    double from = stod(range_info[0]);
    double to   = stod(range_info[1]);
    double step = stod(range_info[2]);

    while( from <= to ){
      vals.push_back( from );
      from += step;
    }
  }
}

void BDTOptions::CalcEfficiency(TMVA::IMethod * method , TDirectory* dir , TString dsname, TString ext ){
  ROCIntegral = ReadROC( method , ext );

  Signal_TrainTest_Kolmo = ReadOverTrainingParam(dir , method->GetName() , dsname , true , true );
  Signal_TrainTest_Chi2 = ReadOverTrainingParam(dir , method->GetName() , dsname , true , false );

  Bkg_TrainTest_Kolmo = ReadOverTrainingParam(dir , method->GetName() , dsname , false , true );
  Bkg_TrainTest_Chi2 = ReadOverTrainingParam(dir , method->GetName() , dsname , false , false );
}

double BDTOptions::ReadROC(TMVA::IMethod * method , TString ext ){
  TString title = method->GetName() ;
  MethodBase* theMethod = dynamic_cast<MethodBase*>(method);
  if(theMethod==0) return -1.0;
  TMVA::Results *results= theMethod->Data()->GetResults(name+ext,Types::kTesting,Types::kClassification);
  std::vector<Float_t> *mvaRes = dynamic_cast<ResultsClassification *>(results)->GetValueVector();
  std::vector<Bool_t>  *mvaResType = dynamic_cast<ResultsClassification *>(results)->GetValueVectorTypes();
  Double_t fROCalcValue = 0;
  TMVA::ROCCurve *fROCCurve = nullptr;
  if (mvaResType->size() != 0) { 
    fROCCurve = new TMVA::ROCCurve(*mvaRes, *mvaResType);
    fROCalcValue = fROCCurve->GetROCIntegral();
  }
  return fROCalcValue ;           
}

double BDTOptions::ReadOverTrainingParam(TDirectory* dir , TString methodTitle , TString dsname , bool signal , bool kolmo ){
  TString hname = dsname + "/Method_" + methodTitle + "/" + methodTitle + "/MVA_" + methodTitle;
  TString appendix = signal ? "_S" : "_B" ;
  TH1* Test  = (TH1*) ( dir->Get( hname + appendix ) );
  TH1* Train = (TH1*) ( dir->Get( hname + "_Train" + appendix ) );
    
  TMVAGlob::NormalizeHists( Test , Train );

  for(int i = 0 ; i < Test->GetNbinsX() + 1 ; i ++){
    if( Test->GetBinContent(i) < 0 )
      Test->SetBinContent( i , 0 );
    if( Train->GetBinContent(i) < 0 )
      Train->SetBinContent( i , 0 );
  }
  
  Double_t ret = 0;
  if(kolmo )
    ret = Test->KolmogorovTest( Train, "X" ) ;
  else
    ret = Test->Chi2Test( Train , "WW CHI2/NDF" ); 
  return ret;
}

TMVA::MethodBDT* BDTOptions::BookMethod( Factory* factory , DataLoader* dataloader, TString extension ){
  TMVA::MethodBase* method = factory->BookMethod( dataloader, TMVA::Types::kBDT, name+extension,
						  "NegWeightTreatment=InverseBoostNegWeights:!H:!V:BoostType=AdaBoost:SeparationType=GiniIndex:CreateMVAPdfs:nCuts="+to_string(ncuts) );
  MethodBDT* ret = dynamic_cast<TMVA::MethodBDT*>( method );
  
  ret->SetNTrees(ntrees);
  ret->SetMinNodeSize(minnodesize);
  ret->SetMaxDepth(maxdepth);
  ret->SetAdaBoostBeta(adaboostbeta);
  cout << ret->GetName() << endl;
  return ret;
}
