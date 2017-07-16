#include <stdio.h>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TROOT.h"

#include "TMVA/Tools.h"

void FourTopsTMVA_606()
{
    TMVA::Tools::Instance();
    
    (TMVA::gConfig().GetVariablePlotting()).fMaxNumOfAllowedVariablesForScatterPlots = 1000;
    
    TChain *FourTopsChain  = new TChain("TMVAanalysis");
    TChain *TTChain        = new TChain("TMVAanalysis");
    TChain *SingleTopChain = new TChain("TMVAanalysis");
    TChain *WChain         = new TChain("TMVAanalysis");
    TChain *DYChain        = new TChain("TMVAanalysis");
    TChain *BosonChain     = new TChain("TMVAanalysis");
    TChain *ToNuChain      = new TChain("TMVAanalysis");

    string dir("/afs/cern.ch/user/v/vwachira/CMSSW_8_0_28/src/TopLJets2015/TopAnalysis/test/summer2017/");
    FourTopsChain-> Add((dir+string("MC13TeV_TTTT.root")).c_str());
    TTChain->       Add((dir+string("MC13TeV_TTJets.root")).c_str());
    SingleTopChain->Add((dir+string("MC13TeV_SingleTbar_tW.root")).c_str());
    SingleTopChain->Add((dir+string("MC13TeV_SingleT_tW.root")).c_str());
    SingleTopChain->Add((dir+string("MC13TeV_SingleTbar_t.root")).c_str());
    SingleTopChain->Add((dir+string("MC13TeV_SingleT_t.root")).c_str());
    WChain->        Add((dir+string("MC13TeV_W0Jets.root")).c_str());
    WChain->        Add((dir+string("MC13TeV_W1Jets.root")).c_str());
    WChain->        Add((dir+string("MC13TeV_W2Jets.root")).c_str());
    DYChain->       Add((dir+string("MC13TeV_DY50toInf_mlm.root")).c_str());
    DYChain->       Add((dir+string("MC13TeV_DY10to50.root")).c_str());
    BosonChain->    Add((dir+string("MC13TeV_ZZTo2L2Nu.root")).c_str());
    BosonChain->    Add((dir+string("MC13TeV_ZZTo2L2Q.root)")).c_str());
    BosonChain->    Add((dir+string("MC13TeV_WWToLNuQQ.root")).c_str());
    BosonChain->    Add((dir+string("MC13TeV_WWTo2L2Nu.root")).c_str());
    BosonChain->    Add((dir+string("MC13TeV_WZTo3LNu.root")).c_str());
    ToNuChain->     Add((dir+string("MC13TeV_TTWToLNu.root")).c_str());
    ToNuChain->     Add((dir+string("MC13TeV_TTWToQQ.root")).c_str());
    ToNuChain->     Add((dir+string("MC13TeV_TTZToQQ.root")).c_str());
    ToNuChain->     Add((dir+string("MC13TeV_TTZToLLNuNu.root")).c_str());
    
    TString outFileName("FourTopsTMVA_out.root");
    TFile *outFile = new TFile(outFileName, "RECREATE");

    TMVA::Factory *factory = new TMVA::Factory("FourTopsML", outFile, "V:Color");

    factory->AddSignalTree(FourTopsChain, 1.0);
    factory->AddBackgroundTree(TTChain, 1.0);
    factory->AddBackgroundTree(SingleTopChain, 1.0);
    factory->AddBackgroundTree(WChain, 1.0);
    factory->AddBackgroundTree(DYChain, 1.0);
    factory->AddBackgroundTree(BosonChain, 1.0);
    factory->AddBackgroundTree(ToNuChain, 1.0);
    
    factory->AddVariable("jet_1_highest_csv",'F');
    factory->AddVariable("jet_2_highest_csv",'F');
    factory->AddVariable("jet_3_highest_csv",'F');
    factory->AddVariable("jet_4_highest_csv",'F');
    // TODO: Add number of jets, b-jets, leptons
    factory->AddVariable("jet_ht",'F');
    factory->AddVariable("jet_btag_ht",'F');
    factory->AddVariable("jet_non_btag_ht",'F');
    factory->AddVariable("jet_smallest_angle",'F');
    factory->AddVariable("jet_smallest_angle_2b",'F');

    factory->AddVariable("lepton_ht",'F');
    factory->AddVariable("lepton_3_chIso",'F');
    factory->AddVariable("lepton_4_chIso",'F');
    factory->AddVariable("lepton_3_minIso",'F');
    factory->AddVariable("lepton_4_minIso",'F');

    factory->AddVariable("met_pt",'F');

    factory->AddVariable("n_jets",'F');
    factory->AddVariable("n_bjets",'F');
    factory->AddVariable("n_leptons",'F');

    factory->AddVariable("n_mu_p",'F');
    factory->AddVariable("n_mu_m",'F');
    factory->AddVariable("n_ele_p",'F');
    factory->AddVariable("n_ele_m",'F');

    TCut mycuts = ""; // Signal cut
    TCut mycutb = ""; // Background cut

    factory->PrepareTrainingAndTestTree(mycuts, mycutb, "SplitMode=Random:V");
    
    //factory->BookMethod(TMVA::Types::kLikelihood,"LikelihoodD","!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:VarTransform=Decorrelate");	
    //factory->BookMethod(TMVA::Types::kCuts,"Cuts","V:VarTransform=N"); // Cuts
    //factory->BookMethod(TMVA::Types::kLikelihood,"Likelihood","V:VarTransform=N"); // Projective likelihood estimator
    //factory->BookMethod(TMVA::Types::kPDERS,"PDERS","V:VarTransform=N"); // PDE range-search approach
    //factory->BookMethod(TMVA::Types::kPDEFoam,"PDEFoam","V"); // PDE-Foam
    //factory->BookMethod(TMVA::Types::kKNN,"kNN","V:VarTransform=N"); // k-Nearest Neighbour classifier
    //factory->BookMethod(TMVA::Types::kHMatrix,"HMatrix","V"); // H-Matrix discriminant
    //factory->BookMethod(TMVA::Types::kFisher,"Fisher","V"); // Fisher discriminants (linear discriminant analysis)
    //factory->BookMethod(TMVA::Types::kFDA,"FDA","V"); // Function discriminant analysis (FDA)
    //factory->BookMethod(TMVA::Types::kMLP,"MLP_ANN","V:NeuronType=tanh:VarTransform=N:NCycles=1000:HiddenLayers=N+10,N"); // ROOT neural network
    factory->BookMethod(TMVA::Types::kMLP,"MLP_ANN","V:VarTransform=N"); // ROOT neural network
    //factory->BookMethod(TMVA::Types::kSVM,"SVM","V:VarTransform=N"); // Support Vector Machine
    factory->BookMethod(TMVA::Types::kBDT,"BDT","V"); // Boosted Decision Tree
    //factory->BookMethod(TMVA::Types::kRuleFit,"RuleFit","V"); // Predictive learning via rule ensembles
    
    factory->OptimizeAllMethodsForClassification(); // Takes too much time!
    factory->TrainAllMethods();
    factory->TestAllMethods();
    
    factory->EvaluateAllMethods();

    outFile->Close();

    TMVA::TMVAGui(outFileName);
}
