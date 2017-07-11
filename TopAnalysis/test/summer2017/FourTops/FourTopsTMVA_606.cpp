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
    
    TChain *TTChain       = new TChain("TMVAanalysis");
    TChain *FourTopsChain = new TChain("TMVAanalysis");

    TTChain->Add("MC13TeV_TTJets2l2nu_amcatnlo.root");
    FourTopsChain->Add("MC13TeV_TTTT.root");
    
    printf("TTChain = %lld events\n",TTChain->GetEntries());
    printf("FourTopsChain = %lld events\n",FourTopsChain->GetEntries());
    getchar();
    
    TString outFileName("FourTopsTMVA_out.root");
    TFile *outFile = new TFile(outFileName, "RECREATE");

    TMVA::Factory *factory = new TMVA::Factory("FourTopsML", outFile, "V:Color");

    factory->AddSignalTree(FourTopsChain, 1.0);
    factory->AddBackgroundTree(TTChain, 1.0);
    
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
    //factory->BookMethod(TMVA::Types::kMLP,"MLP_ANN","V:VarTransform=N"); // ROOT neural network
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
