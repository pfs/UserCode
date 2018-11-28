#ifndef _VBFVectorBoson_h_
#define _VBFVectorBoson_h_

#include "TLorentzVector.h"
#include "TopLJets2015/TopAnalysis/interface/ObjectTools.h"
#include "TopLJets2015/TopAnalysis/interface/SelectionTools.h"
#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>

#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"
#include "TopLJets2015/TopAnalysis/interface/VBFVectorBoson.h"
#include "TopLJets2015/TopAnalysis/interface/EfficiencyScaleFactorsWrapper.h"
#include "TopLJets2015/TopAnalysis/interface/L1PrefireEfficiencyWrapper.h"
#include "TopLJets2015/TopAnalysis/interface/GeneratorTools.h"
#include "TopLJets2015/TopAnalysis/interface/VBFDiscriminatorInputs.h"
#include "TopLJets2015/TopAnalysis/interface/VBFAnalysisCategories.h"

#include <vector>
#include <set>
#include <iostream>
#include <algorithm>
#include "TRandom3.h"
#include "TMath.h"

class VBFVectorBoson{
 public:

  /**
     @short CTOR for VBF V analysis
  */
 VBFVectorBoson(TString filename,
                TString outname,
                TH1F *normH, 
                TH1F *genPU,
                TString era,
                Float_t xsec,
                Bool_t debug=false, 
                Bool_t CR=false, 
                Bool_t QCDTemp=true, 
                Bool_t SRfake= false, 
                Bool_t skimtree=false, 
                bool runSysts=false,
                bool blind =true)
   : filename_(filename), outname_(outname), normH_(normH), genPU_(genPU), era_(era), 
    debug_(debug), CR_(CR), QCDTemp_(QCDTemp), SRfake_(SRfake), skimtree_(skimtree), 
    runSysts_(runSysts), doBlindAnalysis_(blind),
    xsec_(xsec)
  {
    fMVATree_ = NULL;
    newTree_ = NULL;
    init();
    setSelectionCuts();
    rnd_.SetSeed(123456789);
  };

  /**
     @short DTOR
  */
  ~VBFVectorBoson()
    {
    }

  /**
     @short initialize all variables
   */
  inline void init()
  {
    this->readTree();
    std::cout << "...producing " << outname_ << " from " << nentries_ << " events" << std::endl;
    this->prepareOutput();
    this->bookHistograms();
    this->loadCorrections();
    if(skimtree_) this->addMVAvars();
    selector_ = new SelectionTool(filename_, debug_, triggerList_,SelectionTool::VBF);
    std::cout << "init done" << std::endl;
  }

  /**
     @short define baseline selection cuts
  */
  inline void setSelectionCuts()
  {
    zMassWindow_  = 15.;

    leadJetPt_    = 50.;
    subLeadJetPt_ = 50.;
    jetPuId_      = 1;
    cleanEENoise_ = true;

    lowVPtCut_        = 75.;
    lowVPtDetaJJCut_  = 3.0;
    lowVPtMaxRapCut_  = 1.4442;
    lowVPtPhotonTrigs_.push_back("HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_v");

    highVPtCut_       = 200.;
    highVPtPhotonTrigs_.push_back("HLT_Photon200_v");

    lowMJJCut_  = 500.;
    highMJJCut_ = 1000.;        

    applyTrigSafePhoton_=true;
  } 

  void saveHistos();
  void readTree();
  void prepareOutput();
  void bookHistograms();  
  void loadCorrections();
  void addMVAvars();

  /**
     @short the working horse method
   */
  void runAnalysis();

  /**
     @short histogram filler
   */
  void fillControlHistos(TLorentzVector boson, 
                         std::vector<Jet> jets, 
                         float wgt,
                         TString c,
                         std::map<TString, int> mults, 
                         std::vector<Particle> fakePhotonCR ={}, 
                         std::vector<Particle> tightPhotonCR={});


private:

  /**
     @short dice a random number and select for training
  */
  inline int useForTraining(float selFrac=0.5)
  {
    double myRnd = rnd_.Rndm();
    return myRnd>=selFrac ? 1 : 0;
  }
  
  
  //class variables
  TString filename_, outname_, baseName_;
  Int_t nentries_;
  TH1F *normH_, *genPU_, *triggerList_;
  TString era_;
  Bool_t debug_, CR_, QCDTemp_, SRfake_, skimtree_, vetoPromptPhotons_;
  TFile *f_ /*inFile*/, *fMVATree_, *fOut_;
  TTree *t_ /*inTree*/, *newTree_ /*MVA*/;
  HistTool *ht_;
  MiniEvent_t ev_;

  //list of systematics
  std::vector<TString> expSysts_;
  std::vector<WeightSysts_t> weightSysts_;


  TRandom3 rnd_;
  
  //corrections
  LumiTools * lumi_;
  EfficiencyScaleFactorsWrapper * gammaEffWR_;
  L1PrefireEfficiencyWrapper *l1PrefireWR_;
  FakeRateTool * fr_;
  std::vector<Particle> photons_,relaxedTightPhotons_,tmpPhotons_; 


  //Variables to be added to the MVA Tree and additional variables
  vbf::DiscriminatorInputs vbfVars_;
  float vbfmva_;
  float evtWeight_;
  float sihih_,chiso_,r9_,hoe_,mindrl_,mindrj_;
  int training_;
  bool runSysts_;
  bool doBlindAnalysis_;  

  float xsec_;

  //selection configuration
  float leadJetPt_, subLeadJetPt_;
  int jetPuId_;
  bool cleanEENoise_;
  float zMassWindow_;
  float lowVPtCut_,  lowVPtDetaJJCut_, lowVPtMaxRapCut_;
  float highVPtCut_;
  float lowMJJCut_, highMJJCut_;
  std::vector<TString> lowVPtPhotonTrigs_, highVPtPhotonTrigs_;
  bool applyTrigSafePhoton_;
  
  //categorizer and tree reader
  vbf::Category category_;  
  SelectionTool *selector_;	
};

#endif
