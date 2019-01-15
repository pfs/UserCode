#ifndef _TOP17010_h_
#define _TOP17010_h_

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
#include "TopLJets2015/TopAnalysis/interface/EfficiencyScaleFactorsWrapper.h"
#include "TopLJets2015/TopAnalysis/interface/L1PrefireEfficiencyWrapper.h"
#include "TopLJets2015/TopAnalysis/interface/GeneratorTools.h"
#include "TopLJets2015/TopAnalysis/interface/TopWidthEvent.h"

#include <vector>
#include <set>
#include <iostream>
#include <algorithm>
#include "TMath.h"


class TOP17010 {
  
 public:

  /**
     @short CTOR for TOP-17-010 analysis
     scenario is an unsigned 32b integer coded as (imass&0xffff<<16)+(iwidth&0xffff)
     the indices will be used to set the target mass and width of the reweighting
     target mass  = 169 + imass*0.25
     target width = 0.7 + iwidth*0.01
  */
 TOP17010(TString filename, TString outname, TH1F *normH,  TH1F *genPU, TString era, UInt_t scenario=0, Bool_t debug=false)        
   : filename_(filename), outname_(outname), normH_(normH), genPU_(genPU), era_(era), debug_(debug),
    targetGt_(-1), targetMt_(-1)
  {
    init(scenario);
    bookHistograms();
  };

  /** 
      @short the the target mass and width to reweight to
  */
  void setTarget(float mass,float width) { targetMt_=mass; targetGt_=width; }

  /**
     @short the working horse method
   */
  void runAnalysis();


  /**
     @short DTOR
  */
  ~TOP17010()
    {
    }


private:

  void init(UInt_t scenario);
  void readTree();
  void bookHistograms();
  void fillControlHistograms(TopWidthEvent &twe,float &wgt);
  void applyMC2MC(std::vector<Jet> &jetColl);
  float getJERSFBreakdown(TString key,float abseta);

  //class variables
  TString filename_, outname_;
  Int_t nentries_;
  TH1F *normH_, *genPU_, *triggerList_;
  TString era_;
  Bool_t debug_, isSignal_;
  TFile *f_,*fOut_;
  TTree *t_;
  HistTool *ht_;
  MiniEvent_t ev_;

  //list of systematics
  std::vector<TString> expSysts_;
  std::vector<WeightSysts_t> weightSysts_;

  //corrections
  LumiTools * lumi_;
  EfficiencyScaleFactorsWrapper * gammaEffWR_;
  L1PrefireEfficiencyWrapper *l1PrefireWR_;
  BTagSFUtil *btvSF_;
  std::map<TString, TGraph*> fragWeights_;
  std::map<TString, std::map<int, float> > semilepBRwgts_;
  std::map<TString, TGraphErrors *> mc2mcCorr_;
  std::map<TString, TH1 *> jerSFBreakdown_;

  //mass,width and relativistic breit-wigner parameterization
  float targetGt_, origGt_;
  float targetMt_, origMt_;
  TF1 *rbwigner_;

  SelectionTool *selector_;
};

#endif
