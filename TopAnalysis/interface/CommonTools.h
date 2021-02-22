#ifndef _common_tools_h_
#define _common_tools_h_

#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TopLJets2015/TopAnalysis/interface/ObjectTools.h"
#include "TopLJets2015/TopAnalysis/interface/SelectionTools.h"
#include "TopLJets2015/TopAnalysis/interface/EfficiencyScaleFactorsWrapper.h"
#include "TopLJets2015/TopAnalysis/interface/JECTools.h"
#include "TopLJets2015/TopAnalysis/interface/BtagUncertaintyComputer.h"
#include "TopLJets2015/TopAnalysis/interface/LumiTools.h"
#include "TopLJets2015/TopAnalysis/interface/KinematicsTools.h"
#include "TopLJets2015/TopAnalysis/interface/Davismt2.h"
#include "TopLJets2015/TopAnalysis/interface/HistTool.h"
#include "TopLJets2015/TopAnalysis/interface/TemplatedFitTools.h"
#include "TopLJets2015/TopAnalysis/interface/FakeRateTool.h"

Float_t computeMT_old(TLorentzVector &a, TLorentzVector &b);
FactorizedJetCorrector * getFactorizedJetEnergyCorrector_old(TString, bool);
std::vector<float> getJetResolutionScales_old(float pt, float eta, float genjpt);
float getLeptonEnergyScaleUncertainty_old(int l_id, float l_pt, float l_eta);


inline void fill_selection_histo(std::map<TString, TH1 *> & plots, TString ch, const char * level, const char * bin, double  wgt)
{
  /* if (TString(level) == "reco") */
  /*   { */
  /*     printf("fill sel histo %s\n", bin); */
  /*   } */
  {
    const TString key(ch + "_" + level + "_selection");
    plots[key] -> Fill(bin, wgt);
  }
  {
    const TString key(TString("L") + "_" + level + "_selection");
    plots[key] -> Fill(bin, wgt);
  }
}


#endif
