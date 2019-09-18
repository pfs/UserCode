#ifndef _GENERATORTOOLS_H_
#define _GENERATORTOOLS_H_

#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"

#include <vector>
#include <string>
#include <map>

#include "TSystem.h"
#include "TFile.h"
#include "TString.h"
#include "TGraph.h"
#include "TF1.h"

//available theory systs
typedef std::pair<TString,int> WeightSysts_t;
std::vector< WeightSysts_t > getWeightSysts(TFile *,TString proc="EWKAJJ2017");
std::vector< WeightSysts_t > getPartonShowerWeightSysts(TFile *f);

//b-fragmentation
std::map<TString, TGraph*> getBFragmentationWeights(TString era);
float computeBFragmentationWeight(MiniEvent_t &ev, TGraph* wgtGr);

//semi-leptonic B branching ratios
std::map<TString, std::map<int, float> > getSemilepBRWeights(TString era);
float computeSemilepBRWeight(MiniEvent_t &ev, std::map<int, float> corr, int pid=0, bool useabs=true);

//resonance tools
TF1 *getRBW(float m,float g);
float weightBW(TF1 *bwigner,std::vector<float> &obsm,float g,float m,float gini,float mini);

#endif
