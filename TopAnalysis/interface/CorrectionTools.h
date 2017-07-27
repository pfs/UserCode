#ifndef _correction_tools_h_
#define _correction_tools_h_

#include <vector>
#include <string>

#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"
#include "TopLJets2015/TopAnalysis/interface/BtagUncertaintyComputer.h"

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"

#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include <CondFormats/JetMETObjects/interface/JetResolutionObject.h>

#include "TGraphAsymmErrors.h"

typedef std::pair<TString,float> RunPeriod_t;
std::vector<RunPeriod_t> getRunPeriods(TString era);
TString assignRunPeriod(std::vector<RunPeriod_t> &runPeriods,TRandom *rand=0);

std::pair<std::map<Int_t,Float_t>, TH1F *> parseLumiInfo(TString era);
std::map<Int_t,Float_t> lumiPerRun(TString era="era2016");

//pileup weighting
std::vector<TGraph *> getPileupWeights(TString era, TH1 *genPU, TString period = "");
std::map<TString, std::vector<TGraph *> > getPileupWeightsMap(TString era, TH1 *genPU);

//apply jec uncertainty
void applyJetCorrectionUncertainty(MiniEvent_t &ev, JetCorrectionUncertainty *jecUnc, TString jecVar, TString direction);
void applyJetCorrectionUncertainty(TLorentzVector &jp4,JetCorrectionUncertainty *jecUnc,TString direction);

//apply jet energy resolutions
void smearJetEnergies(MiniEvent_t &ev, std::string option = "central");
void smearJetEnergies(MiniEvent_t &ev, JME::JetResolution* jer, std::string option = "central");
void smearJetEnergy(TLorentzVector &jp4, float genJet_pt,std::string option = "central");
void smearJetEnergyStochastic(TLorentzVector &jp4, TRandom* random, double resolution, std::string option = "central");

//see working points in https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation80XReReco
void addBTagDecisions(MiniEvent_t &ev,float wp=0.8484,float wpl=0.8484);

//details in https://twiki.cern.ch/twiki/bin/view/CMS/BTagCalibration
void updateBTagDecisions(MiniEvent_t &ev, 
				std::map<BTagEntry::JetFlavor,BTagCalibrationReader *> &btvsfReaders,
				std::map<BTagEntry::JetFlavor, TGraphAsymmErrors*> &expBtagEff, 
				std::map<BTagEntry::JetFlavor, TGraphAsymmErrors*> &expBtagEffPy8, 
				BTagSFUtil *myBTagSFUtil, 
				std::string optionbc = "central", std::string optionlight = "central");

//details in https://twiki.cern.ch/twiki/bin/view/CMS/BTagCalibration
std::map<BTagEntry::JetFlavor,BTagCalibrationReader *> getBTVcalibrationReaders(TString era,
										BTagEntry::OperatingPoint btagOP=BTagEntry::OP_MEDIUM, 
										TString period = "");
std::map<TString, std::map<BTagEntry::JetFlavor,BTagCalibrationReader *> > getBTVcalibrationReadersMap(TString era,
												       BTagEntry::OperatingPoint btagOP=BTagEntry::OP_MEDIUM); 

//the expections are created with the script scripts/saveExpectedBtagEff.py (cf README)
std::map<BTagEntry::JetFlavor, TGraphAsymmErrors *> readExpectedBtagEff(TString era,TString btagExpPostFix="");

//b fragmentation, see https://gitlab.cern.ch/CMS-TOPPAG/BFragmentationAnalyzer
std::map<TString, TGraph*> getBFragmentationWeights(TString era);
double computeBFragmentationWeight(MiniEvent_t &ev, TGraph* wgtGr);
std::map<TString, std::map<int, double> > getSemilepBRWeights(TString era);
double computeSemilepBRWeight(MiniEvent_t &ev, std::map<int, double> corr, int pid = 0, bool useabs = true);

//apply tracking efficiency
void applyTrackingEfficiencySF(MiniEvent_t &ev, double sf);

#endif
