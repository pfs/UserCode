#ifndef BTagSFUtil_lite_h
#define BTagSFUtil_lite_h

#include <Riostream.h>
#include "TRandom3.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"

//details about calibration and working points in
// https://twiki.cern.ch/twiki/bin/view/CMS/BTagCalibration
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X

class BTagSFUtil{

 public:
    
  BTagSFUtil(TString era="era2017", BTagEntry::OperatingPoint btagOp=BTagEntry::OperatingPoint::OP_MEDIUM, TString btagExp="", int seed=0 );
  ~BTagSFUtil();

  void addBTagDecisions(MiniEvent_t &ev,float wp=0.4941,float wpl=0.4941);
  void updateBTagDecisions(MiniEvent_t &ev, std::string optionbc = "central", std::string optionlight = "central");
  void modifyBTagsWithSF( bool& isBTagged, float Btag_SF = 0.98, float Btag_eff = 1.0);  
  bool applySF(bool& isBTagged, float Btag_SF = 0.98, float Btag_eff = 1.0);
  
 private:
  
  void startBTVcalibrationReaders(TString era,BTagEntry::OperatingPoint btagOp);
  void readExpectedBtagEff(TString era,BTagEntry::OperatingPoint btagOp,TString btagExp);
  
  TRandom3* rand_;
  std::map<BTagEntry::JetFlavor, TGraphAsymmErrors *> expBtagEff_;
  std::map<BTagEntry::JetFlavor,BTagCalibrationReader *> btvCalibReaders_;
};

#endif
