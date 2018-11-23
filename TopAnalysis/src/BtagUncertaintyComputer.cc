#include "TopLJets2015/TopAnalysis/interface/BtagUncertaintyComputer.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TFile.h"

//
BTagSFUtil::BTagSFUtil(TString era,BTagEntry::OperatingPoint btagOp,TString btagExp, int seed) {

  readExpectedBtagEff(era,btagOp,btagExp);
  startBTVcalibrationReaders(era,btagOp);
  rand_ = new TRandom3(seed);

}

BTagSFUtil::~BTagSFUtil() {

  delete rand_;

}

//
void BTagSFUtil::addBTagDecisions(MiniEvent_t &ev,float wp,float wpl) {
  for (int k = 0; k < ev.nj; k++) {
    if (ev.j_hadflav[k] >= 4) ev.j_btag[k] = (ev.j_deepcsv[k] > wp);
    else                      ev.j_btag[k] = (ev.j_deepcsv[k] > wpl);
  }
}

//
void BTagSFUtil::updateBTagDecisions(MiniEvent_t &ev,std::string optionbc, std::string optionlight) {
  for (int k = 0; k < ev.nj; k++) {
    TLorentzVector jp4;
    jp4.SetPtEtaPhiM(ev.j_pt[k],ev.j_eta[k],ev.j_phi[k],ev.j_mass[k]);
    
    bool isBTagged(ev.j_btag[k]);
    if(!ev.isData) {
      float jptForBtag(jp4.Pt()>1000. ? 999. : jp4.Pt()), jetaForBtag(fabs(jp4.Eta()));
      float expEff(1.0), jetBtagSF(1.0);
      
      BTagEntry::JetFlavor hadFlav=BTagEntry::FLAV_UDSG;
      std::string option = optionlight;
      if(abs(ev.j_hadflav[k])==4) { hadFlav=BTagEntry::FLAV_C; option = optionbc; }
      if(abs(ev.j_hadflav[k])==5) { hadFlav=BTagEntry::FLAV_B; option = optionbc; }

      expEff    = expBtagEff_[hadFlav]->Eval(jptForBtag); 
      jetBtagSF = btvCalibReaders_[hadFlav]->eval_auto_bounds( option, hadFlav, jetaForBtag, jptForBtag);
      
      //updated b-tagging decision with the data/MC scale factor
      modifyBTagsWithSF(isBTagged, jetBtagSF, expEff);
      ev.j_btag[k] = isBTagged;
    }
  }
}


//
void BTagSFUtil::modifyBTagsWithSF(bool& isBTagged, float tag_SF, float tag_Eff){
  bool newBTag = isBTagged;
  newBTag = applySF(isBTagged, tag_SF, tag_Eff);
  isBTagged = newBTag;
}


bool BTagSFUtil::applySF(bool& isBTagged, float Btag_SF, float Btag_eff){
  
  bool newBTag = isBTagged;

  if (Btag_SF == 1) return newBTag; //no correction needed 

  //throw die
  float coin = rand_->Uniform(1.);    
  
  if(Btag_SF > 1){  // use this if SF>1

    if( !isBTagged ) {

      //fraction of jets that need to be upgraded
      float mistagPercent = (1.0 - Btag_SF) / (1.0 - (Btag_SF/Btag_eff) );

      //upgrade to tagged
      if( coin < mistagPercent ) {newBTag = true;}
    }

  }else{  // use this if SF<1
      
    //downgrade tagged to untagged
    if( isBTagged && coin > Btag_SF ) {newBTag = false;}

  }

  return newBTag;
}

//
void BTagSFUtil::startBTVcalibrationReaders(TString era,BTagEntry::OperatingPoint btagOP)
{
  //start the btag calibration
  TString btagUncUrl( era+"/DeepCSV_94XSF_V3_B_F.csv");
  gSystem->ExpandPathName(btagUncUrl);
  BTagCalibration btvcalib("DeepCSV",btagUncUrl.Data());

  //start calibration readers for b,c and udsg separately including the up/down variations
  btvCalibReaders_[BTagEntry::FLAV_B]=new BTagCalibrationReader(btagOP, "central", {"up", "down"});
  btvCalibReaders_[BTagEntry::FLAV_B]->load(btvcalib,BTagEntry::FLAV_B,"mujets");
  btvCalibReaders_[BTagEntry::FLAV_C]=new BTagCalibrationReader(btagOP, "central", {"up", "down"});
  btvCalibReaders_[BTagEntry::FLAV_C]->load(btvcalib,BTagEntry::FLAV_C,"mujets");
  btvCalibReaders_[BTagEntry::FLAV_UDSG]=new BTagCalibrationReader(btagOP, "central", {"up", "down"});
  btvCalibReaders_[BTagEntry::FLAV_UDSG]->load(btvcalib,BTagEntry::FLAV_UDSG,"incl");
}

//
void BTagSFUtil::readExpectedBtagEff(TString era,BTagEntry::OperatingPoint btagOp,TString btagExp)
{
  //open up the ROOT file with the expected efficiencies
  TString btagEffExpUrl(Form("%s/expectedBtagEff.root",era.Data()));
  gSystem->ExpandPathName(btagEffExpUrl);

  TFile *beffIn=TFile::Open(btagEffExpUrl);
  TString baseDir("DeepCSV_");
  if(btagOp==BTagEntry::OP_LOOSE) baseDir+="loose";
  if(btagOp==BTagEntry::OP_MEDIUM) baseDir+="medium";
  if(btagOp==BTagEntry::OP_TIGHT) baseDir+="tight";
  expBtagEff_[BTagEntry::FLAV_B]    = (TGraphAsymmErrors *)beffIn->Get(baseDir+"/b");
  expBtagEff_[BTagEntry::FLAV_C]    = (TGraphAsymmErrors *)beffIn->Get(baseDir+"/c");
  expBtagEff_[BTagEntry::FLAV_UDSG] = (TGraphAsymmErrors *)beffIn->Get(baseDir+"/udsg");
  beffIn->Close();
}
