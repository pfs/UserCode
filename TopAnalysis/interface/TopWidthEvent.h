#ifndef _topwidthevent_h_
#define _topwidthevent_h_

#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"

/**
   @short an helper class to store the information about a lepton-b pair candidate
 */
class LeptonBJetPair : public TLorentzVector {
 public:
 LeptonBJetPair(TLorentzVector p4,float drlb,int lid,int bid) : TLorentzVector(p4) {
    drlb_=drlb;
    lid_=lid;
    bid_=bid;
  }
 LeptonBJetPair( const LeptonBJetPair &lb) 
   : TLorentzVector(lb.Px(),lb.Py(),lb.Pz(),lb.E()), drlb_(lb.drlb_), lid_(lb.lid_), bid_(lb.bid_) {}
  float getDR( ) { return drlb_; }
  int leptonId() { return lid_; }
  int bId()      { return bid_; }
 private:
  float drlb_;
  int lid_,bid_;
};

class TopWidthEvent {

 public:

  TopWidthEvent(std::vector<Particle> &leptons,std::vector<Jet> &jets);
  std::vector<LeptonBJetPair> getPairs() { return lbPairs_; }
  ~TopWidthEvent() { }

  //access directly
  TString cat;
  unsigned int dilcode;
  float mll,ptll,l1pt,l1eta,l2pt,l2eta,njets,nbjets,j1pt,j1eta,j2pt,j2eta;
  std::vector<size_t> selJetsIdx;
 private:
  std::vector<LeptonBJetPair> lbPairs_;
  void initSelectionCuts();
  float leadLeptonPt_,subLeadLeptonPt_,maxLeptonEta_,zMassWindow_,minMll_,jetPt_,maxJetEta_;

};


#endif
