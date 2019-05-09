#include "TopLJets2015/TopAnalysis/interface/TopWidthEvent.h"
#include <algorithm>
#include <iostream>

using namespace std;

//
TopWidthEvent::TopWidthEvent(std::vector<Particle> &leptons,std::vector<Jet> &jets) : cat(""), dilcode(0) {

  initSelectionCuts();

  if(leptons.size()<2 || jets.size()<2) return;
  
  //lepton selection
  if(leptons[0].Pt()<leadLeptonPt_    || fabs(leptons[0].Eta())>maxLeptonEta_) return;
  if(leptons[1].Pt()<subLeadLeptonPt_ || fabs(leptons[1].Eta())>maxLeptonEta_) return;
  l1pt  = leptons[0].Pt();
  l1eta = leptons[0].Eta();
  l2pt  = leptons[1].Pt();
  l2eta = leptons[1].Eta();
  
  //dilepton selection
  if(leptons[0].charge()*leptons[1].charge()>0) return;
  mll   = ( (leptons[0]+leptons[1]).M() );
  ptll  = ( (leptons[0]+leptons[1]).Pt() );
  if(mll<20) return;
  bool isZ( fabs(mll-91)<15 );
  
  //jet selection
  std::vector<Jet> bJets,lJets;
  selJetsIdx.clear();
  size_t i(0);
  for(auto j : jets) {
    i++;
    if(j.Pt()<jetPt_) continue;
    if(fabs(j.Eta())>maxJetEta_) continue;    
    bool isBTagged(j.flavor()==5);
    if(isBTagged) { bJets.push_back(j); selJetsIdx.push_back(i-1); }
    else          { lJets.push_back(j); }
  }
  nbjets=bJets.size();
  njets=nbjets+lJets.size();
  if(njets<2) return;
  if(nbjets<1) return;

  j1pt  = bJets[0].Pt();
  j1eta = bJets[0].Eta();
  j2pt  = bJets.size()>1 ? bJets[1].Pt()  : lJets[0].Pt();
  j2eta = bJets.size()>1 ? bJets[1].Eta() : lJets[0].Eta();
  
  //lepton-b pairing
  lbPairs_.clear();
  if(bJets.size()==1) {
    
    //choose pair with leading ptlb
    LeptonBJetPair lb1(leptons[0]+bJets[0],leptons[0].DeltaR(bJets[0]),leptons[0].originalReference(),bJets[0].getJetIndex());
    LeptonBJetPair lb2(leptons[1]+bJets[0],leptons[1].DeltaR(bJets[0]),leptons[1].originalReference(),bJets[1].getJetIndex());
    lbPairs_.push_back( lb1.Pt()>lb2.Pt() ? lb1 : lb2 );
  }
  else {
    
    //sort pairs by decreasing ptlb
    std::vector< LeptonBJetPair > allPairs;
    for(size_t il=0; il<2; il++)
      for(size_t ib=0; ib<2; ib++) {
        LeptonBJetPair lb(leptons[il]+bJets[ib],leptons[il].DeltaR(bJets[ib]),leptons[il].originalReference(),bJets[ib].getJetIndex());
        allPairs.push_back( lb );
      }
    std::sort(allPairs.begin(),allPairs.end(),
              [](const LeptonBJetPair& a, const LeptonBJetPair& b) -> bool
              {
                return a.Pt() > b.Pt();
              });
    
    lbPairs_.push_back( allPairs[0] );
    for(size_t ipair=1; ipair<4; ipair++) {
      if( allPairs[0].leptonId() ==  allPairs[ipair].leptonId() ) continue;
      if( allPairs[0].bId()      ==  allPairs[ipair].bId() ) continue;
      lbPairs_.push_back(allPairs[ipair]);
      break;
    }
  }

  //assign category
  dilcode=(abs(leptons[0].id()*leptons[1].id()));
  if(dilcode==11*11) cat=isZ ? "zee":"ee";
  if(dilcode==13*13) cat=isZ ? "zmm":"mm";
  if(dilcode==11*13) cat="em";
}


//
void TopWidthEvent::initSelectionCuts() {
  leadLeptonPt_    = 30.;
  subLeadLeptonPt_ = 20.;
  maxLeptonEta_    = 2.4;
  zMassWindow_     = 15.;
  minMll_          = 20.;
  jetPt_           = 30.;
  maxJetEta_       = 2.4;
} 
