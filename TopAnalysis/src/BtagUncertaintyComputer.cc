#include "TopLJets2015/TopAnalysis/interface/BtagUncertaintyComputer.h"

BTagSFUtil::BTagSFUtil( int seed ) {

  rand_ = new TRandom3(seed);

}

BTagSFUtil::~BTagSFUtil() {

  delete rand_;

}


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

