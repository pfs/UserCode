#ifndef _correction_tools_h_
#define _correction_tools_h_

#include <vector>

#include "TString.h"

#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TopLJets2015/TopAnalysis/interface/BtagUncertaintyComputer.h"

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"

MiniEvent_t smearJetEnergies(MiniEvent_t ev, TString option = "") {
  for (int k = 0; k < ev.nj; k++) {
    TLorentzVector jp4;
    jp4.SetPtEtaPhiM(ev.j_pt[k],ev.j_eta[k],ev.j_phi[k],ev.j_mass[k]);

    //smear jet energy resolution for MC
    float genJet_pt(0);
    if(ev.j_g[k]>-1) genJet_pt = ev.g_pt[ ev.j_g[k] ];
    if(!ev.isData && genJet_pt>0) {
      float jerSmear=getJetResolutionScales(jp4.Pt(),jp4.Eta(),genJet_pt)[0];
      jp4 *= jerSmear;
      ev.j_pt[k]   = jp4.Pt();
      ev.j_eta[k]  = jp4.Eta();
      ev.j_phi[k]  = jp4.Phi();
      ev.j_mass[k] = jp4.M();
    }
  }
  
  return ev;
}

MiniEvent_t addBTagDecisions(MiniEvent_t ev) {
  for (int k = 0; k < ev.nj; k++) {
    ev.j_btag[k] = (ev.j_csv[k] > 0.800);
  }
  
  return ev;
}

MiniEvent_t updateBTagDecisions(MiniEvent_t ev, BTagCalibrationReader* sfbReader, BTagCalibrationReader* sflReader, std::map<TString, TGraphAsymmErrors*> expBtagEff, std::map<TString, TGraphAsymmErrors*> expBtagEffPy8, BTagSFUtil* myBTagSFUtil, TString option = "") {
  for (int k = 0; k < ev.nj; k++) {
    TLorentzVector jp4;
    jp4.SetPtEtaPhiM(ev.j_pt[k],ev.j_eta[k],ev.j_phi[k],ev.j_mass[k]);

    bool isBTagged(ev.j_btag[k]);
    if(!ev.isData) {
      float jptForBtag(jp4.Pt()>1000. ? 999. : jp4.Pt()), jetaForBtag(fabs(jp4.Eta()));
      float expEff(1.0), jetBtagSF(1.0);
      if(abs(ev.j_hadflav[k])==4) {         
        expEff    = expBtagEff["c"]->Eval(jptForBtag); 
        jetBtagSF = sfbReader->eval( BTagEntry::FLAV_C, jetaForBtag, jptForBtag);
        jetBtagSF *= expEff>0 ? expBtagEffPy8["c"]->Eval(jptForBtag)/expBtagEff["c"]->Eval(jptForBtag) : 0.;
      }
      else if(abs(ev.j_hadflav[k])==5) { 
        expEff    = expBtagEff["b"]->Eval(jptForBtag); 
        jetBtagSF = sfbReader->eval( BTagEntry::FLAV_B, jetaForBtag, jptForBtag);
        jetBtagSF *= expEff>0 ? expBtagEffPy8["b"]->Eval(jptForBtag)/expBtagEff["b"]->Eval(jptForBtag) : 0.;
      }
      else {
        expEff    = expBtagEff["udsg"]->Eval(jptForBtag);
        jetBtagSF = sflReader->eval( BTagEntry::FLAV_UDSG, jetaForBtag, jptForBtag);
        jetBtagSF *= expEff> 0 ? expBtagEffPy8["udsg"]->Eval(jptForBtag)/expBtagEff["udsg"]->Eval(jptForBtag) : 0.;
      }
      
      //updated b-tagging decision with the data/MC scale factor
      ev.j_btag[k] = myBTagSFUtil->modifyBTagsWithSF(isBTagged, jetBtagSF, expEff);
    }
  }
  
  return ev;
}

#endif
