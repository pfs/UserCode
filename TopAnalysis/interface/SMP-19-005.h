#ifndef _SMP19005_H_
#define _SMP19005_H_

#include "PhysicsTools/CandUtils/interface/EventShapeVariables.h"
#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TopLJets2015/TopAnalysis/interface/SelectionTools.h"
#include "TTree.h"


/**
   @short steers the event selection and corrections
 */
void RunSMP19005(TString filename, TString outname,TH1F *normH, TH1F *genPU,TString era,float xsec,bool debug,bool skimTree);


/**
   @short wrap up the build of variables needed for the analysis and categorization in a class
 */

class VBFDiscriminatorInputs{
 public:
  
  int evcats[2];
  float wgt;
  float a_pt,        a_eta,        a_phi,       a_ystar;
  float lead_pt,     lead_eta,     lead_phi,    lead_qg,    lead_dr2a;
  float sublead_pt,  sublead_eta,  sublead_phi, sublead_qg, sublead_dr2a;
  int lead_flav, sublead_flav;
  float jj_pt,       jj_eta,       jj_phi,      jj_m,       jj_dr2a, jj_scalarht, jj_deta, jj_dphi, jj_seta;
  float ajj_pt,      ajj_eta,      ajj_phi,     ajj_m;
  float scalarht,    mht, isotropy,circularity,sphericity,aplanarity,C,D;
  float centj_pt,     centj_eta,  centj_phi, centj_ystar;
  int ncentj;
  
  VBFDiscriminatorInputs() { reset(); }

 VBFDiscriminatorInputs(const VBFDiscriminatorInputs &o) : wgt(o.wgt),
    a_pt(o.a_pt),a_eta(o.a_eta),a_phi(o.a_phi),a_ystar(o.a_ystar),
    lead_pt(o.lead_pt),lead_eta(o.lead_eta),lead_phi(o.lead_phi),lead_qg(o.lead_qg),lead_dr2a(o.lead_dr2a),
    sublead_pt(o.sublead_pt),sublead_eta(o.sublead_eta),sublead_phi(o.sublead_phi),sublead_qg(o.sublead_qg),sublead_dr2a(o.sublead_dr2a),
    lead_flav(o.lead_flav), sublead_flav(o.sublead_flav),
    jj_pt(o.jj_pt),jj_eta(o.jj_eta),jj_phi(o.jj_phi),jj_m(o.jj_m),jj_dr2a(o.jj_dr2a),jj_scalarht(o.jj_scalarht),jj_deta(o.jj_deta), jj_dphi(o.jj_dphi), jj_seta(o.jj_seta),
    ajj_pt(o.ajj_pt),ajj_eta(o.ajj_eta),ajj_phi(o.ajj_phi),ajj_m(o.ajj_m),
    scalarht(o.scalarht),isotropy(o.isotropy),circularity(o.circularity),sphericity(o.sphericity),aplanarity(o.aplanarity),C(o.C),D(o.D),
    centj_pt(o.centj_pt),centj_eta(o.centj_eta),centj_phi(o.centj_phi),centj_ystar(o.centj_ystar),
    ncentj(o.ncentj)
   {
   }
    
 VBFDiscriminatorInputs& operator=(VBFDiscriminatorInputs &o) 
   {
     if (this == &o) return *this;
     wgt = o.wgt; 
     for(size_t i=0; i<2; i++) evcats[i]=o.evcats[i];
     a_pt = o.a_pt; a_eta = o.a_eta; a_phi = o.a_phi; a_ystar = o.a_ystar; 
     lead_pt = o.lead_pt; lead_eta = o.lead_eta; lead_phi = o.lead_phi; lead_qg = o.lead_qg; lead_dr2a = o.lead_dr2a; 
     sublead_pt = o.sublead_pt; sublead_eta = o.sublead_eta; sublead_phi = o.sublead_phi; sublead_qg = o.sublead_qg; sublead_dr2a = o.sublead_dr2a; 
     lead_flav = o.lead_flav;  sublead_flav = o.sublead_flav;
     jj_deta = o.jj_deta;  jj_dphi = o.jj_dphi;
     jj_pt = o.jj_pt; jj_eta = o.jj_eta; jj_phi = o.jj_phi; jj_m = o.jj_m; jj_dr2a = o.jj_dr2a; jj_scalarht = o.jj_scalarht; jj_deta=o.jj_deta; jj_dphi=o.jj_dphi; jj_seta=o.jj_seta;
     ajj_pt = o.ajj_pt; ajj_eta = o.ajj_eta; ajj_phi = o.ajj_phi; ajj_m = o.ajj_m; 
     scalarht = o.scalarht; isotropy = o.isotropy; circularity = o.circularity; sphericity = o.sphericity; aplanarity = o.aplanarity; C = o.C; D = o.D; 
     ncentj = o.ncentj; centj_pt = o.centj_pt; centj_eta = o.centj_eta; centj_phi = o.centj_phi; centj_ystar = o.centj_ystar;
     return *this;
   }
 
 //filler method
 inline bool fill(TLorentzVector &a,
                  std::vector<Jet> &jets,
                  MiniEvent_t &ev,
                  bool passLowPtTrig,
                  bool passHighPtTrig,
                  float w=1.0,
                  bool isGenLevel=false) {
   
   reset();
   wgt=w;
   
   //jet selection
   std::vector<Jet> cleanJets, cleanTagJets;
   for(size_t ij=0; ij<jets.size(); ij++){
     if( jets[ij].DeltaR(a)<0.4) continue;
     cleanJets.push_back( jets[ij] );
     if( jets[ij].Pt()<40 ) continue;
     cleanTagJets.push_back( jets[ij] );
   }
   if(cleanTagJets.size()<2) return false;
   if(cleanTagJets[0].Pt()<50) return false;
   TLorentzVector jj=cleanTagJets[0]+cleanTagJets[1];
   lead_pt=cleanTagJets[0].Pt();       lead_eta=cleanTagJets[0].Eta();        lead_phi=cleanTagJets[0].Phi();     lead_dr2a=cleanTagJets[0].DeltaR(a);
   sublead_pt=cleanTagJets[1].Pt();    sublead_eta=cleanTagJets[1].Eta();     sublead_phi=cleanTagJets[1].Phi();  sublead_dr2a=cleanTagJets[1].DeltaR(a);

   int lead_idx=cleanTagJets[0].getJetIndex();
   lead_flav = abs(ev.j_flav[lead_idx]);
   int sublead_idx=cleanTagJets[1].getJetIndex();
   sublead_flav = abs(ev.j_flav[sublead_idx]);

   jj_pt=jj.Pt();                      jj_eta=jj.Eta();                       jj_phi=jj.Phi();                    jj_m=jj.M();
   jj_dr2a=jj.DeltaR(a);               jj_scalarht=cleanTagJets[0].Pt()+cleanTagJets[1].Pt();
   float maxEta=TMath::Max(lead_eta,sublead_eta);
   float minEta=TMath::Min(lead_eta,sublead_eta);
   jj_deta=(maxEta-minEta);
   jj_seta=(maxEta+minEta); 
   jj_dphi=cleanTagJets[0].DeltaPhi( cleanTagJets[1] );      
   if(!isGenLevel) {
     lead_qg    = ev.j_qg[cleanTagJets[0].getJetIndex()];
     sublead_qg = ev.j_qg[cleanTagJets[1].getJetIndex()];
   }else {
     lead_qg    = 0;
     sublead_qg = 0;
   }      
   if(jj_m<200) return false;
   
   //boson selection
   a_pt=a.Pt();
   a_eta=a.Eta();
   a_phi=a.Phi();
   a_ystar=a_eta-0.5*jj_seta;
   if(isGenLevel || passLowPtTrig)      evcats[0] |= 1;
   if(a.Pt()>75 && fabs(a.Eta())<1.442) evcats[0] |= 2;
   if(isGenLevel || passHighPtTrig)     evcats[1] |= 1;
   if(a.Pt()>200)                       evcats[1] |= 2;
   if(evcats[0]==0 && evcats[1]==0) return false;
   
   //hard process object candidates
   TLorentzVector ajj=a+jj;
   ajj_pt=ajj.Pt(); ajj_eta=ajj.Eta(); ajj_phi=ajj.Phi();  ajj_m=ajj.M();
   scalarht=a.Pt()+jj_scalarht;
   std::vector<math::XYZVector> inputVectors;
   inputVectors.push_back( math::XYZVector(a.Px(),a.Py(),a.Pz()) );
   for(size_t ij=0; ij<2; ij++) {
     inputVectors.push_back( math::XYZVector(cleanTagJets[ij].Px(),cleanTagJets[ij].Py(),cleanTagJets[ij].Pz()) );
   }
   EventShapeVariables esv(inputVectors);
   isotropy    = esv.isotropy();
   circularity = esv.circularity();
   sphericity  = ajj.Pt()/scalarht;
   aplanarity  = esv.aplanarity();
   C           = esv.C();
   D           = esv.D();
   
   //extra radiation activity
   ncentj=0;
   centj_pt=0; centj_eta=0; centj_phi=0; centj_ystar=0;
   for(size_t ij=0; ij<cleanJets.size(); ij++) {
     if(cleanJets[ij].Pt()<20 ) continue;
     if(cleanJets[ij].Eta()<minEta+0.2) continue;
     if(cleanJets[ij].Eta()>maxEta-0.2) continue;
     ncentj+=1;
     if(ncentj>1) continue;
     centj_pt=cleanJets[ij].Pt();
     centj_eta=cleanJets[ij].Eta();
     centj_phi=cleanJets[ij].Phi();
     centj_ystar=centj_eta-0.5*jj_seta;
   }
   
   return true;
 }

 void assignTo(TTree *t) {
   t->Branch("evcats",(this->evcats),"evcats[2]/I");
   t->Branch("wgt",&(this->wgt),"wgt/F");
   t->Branch("a_pt",&(this->a_pt),"a_pt/F");
   t->Branch("a_eta",&(this->a_eta),"a_eta/F");
   t->Branch("a_phi",&(this->a_phi),"a_phi/F");
   t->Branch("a_ystar",&(this->a_ystar),"a_ystar/F");
   t->Branch("lead_pt",&(this->lead_pt),"lead_pt/F");
   t->Branch("lead_eta",&(this->lead_eta),"lead_eta/F");
   t->Branch("lead_phi",&(this->lead_phi),"lead_phi/F");
   t->Branch("lead_qg",&(this->lead_qg),"lead_qg/F");
   t->Branch("lead_dr2a",&(this->lead_dr2a),"lead_dr2a/F");
   t->Branch("lead_flav",&(this->lead_flav),"lead_flav/I");
   t->Branch("sublead_pt",&(this->sublead_pt),"sublead_pt/F");
   t->Branch("sublead_eta",&(this->sublead_eta),"sublead_eta/F");
   t->Branch("sublead_phi",&(this->sublead_phi),"sublead_phi/F");
   t->Branch("sublead_qg",&(this->sublead_qg),"sublead_qg/F");
   t->Branch("sublead_dr2a",&(this->sublead_dr2a),"sublead_dr2a/F");
   t->Branch("sublead_flav",&(this->sublead_flav),"sublead_flav/I");
   t->Branch("jj_pt",&(this->jj_pt),"jj_pt/F");
   t->Branch("jj_eta",&(this->jj_eta),"jj_eta/F");
   t->Branch("jj_phi",&(this->jj_phi),"jj_phi/F");
   t->Branch("jj_m",&(this->jj_m),"jj_m/F");
   t->Branch("jj_dr2a",&(this->jj_dr2a),"jj_dr2a/F");
   t->Branch("jj_scalarht",&(this->jj_scalarht),"jj_scalarht/F");
   t->Branch("jj_deta",&(this->jj_deta),"jj_deta/F");
   t->Branch("jj_dphi",&(this->jj_dphi),"jj_dphi/F");
   t->Branch("jj_seta",&(this->jj_seta),"jj_seta/F");
   t->Branch("ajj_pt",&(this->ajj_pt),"ajj_pt/F");
   t->Branch("ajj_eta",&(this->ajj_eta),"ajj_eta/F");
   t->Branch("ajj_phi",&(this->ajj_phi),"ajj_phi/F");
   t->Branch("ajj_m",&(this->ajj_m),"ajj_m/F");
   t->Branch("scalarht",&(this->scalarht),"scalarht/F");
   t->Branch("mht",&(this->mht),"mht/F");
   t->Branch("isotropy",&(this->isotropy),"isotropy/F");
   t->Branch("circularity",&(this->circularity),"circularity/F");
   t->Branch("sphericity",&(this->sphericity),"sphericity/F");
   t->Branch("aplanarity",&(this->aplanarity),"aplanarity/F");
   t->Branch("C",&(this->C),"C/F");
   t->Branch("D",&(this->D),"D/F");
   t->Branch("ncentj",&(this->ncentj),"ncentj/I");
   t->Branch("centj_pt",&(this->centj_pt),"centj_pt/F");
   t->Branch("centj_eta",&(this->centj_eta),"centj_eta/F");
   t->Branch("centj_phi",&(this->centj_phi),"centj_phi/F");
   t->Branch("centj_ystar",&(this->centj_ystar),"centj_ystar/F");
 }

 private:

 inline void reset() 
 {
   wgt=0;
   for(size_t i=0; i<2; i++) evcats[i]=0;
   a_pt = 0; a_eta = 0; a_phi = 0; a_ystar = 0;
   lead_pt = 0; lead_eta = 0; lead_phi = 0; lead_qg = 0; lead_dr2a = 0; lead_flav=0;
   sublead_pt = 0; sublead_eta = 0; sublead_phi = 0; sublead_qg = 0; sublead_dr2a = 0; sublead_flav=0;
   jj_deta = 0;  jj_dphi = 0;
   jj_pt = 0; jj_eta = 0; jj_phi = 0; jj_m = 0; jj_dr2a = 0; jj_scalarht = 0; jj_deta=0; jj_dphi=0; jj_seta=0;
   ajj_pt = 0; ajj_eta = 0; ajj_phi = 0; ajj_m = 0;
   scalarht = 0; isotropy = 0; circularity = 0; sphericity = 0; aplanarity = 0; C = 0; D = 0; 
   ncentj = 0; centj_pt = 0; centj_eta = 0; centj_phi = 0; centj_ystar = 0;  
 }

};

#endif
