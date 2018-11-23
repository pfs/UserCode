#ifndef _VBFDISCRIMINATORINPUTS_H_
#define _VBFDISCRIMINATORINPUTS_H_

#include "PhysicsTools/CandUtils/interface/EventShapeVariables.h"
#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TopLJets2015/TopAnalysis/interface/SelectionTools.h"

namespace vbf{

  class DiscriminatorInputs{
  public:
    
    //variables to be used in a small tree/MVA
    float leadj_pt,leadj_gawidth,leadj_c2_02,leadj_qg;
    float subleadj_pt,subleadj_gawidth,subleadj_c2_02,subleadj_qg;
    float centraleta,forwardeta;
    float mjj,detajj,dphijj,jjpt,jjetas,ystar,balance,relbpt,dphibjj,dphivj0,dphivj1;
    float isotropy,circularity,sphericity,aplanarity,C,D;
    float scalarht,mht;
    float ncentj;
    std::vector<float> dphivcentj,centjy;
    
    //CTOR
  DiscriminatorInputs(): 
    leadj_pt(0),leadj_gawidth(-99),leadj_c2_02(-99),leadj_qg(-99),
      subleadj_pt(0),subleadj_gawidth(-99),subleadj_c2_02(-99),subleadj_qg(-99),
      centraleta(-99),forwardeta(-99),
      mjj(0),detajj(-99),dphijj(-99),jjpt(0),jjetas(-99),ystar(-99),balance(-99),relbpt(-99),dphibjj(-99),dphivj0(-99),dphivj1(-99),
      isotropy(-99),circularity(-99),sphericity(-99),aplanarity(-99),C(-99),D(-99),scalarht(0),mht(0),
      ncentj(0), dphivcentj(10,-99.), centjy(10,-99.){
    }

    //copy constructor
  DiscriminatorInputs(const DiscriminatorInputs &o): 
    leadj_pt(o.leadj_pt),leadj_gawidth(o.leadj_gawidth),leadj_c2_02(o.leadj_c2_02),leadj_qg(o.leadj_qg),
      subleadj_pt(o.subleadj_pt),subleadj_gawidth(o.subleadj_gawidth),subleadj_c2_02(o.subleadj_c2_02),subleadj_qg(o.subleadj_qg),
      centraleta(o.centraleta),forwardeta(o.forwardeta),
      mjj(o.mjj),detajj(o.detajj),dphijj(o.dphijj),jjpt(o.jjpt),jjetas(o.jjetas),ystar(o.ystar),balance(o.balance),relbpt(o.relbpt),dphibjj(o.dphibjj),dphivj0(o.dphivj0),dphivj1(o.dphivj1),
      isotropy(o.isotropy),circularity(o.circularity),sphericity(o.sphericity),aplanarity(o.aplanarity),C(o.C),D(o.D),
      scalarht(o.scalarht), mht(o.mht),
      ncentj(o.ncentj) {
      dphivcentj.resize(o.dphivcentj.size());
      centjy.resize(o.centjy.size());
      for(size_t ij=0; ij<dphivcentj.size(); ij++){
        dphivcentj[ij] = o.dphivcentj[ij];
        centjy[ij]      = o.centjy[ij];
      }
    }

    //assignment operator
    void assignValuesFrom(const DiscriminatorInputs& o)
      {
        leadj_pt         = o.leadj_pt;  
        leadj_gawidth    = o.leadj_gawidth;  
        leadj_c2_02      = o.leadj_c2_02;  
        leadj_qg         = o.leadj_qg;  
        subleadj_pt      = o.subleadj_pt;  
        subleadj_gawidth = o.subleadj_gawidth;  
        subleadj_c2_02   = o.subleadj_c2_02;  
        subleadj_qg      = o.subleadj_qg;  
        centraleta       = o.centraleta; 
        forwardeta       = o.forwardeta;
        mjj              = o.mjj;  
        detajj           = o.detajj;  
        dphijj           = o.dphijj;  
        jjpt             = o.jjpt;  
        jjetas           = o.jjetas;  
        ystar            = o.ystar;  
        balance          = o.balance;  
        relbpt           = o.relbpt;  
        dphibjj          = o.dphibjj;  
        dphivj0          = o.dphivj0;  
        dphivj1          = o.dphivj1;  
        isotropy         = o.isotropy;  
        circularity      = o.circularity;  
        sphericity       = o.sphericity;  
        aplanarity       = o.aplanarity;  
        C                = o.C;  
        D                = o.D;
        scalarht         = o.scalarht; 
        mht              = o.mht;
        ncentj           = o.ncentj;
        for(size_t ij=0; ij<dphivcentj.size(); ij++){
          dphivcentj[ij] = o.dphivcentj[ij];
          centjy[ij]      = o.centjy[ij];
        } 
      }

    //filler method
    inline void fillDiscriminatorVariables(TLorentzVector &boson,std::vector<Jet> &jets,MiniEvent_t &ev) {

      if(jets.size()>0){
        leadj_pt      = jets[0].Pt();
        leadj_gawidth = ev.j_gawidth[jets[0].getJetIndex()];
        leadj_c2_02   = ev.j_c2_02[jets[0].getJetIndex()];
        leadj_qg      = ev.j_qg[jets[0].getJetIndex()];
        dphivj0       = fabs(jets[0].DeltaPhi(boson));

        centraleta    = jets[0].Eta();
      }
    
      if(jets.size()>=2){
        subleadj_pt      = jets[1].Pt();
        subleadj_gawidth = ev.j_gawidth[jets[1].getJetIndex()];
        subleadj_c2_02   = ev.j_c2_02[jets[1].getJetIndex()];
        subleadj_qg      = ev.j_qg[jets[1].getJetIndex()];
        dphivj1          = fabs(jets[1].DeltaPhi(boson));

        centraleta       = min(fabs(jets[0].Eta()),fabs(jets[1].Eta()));
        forwardeta       = max(fabs(jets[0].Eta()),fabs(jets[1].Eta()));

        TLorentzVector jj(jets[0]+jets[1]);
        mjj     = jj.M();
        detajj  = fabs(jets[0].Eta()-jets[1].Eta());
        dphijj  = jets[0].DeltaPhi(jets[1]);
        jjpt    = jj.Pt();
        jjetas  = jets[0].Eta()*jets[1].Eta();
        ystar   = boson.Rapidity()-0.5*(jets[0].Rapidity()+jets[1].Rapidity());
        balance = (boson+jj).Pt();
        relbpt  = (jets[0].Pt()+jets[1].Pt())/boson.Pt();
        dphibjj = boson.DeltaPhi( jets[0]+jets[1] );    
      }


      //tag jet + addition jet activity
      scalarht = 0.;
      TLorentzVector mhtP4(0,0,0,0);
      ncentj=0;
      std::fill(dphivcentj.begin(),dphivcentj.end(),-99);
      std::fill(centjy.begin(),centjy.end(),-99);
      for(size_t ij=0; ij<jets.size(); ij++){
        scalarht += jets[ij].Pt();
        mhtP4 += jets[ij];
        if(ij<2) continue;

        float dy = fabs(jets[0].Rapidity() - jets[1].Rapidity())/2;
        float sumy = (jets[0].Rapidity() + jets[1].Rapidity())/2;
        if(fabs(jets[ij].Rapidity() - sumy) > dy) continue;
        dphivcentj[ij-2]=fabs(jets[ij].DeltaPhi(boson));
        centjy[ij-2]=jets[ij].Rapidity();
        ncentj++;
      }
      mht = mhtP4.Pt();

      //event shapes
      std::vector<math::XYZVector> inputVectors;
      inputVectors.push_back( math::XYZVector(boson.Px(),boson.Py(),boson.Pz()) );
      for(size_t ij=0; ij<min(size_t(2),jets.size());ij++) {
        inputVectors.push_back( math::XYZVector(jets[ij].Px(),jets[ij].Py(),jets[ij].Pz()) );
      }
      EventShapeVariables esv(inputVectors);
      isotropy    = esv.isotropy();
      circularity = esv.circularity();
      sphericity  = esv.sphericity(1.);
      aplanarity  = esv.aplanarity(1.);
      C           = esv.C(1.);
      D           = esv.D(1.);
    }
  };

}

#endif
