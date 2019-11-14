#ifndef _VBFDISCRIMINATORINPUTS_H_
#define _VBFDISCRIMINATORINPUTS_H_

#include "PhysicsTools/CandUtils/interface/EventShapeVariables.h"
#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TopLJets2015/TopAnalysis/interface/SelectionTools.h"

namespace vbf{

  class DiscriminatorInputs{
  public:
    
    //variables to be used in a small tree/MVA
    float leadj_pt,leadj_eta,leadj_gawidth,leadj_c2_02,leadj_qg;
    float subleadj_pt,subleadj_eta,subleadj_gawidth,subleadj_c2_02,subleadj_qg;
    float centraleta,forwardeta;
    float mjj,detajj,dphijj,jjpt,jjetas,ystar,balance,relbpt,dphibjj,dphivj0,dphivj1,dphivj2,dphivj3;
    float isotropy,circularity,sphericity,aplanarity,C,D;

    float scalarht,mht,dijetht,syspt;
    float ncentj;
   
    //addition production and colour flow
    float cosqj1, cosqjj;
    float beta_v_j2,  beta_j1_j2, beta_v_j3, beta_closej_j3;

    std::vector<float> dphivcentj,centjy;
    
    //CTOR
  DiscriminatorInputs(): 
    leadj_pt(0),leadj_eta(0),leadj_gawidth(-99),leadj_c2_02(-99),leadj_qg(-99),
      subleadj_pt(0),subleadj_eta(0),subleadj_gawidth(-99),subleadj_c2_02(-99),subleadj_qg(-99),
      centraleta(-99),forwardeta(-99),
      mjj(0),detajj(-99),dphijj(-99),jjpt(0),jjetas(-99),ystar(-99),balance(-99),relbpt(-99),dphibjj(-99),dphivj0(-99),dphivj1(-99),dphivj2(-99),dphivj3(-99),

      isotropy(-99),circularity(-99),sphericity(-99),aplanarity(-99),C(-99),D(-99),scalarht(0),mht(0),dijetht(0),syspt(0),

      ncentj(0),
      cosqj1(-99),  cosqjj(-99),
      beta_v_j2(-99),  beta_j1_j2(-99), beta_v_j3(-99), beta_closej_j3(-99),
      dphivcentj(10,-99.), centjy(10,-99.)
      {
      }

    //copy constructor
  DiscriminatorInputs(const DiscriminatorInputs &o): 
    leadj_pt(o.leadj_pt),leadj_gawidth(o.leadj_gawidth),leadj_c2_02(o.leadj_c2_02),leadj_qg(o.leadj_qg),
      subleadj_pt(o.subleadj_pt),subleadj_gawidth(o.subleadj_gawidth),subleadj_c2_02(o.subleadj_c2_02),subleadj_qg(o.subleadj_qg),
      centraleta(o.centraleta),forwardeta(o.forwardeta),
      mjj(o.mjj),detajj(o.detajj),dphijj(o.dphijj),jjpt(o.jjpt),jjetas(o.jjetas),ystar(o.ystar),balance(o.balance),relbpt(o.relbpt),dphibjj(o.dphibjj),dphivj0(o.dphivj0),dphivj1(o.dphivj1),dphivj2(o.dphivj2),dphivj3(o.dphivj3),
      isotropy(o.isotropy),circularity(o.circularity),sphericity(o.sphericity),aplanarity(o.aplanarity),C(o.C),D(o.D),
      scalarht(o.scalarht), mht(o.mht),dijetht(o.dijetht),syspt(o.syspt),
      ncentj(o.ncentj),
      cosqj1(o.cosqj1),  cosqjj(o.cosqjj), 
      beta_v_j2(o.beta_v_j2),  beta_j1_j2(o.beta_j1_j2), beta_v_j3(o.beta_v_j3), beta_closej_j3(o.beta_closej_j3)      
      {
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
	dphivj2          = o.dphivj2;  
        dphivj3          = o.dphivj3;  
        isotropy         = o.isotropy;  
        circularity      = o.circularity;  
        sphericity       = o.sphericity;  
        aplanarity       = o.aplanarity;  
        C                = o.C;  
        D                = o.D;
        scalarht         = o.scalarht;
	//scalarGht         = o.scalarGht; 
        mht              = o.mht;
	dijetht          = o.dijetht;
	syspt            = o.syspt;
        ncentj           = o.ncentj;
        cosqj1           = o.cosqj1;
        cosqjj           = o.cosqjj;        
        beta_v_j2        = o.beta_v_j2;
        beta_j1_j2       = o.beta_j1_j2;
        beta_v_j3        = o.beta_v_j3;
        beta_closej_j3   = o.beta_closej_j3;      
        for(size_t ij=0; ij<dphivcentj.size(); ij++){
          dphivcentj[ij] = o.dphivcentj[ij];
          centjy[ij]      = o.centjy[ij];
        } 
      }

    //filler method
    inline void fillDiscriminatorVariables(TLorentzVector &boson,std::vector<Jet> &jets,MiniEvent_t &ev) {
     
      if(jets.size()>0){
        leadj_pt      = jets[0].Pt();
	leadj_eta      = fabs(jets[0].Eta());
        leadj_gawidth = ev.j_gawidth[jets[0].getJetIndex()];
        leadj_c2_02   = ev.j_c2_02[jets[0].getJetIndex()];
        leadj_qg      = ev.j_qg[jets[0].getJetIndex()];
        dphivj0       = (jets[0].DeltaPhi(boson));

        centraleta    = jets[0].Eta();
      }
    
      if(jets.size()>=2){
        subleadj_pt      = jets[1].Pt();
	subleadj_eta = fabs(jets[1].Eta());
        subleadj_gawidth = ev.j_gawidth[jets[1].getJetIndex()];
        subleadj_c2_02   = ev.j_c2_02[jets[1].getJetIndex()];
        subleadj_qg      = ev.j_qg[jets[1].getJetIndex()];
        dphivj1          = (jets[1].DeltaPhi(boson));
	if (jets.size()>2) dphivj2 = fabs(jets[2].DeltaPhi(boson));
	if (jets.size()>3) dphivj3 = fabs(jets[3].DeltaPhi(boson));
        centraleta       = min(fabs(jets[0].Eta()),fabs(jets[1].Eta()));
        forwardeta       = max(fabs(jets[0].Eta()),fabs(jets[1].Eta()));
	dijetht  = jets[0].Pt()+jets[1].Pt();
	syspt    = dijetht+boson.Pt();
        TLorentzVector jj(jets[0]+jets[1]);
        mjj     = jj.M();
        detajj  = fabs(jets[0].Eta()-jets[1].Eta());
	dphijj  = (jets[0].DeltaPhi(jets[1]));
        jjpt    = jj.Pt();
        jjetas  = jets[0].Eta()*jets[1].Eta();
        ystar   = boson.Rapidity()-0.5*(jets[0].Rapidity()+jets[1].Rapidity());
        balance = (boson+jj).Pt();
        relbpt  = (jets[0].Pt()+jets[1].Pt())/boson.Pt();
	dphibjj = (boson.DeltaPhi( jets[0]+jets[1] ));  
        cosqj1      = TMath::TanH( 0.5*(boson.Rapidity()-jets[0].Rapidity()) );
        cosqjj      = TMath::TanH( 0.5*(boson.Rapidity()-jj.Rapidity()) );
        beta_v_j2   = fabs(boson.DeltaPhi(jets[1]))/((boson.Eta()<0 ? -1 : 1)*(jets[1].Eta()-boson.Eta()+1e-6));
        beta_j1_j2  = fabs(jets[0].DeltaPhi(jets[1]))/((jets[0].Eta()<0 ? -1 : 1)*(jets[1].Eta()-jets[0].Eta()+1e-6));
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
	//	scalarGht += scalarht+boson.Pt();
        if(ij<2) continue;

        float dy = fabs(jets[0].Rapidity() - jets[1].Rapidity())/2;
        float sumy = (jets[0].Rapidity() + jets[1].Rapidity())/2;
        if(fabs(jets[ij].Rapidity() - sumy) > dy) continue;
        dphivcentj[ij-2]=fabs(jets[ij].DeltaPhi(boson));
        centjy[ij-2]=jets[ij].Rapidity();
        ncentj++;

        if(ncentj>1) continue;
        
        beta_v_j3       = fabs(boson.DeltaPhi(jets[ij]))/((boson.Eta()<0 ? -1 : 1)*(jets[ij].Eta()-boson.Eta()+1e-6));        
        int closeJ( jets[0].DeltaR(jets[ij])<jets[1].DeltaR(jets[ij]) ? 0 : 1 );
        beta_closej_j3  = fabs(jets[closeJ].DeltaPhi(jets[ij]))/((jets[closeJ].Eta()<0 ? -1 : 1)*(jets[ij].Eta()-jets[closeJ].Eta()+1e-6));

	
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
