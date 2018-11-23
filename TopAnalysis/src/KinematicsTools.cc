#include "TopLJets2015/TopAnalysis/interface/KinematicsTools.h"
#include <iostream>


//
JetPullInfo_t getPullVector( MiniEvent_t &ev, int ijet)
{
  JetPullInfo_t result;
  result.n=0; result.nch=0;
  result.pull=TVector2(0,0);
  result.chPull=TVector2(0,0);
  
  //re-reconstruct the jet direction with the charged tracks
  TLorentzVector jet(0,0,0,0);
  jet.SetPtEtaPhiM(ev.j_pt[ijet], ev.j_eta[ijet], ev.j_phi[ijet], ev.j_mass[ijet]);
  TLorentzVector chargedJet(0,0,0,0);
  TLorentzVector constituent(0,0,0,0);
  std::vector<std::pair<TLorentzVector,bool> > allConstituents;
  unsigned int nCharged = 0;
  for(Int_t idx = 0; idx<ev.npf; idx++)
    {
      if(ev.pf_j[idx]!=ijet) continue;
      constituent.SetPtEtaPhiM( ev.pf_pt[idx], ev.pf_eta[idx], ev.pf_phi[idx], ev.pf_m[idx]);
      bool isCharged(abs(ev.pf_id[idx])==11 ||
		     abs(ev.pf_id[idx])==13 ||
		     abs(ev.pf_id[idx])==211 );
      allConstituents.push_back(std::make_pair(constituent,isCharged) );
      if(isCharged)
	{
	  chargedJet += constituent;
	  ++nCharged;      
	}
    }
  std::cout << std::endl;
  result.n=(Int_t) allConstituents.size();
  result.nch=nCharged;

  //stop here if <2 charged
  if( nCharged < 2 ) return result;

  //compute the pull
  double jetPt        = jet.Pt(),        jetPhi=jet.Phi(),                 jetRapidity=jet.Rapidity();
  double jetPtCharged = chargedJet.Pt(), jetPhiCharged = chargedJet.Phi(), jetRapidityCharged = chargedJet.Rapidity();
  TVector2 r(0,0);
  TVector2 pullAll(0,0);
  TVector2 pullCharged(0,0);
  for(size_t idx = 0; idx<allConstituents.size(); ++idx)
    {
      TLorentzVector &cp4=allConstituents[idx].first;
      bool &isCharged=allConstituents[idx].second;
      double constituentPt       = cp4.Pt();
      double constituentPhi      = cp4.Phi();
      double constituentRapidity = cp4.Rapidity();
      r.Set( constituentRapidity - jetRapidity, TVector2::Phi_mpi_pi( constituentPhi - jetPhi ) );
      pullAll += ( constituentPt / jetPt ) * r.Mod() * r;
      //calculate TVector using only charged tracks
      if( isCharged )
	r.Set( constituentRapidity - jetRapidityCharged, TVector2::Phi_mpi_pi( constituentPhi - jetPhiCharged ) );
      pullCharged += ( constituentPt / jetPtCharged ) * r.Mod() * r;
    }
  
  result.pull=pullAll;
  result.chPull=pullCharged;
  return result;
}

//
Float_t computeRecoilMT(TLorentzVector &a, TVector2 &h){
  float mt2_val = pow(a.Pt(),2);
  mt2_val += a.Pt()*((a.Vect().XYvector()+h).Mod());
  mt2_val += a.Vect().XYvector()*h;
  mt2_val *= 2;
  return TMath::Sqrt(mt2_val);
}

//
Float_t computeMT(TLorentzVector &a, TLorentzVector &b)
{
  return TMath::Sqrt(2*a.Pt()*b.Pt()*(1-TMath::Cos(a.DeltaPhi(b))));
}

Float_t computeMR(TLorentzVector &hem1, TLorentzVector &hem2){
  return sqrt(pow(hem1.P() + hem2.P(), 2) - pow(hem1.Pz() + hem2.Pz(), 2));
}

Float_t computeRsq(TLorentzVector &hem1, TLorentzVector &hem2, TLorentzVector &pfMet){
  Float_t mR = computeMR(hem1, hem2);
  Float_t term1 = pfMet.Pt()/2*(hem1.Pt() + hem2.Pt());
  Float_t term2 = pfMet.Px()/2*(hem1.Px() + hem2.Px()) + pfMet.Py()/2*(hem1.Py() + hem2.Py()); //dot product of MET with (p1T + p2T)
  Float_t mTR = sqrt(term1 - term2);
  return (mTR / mR) * (mTR / mR);
}

Float_t computePhiStar(TLorentzVector &lm,TLorentzVector &lp){
  Float_t dphill        = fabs(lm.DeltaPhi(lp));
  Float_t phi_acop      = TMath::Pi()-dphill;
  Float_t costhetastar  = TMath::TanH(0.5*(lm.Eta()-lp.Eta()));
  Float_t sin2thetastar = costhetastar>1 ? 0.0 : (1.0 - TMath::Sqrt(costhetastar));
  Float_t phistar       = TMath::Tan(0.5*phi_acop) * TMath::Sqrt( sin2thetastar );
  return phistar;
}

Float_t computeCosThetaStar(TLorentzVector &lm,TLorentzVector &lp){
  TLorentzVector dil(lm+lp);
  Float_t costhetaCS((dil.Pz()>0 ? 1.0 : -1.0)/dil.M());
  costhetaCS *= (lm.E() + lm.Pz()) * (lp.E() - lp.Pz()) - (lm.E() - lm.Pz()) * (lp.E() + lp.Pz());
  costhetaCS /= TMath::Sqrt( pow(dil.M(),2) + pow(dil.Pt(),2) );
  return costhetaCS;
}

Float_t computeAcoplanarity(TLorentzVector &lm,TLorentzVector &lp){
  return 1.0-fabs(lm.DeltaPhi(lp)/TMath::Pi());
}

//
ValueCollection_t calcCsi(TLorentzVector &a, Float_t &aUnc, TLorentzVector &b, Float_t &bUnc, Float_t bEtaUnc, Float_t sqrts){
  
  Float_t csip(a.Pt()*exp(a.Eta())+b.Pt()*exp(b.Eta()));
  Float_t csipUnc( TMath::Sqrt( pow(aUnc*a.Pt()*exp(a.Eta()),2)
                                +pow(bUnc*b.Pt()*exp(b.Eta()),2)
                                +pow(b.Pt()*exp(b.Eta())*bEtaUnc,2)
                                ) );
  Float_t csim(a.Pt()*exp(-a.Eta())+b.Pt()*exp(-b.Eta()));
  Float_t csimUnc( TMath::Sqrt( pow(aUnc*a.Pt()*exp(-a.Eta()),2)
                                +pow(bUnc*b.Pt()*exp(-b.Eta()),2)
                                +pow(b.Pt()*exp(-b.Eta())*bEtaUnc,2) ) );
  ValueCollection_t toRet;
  toRet.push_back( Value_t(csip/sqrts,csipUnc/sqrts) );
  toRet.push_back( Value_t(csim/sqrts,csimUnc/sqrts) );
  return toRet;
}
