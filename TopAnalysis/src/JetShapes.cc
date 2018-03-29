#include "TopLJets2015/TopAnalysis/interface/JetShapes.h"
#include "TopLJets2015/TopAnalysis/interface/EnergyCorrelations.h"

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>

#include "fastjet/tools/Recluster.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"

#include "Rivet/Math/MatrixN.hh"
#include "Rivet/Math/MatrixDiag.hh"
using Rivet::Matrix;
using Rivet::EigenSystem;

#include "TVector2.h"

using namespace fastjet;
using namespace fastjet::contrib;
using namespace std;

typedef math::XYZTLorentzVector LorentzVector;

//Calculate generalized angularities
/*
  Definition from arXiv:1408.3122
  (beta,kappa) =
    (0,0) -> multiplicity
    (0,2) -> pt dispersion
    (1,1) -> broadening/width/girth
    (2,1) -> thrust (m^2/e)
*/
double calcGA(double beta, double kappa, const pat::Jet *jet, bool includeNeutrals, double ptcut) {
  int mult = 0;
  double sumpt = 0.;
  LorentzVector axis(0., 0., 0., 0.);
  for(size_t ipf=0; ipf<jet->numberOfDaughters(); ipf++) {
    const reco::Candidate *pf=jet->daughter(ipf);
    if (not includeNeutrals and pf->charge() == 0) continue;
    if (ptcut > pf->pt()) continue;
    mult += 1;
    sumpt += pf->pt();
    axis  += pf->p4();
  }
  if (mult < 2) return -1.;
  //std::cout << "sumpt=" << sumpt << " mult=" << mult << std::endl;
  
  double ga = 0.;
  for(size_t ipf=0; ipf<jet->numberOfDaughters(); ipf++) {
    const reco::Candidate *pf=jet->daughter(ipf);
    if (not includeNeutrals and pf->charge() == 0) continue;
    if (ptcut > pf->pt()) continue;
    ga += pow(pf->p4().Pt()/sumpt, kappa) * pow(deltaR(axis, pf->p4())/0.4, beta);
  }
  //std::cout << "ga(" << beta <<","<< kappa << ") ga=" << ga << std::endl;
  
  /*
  if (ga > 1. && beta == 1 && kappa == 1 && iptcut == 0 && icharge == 0) {
    std::cout << "ga(1,1) = " << ga << std::endl;
    std::cout << "pt/sum(pt) deltaR" << std::endl;
    double sumpt_test = 0;
    for(size_t ipf=0; ipf<jet->numberOfDaughters(); ipf++) {
      const reco::Candidate *pf=jet->daughter(ipf);
      if (pf->p4().Pt() < (iptcut+1)*0.500) continue;
      if(pf->charge()==0) continue;
      std::cout << pf->p4().Pt()/sumpt << " " << jet.p4().DeltaR(pf->p4()) << std::endl;
      sumpt_test += pf->p4().Pt()/sumpt;
    }
    std::cout << "sumpt_test = " << sumpt_test << "\n" << std::endl;
  }
  */
  
  return ga;
}

double getMult(const pat::Jet *jet, bool includeNeutrals, double ptcut) {
  //std::cout << "getMult()" << std::endl;
  double sumWeight = 0.;
  for(size_t ipf=0; ipf<jet->numberOfDaughters(); ipf++) {
    const reco::Candidate *pf=jet->daughter(ipf);
    if (not includeNeutrals and pf->charge() == 0) continue;
    if (ptcut > pf->pt()) continue;
    sumWeight += 1.0;
  }
  return sumWeight;
}

double getPtD(const pat::Jet *jet, bool includeNeutrals,  double ptcut) {
  //std::cout << "getPtD()" << std::endl;
  int mult = 0;
  double sumpt  = 0.;
  double sumpt2 = 0.;
  for(size_t ipf=0; ipf<jet->numberOfDaughters(); ipf++) {
    const reco::Candidate *pf=jet->daughter(ipf);
    if (not includeNeutrals and pf->charge() == 0) continue;
    if (ptcut > pf->pt()) continue;
    ++mult;
    sumpt  += pf->pt();
    sumpt2 += pow(pf->pt(), 2);
  }
  if (mult < 2) return -1.;
  double ptd = sumpt2/pow(sumpt,2);
  return ptd;
}

double getPtDs(const pat::Jet *jet, bool includeNeutrals, double ptcut) {
  //std::cout << "getPtDs()" << std::endl;
  double mult   = 0.;
  double sumpt  = 0.;
  double sumpt2 = 0.;
  for(size_t ipf=0; ipf<jet->numberOfDaughters(); ipf++) {
    const reco::Candidate *pf=jet->daughter(ipf);
    if (not includeNeutrals and pf->charge() == 0) continue;
    if (ptcut > pf->pt()) continue;    
    sumpt  += pf->pt();
    sumpt2 += pow(pf->pt(), 2);
  }
  if (mult < 2.) return -1.;
  double ptd = sumpt2/pow(sumpt,2);
  return max(0., sqrt((ptd-1./mult) * mult/(mult-1.)));
}

double getWidth(const pat::Jet *jet, bool includeNeutrals, double ptcut) {
  //std::cout << "getWidth()" << std::endl;
  int mult = 0;
  double sumpt   = 0.;
  double sumptdr = 0.;
  LorentzVector axis(0., 0., 0., 0.);
  for(size_t ipf=0; ipf<jet->numberOfDaughters(); ipf++) {
    const reco::Candidate *pf=jet->daughter(ipf);
    if (not includeNeutrals and pf->charge() == 0) continue;
    if (ptcut > pf->pt()) continue;
    ++mult;
    axis += pf->p4();
  }
  if (mult < 2) return -1.;
  for(size_t ipf=0; ipf<jet->numberOfDaughters(); ipf++) {
    const reco::Candidate *pf=jet->daughter(ipf);
    if (not includeNeutrals and pf->charge() == 0) continue;
    if (ptcut > pf->pt()) continue;    
    sumpt   += pf->pt();
    sumptdr += pf->pt() * deltaR(axis, pf->p4());
  }
  double width = sumptdr/sumpt;
  return width;
}

double getEcc(const pat::Jet *jet, bool includeNeutrals, double ptcut) {
  //std::cout << "getEcc()" << std::endl;
  // Get mean axis
  int mult = 0;
  LorentzVector axis(0., 0., 0., 0.);
  for(size_t ipf=0; ipf<jet->numberOfDaughters(); ipf++) {
    const reco::Candidate *pf=jet->daughter(ipf);
    if (not includeNeutrals and pf->charge() == 0) continue;
    if (ptcut > pf->pt()) continue;
    ++mult;
    axis += pf->p4();
  }
  if (mult < 4) return -1.;
  // Covariance matrix
  Matrix<2> M;
  for(size_t ipf=0; ipf<jet->numberOfDaughters(); ipf++) {
    const reco::Candidate *pf=jet->daughter(ipf);
    if (not includeNeutrals and pf->charge() == 0) continue;
    if (ptcut > pf->pt()) continue;
    Matrix<2> MPart;
    MPart.set(0, 0, (pf->eta() - axis.Eta()) * (pf->eta() - axis.Eta()));
    MPart.set(0, 1, (pf->eta() - axis.Eta()) * mapAngleMPiToPi(pf->phi() - axis.Phi()));
    MPart.set(1, 0, mapAngleMPiToPi(pf->phi() - axis.Phi()) * (pf->eta() - axis.Eta()));
    MPart.set(1, 1, mapAngleMPiToPi(pf->phi() - axis.Phi()) * mapAngleMPiToPi(pf->phi() - axis.Phi()));
    M += MPart * pf->energy();
  }
  // Calculate eccentricity from eigenvalues
  const EigenSystem<2> eigen = diagonalize(M);
  double ecc = 1. - eigen.getEigenValues()[1]/eigen.getEigenValues()[0];
  
  return ecc;
}

double getTau(int N, int M, const pat::Jet *jet, bool includeNeutrals, double ptcut) {
  // Recluster constituents with CA
  int mult = 0.;
  vector<PseudoJet> particles;
  for(size_t ipf=0; ipf<jet->numberOfDaughters(); ipf++) {
    const reco::Candidate *pf=jet->daughter(ipf);
    if (not includeNeutrals and pf->charge() == 0) continue;
    if (ptcut > pf->pt()) continue;    
    ++mult;
    particles.push_back( PseudoJet(pf->px(), pf->py(), pf->pz(), pf->energy()) );
  }
  if (mult < N+1) return -1.;
  
  JetDefinition jet_def(fastjet::cambridge_algorithm, fastjet::JetDefinition::max_allowable_R);
  
  ClusterSequence cs(particles, jet_def);
  vector<PseudoJet> jets = sorted_by_pt(cs.exclusive_jets(1));
  
  PseudoJet cajet = jets[0];
  
  NsubjettinessRatio tau_ratio(N, M, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(1.0));
  
  return tau_ratio(cajet);
}

double getC(int N, double beta, const pat::Jet *jet, bool includeNeutrals, double ptcut) {
  // Recluster constituents with CA
  int mult = 0.;
  vector<PseudoJet> particles;
  for(size_t ipf=0; ipf<jet->numberOfDaughters(); ipf++) {
    const reco::Candidate *pf=jet->daughter(ipf);
    if (not includeNeutrals and pf->charge() == 0) continue;
    if (ptcut > pf->pt()) continue;   
    ++mult;
    particles.push_back( PseudoJet(pf->px(), pf->py(), pf->pz(), pf->energy()) );
  }
  if (mult < N+1) return -1.;
  
  JetDefinition jet_def(fastjet::cambridge_algorithm, fastjet::JetDefinition::max_allowable_R);
  
  ClusterSequence cs(particles, jet_def);
  vector<PseudoJet> jets = sorted_by_pt(cs.exclusive_jets(1));
  
  PseudoJet cajet = jets[0];
  
  EnergyCorrelatorDoubleRatio C(N, beta);
  
  return C(cajet);
}

std::vector<double> getZg(const pat::Jet *jet, bool includeNeutrals, double ptcut) {
  //std::cout << "getZg()" << std::endl;
  // Recluster constituents with CA
  int mult = 0.;
  vector<PseudoJet> particles;
  for(size_t ipf=0; ipf<jet->numberOfDaughters(); ipf++) {
    const reco::Candidate *pf=jet->daughter(ipf);
    if (not includeNeutrals and pf->charge() == 0) continue;
    if (ptcut > pf->pt()) continue;
    ++mult;
    particles.push_back( PseudoJet(pf->px(), pf->py(), pf->pz(), pf->energy()) );
  }
  if (mult < 2) return {-1., -1., -1.};
  
  JetDefinition jet_def(fastjet::cambridge_algorithm, fastjet::JetDefinition::max_allowable_R);
  
  ClusterSequence cs(particles, jet_def);
  vector<PseudoJet> jets = sorted_by_pt(cs.exclusive_jets(1));
  
  PseudoJet cajet = jets[0];
  PseudoJet cajetp1, cajetp2;
  double zg = 0.;
  while (zg < 0.1 and cajet.has_parents(cajetp1, cajetp2)) {
    zg    = cajetp2.pt()/cajet.pt();
    cajet = cajetp1;
  }
  if (zg < 0.1) return {-1., -1., -1.};
  //std::cout << "zg = " << zg << std::endl;
  std::vector<double> results;
  results.push_back(zg);
  results.push_back(zg*cajetp1.delta_R(cajetp2));
  results.push_back(cajetp1.delta_R(cajetp2));
  return results;
}

double mapAngleMPiToPi(double phi) { return TVector2::Phi_mpi_pi(phi); }

double getNSD(double zcut, double beta, const pat::Jet *jet, bool includeNeutrals, double ptcut) {
  // Recluster constituents with CA
  int mult = 0.;
  vector<PseudoJet> particles;
  for(size_t ipf=0; ipf<jet->numberOfDaughters(); ipf++) {
    const reco::Candidate *pf=jet->daughter(ipf);
    if (not includeNeutrals and pf->charge() == 0) continue;
    if (ptcut > pf->pt()) continue;
    ++mult;
    particles.push_back( PseudoJet(pf->px(), pf->py(), pf->pz(), pf->energy()) );
  }
  if (mult < 1) return 0;
  
  JetDefinition jet_def(fastjet::cambridge_algorithm, fastjet::JetDefinition::max_allowable_R);
  
  ClusterSequence cs(particles, jet_def);
  vector<PseudoJet> jets = sorted_by_pt(cs.exclusive_jets(1));
  
  PseudoJet cajet = jets[0];
  PseudoJet cajetp1, cajetp2;
  double nsd = 0.;
  double zg = 0.;
  while (cajet.has_parents(cajetp1, cajetp2)) {
    zg    = cajetp2.pt()/cajet.pt();
    if (zg > zcut * pow(cajetp1.delta_R(cajetp2)/0.4, beta))
      nsd += 1.;
    cajet = cajetp1;
  }
  
  return nsd;
}

std::map<TString,double> getECF(const pat::Jet *jet, bool includeNeutrals, double ptcut) {
  int mult = 0.;
  vector<PseudoJet> particles;
  for(size_t ipf=0; ipf<jet->numberOfDaughters(); ipf++) {
    const reco::Candidate *pf=jet->daughter(ipf);
    if (not includeNeutrals and pf->charge() == 0) continue;
    if (ptcut > pf->pt()) continue;
    ++mult;
    particles.push_back( PseudoJet(pf->px(), pf->py(), pf->pz(), pf->energy()) );
  }
  
  std::map<TString,double> results;
  
  EnergyCorrelations* fECF = new EnergyCorrelations();
  
  fECF->calcECFN(1., particles, true);
  std::map<TString,double> ecfns = fECF->manager->ecfns;
  
  if (mult >= 3) results["m2_b1"] = ecfns["3_1"]/ecfns["2_1"];
  else results["m2_b1"] = -1.;
  
  if (mult >= 3) results["n2_b1"] = ecfns["3_2"]/pow(ecfns["2_1"], 2);
  else results["n2_b1"] = -1.;
  
  if (mult >= 4) results["n3_b1"] = ecfns["4_2"]/pow(ecfns["3_1"], 2);
  else results["n3_b1"] = -1.;
  
  fECF->calcECFN(2., particles);
  std::map<TString,double> ecfns2 = fECF->manager->ecfns;
  
  if (mult >= 3) results["m2_b2"] = ecfns2["3_1"]/ecfns2["2_1"];
  else results["m2_b2"] = -1.;
  
  if (mult >= 3) results["n2_b2"] = ecfns2["3_2"]/pow(ecfns2["2_1"], 2);
  else results["n2_b2"] = -1.;
  
  if (mult >= 4) results["n3_b2"] = ecfns2["4_2"]/pow(ecfns2["3_1"], 2);
  else results["n3_b2"] = -1.;
  
  delete fECF;
  
  return results;
}

double getPFFraction(std::vector<int> pids, const pat::Jet *jet) {
  double sumpt_pid = 0.;
  double sumpt_all = 0.;
  for(size_t ipf=0; ipf<jet->numberOfDaughters(); ipf++) {
    const reco::Candidate *pf=jet->daughter(ipf);
    sumpt_all += pf->pt();
    if (std::find(pids.begin(), pids.end(), abs(pf->pdgId())) != pids.end()) {
      sumpt_pid += pf->pt();
    };
  }
  return sumpt_pid/sumpt_all;
}

