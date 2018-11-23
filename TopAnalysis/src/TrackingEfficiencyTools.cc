#include "TopLJets2015/TopAnalysis/interface/TrackingEfficiencyTools.h"
#include "TH2F.h"
#include "TRandom3.h"


std::map<TString, std::map<TString, std::vector<double> > > getTrackingEfficiencyMap(TString era) {
  std::map<TString, std::map<TString, std::vector<double> > > trackEffMap;
  
  if(era.Contains("era2016")) {
    trackEffMap["BCDEF"]["binning"] = {-2.4, -1.5, -0.8, 0.8, 1.5, 2.4};
    trackEffMap["BCDEF"]["nominal"] = {0.93, 1.08, 1.01, 1.08, 0.93};
    trackEffMap["BCDEF"]["unc"]     = {0.04, 0.04, 0.03, 0.04, 0.04};
    for (unsigned int i = 0; i < trackEffMap["BCDEF"]["nominal"].size(); i++) {
      trackEffMap["BCDEF"]["up"].push_back(trackEffMap["BCDEF"]["nominal"][i]+trackEffMap["BCDEF"]["unc"][i]);
      trackEffMap["BCDEF"]["down"].push_back(trackEffMap["BCDEF"]["nominal"][i]-trackEffMap["BCDEF"]["unc"][i]);
    }
    trackEffMap["GH"]["binning"] = {-2.4, -1.5, -0.8, 0.8, 1.5, 2.4};
    trackEffMap["GH"]["nominal"] = {1.12, 1.07, 1.04, 1.07, 1.12};
    trackEffMap["GH"]["unc"]     = {0.05, 0.06, 0.03, 0.06, 0.05};
    for (unsigned int i = 0; i < trackEffMap["GH"]["nominal"].size(); i++) {
      trackEffMap["GH"]["up"].push_back(trackEffMap["GH"]["nominal"][i]+trackEffMap["GH"]["unc"][i]);
      trackEffMap["GH"]["down"].push_back(trackEffMap["GH"]["nominal"][i]-trackEffMap["GH"]["unc"][i]);
    }
  }
  
  return trackEffMap;
}

void applyEtaDepTrackingEfficiencySF(MiniEvent_t &ev, std::vector<double> sfs, std::vector<double> etas) {
  if (sfs.size() != (etas.size() - 1)) std::cout << "applyEtaDepTrackingEfficiencySF error: need one more bin boundary than scale factors: " << sfs.size() << "," << etas.size() << std::endl;
  for (unsigned int i = 0; i < sfs.size(); i++) {
    applyTrackingEfficiencySF(ev, sfs[i], etas[i], etas[i+1]);
  }
}

void applyTrackingEfficiencySF(MiniEvent_t &ev, double sf, double minEta, double maxEta) {
  if(ev.isData) return;
  
  TRandom* random = new TRandom3(0); // random seed

  if (sf <= 1) {
    for (int k = 0; k < ev.npf; k++) {
      if (abs(ev.pf_id[k]) != 211) continue;
      if (ev.pf_eta[k] < minEta) continue;
      if (ev.pf_eta[k] > maxEta) continue;
      if (random->Rndm() > sf) {
        //make sure that particle does not pass any cuts
        ev.pf_pt[k]  = 1e-20;
        ev.pf_m[k]   = 1e-20;
        ev.pf_eta[k] = 999.;
        ev.pf_c[k]   = 0;
      }
    }
  }
  else { // sf > 1
    // find charged hadrons that were not reconstructed
    double dRcut = 0.01;
    std::vector<int> chGenNonRecoHadrons;
    int NchGenHadrons = 0;
    for (int g = 0; g < ev.ngpf; g++) {
      if (ev.gpf_pt[g] < 0.9) continue;
      if (ev.gpf_eta[g] < minEta) continue;
      if (ev.gpf_eta[g] > maxEta) continue;
      if (ev.gpf_c[g] == 0) continue;
      if (abs(ev.gpf_id[g]) < 100) continue;
      NchGenHadrons++;
      bool matched = false;
      for (int k = 0; k < ev.npf; k++) {
        if (ev.pf_pt[k] < 0.8) continue;
        if (abs(ev.pf_id[k]) != 211) continue;
        double dEta = ev.gpf_eta[g] - ev.pf_eta[k];
        double dPhi = TVector2::Phi_mpi_pi(ev.gpf_phi[g] - ev.pf_phi[k]);
        double dR = sqrt(pow(dEta, 2) + pow(dPhi, 2));
        if (dR < dRcut) {
          matched = true;
          break;
        }
      }
      if (!matched) chGenNonRecoHadrons.push_back(g);
    }
    if (chGenNonRecoHadrons.size() == 0) return;
    double promotionProb = TMath::Min(1., NchGenHadrons*(sf-1.)/chGenNonRecoHadrons.size());
    std::vector<int> chGenNonRecoHadronsToPromote;
    for (const int g : chGenNonRecoHadrons) {
      if (random->Rndm() < promotionProb) {
        chGenNonRecoHadronsToPromote.push_back(g);
      }
    }
    for (unsigned int i = 0; i < chGenNonRecoHadronsToPromote.size(); i++) {
      int k = ev.npf + i;
      int g = chGenNonRecoHadronsToPromote[i];
      // jet association
      int j = -1;
      double jetR = 0.4;
      for (int ij = 0; ij < ev.nj; ij++) {
        double dEta = ev.gpf_eta[g] - ev.j_eta[ij];
        double dPhi = TVector2::Phi_mpi_pi(ev.gpf_phi[g] - ev.j_phi[ij]);
        double dR = sqrt(pow(dEta, 2) + pow(dPhi, 2));
        if (dR < jetR) {
          j = ij;
          break;
        }
      }
      ev.pf_j[k]   = j;
      ev.pf_id[k]  = ev.gpf_id[g];
      ev.pf_c[k]   = ev.gpf_c[g];
      ev.pf_pt[k]  = ev.gpf_pt[g];
      ev.pf_eta[k] = ev.gpf_eta[g];
      ev.pf_phi[k] = ev.gpf_phi[g];
      ev.pf_m[k]   = ev.gpf_m[g];
      ev.pf_dxy[k] = 0.;
      ev.pf_dz[k]  = 0.;
    }
    ev.npf    = ev.npf + chGenNonRecoHadronsToPromote.size();
  }
  
  delete random;
}

