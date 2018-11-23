#ifndef _jetshapes_h_
#define _jetshapes_h_

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

double calcGA(double beta, double kappa, const pat::Jet *jet, bool includeNeutrals, double ptcut);
double getMult(const pat::Jet *jet, bool includeNeutrals, double ptcut);
double getPtD(const pat::Jet *jet, bool includeNeutrals, double ptcut);
double getPtDs(const pat::Jet *jet, bool includeNeutrals, double ptcut);
double getWidth(const pat::Jet *jet, bool includeNeutrals, double ptcut);
double getEcc(const pat::Jet *jet, bool includeNeutrals, double ptcut);
double getTau(int N, int M, const pat::Jet *jet, bool includeNeutrals, double ptcut);
double getC(int N, double beta, const pat::Jet *jet, bool includeNeutrals, double ptcut);
std::vector<double> getZg(const pat::Jet *jet, bool includeNeutrals, double ptcut);
double mapAngleMPiToPi(double phi);
double getNSD(double zcut, double beta, const pat::Jet *jet, bool includeNeutrals, double ptcut);
std::map<TString,double> getECF(const pat::Jet *jet, bool includeNeutrals, double ptcut);
double getPFFraction(std::vector<int> pids, const pat::Jet *jet);
double mapAngleMPiToPi(double phi);

#endif
