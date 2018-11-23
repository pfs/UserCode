#ifndef __FASTJET_CONTRIB_ENERGYCORRELATIONS_HH__
#define __FASTJET_CONTRIB_ENERGYCORRELATIONS_HH__
#include "fastjet/PseudoJet.hh"
#include <vector>
#include <map>
#include "TMath.h"
#include "TString.h"


class ECFNManager {
public:
  // just a bunch of floats and bools to hold different values of normalized ECFs
  ECFNManager() {
    flags["3_1"]=true; 
    flags["3_2"]=true; 
    flags["3_3"]=true; 
    flags["4_1"]=true; 
    flags["4_2"]=true; 
    flags["4_3"]=false; 
  }
  ~ECFNManager() {}

  std::map<TString,double> ecfns; // maps "N_I" to ECFN
  std::map<TString,bool>   flags; // maps "N_I" to flag

  bool doN1=true, doN2=true, doN3=true, doN4=true;
  inline void clear() { for (std::map<TString,double>::iterator it=ecfns.begin(); it!=ecfns.end(); ++it) it->second = -999;}
};


class EnergyCorrelations { 
 public:
  EnergyCorrelations();
  ~EnergyCorrelations(){}
  double DeltaR2(fastjet::PseudoJet j1, fastjet::PseudoJet j2);
  double DeltaR2(double iEta1,double iPhi1,double iEta2,double iPhi2);
  void calcECF(double beta, std::vector<fastjet::PseudoJet> &constituents, double *n1=0, double *n2=0, double *n3=0, double *n4=0);
  void calcECFN(double beta, std::vector<fastjet::PseudoJet> &constituents, bool iClear=false, bool useMin=true);
  ECFNManager *manager;
};

#endif
