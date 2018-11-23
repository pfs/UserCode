#include "../interface/EnergyCorrelations.h"

EnergyCorrelations::EnergyCorrelations() { 
  manager = new ECFNManager();
}
double EnergyCorrelations::DeltaR2(fastjet::PseudoJet j1, fastjet::PseudoJet j2) {
  return DeltaR2(j1.eta(),j1.phi(),j2.eta(),j2.phi());
}
double EnergyCorrelations::DeltaR2(double iEta1,double iPhi1,double iEta2,double iPhi2) { 
  double pDPhi = fabs(iPhi1-iPhi2);
  if(2.*TMath::Pi()-pDPhi < pDPhi) pDPhi = 2.*TMath::Pi()-pDPhi;
  return pDPhi*pDPhi+fabs(iEta1-iEta2)*fabs(iEta1-iEta2);
}

void EnergyCorrelations::calcECF(double beta, std::vector<fastjet::PseudoJet> &constituents, double *n1, double *n2, double *n3, double *n4) {
  unsigned int nC = constituents.size();
  double halfBeta = beta/2.;

  // if only N=1,2, do not bother caching kinematics
  if (!n3 && !n4) {
    if (n1) { // N=1
      double val=0;
      for (unsigned int iC=0; iC!=nC; ++iC) {    
        val += constituents[iC].perp();
      }
      *n1 = val;
    }
    if (n2) { // N=2
      double val=0;
      for (unsigned int iC=0; iC!=nC; ++iC) {
        fastjet::PseudoJet iconst = constituents[iC];
        for (unsigned int jC=0; jC!=iC; ++jC) {
          fastjet::PseudoJet jconst = constituents[jC];
          val += iconst.perp() * jconst.perp() * pow(DeltaR2(iconst,jconst),halfBeta);
        }
      }
      *n2 = val;
    }
    return;
  }

  // cache kinematics
  double *pTs = new double[nC];
  double **dRs = new double*[nC];
  for (unsigned int iC=0; iC!=nC; ++iC) {
    dRs[iC] = new double[iC];
  }
  for (unsigned int iC=0; iC!=nC; ++iC) {
    fastjet::PseudoJet iconst = constituents[iC];
    pTs[iC] = iconst.perp();
    for (unsigned int jC=0; jC!=iC; ++jC) {
      fastjet::PseudoJet jconst = constituents[jC];
      dRs[iC][jC] = pow(DeltaR2(iconst,jconst),halfBeta);
    }
  }
  
  // now we calculate the real ECFs
  if (n1) { // N=1
    double val=0;
    for (unsigned int iC=0; iC!=nC; ++iC) {    
      val += pTs[iC];
    } // iC
    *n1 = val;
  }
  if (n2) { // N=2
    double val=0;
    for (unsigned int iC=0; iC!=nC; ++iC) {
      for (unsigned int jC=0; jC!=iC; ++jC) {
        val += pTs[iC] * pTs[jC] * dRs[iC][jC]; 
      } // jC
    } // iC
    *n2 = val;
  }
  if (n3) {
    double val=0;
    for (unsigned int iC=0; iC!=nC; ++iC) {
      for (unsigned int jC=0; jC!=iC; ++jC) {
        double val_ij = pTs[iC]*pTs[jC]*dRs[iC][jC];
        for (unsigned int kC=0; kC!=jC; ++kC) {
          val += val_ij * pTs[kC] * dRs[iC][kC] * dRs[jC][kC];
        } // kC
      } // jC
    } // iC
    *n3 = val;
  }
  if (n4) {
    double val=0;
    for (unsigned int iC=0; iC!=nC; ++iC) {
      for (unsigned int jC=0; jC!=iC; ++jC) {
        double val_ij = pTs[iC]*pTs[jC]*dRs[iC][jC];
        for (unsigned int kC=0; kC!=jC; ++kC) {
          double val_ijk = val_ij * pTs[kC] * dRs[iC][kC] * dRs[jC][kC];
          for (unsigned int lC=0; lC!=kC; ++lC) {
            val += val_ijk * pTs[lC] * dRs[iC][lC] * dRs[jC][lC] * dRs[kC][lC];
          } // lC
        } // kC
      } // jC
    } // iC
    *n4 = val;
  }

  // cleanup
  delete[] pTs;
  for (unsigned int iC=0; iC!=nC; ++iC) {
    delete[] dRs[iC];
  }
  delete[] dRs;
}
void EnergyCorrelations::calcECFN(double beta, std::vector<fastjet::PseudoJet> &constituents, bool iClear, bool useMin) {
  unsigned int nC = constituents.size();
  double halfBeta = beta/2.;
  if(iClear) manager->clear();
  // get the normalization factor
  double baseNorm=0; 
  calcECF(beta,constituents,&baseNorm,0,0,0);

  // cache kinematics
  double *pTs = new double[nC];
  double **dRs = new double*[nC];
  for (unsigned int iC=0; iC!=nC; ++iC) {
    dRs[iC] = new double[iC];
  }
  for (unsigned int iC=0; iC!=nC; ++iC) {
    fastjet::PseudoJet iconst = constituents[iC];
    pTs[iC] = iconst.perp();
    for (unsigned int jC=0; jC!=iC; ++jC) {
      fastjet::PseudoJet jconst = constituents[jC];
      dRs[iC][jC] = pow(DeltaR2(iconst,jconst),halfBeta);
    }
  }
  
  // now we calculate the ECFNs
  if (manager->doN1) { // N=1
    manager->ecfns["1_1"] = 1;
    manager->ecfns["1_2"] = 1;
    manager->ecfns["1_3"] = 1;
  }
  if (manager->doN2) { // N=2
    double norm = pow(baseNorm,2);
    double val=0;
    for (unsigned int iC=0; iC!=nC; ++iC) {
      for (unsigned int jC=0; jC!=iC; ++jC) {
        val += pTs[iC] * pTs[jC] * dRs[iC][jC] / norm; 
      } // jC
    } // iC
    manager->ecfns["2_1"] = val;
    manager->ecfns["2_2"] = val;
    manager->ecfns["2_3"] = val;
  }

  bool doI1=manager->flags["3_1"];
  bool doI2=manager->flags["3_2"];
  bool doI3=manager->flags["3_3"];
  if (manager->doN3 && (doI1||doI2||doI3)) {
    double norm = pow(baseNorm,3);
    double val1=0,val2=0,val3=0;
    unsigned int nAngles=3;
    std::vector<double> angles(nAngles);

    for (unsigned int iC=0; iC!=nC; ++iC) {
      for (unsigned int jC=0; jC!=iC; ++jC) {
        double val_ij = pTs[iC]*pTs[jC];
        angles[0] = dRs[iC][jC];

        for (unsigned int kC=0; kC!=jC; ++kC) {
          angles[1] = dRs[iC][kC];
          angles[2] = dRs[jC][kC];

          if (doI1||doI2||doI3) {
            double angle_1=999; unsigned int index_1=999;
            for (unsigned int iA=0; iA!=nAngles; ++iA) {
              if ((useMin && angles[iA]<angle_1)  || 
                  (!useMin && angles[iA]>angle_1)) {
                angle_1 = angles[iA];
                index_1 = iA;
              }
            }
            if (doI2||doI3) {
              double angle_2=999; 
              for (unsigned int jA=0; jA!=nAngles; ++jA) {
                if (jA==index_1) continue;
                if ((useMin && angles[jA]<angle_2)  || 
                    (!useMin && angles[jA]>angle_2)) {
                  angle_2 = angles[jA];
                }
              }
              if (doI3) {
                val3 += val_ij * pTs[kC] * angles[0] * angles[1] * angles[2] / norm;
              }
              if (doI2)
                val2 += val_ij * pTs[kC] * angle_1 * angle_2 / norm;
            }
            if (doI1)
              val1 += val_ij * pTs[kC] * angle_1 / norm;
          }

        } // kC
      } // jC
    } // iC
    manager->ecfns["3_1"] = val1;
    manager->ecfns["3_2"] = val2;
    manager->ecfns["3_3"] = val3;
  }

  doI1=manager->flags["4_1"];
  doI2=manager->flags["4_2"];
  if (manager->doN4 && (doI1||doI2)) {
    double norm = pow(baseNorm,4);
    double val1=0,val2=0;
    unsigned int nAngles=6;
    std::vector<double> angles(nAngles);

    for (unsigned int iC=0; iC!=nC; ++iC) {
      for (unsigned int jC=0; jC!=iC; ++jC) {
        double val_ij = pTs[iC]*pTs[jC];
        angles[0] = dRs[iC][jC];

        for (unsigned int kC=0; kC!=jC; ++kC) {
          double val_ijk = val_ij * pTs[kC];
          angles[1] = dRs[iC][kC];
          angles[2] = dRs[jC][kC];

          for (unsigned int lC=0; lC!=kC; ++lC) {
            angles[3] = dRs[iC][lC];
            angles[4] = dRs[jC][lC];
            angles[5] = dRs[kC][lC];

            if (doI1||doI2) {
              double angle_1=999; unsigned int index_1=999;
              for (unsigned int iA=0; iA!=nAngles; ++iA) {
                if ((useMin && angles[iA]<angle_1)  || 
                    (!useMin && angles[iA]>angle_1)) {
                  angle_1 = angles[iA];
                  index_1 = iA;
                }
              }
              if (doI2) {
                double angle_2=999; 
                for (unsigned int jA=0; jA!=nAngles; ++jA) {
                  if (jA==index_1) continue;
                  if ((useMin && angles[jA]<angle_2)  || 
                      (!useMin && angles[jA]>angle_2)) {
                    angle_2 = angles[jA];
                  }
                }
                val2 += val_ijk * pTs[lC] * angle_1 * angle_2 / norm;
              }
              if (doI1)
                val1 += val_ijk * pTs[lC] * angle_1 / norm;
            }
          } // lC
        } // kC
      } // jC
    } // iC
    manager->ecfns["4_1"] = val1;
    manager->ecfns["4_2"] = val2;
    manager->ecfns["4_3"] = 0;
  }
  // cleanup
  delete[] pTs;
  for (unsigned int iC=0; iC!=nC; ++iC) {
    delete[] dRs[iC];
  }
  delete[] dRs;
}

