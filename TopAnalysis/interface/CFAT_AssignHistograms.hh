#ifndef _CFAT_AssignHistograms_hh_
#define _CFAT_AssignHistograms_hh_
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include <map>
#include <utility>
using namespace std;
typedef map<TString, pair<TH1*, TH2*>* > map_hist_1D;
void AssignHistograms(map_hist_1D &);
void AssignSpecificHistograms2D(map<TString, TH2*> &);
void CreateHistogram1D(map_hist_1D &, const char * title, const char * axes, unsigned int nbins, double min , double max);
void CreateHistogram1Darray(map_hist_1D &, const char * title, const char * axes, unsigned int nbins, double *bins);
void CreateHistogram2D(map<TString, TH2*> &, const char * title, const char * axes, unsigned int nbins_x, double min_x, double max_x, unsigned int nbins_y, double min_y, double max_y);
//void CreateHistogram2D(map<TString, TH2*> &, const char * title, const char * axes, unsigned int nbins_x, double * bins, unsigned int nbins_y, double min_y, double max_y);
#endif
