//compile with comp.sh
//./analyse_masses 0 lx
#include "ERRORS.h"

#include "TFile.h"
#include "TH1F.h"
#include "TApplication.h"
#include "Definitions.hh"
#include "interface/Definitions_cmssw.hh"
#include "TMath.h"
#include "TF1.h"
#include <map>
#include "TCanvas.h"
#include <string>
#include "TMarker.h"
#include "TROOT.h"
#include "TFitResult.h"
using namespace std;
using namespace Definitions;

static const unsigned short N_MC_histo_names = 1;
static const char * MC_histo_names[N_MC_histo_names] = {"t#bar{t}"};
TH1F * sum_MC(const char *);
enum model_t {SM, CFLIP};
unsigned char model(SM);
TFile * plotter_file(nullptr);
const char * model_tag[] = {"_nominal", "_cflip"};
const unsigned char N_interesting_jets = 4;
const unsigned char interesting_jets[N_interesting_jets] = {HAD_W, HAD_T, LEPT_W, LEPT_T};
double CalculateMass(TH1 *, void *params);
int main()
{
  //gROOT -> SetBatch(kTRUE);
  plotter_file = TFile::Open("$EOS/analysis_MC13TeV_TTJets/plots/plotter.root");
  //TApplication app("myapp", 0, 0);
  unsigned short ch_ind = L; 
  unsigned short level_ind = RECO; 
  unsigned short jet_ind = 0;
  const unsigned char interesting_jet_ind = interesting_jets[jet_ind];

  const TString dir = TString(tag_channels_types_[ch_ind]) + "_jet_mass_" + tag_levels_types_[level_ind] + "_" + tag_jet_types_[interesting_jet_ind];
  TH1F * h_chi =  sum_MC(dir) ;
  TH1F * h_chi_syst =  (TH1F *) plotter_file -> Get(dir + "/totalmcunc");
  double mass(calculate_result(h_chi, h_chi_syst, CalculateMass, nullptr));
  plotter_file -> Close();

}

TH1F * sum_MC(const char * dir)
{
  TH1F * h = (TH1F*) (plotter_file -> GetDirectory(dir) -> Get(TString(dir) + "_" +  MC_histo_names[0]) -> Clone(TString(dir) + "_MCsum"));
  return h;
	       
}

double CalculateMass(TH1 * h, void *params)
{
 
  return 0.0;
}



