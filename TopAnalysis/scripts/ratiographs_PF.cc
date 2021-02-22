// ./ratiographs_PF merged SM
#include "TFile.h"
#include "TH1F.h"
#include "TApplication.h"
#include "Definitions.hh"
#include "interface/Definitions_cmssw.hh"
#include "TMath.h"
#include <map>
#include "TCanvas.h"
#include "TGraph.h"
#include "TROOT.h"
#include "TLegend.h"
#include <assert.h>
#include "TStyle.h"
#include "THStack.h"
static const unsigned short N_MC_histo_names = 7;
static const char * MC_histo_names[N_MC_histo_names] = {"t#bar{t}", "Single top", "W", "DY", "Multiboson", "t#bar{t}+V", "QCD"};
TH1F * sum_MC(const char *);
enum model_t {SM, CFLIP};
enum flow_t {N, E, PT};
const char* tag_flow[] = {"N", "E", "Pt"};
const char * title_flow[] = {"particle", "energy", "Pt"};
//unsigned char model = SM;
Color_t color[2] = {kRed, kBlue};

const char * model_tag[] = {"SM", "cflip"};
const char * model_title[] = {"", " for the colour octet $W$ model"};
const unsigned char npairs = 4;
const char * jetpairs[npairs] = {"blb2l", "qfhb", "hbqc", "qlq2l"};
const char * jetpairstitle[npairs] = {"j_{1}^{b},j_{2}^{b}", 
				      "j_{f}^{W},j_{h}^{b}", 
				      "j_{h}^{b},j_{c}^{W}",
				      "j_{1}^{W},j_{2}^{W}", 
};
const unsigned char nsources = 2;
const char * tag_source_type[nsources] = {"MC", "data"};
double R(TH1F *);
using namespace std;
using namespace Definitions;
TFile * plotter_file[2] = {
  TFile::Open("$EOS/analysis_MC13TeV_TTJets/plots/plotter.root"), 
  TFile::Open("$EOS/analysis_MC13TeV_TTJets_cflip/plots/plotter.root")};

int main(int argc, char * argv[])
{
  assert(argc == 3);

  const string method(argv[1]); 
  const string rel(argv[2]);
  const string dir(string("ratiographs_") + method + "_" + rel); 
  system(TString("mkdir -p ") + dir);
  system(TString("rm ") + dir + "/*");
  gStyle -> SetOptStat(0);
  gStyle -> SetOptTitle(0);
  gROOT -> SetBatch(kTRUE);
  gStyle -> SetHatchesLineWidth(1);

  unsigned char modelstart = 0;
  unsigned char modelend = 0;
  unsigned char nmodels = 0;
  
  if (method.compare("SM") == 0)
    {
      modelstart = 0;
      modelend = 0;
      nmodels = 1;
    }
  if (method.compare("cflip") == 0)
    {
      modelstart = 1;
      modelend = 1;
      nmodels = 1;
    }
  if (method.compare("merged") == 0)
    {
      modelstart = 0;
      modelend = 1;
      nmodels = 2;
    }

  for (unsigned short ch_ind = 0; ch_ind < N_channels_types_; ch_ind ++)
    {
      for (unsigned short level_ind = 0; level_ind < N_levels_types_; level_ind ++)
	{
	  const unsigned char ntypes = level_ind == RECO ? 2 : 1;
	  for (unsigned short flow_ind = 0; flow_ind <= PT; flow_ind ++)
	    {
	      for (unsigned short charge_ind = 0; charge_ind < N_charge_types_; charge_ind ++)
		{
		  TH1F * h[npairs][ntypes][nmodels];
		  TH1F * hunc[npairs];
		  for (unsigned short pair_ind = 0; pair_ind < npairs; pair_ind ++)
		    {
		      hunc[pair_ind] = nullptr;
		      for (unsigned short type_ind = 0; type_ind < ntypes; type_ind ++)
			{
			  for (unsigned short model_ind = modelstart; model_ind <= modelend; model_ind ++)
			    {
			      h[pair_ind][type_ind][model_ind] = nullptr;
			    }
			}
		    }
		  for (unsigned short type_ind = 0; type_ind < ntypes; type_ind ++)
		    {
		      for (unsigned short pair_ind = 0; pair_ind < npairs; pair_ind ++)
			{
			  const TString dir = TString(tag_channels_types_[ch_ind]) + 
			    "_chi" + jetpairs[pair_ind] + 
			    "_" + tag_flow[flow_ind] + "_" + 
			    tag_charge_types_[charge_ind] + "_" + 
			    tag_levels_types_[level_ind] + "_jetprt";
			  for (unsigned short model_ind = modelstart; model_ind <= modelend; model_ind ++)
			    {
			      if (model_ind == SM)
				{
				  h[pair_ind][type_ind][model_ind] = type_ind == 0 ? sum_MC(dir) : (TH1F*) plotter_file[model_ind] -> GetDirectory(dir) -> Get(dir);
				  if (type_ind == MC)
				    {
				      hunc[pair_ind] =  (TH1F *) plotter_file[model_ind] -> Get(dir + "/" + dir + "_totalMCUncShape");
				      // plotter_file[model_ind] -> Get(dir) -> ls();
				      // printf("hunc[pair_ind] %p\n", hunc[pair_ind]);
				      // getchar();
				    }
				}
			      if (model_ind == CFLIP and type_ind == MC)
				{
				  h[pair_ind][type_ind][model_ind] = (TH1F*) plotter_file[model_ind] -> GetDirectory(dir) -> Get(dir + "_t#bar{t} cflip");
				  if (type_ind == MC and model_ind == SM)
				    hunc[pair_ind] = (TH1F *) plotter_file[model_ind] -> Get(dir + "/" + dir + "_totalMCUncShape") ;

				}
			    }
			}
		    }
		  for (unsigned short type_ind = 0; type_ind < ntypes; type_ind ++)
		    {
		      for (unsigned short model_ind = modelstart; model_ind <= modelend; model_ind ++)
			{
			  if (h[3][type_ind][model_ind])
			    h[3][type_ind][model_ind] -> Scale(1.0/h[3][type_ind][model_ind] -> Integral());
			}
		    }
		  char npairscomp = npairs -1;
		  if (rel.compare("SM") == 0)
		    {
		      npairscomp = npairs;
		    }
		  for (unsigned short pair_ind = 0; pair_ind < npairscomp; pair_ind ++)
		    {
		      const unsigned char nbins = h[3][0][0] -> GetNbinsX();
		      char title[128]; 
		      if (rel.compare("self") == 0)
			{
			  sprintf(title, "#frac{#Delta N(%s)}{#Delta N(%s)}#frac{N(%s)}{N(%s)}", 
				  jetpairstitle[pair_ind], 
				  jetpairstitle[3],
				  jetpairstitle[3],
				  jetpairstitle[pair_ind]); 
			}
		      if (rel.compare("SM") == 0)
			{
			  sprintf(title, "#frac{#Delta N_{cflip}(%s)}{#Delta N_{SM}(%s)}#frac{N_{SM}(%s)}{N_{cflip}(%s)}", 
				  jetpairstitle[pair_ind], 
				  jetpairstitle[3],
				  jetpairstitle[3],
				  jetpairstitle[pair_ind]
				  ); 
			}
		      const TString name = TString(tag_channels_types_[ch_ind]) + "_" +
			jetpairs[pair_ind] + "_" + 
			tag_flow[flow_ind] + "_" + 
			tag_charge_types_[charge_ind] + "_" + 
			tag_levels_types_[level_ind];
		      TCanvas c(name, name);
		      c.SetLeftMargin(0.18);
		      THStack stack("stack", TString("stack; #chi; ") + title);
		      TLegend legend(0.7, 0.6, 0.88, 0.88);
		      legend .SetLineWidth(0);
		      for (unsigned short type_ind = 0; type_ind < ntypes; type_ind ++)
			{
			  for (unsigned short model_ind = modelstart; model_ind <= modelend; model_ind ++)
			    {
			      unsigned char modelcomp;
			      if (rel.compare("self") == 0)
				{
				  modelcomp = model_ind;
				}
			      if (rel.compare("SM") == 0)
			
				{
				  modelcomp = SM;
				  if (model_ind == SM)
				    continue;
				}




			      if (h[pair_ind][type_ind][model_ind])
				{
				  h[pair_ind][type_ind][model_ind] -> Scale(1.0 / h[pair_ind][type_ind][model_ind] -> Integral());
				  h[pair_ind][type_ind][model_ind] -> Divide(h[3][type_ind][modelcomp]);
				  /*		  h[pair_ind][type_ind][model_ind] -> SetMinimum(0.0);
						  h[pair_ind][type_ind][model_ind] -> SetMaximum(1.4 * h[pair_ind][type_ind][model_ind] -> GetMaximum());*/
				}
			      else
				continue;
			      if (type_ind == MC)
				{
				  h[pair_ind][type_ind][model_ind] -> SetLineColor(color[model_ind]);
				  char ytitle[128];
				  sprintf(ytitle,"#frac{#Delta N_{%s}}{#Delta N_{%s}}#frac{N_{%s}}{N_{%s}}", 
					  jetpairstitle[pair_ind], 
					  jetpairstitle[3], 
					  jetpairstitle[3], 
					  jetpairstitle[pair_ind]);
				  stack.Add(h[pair_ind][type_ind][model_ind]);
				  //h[pair_ind][type_ind][model_ind] j_{f}^{W}-> GetYaxis() -> SetTitle(ytitle);
				  //h[pair_ind][type_ind][model_ind] -> GetYaxis() -> SetTitleOffset(1.6);
				  //h[pair_ind][type_ind][model_ind] -> Draw(model_ind == modelstart ? "" : "SAME");
				  legend.AddEntry(h[pair_ind][type_ind][model_ind], TString(jetpairstitle[pair_ind]) + " " + 
						  model_tag[model_ind], "l");
				}
			      else if (model_ind == SM)
				{
				  h[pair_ind][type_ind][model_ind] -> SetMarkerStyle(kFullCircle);
				  stack.Add(h[pair_ind][type_ind][model_ind], "P");
				  //h[pair_ind][type_ind][model_ind] -> Draw("SAME");
				  legend.AddEntry(h[pair_ind][type_ind][model_ind], "data", "p");
				}
			      if (model_ind == SM and type_ind == MC)
				{
				  hunc[pair_ind] -> Scale(1.0/hunc[pair_ind] -> Integral());
				  hunc[pair_ind] -> Divide(h[3][type_ind][model_ind]);
				  hunc[pair_ind] -> SetFillColor(TColor::GetColor("#99d8c9"));
				  hunc[pair_ind] -> SetFillStyle(3254);
				  stack.Add(hunc[pair_ind], "E2");
				  //hunc[pair_ind] -> Draw("E2");
				  /*				  if (name == "L_blb2l_N_allconst_reco")
								  hunc[pair_ind] -> Print("all");*/
				}
			  
			    }      
			}
		      stack.Draw("nostack");
		      stack.SetMinimum(0.0);
		      stack.SetMaximum(1.6 * stack.GetMaximum("nostack"));
		      stack. GetYaxis() -> SetTitleOffset(1.6);
		      legend.Draw("SAME");
		      c.SaveAs((dir + "/" + c.GetName() + ".png").c_str());
		      c.SaveAs((dir + "/" + c.GetName() + ".eps").c_str());
		    }
		}
	    }
	}
	
    }
  plotter_file[SM] -> Close();
  plotter_file[CFLIP] -> Close();
}

TH1F * sum_MC(const char * dir)
{
  //printf("dir %s\n", dir);
  
  TH1F * h = (TH1F*) (plotter_file[SM] -> GetDirectory(dir) -> Get(TString(dir) + "_" +  MC_histo_names[0]) -> Clone(TString(dir) + "_MCsum"));
  //printf("hname %s h int %f\n", h-> GetName(), h -> Integral()); 
  for (unsigned short MC_ind = 1; MC_ind < N_MC_histo_names; MC_ind ++)
    {
      h -> Add((TH1F*)  (plotter_file[SM] -> GetDirectory(dir) -> Get(TString(dir) + "_" + MC_histo_names[MC_ind])));

    }
  return h;
	       
}
