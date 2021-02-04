// ./event_yields_script lx
#include "TFile.h"
#include "TApplication.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TMath.h"
#include <map>
#include <assert.h>
//#include <string>
using namespace std;
const unsigned short N_bins = 4;
const unsigned short N_MC_histo_names = 7;
const char * MC_syst_labels[N_MC_histo_names + 1] = {"\\ttbar", "Single top", "\\PW", "DY", "Multiboson", "\\ttbar+V", "QCD", "statistics"};
static unsigned char const subMC_N[N_MC_histo_names + 1] = {1, 4, 3, 2, 5, 4, 17, 1};

TString insertnumber(float, const char *);

void insert_header(FILE *, unsigned char level_ind);
void insert_channel_header(FILE *, const char *);
void insert_MC(FILE *, double (*)[N_bins], map<TString, float *>);
void insert_MC_sum(FILE *, double *, TH1 * h_MC_Unc, TH1 * hstat);
void insert_MC_unc(FILE *, TH1 *);
void insert_data(FILE *, TH1 *);
void insert_footer(FILE *);
const char * pm = nullptr;
const char * pmlx = "$\\pm$";
const char * pmnonlx = "+-";
TH1 * hstat    = nullptr;
const unsigned short N_levels                 = 2;
const char * level_names[N_levels]            = {"reco", "gen"};
string environment;
int main(int argc, char * argv[])
{
  assert(argc == 2);
  environment = string(argv[1]);
  if (environment.compare("lx") == 0)
    pm = pmlx;
  else
    pm = pmnonlx;

  const char * subMC_ttbar[subMC_N[0]] = {
    "MC13TeV_TTJets"
  };

  const char * subMC_SingleT[subMC_N[1]] = {
    "MC13TeV_SingleTbar_tW",
    "MC13TeV_SingleTbar_t",
    "MC13TeV_SingleT_tW",
    "MC13TeV_SingleT_t"
  };

  const char * subMC_W[subMC_N[2]] = {
    "MC13TeV_W0Jets",
    "MC13TeV_W1Jets",
    "MC13TeV_W2Jets"
  };

  const char * subMC_DY[subMC_N[3]] = {
    "MC13TeV_DY50toInf_mlm",
    "MC13TeV_DY10to50"
  };

  const char * subMC_Multiboson[subMC_N[4]] = {
    "MC13TeV_ZZTo2L2Nu",
    "MC13TeV_ZZTo2L2Q",
    "MC13TeV_WWToLNuQQ",
    "MC13TeV_WWTo2L2Nu",
    "MC13TeV_WZTo3LNu"
  };

  const char * subMC_ttbarpV[subMC_N[5]] = {
    "MC13TeV_TTWToLNu",
    "MC13TeV_TTWToQQ",
    "MC13TeV_TTZToQQ",
    "MC13TeV_TTZToLLNuNu"
  };

  const char * subMC_QCD[subMC_N[6]] = {
    "MC13TeV_QCDEMEnriched120to170",
    "MC13TeV_QCDEMEnriched170to300",
    "MC13TeV_QCDEMEnriched300toInf",
    "MC13TeV_QCDEMEnriched30to50",
    "MC13TeV_QCDEMEnriched50to80",
    "MC13TeV_QCDEMEnriched80to120",
    "MC13TeV_QCDMuEnriched1000toInf",
    "MC13TeV_QCDMuEnriched120to170",
    "MC13TeV_QCDMuEnriched170to300",
    "MC13TeV_QCDMuEnriched300to470",
    "MC13TeV_QCDMuEnriched30to50",
    "MC13TeV_QCDMuEnriched470to600",
    "MC13TeV_QCDMuEnriched50to80",
    "MC13TeV_QCDMuEnriched600to700",
    "MC13TeV_QCDMuEnriched600to800",
    "MC13TeV_QCDMuEnriched800to1000",
    "MC13TeV_QCDMuEnriched80to120"
  };

  const char * subMC_stat[subMC_N[7]] = {
    "stat"
  };

  const char ** subMC[N_MC_histo_names + 1] = {
    subMC_ttbar,
    subMC_SingleT,
    subMC_W,
    subMC_DY,
    subMC_Multiboson,
    subMC_ttbarpV,
    subMC_QCD,
    subMC_stat
  };

  

  //  
  printf("pm %s\n", pm);
  //getchar();
  TFile * plotter                               = TFile::Open("$EOS/analysis_MC13TeV_TTJets/plots/plotter.root");
  // TFile * plotter_cflip                         = TFile::Open("$EOS/analysis_MC13TeV_TTJets_cflip/plots/plotter.root");
  const unsigned short N_ch                     = 3;
  const char * ch_names[N_ch]                   = {"E", "M", "L"};
  const char * ch_titles[N_ch]                  = {"$e$", "$\\mu$", "Combined $\\ell$"};
  const unsigned short N_histo_names            = 8;
  const char * histo_names[N_histo_names]       = {"", "t#bar{t}", "Single top", "W", "DY", "Multiboson", "t#bar{t}+V", "QCD"};
  const char * MC_histo_names[N_MC_histo_names] = {"t#bar{t}", "Single top", "W", "DY", "Multiboson", "t#bar{t}+V", "QCD"};
  system("mkdir -p event_yields_tables");
  for (unsigned short level_ind  = 0; level_ind < N_levels; level_ind ++)
    {
      FILE * file                                   = fopen(TString("event_yields_tables/event_yields_table")  + "_" + level_names[level_ind] + ".txt", "w");
      insert_header(file, level_ind);
      for (unsigned short ch_ind = 0; ch_ind < N_ch; ch_ind ++)
	{
	  insert_channel_header(file, ch_titles[ch_ind]);
	  double MC_sum[N_bins];
	  double MC_results[N_MC_histo_names][N_bins];
	  double MC_results_cflip[N_bins];
	  map<TString, float *> mmcsyst;
	  for (unsigned char sind = 0; sind < N_MC_histo_names + 1; sind ++)
	    {
	      mmcsyst[MC_syst_labels[sind]] = new float[4];
	      for (unsigned char bind = 0; bind < 4; bind ++)
		{
		  mmcsyst[MC_syst_labels[sind]][bind] = 0.0;
		}
	    }
	  for (unsigned short bin_ind = 0; bin_ind < N_bins; bin_ind ++)
	    {
	      MC_sum[bin_ind] = 0.0;
	      MC_results_cflip[bin_ind] = 0.0;
	      for (unsigned short h_ind = 0; h_ind < N_MC_histo_names; h_ind ++)
		{
	      
		  MC_results[h_ind][bin_ind] = 0.0;
		}
	    }
	  TH1 * h_data = NULL;
	  TH1 * h_MC_Unc = NULL;
	  const TString dir(TString(ch_names[ch_ind]) + "_" + level_names[level_ind] + "_selection");
	  TList * list = plotter -> GetDirectory(dir) -> GetListOfKeys();
	  TObject * it = list -> First();
	  do
	    {
	      TObject * obj = plotter  -> Get(dir + "/" + it -> GetName()); 
	      // printf("[%s] %s\n", it -> GetName(), obj -> ClassName());
	      
	      if (obj -> InheritsFrom("TH1"))
		{
		  TH1 * h = (TH1 *) obj;
		  bool isMC = false;
		  bool isMCsyst = false;
		   if (h -> GetName() == dir)
		     {
		       printf("data histo found \n");
		       h_data = h;
		     }
		   if (h -> GetName() == TString(dir + "_totalMCUnc"))
		     {
		       printf("MC unc histo found\n");

		       h_MC_Unc = h;
		     }
		   if (not TString(h->GetName()).Contains("_unc_"))
		     {
		       unsigned short MC_sample_ind = 0;
		       do
			 {
			   const TString comp = dir + "_" + MC_histo_names[MC_sample_ind];
			   if (comp == h -> GetName())
			     {
			       isMC = true;
			       for (unsigned short bin_ind = 1; bin_ind <= h -> GetNbinsX(); bin_ind ++)
				 {
				   MC_sum[bin_ind - 1] += h -> GetBinContent(bin_ind);
				   MC_results[MC_sample_ind][bin_ind - 1] = h -> GetBinContent(bin_ind);
				 }
			     }
			   MC_sample_ind ++;
			 } while (not isMC and MC_sample_ind < N_MC_histo_names); 
		     }
		   if (not isMC)
		    {
		      unsigned short MC_sample_ind_syst      = 0;
		      do
			{
			  unsigned char MC_subsample_ind_syst = 0;
			  do
			    {
			      // printf("%u %u\n", MC_sample_ind_syst,  MC_subsample_ind_syst);
			      // printf("%s\n", subMC[MC_sample_ind_syst][MC_subsample_ind_syst]);
			      const TString comp = dir + "_unc_" + subMC[MC_sample_ind_syst][MC_subsample_ind_syst];
			      if (comp == h -> GetName())
				{
				  isMCsyst = true;
				  for (unsigned char bind = 1; bind < h -> GetNbinsX() + 1; bind ++)
				    {
				      mmcsyst[MC_syst_labels[MC_sample_ind_syst]][bind - 1] = 
					TMath::Sqrt(
						    TMath::Power(mmcsyst[MC_syst_labels[MC_sample_ind_syst]][bind - 1], 2) + 
						    TMath::Power(h -> GetBinError(bind), 2)
						    );
				    }
			  
				} 
			      if (isMCsyst and TString(subMC[MC_sample_ind_syst][MC_subsample_ind_syst]) == "stat")
				{
				  hstat = h;
				}
			      if (isMCsyst and TString(subMC[MC_sample_ind_syst][MC_subsample_ind_syst]) == "stat" and ch_ind == 2 and level_ind == 0)
				{
				  printf("checking statistics\n");
				  for (unsigned char bind = 1; bind < h -> GetNbinsX() + 1; bind ++)
				    {
				      printf("%f %f\n", h -> GetBinContent(bind), h -> GetBinError(bind));
				    }				  
				}
			      MC_subsample_ind_syst ++;
			    } while (not isMCsyst and MC_subsample_ind_syst < subMC_N[MC_sample_ind_syst]);
			  MC_sample_ind_syst ++;
			} while (not isMCsyst and MC_sample_ind_syst < N_MC_histo_names + 1); 
		    }
		  // printf("\n errors \n");
		  // for (unsigned short bin_ind = 0; bin_ind <= h -> GetNbinsX() + 1; bin_ind ++)
		  //   {
		  //     printf("bin %u content %f\n", bin_ind, h -> GetBinError(bin_ind));
		  //   }

		  // printf("\nhisto listed \n");
		}
	      it = list -> After(it);
	    } while (it != list -> After(list -> Last()));
	  //TH1F * h_cflip = (TH1F*) plotter_cflip -> GetDirectory(dir) -> Get(dir + "_t#bar{t} cflip");
	  // if (not plotter_cflip)
	  //   printf("cflip file not found\n");
	  // for (unsigned short bin_ind = 1; bin_ind < h -> GetNbinsX() + 1; bin_ind ++)
	  //   {
	  //     MC_results_cflip[bin_ind - 1] = 1.0; //h_cflip -> GetBinContent(bin_ind);
	  //   }

	  insert_MC(file, MC_results, mmcsyst);
	  insert_MC_sum(file, MC_sum, h_MC_Unc, hstat);
	  //insert_MC_unc(file, h_MC_Unc);
	  insert_data(file, h_data);
	  for (map<TString, float *>::iterator it = mmcsyst.begin(); it != mmcsyst.end(); it ++)
	    {
	      delete [] it -> second;
	    }
	}
      printf("closing file\n");
      insert_footer(file);
      fclose(file);

    }
  /*  TApplication app("myapp", 0, 0);
  app.Run(kTRUE);
  app.Terminate();*/
  plotter -> Close();
  //  plotter_cflip -> Close();
}

void insert_header(FILE * file, unsigned char level_ind)
{
  if (environment.compare("lx") != 0)
    fprintf(file, "\\begin{longtable}{@{\\extracolsep{\\fill}}p{4cm}S[table-format=9.1]S[table-format=7.1]S[table-format=9.1]S[table-format=7.1]S[table-format=9.1]S[table-format=7.1]S[table-format=9.1]S[table-format=7.1]}\n");
  else
    fprintf(file, "\\begin{longtable}{lrrrrrrrr}\n");
  fprintf(file, "\\caption{Event yields at %s level.}\n", level_names[level_ind]);
  fprintf(file, "\\label{tab:yields_%s}\\\\\n", level_names[level_ind]);
  fprintf(file, "\\noalign{\\global\\arrayrulewidth=0.5mm}\\hline");
  if (environment.compare("lx") != 0)
    fprintf(file, "\\makecell{\\textbf{Process}} & \\multicolumn{2}{c}{\\boldmath$1 \\ell$} & \\multicolumn{2}{c}{\\boldmath$1 \\ell + \\geq 4 j $} & \\multicolumn{2}{c}{\\boldmath$1 \\ell + \\geq 4 j (2 b)$} & \\multicolumn{2}{c}{\\boldmath$1 \\ell + 4 j (2 b, 2 lj)$} \n");
  else
    fprintf(file, "\\textbf{Process} & \\multicolumn{2}{c}{\\boldmath$1 \\ell$} & \\multicolumn{2}{c}{\\boldmath$1 \\ell + \\geq 4 j $} & \\multicolumn{2}{c}{\\boldmath$1 \\ell + \\geq 4 j (2 b)$} & \\multicolumn{2}{c}{\\boldmath$1 \\ell + 4 j (2 b, 2 lj)$} \n");

}

void insert_channel_header(FILE * file, const char * channel)
{
  fprintf(file, "\\\\\\noalign{\\global\\arrayrulewidth=0.4pt}\\hline\n");
  fprintf(file, "\\multicolumn{9}{c}{\\textbf{{\\boldmath %s} + jets channel}}\\\\\n", channel);
  fprintf(file, "\\hline\n");
}

void insert_MC(FILE * file, double MC_results[][N_bins], map<TString, float *> mmcyst)
{
  for (unsigned short MC_samples_ind = 0; MC_samples_ind < N_MC_histo_names; MC_samples_ind ++)
    {
      fprintf(file, "%*s\t", 25, MC_syst_labels[MC_samples_ind]);
       for (unsigned short bin_ind = 0; bin_ind < N_bins; bin_ind ++)
	 {
	   fprintf(file, "& %s & %s %s ", insertnumber(MC_results[MC_samples_ind][bin_ind], "20.1").Data(), pm, insertnumber(mmcyst[MC_syst_labels[MC_samples_ind]][bin_ind], "20.1").Data());
	 }
       fprintf(file, "\\\\\n");			
       /*       if (MC_samples_ind == 0)
	 {
	   fprintf(file, "\\small %*s\t", 25, "t\\overline{t} cflip");
	   for (unsigned short bin_ind = 0; bin_ind < N_bins; bin_ind ++)
	     {
	       fprintf(file, "& \\small %20.1f\t\t", MC_results_cflip[bin_ind]);
	     }
	   fprintf(file, "\\\\\n");	
	   }*/
    }
  fprintf(file, "%*s\t", 25, MC_syst_labels[N_MC_histo_names]);
  for (unsigned short bin_ind = 0; bin_ind < N_bins; bin_ind ++)
    {
      fprintf(file, "&   & %s %s", pm, insertnumber(mmcyst[MC_syst_labels[N_MC_histo_names]][bin_ind] * TMath::Sqrt(2), "20.1").Data());
    }
  fprintf(file, "\\\\\n");			

  //  fprintf(file, "\\hline\n");
}


void insert_MC_sum(FILE * file, double * MC_sum, TH1 * h_MC_Unc, TH1 * hstat)
{
  fprintf(file, "\\hline\n");

  fprintf(file, "%*s", 25, "Total MC\t");
  for (unsigned short bin_ind = 1; bin_ind < N_bins + 1; bin_ind ++)
    {
      double err = TMath::Sqrt(TMath::Power(h_MC_Unc -> GetBinError(bin_ind), 2) + TMath::Power(hstat -> GetBinError(bin_ind), 2));
      
      fprintf(file, "& %s & %s %s", insertnumber(MC_sum[bin_ind - 1], "20.1").Data(), pm, insertnumber(err, "16.1").Data());

    }
  fprintf(file, "\\\\\n");
}

void insert_MC_unc(FILE * file, TH1 * h_MC_Unc)
{
  fprintf(file, "%*s", 25, "(\\\\\\ttbar uncertainty)\t");
  for (unsigned short bin_ind = 1; bin_ind < N_bins + 1; bin_ind ++)
    {
      fprintf(file, "& %s %s\t", pm, insertnumber(h_MC_Unc -> GetBinError(bin_ind), "16.1").Data());
    }
  fprintf(file, "\\\\\n");
}

void insert_data(FILE * file, TH1 * h_data)
{
  fprintf(file, "%*s", 25, "Data\t");
  for (unsigned short bin_ind = 1; bin_ind < N_bins + 1; bin_ind ++)
    {
      fprintf(file, "& %s & %s %s ", insertnumber(h_data -> GetBinContent(bin_ind), "20.1").Data(), pm, insertnumber(h_data -> GetBinError(bin_ind) * TMath::Sqrt(2), "16.1").Data());
    }
  //  fprintf(file, "\\\\\n\\hline\n");
  
}

void insert_footer(FILE * file)
{
  fprintf(file, "\\\\\n\\noalign{\\global\\arrayrulewidth=0.5mm}\\hline\n");
  fprintf(file, "\\end{longtable}\n");
}

TString insertnumber(float n, const char * format)
{
  return Form((string("%") + format + "f").c_str(), n/*/2.0*/); 
  return TString("\\num{") + Form((string("%") + format + "f}").c_str(), n); 
}
