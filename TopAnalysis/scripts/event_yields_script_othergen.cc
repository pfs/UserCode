// ./event_yields_script_othergen lx
#include "TFile.h"
#include "TApplication.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TMath.h"
#include <assert.h>
const float lumi = 35874.8;
const float xsec = 832.0;
const unsigned short N_bins = 4;
//const unsigned short N_MC_histo_names = 8;
const unsigned char Nsamples = 3;
const char * samples[Nsamples] = {"sherpa", "dire2002", "herwig7"};

TString insertnumber(float, const char *);


void insert_header(FILE *);
void insert_channel_header(FILE *, const char *);
void insert_values(FILE *, TH1 *);
void insert_footer(FILE *);
const char *pm = 0;
const char *pmlx = "$\\pm$";
using namespace std;
string environment;
int main(int argc, char * argv[])
{
  assert(argc == 2);
  environment = string(argv[1]);
  if (environment.compare("lx") == 0)
    pm = pmlx;
  FILE * file = fopen("event_yields_tables/event_yields_table_othergen.txt", "w");
  insert_header(file);
  const char * level_name = "gen";
  const unsigned short N_ch = 3;
  const char * ch_names[N_ch] = {"E", "M", "L"};
  const char * ch_titles[N_ch] = {"$e$", "$\\mu$", "Combined $\\ell$"};
  TFile *inputfile[Nsamples];
  TFile * ttbarfile =  TFile::Open("root://eosuser.cern.ch//$EOS/analysis_MC13TeV_TTJets/HADDChunks/MC13TeV_TTJets.root", "READ");
  for (unsigned char sampleind = 0; sampleind < Nsamples; sampleind ++)
    {
      inputfile[sampleind] = TFile::Open(TString("root://eosuser.cern.ch//$EOS/analysis_MC13TeV_TTJets_") + 
				      samples[sampleind] + 
				      "/HADDChunks/MC13TeV_TTJets_" + 
				      samples[sampleind] + 
					 ".root", "READ");
    }
  for (unsigned short ch_ind = 0; ch_ind < N_ch; ch_ind ++)
    {
      insert_channel_header(file, ch_titles[ch_ind]);
      const TString dir(TString(ch_names[ch_ind]) + "_" + level_name + "_selection");
      TH1 * httbar = (TH1*) ttbarfile -> Get(dir);
      for (unsigned char sampleind = 0; sampleind < Nsamples; sampleind ++)
	{
	  TH1 * hselection = (TH1*) inputfile[sampleind] -> Get(dir);
	  const double sumow = hselection -> GetSumOfWeights();
	  printf("sum of w %f\n", hselection -> GetSumOfWeights());
	  hselection -> Scale(lumi * xsec * httbar -> GetSumOfWeights()/ sumow);
	  hselection -> SetTitle(samples[sampleind]);
	  insert_values(file, hselection);
	}
    }


  for (unsigned char sampleind = 0; sampleind < Nsamples; sampleind ++)
    {
      inputfile[sampleind] -> Close();
    }
  ttbarfile -> Close();
  insert_footer(file);
  fclose(file);
}

void insert_header(FILE * file)
{
  if (environment.compare("lx") != 0)
    fprintf(file, "\\begin{longtable}{@{\\extracolsep{\\fill}}p{4cm}S[table-format=9.1]S[table-format=7.1]S[table-format=9.1]S[table-format=7.1]S[table-format=9.1]S[table-format=7.1]S[table-format=9.1]S[table-format=7.1]}\n");
  else
    fprintf(file, "\\begin{longtable}{lrrrrrrrr}\n");
  //  fprintf(file, "\\begin{longtable}{lrrrr}\n");
  fprintf(file, "\\caption{Event yields for the alternative generator samples.}\n");
  fprintf(file, "\\label{tab:yields_alternative}\\\\\n");
  fprintf(file, "\\noalign{\\global\\arrayrulewidth=0.5mm}\\hline");
  if (environment.compare("lx") != 0)
    fprintf(file, "\\makecell{\\textbf{Process}} & \\multicolumn{2}{c}{\\boldmath$1 \\ell$} & \\multicolumn{2}{c}{\\boldmath$1 \\ell + \\geq 4 j $} & \\multicolumn{2}{c}{\\boldmath$1 \\ell + \\geq 4 j (2 b)$} & \\multicolumn{2}{c}{\\boldmath$1 \\ell + 4 j (2 b, 2 lj)$} \n");
  else
    fprintf(file, "\\textbf{Process} & \\multicolumn{2}{c}{\\boldmath$1 \\ell$} & \\multicolumn{2}{c}{\\boldmath$1 \\ell + \\geq 4 j $} & \\multicolumn{2}{c}{\\boldmath$1 \\ell + \\geq 4 j (2 b)$} & \\multicolumn{2}{c}{\\boldmath$1 \\ell + 4 j (2 b, 2 lj)$} \n");


  // fprintf(file, "Process & $1 \\ell$ & $1 \\ell + \\geq 4 j$ & $1 \\ell + \\geq 4 j (2 b)$ & $1 \\ell + 4 j (2 b, 2 lj)$ \\\\\n");
  // fprintf(file, "\\hline\n");

}

void insert_channel_header(FILE * file, const char * channel)
{
  fprintf(file, "\\\\\\noalign{\\global\\arrayrulewidth=0.4pt}\\hline\n");

  fprintf(file, "\\multicolumn{9}{c}{\\textbf{{\\boldmath %s} + jets channel}}\\\\\n", channel);
  //  fprintf(file, "%s + jets channel\\\\\n", channel);
  fprintf(file, "\\hline\n");
}

void insert_values(FILE * file, TH1 *h)
{
  //  fprintf(file, "\\hline\n");

  fprintf(file, "%*s\t", 25, h -> GetTitle());
  for (unsigned short bin_ind = 1; bin_ind <= N_bins; bin_ind ++)
    {
      fprintf(file, "& %s\t\t & %s %s\t\t", insertnumber(h -> GetBinContent(bin_ind), "20.1").Data(), pm, insertnumber(h -> GetBinError(bin_ind), "16.1").Data());
    }
    fprintf(file, "\\\\\n");
}


void insert_footer(FILE * file)
{
  fprintf(file, "\\\\\n\\noalign{\\global\\arrayrulewidth=0.5mm}\\hline\n");
  fprintf(file, "\\end{longtable}\n");
}

TString insertnumber(float n, const char * format)
{
  return Form((string("%") + format + "f").c_str(), n); 
  return TString("\\num{") + Form((string("%") + format + "f}").c_str(), n); 
}
