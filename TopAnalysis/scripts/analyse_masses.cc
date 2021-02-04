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

static const unsigned short N_MC_histo_names = 7;
static const char * MC_histo_names[N_MC_histo_names] = {"t#bar{t}", "Single top", "W", "DY", "Multiboson", "t#bar{t}+V", "QCD"};
TH1F * sum_MC(const char *);
enum model_t {SM, CFLIP};
unsigned char model(0);
TFile * plotter_file(nullptr);// = model == SM ? TFile::Open("$EOS/analysis_24_12_2017/plots/plotter.root") : TFile::Open("$EOS/analysis_24_12_2017/MC13TeV_TTJets_cflip.root");
const char * model_tag[] = {"_nominal", "_cflip"};
const char * model_title[] = {" for the SM model", " for the colour octet $W$ model"};
const unsigned char N_interesting_jets = 4;
const unsigned char interesting_jets[N_interesting_jets] = {HAD_W, HAD_T, LEPT_W, LEPT_T};
double CalculateMass(TH1 *, void *params);
double CalculateWidth(TH1 *, void *params);
double determinant(unsigned int, void *);

void prepare_header(FILE *, const char *, const char *, SourceCode_t);
void prepare_footer(FILE *);
void add_entry(FILE *, Result, Result);
//Result CalculateResult(TH1 *, TH1 *, double (*)(TH1*, void *));
const char * env = nullptr;
int main(int argc, char* argv[])
{
  gROOT -> SetBatch(kTRUE);
  if (argc < 3)
    throw "Included must be at least 1 program argument";
  model = stoi(argv[1]);
  env = argv[2];
  plotter_file = model == SM ? TFile::Open("$EOS/analysis_MC13TeV_TTJets/plots/plotter.root") : TFile::Open("$EOS/analysis_MC13TeV_TTJets_cflip/plots/plotter.root");
  const char * cflip_title = "t#bar{t} cflip";
  if (model == CFLIP)
    {
      MC_histo_names[0] = cflip_title;
    }
  map<TString, FILE * > file_map;
  for (unsigned short ch_ind = L; ch_ind < N_channels_types_; ch_ind ++)
    {
      for (unsigned short level_ind = RECO; level_ind < N_levels_types_; level_ind ++)
	{
	  for (unsigned short type_ind = 0; type_ind < 2; type_ind ++)
	    {
	      const TString file_name = TString("mass_") + tag_channels_types_[ch_ind] + "_" + tag_levels_types_[level_ind] + "_" + tag_sources_types_[type_ind] + model_tag[model] + ".txt";
	      file_map[file_name] = fopen((TString("masses/tables/") + file_name).Data(), "w");
	      char caption[128];
	      sprintf(caption, "Observed masses of objects for the %s channel at %s level for %s%s.", channel_titles_[ch_ind], level_titles_[level_ind], title_sources_types_[type_ind], model_title[model]); 
	      char label[128];
	      sprintf(label, "mass_%s_%s_%s%s", tag_channels_types_[ch_ind], tag_levels_types_[level_ind], tag_sources_types_[type_ind], model_tag[model]); 
	      prepare_header(file_map[file_name], caption, label, type_ind);
	    }
	}
    }
  //TApplication app("myapp", 0, 0);
  for (unsigned short ch_ind = L; ch_ind < N_channels_types_; ch_ind ++)
    {
      for (unsigned short level_ind = RECO; level_ind < N_levels_types_; level_ind ++)
	{
	  for (unsigned short type_ind = 0; type_ind < (model == SM ? 2 : 1); type_ind ++) //interface/Definitions_cmssw.hh
	    {
	      for (unsigned short jet_ind = 0; jet_ind < N_interesting_jets; jet_ind ++)
		{	 
		  const unsigned char interesting_jet_ind = interesting_jets[jet_ind];
		  const TString dir = TString(tag_channels_types_[ch_ind]) + "_jet_mass_" + tag_levels_types_[level_ind] + "_" + tag_jet_types_[interesting_jet_ind];
		  // printf("dir %s\n", dir.Data());
		  // printf("ch_ind %u E %u\n", ch_ind, E);
		  // printf("level_ind %u RECO %u\n", level_ind, RECO);
		  // printf("jet_ind %u HAD_W %u\n", level_ind, RECO);
		  // if (ch_ind != L or level_ind != GEN or interesting_jet_ind !=HAD_T )
		  //   continue;

		  const TString file_name = TString("mass_") + tag_channels_types_[ch_ind] + "_" + tag_levels_types_[level_ind] + "_" + tag_sources_types_[type_ind] + model_tag[model] + ".txt";
		  FILE * table_file = file_map[file_name];
		  fprintf(table_file, "\t\t%s ", jet_titles_[interesting_jet_ind]);
		  if (level_ind == GEN and type_ind == DATA)
		    {
		      fprintf(table_file, "\t& x\\\\\n");
		      continue;
		    }
		  printf("analysing type %s channel %s level %s source %s jet %s\n", type_ind == 0 ? "MC" : "DATA", tag_channels_types_[ch_ind], tag_levels_types_[level_ind], tag_sources_types_[type_ind], jet_titles_[interesting_jet_ind]);
		  TH1F * h_chi(nullptr);
		  TH1F * h_chi_syst(nullptr);
		  h_chi = type_ind == 0 ? sum_MC(dir) : (TH1F*) plotter_file -> GetDirectory(dir) -> Get(dir);
		  //		  plotter_file -> GetDirectory(dir) -> ls();
		  // printf("h_chi %p type_ind %u\n", h_chi, type_ind);
		  // plotter_file -> Get(dir) -> ls();

		  h_chi_syst = type_ind == 0 ? (TH1F *) plotter_file -> Get(dir + "/" + dir + "_" + "totalMCUncShape") : nullptr;
		  Result mass(calculate_result(h_chi, h_chi_syst, CalculateMass, nullptr, 3.0, 1.0));
		  mass.SetName("mass");
		  //		  printf("calculating width\n"); getchar();
		  Result width(calculate_result(h_chi, h_chi_syst, CalculateWidth, nullptr, 3.0, 1.0));
		  width.SetName("width");
		  add_entry(table_file, mass, width);

		  mass.ls();
		  width.ls();
		}
	    }
	}
    }
  for (unsigned short ch_ind = L; ch_ind < N_channels_types_; ch_ind ++)
    {
      for (unsigned short level_ind = RECO; level_ind < N_levels_types_; level_ind ++)
	{
	  for (unsigned short type_ind = 0; type_ind < 2; type_ind ++)
	    {

	      const TString file_name = TString("mass_") + tag_channels_types_[ch_ind] + "_" + tag_levels_types_[level_ind] + "_" + tag_sources_types_[type_ind] + model_tag[model] + ".txt";
	      prepare_footer(file_map[file_name]);
	      fclose(file_map[file_name]);
	    }
	}
    }

  printf("Done\n");
  //app.Run(kTRUE);
  //app.Terminate();

  plotter_file -> Close();

}

TH1F * sum_MC(const char * dir)
{
  TH1F * h = (TH1F*) (plotter_file -> GetDirectory(dir) -> Get(TString(dir) + "_" +  MC_histo_names[0]) -> Clone(TString(dir) + "_MCsum"));
  for (unsigned short MC_ind = 1; MC_ind < N_MC_histo_names; MC_ind ++)
    {
      h -> Add((TH1F*)  (plotter_file -> GetDirectory(dir) -> Get(TString(dir) + "_" + MC_histo_names[MC_ind])));

    }
  return h;
	       
}

double CalculateMass(TH1 * h, void *params)
{
  for (unsigned char bind = 1; bind <= h-> GetNbinsX(); bind ++)
    {
      //h -> SetBinError(bind, h -> GetBinError(bind) * 0.001);
    }
  const double y_max = h -> GetMaximum();
  int bin_max = 0;
  h -> GetBinWithContent(y_max, bin_max);
  const double x_max = h -> GetBinCenter(bin_max);
  //  printf("x_max %f y_max %f underflow %f\n", x_max, y_max, h -> GetBinContent(0));
  const double a = 1.0;
  const double b = -2*x_max;
  const double c = y_max + 2 * x_max * x_max;

  TF1 polynomial("polynomial", "[0]*x*x + [1]*x + [2]");
  polynomial.SetParameters(a, b, c);
  polynomial.SetLineColor(kGreen);
  gROOT -> SetBatch(kTRUE);
  // TCanvas canvas(TString("mass_histo") + h -> GetName(), TString("mass_histo") + h -> GetTitle());
  // h -> Draw();
  const double interval = TString(h -> GetName()).Contains("lept_w") ? 1.0: 20.0;
  TFitResultPtr r =h -> Fit(& polynomial, "QS", "", x_max - interval, x_max + interval);
  double results[3];
  polynomial.GetParameters(results);
  const double mass = - results[1]/(2.0 * results[0]);
  const double peak = results[2] - results[1] * results[1]/(4 * results[0]);
  //  canvas.SaveAs(TString("masses/plots/") + h -> GetName() + ((Long_t) mass) + ".png");
//  printf("chi2 %f mass %f\n", r -> Chi2(), mass);
  //delete r.Get();
  return mass; 
}


double CalculateWidth(TH1 * h, void *params)
{
  // printf("!!! calculating width\n");
  // for (unsigned char bind = 1; bind <= h-> GetNbinsX(); bind ++)
  //   {
  //     h -> SetBinError(bind, h -> GetBinError(bind) * 0.001);
  //   }
  auto avg = [](TH1 * h, unsigned char bind) -> double
    {
      return (h -> GetBinContent(bind - 1) + h -> GetBinContent(bind) + h -> GetBinContent(bind + 1))/3.0;
      return (h -> GetBinContent(bind - 2) + h -> GetBinContent(bind - 1) + h -> GetBinContent(bind) + h -> GetBinContent(bind + 1) + h -> GetBinContent(bind + 2))/5.0;
    };
  enum side {LEFT, RIGHT};
  const unsigned int nbins(h -> GetNbinsX());
  double x_FWHM[2];

  const double y_max = h -> GetMaximum();
  int bin_max = 0;
  h -> GetBinWithContent(y_max, bin_max);
  const double x_max = h -> GetBinCenter(bin_max);
  //printf("%s\n", h -> GetName());
  if (TString(h -> GetName()).Contains(tag_jet_types_[LEPT_W]))
    {
      for (unsigned int side_ind = LEFT; side_ind <= RIGHT; side_ind ++)
	{
	  unsigned short bin_ind = side_ind == LEFT ? 2 : h-> GetNbinsX() - 1; 
	  do
	    {
	      if (side_ind == LEFT)
		{
		  bin_ind ++;
		  if (bin_ind == nbins + 1)
		    throw "bin_ind = nbins  leptonic\n";
		}
	      else
		{ 
		  bin_ind --;
		  if (bin_ind == 0)
		    throw "bin_ind = 1 leptonic\n";
		}
	    }while (h -> GetBinContent(side_ind == LEFT ? bin_ind + 1 : bin_ind -1) < y_max/2.0);
	  //	  printf("side_ind %u bin_ind %u low edge %f width %f\n", side_ind, bin_ind, h -> GetBinLowEdge(bin_ind), h-> GetBinWidth(bin_ind)); 
	  x_FWHM[side_ind] = side_ind == LEFT ? h -> GetBinLowEdge(bin_ind) + h -> GetBinWidth(bin_ind) : h -> GetBinLowEdge(bin_ind); 
	}
      const double width = x_FWHM[RIGHT] - x_FWHM[LEFT];
      //printf("leptonic W width %f\n", width);
      
      return width;
    }
  //getchar();
  // printf("x_max %f y_max %f underflow %f\n", x_max, y_max, h -> GetBinContent(0));
  const double a = 1.0;
  const double b = -2*x_max;
  const double c = y_max + 2 * x_max * x_max;

  TF1 polynomial("polynomial", "[0]*x*x + [1]*x + [2]");
  polynomial.SetParameters(a, b, c);
  //TCanvas canvas("canvas_width", "canvas_width");
  //h -> Draw();
  polynomial.SetLineColor(kBlue);
  // h -> Print("all");
  //   printf("x_max %f \n", x_max);
  h -> Fit(& polynomial, "Q", "", x_max - 12.0, x_max + 12.0);

  //canvas.Update();
  //canvas.Modified();
  double results[3];
  polynomial.GetParameters(results);
  
  const double mass = - results[1]/(2.0 * results[0]);
  const double peak = results[2] - results[1] * results[1]/(4.0 * results[0]);
  //  printf("y_max %f peak %f mass %f\n", y_max, peak, mass);
  double y[3] = {peak*1.6/2.0, peak/2.0, peak*0.8/2.0};
  
  const char sign[] = {-1, 1};
  //TCanvas canvas;
  for (unsigned char side_ind = LEFT; side_ind <= RIGHT; side_ind ++)
    {
      double bin[3]; 
      double A[3][3];
      double coeff[3]; 
      double x[3];
      for (unsigned char ind = 0; ind < 3; ind ++)
	{
	  unsigned short bin_ind = side_ind == LEFT ? 1 : h -> GetNbinsX();
	  do
	    {
	      if (side_ind == LEFT)
		{
		  bin_ind ++;
		  // printf("ind %u bin_ind %u avg %f y[ind] %f\n", ind, bin_ind, avg(h, bin_ind -1), y[ind]);
		  if (bin_ind == nbins - 1 )
		    {
		      //		      h -> Print("all");
		      throw "bin_ind = nbins -1\n";
		    }
		}
	      else
		{
		  bin_ind --;
		  //printf("ind %u bin_ind %u avg %f y[ind] %f\n", ind, bin_ind, avg(h, bin_ind -1), y[ind]);
		  if (bin_ind == 2)
		    {
		      //h -> Print("all");
		      
		      throw "bin_ind = 2\n";
		    }
		}
	     
	    } while(avg(h, side_ind == LEFT ? bin_ind + 1 : bin_ind - 1) < y[ind]);
	  //	  h -> GetBinWithContent(y[ind], bin, side_ind == LEFT ? 0 : bin_max);
	  x[ind] = side_ind == LEFT ?  h -> GetBinLowEdge(bin_ind) + h -> GetBinWidth(bin_ind) : h -> GetBinLowEdge(bin_ind);
	  if (ind == 2)
	    x[ind] += sign[side_ind]*5;
	  //	  printf("side %u x[ind] %f y[ind] %f \n", side_ind, x[ind], y[ind]);
	  A[ind][0] = x[ind] * x[ind];
	  A[ind][1] = x[ind];
	  A[ind][2] = 1.0;
	}
      double x_min = x[0] < x[2] ? x[0] : x[2];
      double x_max = x[0] < x[2] ? x[2] : x[0];
      for (unsigned char coeff_ind = 0; coeff_ind < 3; coeff_ind ++)
	{
	  double M[3][3];
	  for (unsigned char ind_M_row = 0; ind_M_row < 3; ind_M_row ++)
	    {
	      for (unsigned char ind_M_col = 0; ind_M_col < 3; ind_M_col ++)
		{
		  M[ind_M_row][ind_M_col] = (ind_M_col == coeff_ind) ? y[ind_M_row] : A[ind_M_row][ind_M_col];
		}
	    }
	  coeff[coeff_ind] = determinant(3, M)/determinant(3, A);
	}
      
      
      TF1 FWHM_polynomial("FWHM_polynomial", "[0]*x*x + [1]*x + [2]", x[0], x[2]);
      //  printf("before fit %f %f %f\n", coeff[0], coeff[1], coeff[2]);
      FWHM_polynomial.SetParameters(coeff[0], coeff[1], coeff[2]);
      FWHM_polynomial.SetLineColor(kGreen);
      //printf("check %f\n", coeff[0] * x[1]*x[1] + coeff[1]*x[1] + coeff[2]);

      /*if (side_ind == LEFT) 
       */

      //h -> Draw();
      
      FWHM_polynomial.Draw("SAME");
      //  printf("x_min %f x_max %f\n", x_min, x_max);
      TFitResultPtr r(h -> Fit(& FWHM_polynomial, "QS", "", x_min, x_max));
      if (not r.Get())
	throw "failed fit\n";
      /* canvas.Update();
	 canvas.Modified();*/
      double results[3];
      FWHM_polynomial.GetParameters(results);
      //printf("after fit %f %f %f\n", results[0], results[1], results[2]);
      //printf("after fit check %f\n", results[0] * x[1]*x[1] + results[1]*x[1] + results[2]);
      //printf("x_min %f x_max %f\n", x_min, x_max);
      const double a = results[0]; const double b = results[1]; const double c = results[2] - y[1];
      const double discriminant = b*b - 4 * a * c;
      const double x1 = (- b + TMath::Sqrt(discriminant)) / (2 * a);
      const double x2 = (- b - TMath::Sqrt(discriminant)) / (2 * a);
      if (x1 > x_min and x1 < x_max and not (x2 > x_min and x2 < x_max))
	x_FWHM[side_ind] = x1; 
      else if (not (x1 > x_min and x1 < x_max) and x2 > x_min and x2 < x_max)
	x_FWHM[side_ind] = x2; 
      else
	x_FWHM[side_ind] = x[1]; 

	//throw "Not in range";
      /*printf("side %u solution %f\n", side_ind, x_FWHM[side_ind]);
      TMarker marker(x_FWHM[side_ind], y[1], 2);
      marker.SetMarkerColor(kRed);
      marker.Draw("SAME");
      canvas.Update();
      canvas.Modified();
      
      getchar();
      */      
    }
  const double width = x_FWHM[RIGHT] - x_FWHM[LEFT];
  //printf("width %f\n", width);
  //  printf("!!!!! width calculated\n");
  //  getchar();
  return width; 
}


void prepare_header(FILE * file, const char * caption, const char * label, SourceCode_t source)
{
  fprintf(file, "\\begin{table}[h!tbp]\n");
  fprintf(file, "\t\\centering\n");
  fprintf(file, "\t\\caption{%s}\n", caption);
  fprintf(file, "\t\\label{tab:%s}\n", label);
  fprintf(file, "\t\\begin{tabular}{l|rr}\n");
  fprintf(file, "\t\t\\noalign{\\global\\arrayrulewidth=0.5mm}\\hline\\noalign{\\global\\arrayrulewidth=0.4pt}\n");
  char err[128];
  if (source == MC)
    {
      sprintf(err, "%s", "$\\pm$(stat)$\\pm$(syst)");
    }
  else
    {
      sprintf(err, "%s", "$\\pm$(stat)");
    }
  if (string(env).compare("lx") == 0)
    {
      fprintf(file, "\t\t{\\bf Object} & {\\bf Mass%s [GeV]} & {\\bf FWHM%s [GeV]}\\\\\n", err, err);
    } 
  else
    {
      fprintf(file, "\t\t\\makecell[c]{\\bf Object}&\\makecell[c]{\\bf Mass%s [GeV]}&\\makecell[c]{\\bf FWHM%s [GeV]}\\\\\n", err, err);
    }
  fprintf(file, "\t\t\\hline\n");
}

void add_entry(FILE * table_file, Result mass_result, Result width_result)
{
  fprintf(table_file, "\t& %s \t& %s \\\\\n", mass_result.latex(".3e").Data(), width_result.latex(".3e").Data());
}

void prepare_footer(FILE *file)
 {
  fprintf(file, "\t\t\\noalign{\\global\\arrayrulewidth=0.5mm}\\hline\n");
  fprintf(file, "\t\\end{tabular}\n");
  fprintf(file, "\\end{table}\n");

}

double determinant(unsigned int N, void * pointer)
{
  double (* matrix)[N] = (double (*)[N]) pointer;
  double (* test)[N] = matrix;
  /*printf("listing matrix \n");
  for (unsigned char row_ind = 0; row_ind < N; row_ind ++)
    {
      for (unsigned char col_ind = 0; col_ind < N; col_ind ++)
	{
	  printf("%f ", *(*(matrix + row_ind) + col_ind)); 
	}
      printf("\n");
    }*/
  
  if (N == 1)
    return ** matrix;
  double det = 0.0;
  char coeff = 1;

  for (unsigned int col_ind = 0; col_ind < N; col_ind ++)
    {
      double minor[N - 1][N - 1];
      unsigned int minor_col_ind = 0;
      for (unsigned int minor_col_iter = 0; minor_col_iter < N; minor_col_iter ++)
	{
	  if (minor_col_iter == col_ind)
	    continue;

	  for (unsigned int minor_row_iter = 1; minor_row_iter < N; minor_row_iter ++)
	    {
	      const unsigned int minor_row_ind = minor_row_iter - 1;
	      minor[minor_row_ind][minor_col_ind] = *(*(matrix + minor_row_iter) + minor_col_iter); 
	      //printf("adding to minor %f minor_row_ind %u minor_col_ind %u \n", *(*(matrix + minor_row_iter) + minor_col_iter), minor_row_ind, minor_col_ind);
	    }
	  minor_col_ind ++;

	}    
      det += coeff * *(*matrix + col_ind) * determinant(N - 1, minor);
      // printf("coeff %d *(*matrix + ind) %f determinant(N - 1, minor) %f \n", coeff, *(*matrix + col_ind),  determinant(N - 1, minor));
      coeff *= -1;
      
    }
  return det;
}
