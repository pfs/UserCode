#include "CFAT_AssignHistograms_raw.hh"
#include "TopLJets2015/TopAnalysis/interface/CFAT_AssignHistograms.hh"
#include "TopLJets2015/TopAnalysis/interface/Definitions_cmssw.hh"
using namespace Definitions;

void AssignHistograms(map_hist_1D & target)
{
  AssignHistograms_raw<map_hist_1D>(target, CreateHistogram1D);
  AssignHistograms_raw_array<map_hist_1D>(target, CreateHistogram1Darray);
}
void AssignSpecificHistograms2D(map<TString, TH2*> & target)
{
  CreateHistogram2D(target, "Migration", "Migration : rec : gen", 100, -TMath::Pi(), TMath::Pi(), 30, -TMath::Pi(), TMath::Pi());
}
void CreateHistogram1D(map_hist_1D & target, const char * title, const char * axes, unsigned int nbins, double min , double max)
{
  for (unsigned short channel_ind = 0; channel_ind < N_channels_types_; channel_ind ++)
    {
      
      TH1F *h1 = new TH1F(TString(tag_channels_types_[channel_ind]) + "_" + title, TString(title) + axes, nbins, min, max);
      TH2F *h2 = nullptr;
      /*printf("nsyst_ %u %s \n", nsyst_, h1->GetName());
	getchar();*/
      if (nsyst_ > 0)
	{
	  h2 = new TH2F(TString(h1 -> GetName()) + "_syst", h1 -> GetTitle(), h1 -> GetNbinsX(), h1-> GetXaxis() -> GetXmin(), h1 -> GetXaxis() ->  GetXmax(), nsyst_ + 1, -0.5, nsyst_ + 0.5);
	}
      target[TString(tag_channels_types_[channel_ind]) + "_" + title] =  new pair<TH1*, TH2*>;
      target.at(TString(tag_channels_types_[channel_ind]) + "_" + title) -> first  = h1;
      // if ((TString(tag_channels_types_[channel_ind]) + "_" + title) == "L_pull_angle_allconst_reco_leading_jet_scnd_leading_jet_DeltaRTotal")
      // 	{
      // 	  if (h2)
      // 	    printf("creating %u\n", h2 -> GetNbinsY());
      // 	  getchar();
      // 	}
      target.at(TString(tag_channels_types_[channel_ind]) + "_" + title) -> second = h2;
    }
}

void CreateHistogram1Darray(map_hist_1D & target, const char * title, const char * axes, unsigned int nbins, double * bins)
{
  for (unsigned short channel_ind = 0; channel_ind < N_channels_types_; channel_ind ++)
    {
      
      TH1F *h1 = new TH1F(TString(tag_channels_types_[channel_ind]) + "_" + title, TString(title) + axes, nbins, bins);
      TH2F *h2 = nullptr;
      /*printf("nsyst_ %u %s \n", nsyst_, h1->GetName());
	getchar();*/
      if (nsyst_ > 0)
	{
	  h2 = new TH2F(TString(h1 -> GetName()) + "_syst", 
			h1 -> GetTitle(), 
			h1 -> GetNbinsX(), 
			h1-> GetXaxis() -> GetXbins() -> GetArray(), 
			nsyst_ + 1, 
			-0.5, 
			nsyst_ + 0.5);
	}
      target[TString(tag_channels_types_[channel_ind]) + "_" + title] =  new pair<TH1*, TH2*>;
      target.at(TString(tag_channels_types_[channel_ind]) + "_" + title) -> first  = h1;
      // if ((TString(tag_channels_types_[channel_ind]) + "_" + title) == "L_pull_angle_allconst_reco_leading_jet_scnd_leading_jet_DeltaRTotal")
      // 	{
      // 	  if (h2)
      // 	    printf("creating %u\n", h2 -> GetNbinsY());
      // 	  getchar();
      // 	}
      target.at(TString(tag_channels_types_[channel_ind]) + "_" + title) -> second = h2;
    }
}


void CreateHistogram2D(map<TString, TH2*> & target, const char * title, const char * axes, unsigned int nbins_x, double min_x , double max_x, unsigned int nbins_y, double min_y , double max_y)
{
  for (unsigned short channel_ind = 0; channel_ind < N_channels_types_; channel_ind ++)
    {
      target[TString(tag_channels_types_[channel_ind]) + "_" + title] = new TH2F(TString(tag_channels_types_[channel_ind]) + "_" + title, TString(title) + axes, nbins_x, min_x, max_x, nbins_y, min_y, max_y);
    }
}

// void CreateHistogram2D(map<TString, TH2*> & target, const char * title, const char * axes, unsigned int nbins_x, double * bins, unsigned int nbins_y, double min_y , double max_y)
// {
//   for (unsigned short channel_ind = 0; channel_ind < N_channels_types_; channel_ind ++)
//     {
//       target[TString(tag_channels_types_[channel_ind]) + "_" + title] = new TH2F(TString(tag_channels_types_[channel_ind]) + "_" + title, TString(title) + axes, nbins_x, bins, nbins_y, min_y, max_y);
//     }
// }

