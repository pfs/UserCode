#include "TopLJets2015/TopAnalysis/interface/CFAT_cmssw.hh"
#include "TopLJets2015/TopAnalysis/interface/Definitions_cmssw.hh"
#include "Definitions.hh"

#include "TTree.h"
#include "TFile.h"
#include <iostream>
CFAT_cmssw::CFAT_cmssw():ColourFlowAnalysisTool()
{

  const char * observables[2] = {"pull_angle", "pvmag"};
  float  (*obs[2])[2][2][4] = {pull_angle_, pvmag_};
  migration_tree_ = 
    new TTree("migration", "migration");
  migration_tree_ -> SetDirectory(nullptr);
  for (ChannelCode_t channel_code = 1; channel_code < N_channels_types_; channel_code ++)
    {
      for (unsigned short observable_ind = 0; observable_ind < 2; observable_ind ++)
	{
	  for (ChargeCode_t charge_code = 0; charge_code < 2; charge_code ++)
	    {
	      for (unsigned short jet1_iter = 0; jet1_iter < 4; jet1_iter ++)
		{
		  for (unsigned short level_ind = 0; level_ind < 2; level_ind ++)
		    {

		      TString obs_name = TString(tag_channels_types_[channel_code]) + "_" + 
						 observables[observable_ind] + "_" + 
						 jet_names_[jet1_iter] + "_" + 
						 tag_charge_types_[charge_code] + "_" + 
						 level_titles_[level_ind];
		      migration_tree_ -> 
			Branch(obs_name.Data(), &obs[observable_ind][channel_code - 1][level_ind][charge_code][jet1_iter],  (obs_name + "/F").Data());
		    }
		}

	    }
	}
    }
  migration_tree_ -> Branch("weight",     &weights_);
  migration_tree_ -> Branch("period",     &_period);
  migrationvaluesstored = false;
  // migration_tree_ -> Fill();
  // migration_tree_ -> Print("all");
  //   getchar();
}

void CFAT_cmssw::SetPeriod(TString period)
{
  _period = period;
}

void CFAT_cmssw::SetMigrationOutput(const char * migration_output)
{
  file_tag_ = TString(migration_output);
}

void CFAT_cmssw::ResetMigrationValues()
{
  _period = "noperiod";
  migrationvaluesstored = false;
  weights_.clear();
  for (ChannelCode_t channel_code = 1; channel_code < N_channels_types_; channel_code ++)
    {
      for (ChargeCode_t charge_code = 0; charge_code < 2; charge_code ++)
	{
	  for (unsigned short jet1_iter = 0; jet1_iter < 4; jet1_iter ++)
	    {
	      //fill_[charge_code][jet1_iter] = false;
	      for (WorkCode_t level_code = 0; level_code < N_levels_types_; level_code ++)
		{
		  pull_angle_[channel_code - 1][level_code][charge_code][jet1_iter] = - 10.0;
		  pvmag_[channel_code - 1][level_code][charge_code][jet1_iter] = - 1.0;
		}
	    }
	}
    }
}

void CFAT_cmssw::StoreMigrationValues(ChargeCode_t chargecode, VectorCode_t jetcode, double pa, double mag)
{
  weights_ = GetEvent() -> weights_;
  //  printf("storing %lu\n", GetEvent() -> weights_.size());
  //  fill_[chargecode][jetcode] = true;
  //  printf("channel_code_ %u %f\n", channel_code_ - 1, pa);
  pull_angle_[channel_code_ - 1][work_mode_][chargecode][jetcode] = pa;
  // printf("%f\n", pull_angle_[channel_code_ - 1][work_mode_][chargecode][jetcode]);;

  pvmag_[channel_code_ - 1][work_mode_][chargecode][jetcode] = mag;
  migrationvaluesstored = true;
}


void CFAT_cmssw::PlotMigrationValues()
{ 
  if (migrationvaluesstored)
    {
      //      printf("CFAT_cmssw::PlotMigrationValues %lu %lu\n", weights_.size(), GetEvent() -> weights_.size());
      // printf("filling\n");
      migration_tree_ -> Fill();
      // printf("period [%s]\n", _period.Data());
      // getchar();
    } 
  // for (ChargeCode_t charge_code = 0; charge_code < 2; charge_code ++)
  //   {
  //     for (unsigned short jet1_iter = 0; jet1_iter < 4; jet1_iter ++)
  // 	{
	  // if (fill_[charge_code][jet1_iter])
	  //   {
	  
	      // migration_tree_[channel_code_ - 1][charge_code][jet1_iter] -> Fill();
	      // if (charge_code == 1 and jet1_iter == 3)
	      // 	{
	      // 	  //		  printf("filling values wgt %.9f pa %f pvmag %f\n", weight_, pull_angle_[0][1][3], pvmag_[0][1][3]);
	      // 	  //getchar();
	      // 	}
    // 	    }
    // 	}
    // }
}


void CFAT_cmssw::Fill1D(const TString & key, double value, double weight) const                                                                                                                        
{                                                                                                                                                                                                       
  {
    TH1 * p1d = nullptr;
    TH2 * p2d = nullptr;
    try
      {
	p1d = ((map<TString, TH1 *> *) plots_ptr_)   -> at(TString(tag_channels_types_[channel_code_]) + "_" + key);
	// if (TString(tag_channels_types_[channel_code_]) + "_" + key == "L_pull_angle_allconst_reco_leading_jet_scnd_leading_jet_DeltaRTotal")
	//   {
	// 	printf("plotting %lu \n", weights_.size()); 
	//   }
      }
    catch(const  char *e)
      {
	//	printf("CFAT_cmssw::Fill1D: (A) caught an exception filling %s\n", key.Data());
      }
    p1d -> Fill(value, GetEvent() -> weights_[0] * weight);
    //printf("%u\n", Definitions::nsyst_);
    if (Definitions::nsyst_ > 0) 
      {
	p2d = ((map<TString, TH2 *> *) plots2D_ptr_) -> at(TString(tag_channels_types_[channel_code_]) + "_" + key + "_syst");
	if (GetEvent() -> weights_.size() > nsyst_)
	  cout << "WARNING: Size of uncertainty weight vector larger than uncertainty histogram size." << endl;
	p2d -> Fill(value, 0.0, GetEvent() -> weights_[0] * weight);
	for (unsigned int ind = 1; ind < GetEvent () -> weights_.size(); ++ ind) 
	  {
	    p2d -> Fill(value, ind, GetEvent() -> weights_[0] * GetEvent() -> weights_[ind] * weight);
	  }
      }
  }
  {
    TH1 * p1d = nullptr;
    TH2 * p2d = nullptr;
    try
      {
	p1d = ((map<TString, TH1 *> *) plots_ptr_) -> at(TString(tag_channels_types_[L]) + "_" + key);
      }
    catch(const  char *e)
      {
	//	printf("CFAT_cmssw::Fill1D: (B) caught an exception filling %s\n", key.Data());
      }
    
    p1d -> Fill(value, GetEvent() -> weights_[0] * weight);
    
    if (Definitions::nsyst_ > 0) 
      {
	if (GetEvent() -> weights_.size() > nsyst_)
	  cout << "WARNING: Size of uncertainty weight vector larger than uncertainty histogram size." << endl;
	p2d = ((map<TString, TH2 *> *) plots2D_ptr_) -> at(TString(tag_channels_types_[L]) + "_" + key + "_syst");
	p2d -> Fill(value, 0.0, GetEvent() -> weights_[0] * weight);
	for (unsigned int ind = 1; ind < GetEvent() -> weights_.size(); ++ ind) 
	  {
	    p2d -> Fill(value, ind, GetEvent() -> weights_[0] * GetEvent() -> weights_[ind] * weight);
	  }
	// const TString title(p -> first -> GetName());
      	// if (title == "L_pull_angle_allconst_reco_leading_jet_scnd_leading_jet_DeltaRTotal")
	//   {
	//     printf("size %lu\n", weights_.size());
	//     printf("after filling %s\n", title.Data());
	//     TH2 * h2 = p -> second;
	//     for (unsigned char bindy = 1; bindy <= h2 -> GetNbinsY(); bindy ++)
	//       {
	// 	const double integral = h2 -> Integral(0, h2 -> GetNbinsX() + 1, bindy, bindy);
	// 	printf("bindy %u %f\n", bindy, 1000000*integral);
	//       }
	//     getchar();
	//   }

      }
    // if (TString(tag_channels_types_[L]) + "_" + key == "L_pull_angle_allconst_reco_leading_jet_scnd_leading_jet_DeltaRTotal")
    //   {
    // 	printf("GetEvent() -> weights_[0] %.9f weight %.9f\n", GetEvent() -> weights_[0], weight);
    // 	for (unsigned char bind = 1; bind < p -> first -> GetNbinsX(); bind ++)
    // 	  {
    // 	    printf("%u value %f %.9f %.9f %.9f\n", bind, value, p -> first -> GetBinContent(bind), TMath::Sqrt(p -> first -> GetBinContent(bind)), p -> first -> GetBinError(bind));
    // 	  }
    // 	// getchar();
    // 	printf ("*******************\n");
    //   }

  }
}

void CFAT_cmssw::Fill2D(const TString & key, double value_x, double value_y, double weight) const
{                                                                                                                                                                                                         ((map<TString, TH2*>*)plots2D_ptr_) -> operator[](TString(tag_channels_types_[channel_code_]) + "_" + key) -> Fill(value_x, value_y, weight);                                                               ((map<TString, TH2*>*)plots2D_ptr_) -> operator[](TString(tag_channels_types_[L]) + "_" + key) -> Fill(value_x, value_y, weight);                                                                                                             
}


void CFAT_cmssw::WriteMigrationTree()
{
  //migration_tree_ -> SetBranchAddress("pull_angle_gen", &y);
  migration_file_ = TFile::Open(file_tag_, "RECREATE");
  printf("Opened migration file %p %s. It contains %llu events\n", migration_file_, file_tag_.Data(), migration_tree_ -> GetEntries());
  migration_tree_ -> SetDirectory(migration_file_);
  
  migration_tree_ -> Write();
  // for (ChannelCode_t channel_code = 1; channel_code < N_channels_types_; channel_code ++)
  //   {
      //     for (ChargeCode_t charge_code = 0; charge_code < 2; charge_code ++)
      // 	{
      // 	  for (unsigned short jet1_iter = 0; jet1_iter < 4; jet1_iter ++)
      // 	    {
      // 	      TTree * tree = migration_tree_[channel_code - 1][charge_code][jet1_iter];
      // 	      tree -> SetDirectory(migration_file_);
      // 	      printf("%s collected %lld entries\n", tree -> GetName(), tree -> GetEntries());
      // 	      tree -> Write();
      // 	    }
      // 	}
      //    }
  migration_file_ -> Close();
}

