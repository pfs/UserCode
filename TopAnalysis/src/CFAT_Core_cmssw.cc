#include "TopLJets2015/TopAnalysis/interface/CFAT_Core_cmssw.hh"
#include "TopLJets2015/TopAnalysis/interface/Definitions_cmssw.hh"

pf_cmssw::pf_cmssw()
{
  core_ptr_ = NULL;
  index_ = 65535;
}
double pf_cmssw::GetCharge() const
{
  return core_ptr_ -> PF.charge[index_];
}

unsigned short pf_cmssw::GetJetConstituentIndex() const
{
  return index_;
}

unsigned short pf_cmssw::GetJetIndex() const
{
  return core_ptr_ -> PF.jet_index[index_];
}

void  pf_cmssw::SetCore(CFAT_Core_cmssw & other) 
{
  core_ptr_ = & other;
}


TLorentzVector pf_cmssw::GetLorentzVector() const
{
  TLorentzVector ret;
  ret.SetPtEtaPhiM(core_ptr_ -> PF.pt[index_], core_ptr_ -> PF.eta[index_], core_ptr_ -> PF.phi[index_], core_ptr_ -> PF.mass[index_]);
  return ret;
}

CFAT_Core_cmssw::pf_iter_core::pf_iter_core() : CFAT_Core::pf_iter_core(&particle)
{
  
}

CFAT_Core_cmssw::pf_iter_core * const CFAT_Core_cmssw::END = new CFAT_Core_cmssw::pf_iter_core;

void CFAT_Core_cmssw::increment(pf_iter &iter)
{
  pf_iter_core* iter_cmssw = (pf_iter_core*) iter.GetIterCore();
  unsigned short jet_const_index = iter_cmssw -> GetPFCMSSW() -> index_;
  
  CFAT_Core_cmssw * core = iter_cmssw -> GetPFCMSSW() -> core_ptr_; 
  do
    {
      jet_const_index ++;
      //      printf("incrementing %u size %u\n", jet_const_index, PF.size);
    } 
  while (iter.GetTrackedVectorCode() != 255 and 
	 core -> PF.jet_index[jet_const_index] != core -> GetIndex(iter_cmssw -> tracked_vector_) and 
	 jet_const_index <= core -> PF.size);
  if (jet_const_index >= core -> PF.size)
    {
      //printf("setting to END\n");
      iter.SetIterCore(*END);
      return;
    }
  pf_cmssw * particle = iter_cmssw -> GetPFCMSSW();
  particle -> index_ = jet_const_index;
  particle -> core_ptr_ = core;
}

pf      * CFAT_Core_cmssw::pf_iter_core::operator -> ()
{
  return particle_ptr_;
}



pf_cmssw      * CFAT_Core_cmssw::pf_iter_core::GetPFCMSSW ()
{
  return (pf_cmssw*)particle_ptr_;
}

pf_cmssw      * const CFAT_Core_cmssw::pf_iter_core::GetPFCMSSW () const 
{
  return (pf_cmssw * const )particle_ptr_;
}

void CFAT_Core_cmssw::pf_iter_core::WhoAmI()
{
  printf("I am CFAT_Core_cmssw::pf_iter_core %p\n", this);
}

bool      CFAT_Core_cmssw::pf_iter_core::compare (const CFAT_Core::pf_iter_core * other) 
{
  return this != other or GetPFCMSSW() -> index_ != ((const CFAT_Core_cmssw::pf_iter_core* )other) -> GetPFCMSSW() -> index_;
}

VectorCode_t   CFAT_Core_cmssw::pf_iter_core::GetVectorCode() const 
{
  const unsigned char n = 4;
  const unsigned char codes[n] = {LEADING_JET, SCND_LEADING_JET, HAD_B, LEPT_B};
  for (VectorCode_t jet_code = 0; jet_code < n; jet_code ++)
    {
      if (GetPFCMSSW() -> GetJetIndex() == ((CFAT_Core_cmssw*)GetIter() -> GetCore()) -> GetIndex(codes[jet_code]))
	return codes[jet_code];

    }
  return 255;
}

CFAT_Core_cmssw::pf_iter_core::~pf_iter_core()
{
  printf("Deconstructing CFAT_Core_cmssw::pf_iter_core %p\n", this);
}

void CFAT_Core_cmssw::check(const char * intro) const
{
}


CFAT_Core_cmssw::CFAT_Core_cmssw()
{
 
  leading_light_jet_ptr_                 = nullptr;
  leading_light_jet_prtcl_all            = TLorentzVector(0.0, 0.0, 0.0, 0.0);
  leading_light_jet_prtcl_charged        = TLorentzVector(0.0, 0.0, 0.0, 0.0);
  leading_light_jet_index_               = 65535;

  second_leading_light_jet_ptr_          = nullptr;
  second_leading_light_jet_prtcl_all     = TLorentzVector(0.0, 0.0, 0.0, 0.0);
  second_leading_light_jet_prtcl_charged = TLorentzVector(0.0, 0.0, 0.0, 0.0);
  second_leading_light_jet_index_        = 65535;

  had_b_ptr_                             = nullptr;
  had_b_prtcl_all                        = TLorentzVector(0.0, 0.0, 0.0, 0.0);
  had_b_prtcl_charged                    = TLorentzVector(0.0, 0.0, 0.0, 0.0);
  had_b_index_                           = 65535;
  
  lepton_ptr_                            = nullptr;
  neutrino_ptr_                          = nullptr;
  lept_b_ptr_                            = nullptr;
  lept_b_prtcl_all                       = TLorentzVector(0.0, 0.0, 0.0, 0.0);
  lept_b_prtcl_charged                   = TLorentzVector(0.0, 0.0, 0.0, 0.0);
  lept_b_index_                          = 65535;
  leading_b_                             = nullptr;
  leading_b_prtcl_all                    = TLorentzVector(0.0, 0.0, 0.0, 0.0);
  leading_b_prtcl_charged                = TLorentzVector(0.0, 0.0, 0.0, 0.0);

  leading_b_index_                       = 65535;

  scnd_leading_b_                        = nullptr;
  scnd_leading_b_prtcl_all               = TLorentzVector(0.0, 0.0, 0.0, 0.0);
  scnd_leading_b_prtcl_charged           = TLorentzVector(0.0, 0.0, 0.0, 0.0);
  scnd_leading_b_index_                  = 65535;

  PF.size                                = 65535;
  PF.jet_index                           = nullptr;
  PF.id                                  = nullptr; 
  PF.charge                              = nullptr;
  PF.pt                                  = nullptr;
  PF.eta                                 = nullptr; 
  PF.phi                                 = nullptr; 
  PF.mass                                = nullptr;
}



CFAT_Core::pf_iter CFAT_Core_cmssw::begin(VectorCode_t vector_code)
{
  unsigned short jet_const_index = 0;
  while (PF.jet_index[jet_const_index] != GetIndex(vector_code) and jet_const_index < PF.size)
    
      {
	jet_const_index ++;
      }
  pf_cmssw particle;
  particle.index_ = jet_const_index;
  particle.SetCore(*this);
  //pf_iter_core * iter_core = new pf_iter_core;
  static pf_iter_core * const RUN = new pf_iter_core;
  
  RUN -> tracked_vector_ = vector_code;
  RUN -> particle = particle;
  CFAT_Core::pf_iter iter;
  iter.SetCore(*this);
  iter.SetTrackedVectorCode(vector_code);
  iter.SetDeleteOption(false);
  iter.SetIterCore(*RUN);
  return iter;
}


CFAT_Core::pf_iter CFAT_Core_cmssw::end(VectorCode_t vector_code)
{
  unsigned short jet_const_index = PF.size;
  do  
      {
	jet_const_index --;

      }while (PF.jet_index[jet_const_index] != GetIndex(vector_code) and jet_const_index > 0);
  
  pf_cmssw particle;
  particle.index_ = jet_const_index;
  particle.SetCore(*this);
  END -> particle = particle;
  CFAT_Core::pf_iter iter;
  iter.SetTrackedVectorCode(vector_code);
  iter.SetDeleteOption(false);
  iter.SetIterCore(*END);
  return iter;
}

CFAT_Core::pf_iter CFAT_Core_cmssw::begin()
{
  unsigned short jet_const_index = 0;
  
  pf_cmssw particle;
  particle.index_ = jet_const_index;
  particle.SetCore(*this);
  //pf_iter_core * iter_core = new pf_iter_core;
  static pf_iter_core * const RUN = new pf_iter_core;
  
  RUN -> tracked_vector_ = 255;
  RUN -> particle = particle;
  CFAT_Core::pf_iter iter;
  iter.SetCore(*this);
  iter.SetDeleteOption(false);
  iter.SetIterCore(*RUN);
  return iter;
}


CFAT_Core::pf_iter CFAT_Core_cmssw::end()
{
  unsigned short jet_const_index = PF.size;
  
  pf_cmssw particle;
  particle.index_ = jet_const_index;
  particle.SetCore(*this);
  END -> particle = particle;
  CFAT_Core::pf_iter iter;
  iter.SetDeleteOption(false);
  iter.SetIterCore(*END);
  return iter;
}


MiniEvent_t * const CFAT_Core_cmssw::GetEvent() const
{
  return (MiniEvent_t * const) event_ptr_;
}

void CFAT_Core_cmssw::SetEvent(const MiniEvent_t & ev)
{
  event_ptr_ = &ev;
}



void CFAT_Core_cmssw::SetWorkMode(WorkCode_t work_code)
{
  work_mode_ = work_code;
  switch (work_mode_)
    {
    case Definitions::GEN:
      PF.size        = (unsigned short)GetEvent() -> ngpf;
      PF.jet_index   = GetEvent() -> gpf_g;
      PF.id          = GetEvent() -> gpf_id; 
      PF.charge      = GetEvent() -> gpf_c;
      PF.pt          = GetEvent() -> gpf_pt;
      PF.eta         = GetEvent() -> gpf_eta; 
      PF.phi         = GetEvent() -> gpf_phi; 
      PF.mass        = GetEvent() -> gpf_m;
    
      break;
    case Definitions::RECO:
      PF.size        = (unsigned short)GetEvent() -> npf;
      PF.jet_index   = GetEvent() -> pf_j;
      PF.id          = GetEvent() -> pf_id; 
      PF.charge      = GetEvent() -> pf_c;
      PF.pt          = GetEvent() -> pf_pt;
      PF.eta         = GetEvent() -> pf_eta; 
      PF.phi         = GetEvent() -> pf_phi; 
      PF.mass        = GetEvent() -> pf_m;
      break;
    default:
      ;
    }
}


unsigned short CFAT_Core_cmssw::GetIndex(VectorCode_t vector_code) const
{
  
  unsigned int index = 32767;
  switch(vector_code)
    {
    case LEADING_JET:
      return leading_light_jet_index_;
    case SCND_LEADING_JET:
      return second_leading_light_jet_index_;
    case HAD_B:
      return had_b_index_;
    case LEPT_B:
      return lept_b_index_;
    case LEADING_B:
      return leading_b_index_;
    case SCND_LEADING_B:
      return scnd_leading_b_index_;
    default:
      char error[128];
      //      sprintf(error, "unsigned short CFAT_Core::GetIndex(VectorCode_t) : Invalid vector code %u,", vector_code);
      printf("%s", error);
      throw error;
    }
}

const TLorentzVector *& CFAT_Core_cmssw::GetVectorRef(VectorCode_t code)
{
  switch(code)
    {
    case LEADING_JET:
      return leading_light_jet_ptr_;
    case SCND_LEADING_JET:
      return second_leading_light_jet_ptr_;
    case HAD_B:
      return had_b_ptr_;
    case LEPTON:
      return lepton_ptr_;
    case NEUTRINO:
      return neutrino_ptr_;
    case LEPT_B:
      return lept_b_ptr_;
    case LEADING_B:
      return leading_b_;
    case SCND_LEADING_B:
      return scnd_leading_b_;
    default:
      //printf("TLorentzVector *& CFAT_Core_cmssw::GetVectorRef(VectorCode_t) : check vector code %u ", code);
      throw "TLorentzVector *& CFAT_Core_cmssw::GetVectorRef(VectorCode_t) : check vector code";
    }
}

const TLorentzVector * CFAT_Core_cmssw::GetVector(VectorCode_t vector_code, const char * particle, ChargeCode_t charge_code)
{
  //  printf("vectore_code %u particle %s, charge code %u\n", vector_code, particle, charge_code);
  // printf("particle [%p]\n", particle);
  // printf("particle [%s]\n", particle);
  switch(vector_code)
    {
    case LEADING_JET:
      if (particle == nullptr)
	{
	  return leading_light_jet_ptr_;
	}
      else if (string(particle).compare("particle") == 0)
	{
	  switch(charge_code)
	    {
	    case ALLCOMP:
	      if (leading_light_jet_prtcl_all == TLorentzVector(0.0, 0.0, 0.0, 0.0))
		{
		  // printf("probe B\n");
		  return nullptr;
		}
	      else
		{
		  return & leading_light_jet_prtcl_all;
		}
	      
	    case CHARGED:
	      if (leading_light_jet_prtcl_charged == TLorentzVector(0.0, 0.0, 0.0, 0.0))
		{
		  return nullptr;
		}
	      else
		{
		  return & leading_light_jet_prtcl_charged;
		}
	    }
	}
      else
	{
	  char error[256];
	  sprintf(error, "const TLorentzVector * CFAT_Core_cmssw::GetVector(VectorCode_t vector_code): incorrect option %s", particle);
	  throw error;
	}
    case SCND_LEADING_JET:
      if (particle == nullptr)
	{
	  return second_leading_light_jet_ptr_;
	}
      else if (string(particle).compare("particle") == 0)
	{
	  switch(charge_code)
	    {
	    case ALLCOMP:
	      if (second_leading_light_jet_prtcl_all == TLorentzVector(0.0, 0.0, 0.0, 0.0))
		{
		  return nullptr;
		}
	      else
		{
		  return & second_leading_light_jet_prtcl_all;
		}
	      
	    case CHARGED:
	      if (second_leading_light_jet_prtcl_charged == TLorentzVector(0.0, 0.0, 0.0, 0.0))
		{
		  return nullptr;
		}
	      else
		{
		  return & second_leading_light_jet_prtcl_charged;
		}
	    }
	}
      else
	{
	  char error[256];
	  sprintf(error, "const TLorentzVector * CFAT_Core_cmssw::GetVector(VectorCode_t vector_code): incorrect option %s", particle);
	  throw error;
	}
    case HAD_B:
      if (particle == nullptr)
	{
	  return had_b_ptr_;
	}
      else if (string(particle).compare("particle") == 0)
	{
	  switch(charge_code)
	    {
	    case ALLCOMP:
	      if (had_b_prtcl_all == TLorentzVector(0.0, 0.0, 0.0, 0.0))
		{
		  return nullptr;
		}
	      else
		{
		  return & had_b_prtcl_all;
		}
	      
	    case CHARGED:
	      if (had_b_prtcl_charged == TLorentzVector(0.0, 0.0, 0.0, 0.0))
		{
		  return nullptr;
		}
	      else
		{
		  return & had_b_prtcl_charged;
		}
	    }
	}
      else
	{
	  char error[256];
	  sprintf(error, "const TLorentzVector * CFAT_Core_cmssw::GetVector(VectorCode_t vector_code): incorrect option %s", particle);
	  throw error;
	}

    case LEPTON:
      return lepton_ptr_;
    case NEUTRINO:
      return neutrino_ptr_;
    case LEPT_B:
      if (particle == nullptr)
	{
	  return lept_b_ptr_;
	}
      else if (string(particle).compare("particle") == 0)
	{
	  switch(charge_code)
	    {
	    case ALLCOMP:
	      if (lept_b_prtcl_all == TLorentzVector(0.0, 0.0, 0.0, 0.0))
		{
		  return nullptr;
		}
	      else
		{
		  return & lept_b_prtcl_all;
		}
	      
	    case CHARGED:
	      if (lept_b_prtcl_charged == TLorentzVector(0.0, 0.0, 0.0, 0.0))
		{
		  return nullptr;
		}
	      else
		{
		  return & lept_b_prtcl_charged;
		}
	    }
	}
      else
	{
	  char error[256];
	  sprintf(error, "const TLorentzVector * CFAT_Core_cmssw::GetVector(VectorCode_t vector_code): incorrect option %s", particle);
	  throw error;
	}
    case LEADING_B:
      if (particle == nullptr)
	{
	  return leading_b_;
	}
      else if (string(particle).compare("particle") == 0)
	{
	  switch(charge_code)
	    {
	    case ALLCOMP:
	      return & leading_b_prtcl_all;
	    case CHARGED:
	      return & leading_b_prtcl_charged;
	    }
	}
      else
	{
	  char error[256];
	  sprintf(error, "const TLorentzVector * CFAT_Core_cmssw::GetVector(VectorCode_t vector_code): incorrect option %s", particle);
	  throw error;
	}
    case SCND_LEADING_B:
      if (particle == nullptr)
	{
	  return scnd_leading_b_;
	}
      else if (string(particle).compare("particle") == 0)
	{
	  switch(charge_code)
	    {
	    case ALLCOMP:
	      return & scnd_leading_b_prtcl_all;
	    case CHARGED:
	      return & scnd_leading_b_prtcl_charged;
	    }
	}
      else
	{
	  char error[256];
	  sprintf(error, "const TLorentzVector * CFAT_Core_cmssw::GetVector(VectorCode_t vector_code): incorrect option %s", particle);
	  throw error;
	}
    default:
      //printf("TLorentzVector *& CFAT_Core_cmssw::GetVectorRef(VectorCode_t) : check vector code %u ", code);
      throw "TLorentzVector *& CFAT_Core_cmssw::GetVectorRef(VectorCode_t) : check vector code";
    }
} 

void CFAT_Core_cmssw::Reset()
{
  leading_light_jet_ptr_                  = nullptr;
  leading_light_jet_prtcl_all             = TLorentzVector(0.0, 0.0, 0.0, 0.0);
  leading_light_jet_prtcl_charged         = TLorentzVector(0.0, 0.0, 0.0, 0.0);
  leading_light_jet_index_                = 65535;

  second_leading_light_jet_ptr_           = nullptr;
  second_leading_light_jet_prtcl_all      = TLorentzVector(0.0, 0.0, 0.0, 0.0);
  second_leading_light_jet_prtcl_charged  = TLorentzVector(0.0, 0.0, 0.0, 0.0);
  second_leading_light_jet_index_         = 65535;
  
  had_b_ptr_                              = nullptr;
  had_b_prtcl_all                         = TLorentzVector(0.0, 0.0, 0.0, 0.0);
  had_b_prtcl_charged                     = TLorentzVector(0.0, 0.0, 0.0, 0.0);
  had_b_index_                            = 65535;
  
  lepton_ptr_                             = nullptr;
  neutrino_ptr_                           = nullptr;
  lept_b_ptr_                             = nullptr;
  lept_b_prtcl_all                        = TLorentzVector(0.0, 0.0, 0.0, 0.0);
  lept_b_prtcl_charged                    = TLorentzVector(0.0, 0.0, 0.0, 0.0);
  lept_b_index_                           = 65535;

  leading_b_                              = nullptr;
  leading_b_prtcl_all                     = TLorentzVector(0.0, 0.0, 0.0, 0.0);
  leading_b_prtcl_charged                 = TLorentzVector(0.0, 0.0, 0.0, 0.0);
  leading_b_index_                        = 65535;

  scnd_leading_b_                         = nullptr;
  scnd_leading_b_prtcl_all                = TLorentzVector(0.0, 0.0, 0.0, 0.0);
  scnd_leading_b_prtcl_charged            = TLorentzVector(0.0, 0.0, 0.0, 0.0);
  scnd_leading_b_index_                   = 65535;
}

void CFAT_Core_cmssw::RecomputeJetsFromParticles()
{
  leading_light_jet_prtcl_all             = GetJetFromParticles(LEADING_JET, ALLCOMP);
  // return;
  leading_light_jet_prtcl_charged         = GetJetFromParticles(LEADING_JET, CHARGED);
  second_leading_light_jet_prtcl_all      = GetJetFromParticles(SCND_LEADING_JET, ALLCOMP);
  second_leading_light_jet_prtcl_charged  = GetJetFromParticles(SCND_LEADING_JET, CHARGED);
  had_b_prtcl_all                         = GetJetFromParticles(HAD_B, ALLCOMP);
  had_b_prtcl_charged                     = GetJetFromParticles(HAD_B, CHARGED);
  lept_b_prtcl_all                        = GetJetFromParticles(LEPT_B, ALLCOMP);
  lept_b_prtcl_charged                    = GetJetFromParticles(LEPT_B, CHARGED);
  leading_b_prtcl_all                     = GetJetFromParticles(LEADING_B, ALLCOMP);
  leading_b_prtcl_charged                 = GetJetFromParticles(LEADING_B, CHARGED);
  scnd_leading_b_prtcl_all                = GetJetFromParticles(SCND_LEADING_B, ALLCOMP);
  scnd_leading_b_prtcl_charged            = GetJetFromParticles(SCND_LEADING_B, CHARGED);
}

void CFAT_Core_cmssw::AddLightJets(const vector<TLorentzVector> & light_jets, const vector<unsigned short> & light_jets_indices)
{
  /*
  printf("1 leading jet pt %f  phi %f rapdidity %f \n", light_jets[0].Pt(), light_jets[0].Phi(), light_jets[0].Rapidity());
  printf("2 leading jet pt %f  phi %f rapdidity %f \n", light_jets[1].Pt(), light_jets[1].Phi(), light_jets[1].Rapidity());
  */
  const bool IsFirstLeading = light_jets[0].Pt() >= light_jets[1].Pt(); 
  GetVectorRef(LEADING_JET) = IsFirstLeading          ?   &light_jets[0]        : &light_jets[1];
  leading_light_jet_index_ = IsFirstLeading        ?   light_jets_indices[0] : light_jets_indices[1];
  //  printf("Leading light jet index %u\n", leading_light_jet_index_);
  GetVectorRef(SCND_LEADING_JET) = IsFirstLeading   ?   &light_jets[1]        : &light_jets[0];
  second_leading_light_jet_index_ = IsFirstLeading ?   light_jets_indices[1] : light_jets_indices[0];
}

void CFAT_Core_cmssw::AddBJets(const vector<TLorentzVector> & b_jets, const vector<unsigned short> & b_jets_indices)
{
  if (neutrino_ptr_)
    {
      const TLorentzVector had_W = *leading_light_jet_ptr_ + *second_leading_light_jet_ptr_;
      const TLorentzVector lept_W = *lepton_ptr_ + *neutrino_ptr_;
      static const float t_mass = 173.34;
      double mass_dif[2][2] = 
	{
	  {1000, 1000}, 
	  {1000, 1000}
	};
      unsigned char min_index[2] = {2, 2}; 
      for (unsigned char index = 0; index < 2; index ++)
	{
	  mass_dif[index][0] = (b_jets[index] + had_W).M() - t_mass;
	  mass_dif[index][1] = (b_jets[index] + lept_W).M() - t_mass;
	}
      for (unsigned char index = 0; index < 2; index ++)
	{
	  min_index[index] = fabs(mass_dif[index][0]) >= fabs(mass_dif[index][1]) ? 1 : 0;
	}
      unsigned char had_b_local_index = 2;
      unsigned char lept_b_local_index = 2;
      if (min_index[0] != min_index[1])
	{
	  had_b_local_index  = min_index[0] == 0 ? 0 : 1;
	  lept_b_local_index = min_index[0] == 0 ? 1 : 0;
	  GetVectorRef(HAD_B)         = &b_jets.at(had_b_local_index);
	  GetVectorRef(LEPT_B)        = &b_jets.at(lept_b_local_index);
	  had_b_index_       = b_jets_indices.at(had_b_local_index);
	  lept_b_index_      = b_jets_indices.at(lept_b_local_index);
	}
    }
  if (b_jets.at(0) . Pt() >= b_jets.at(1) . Pt())
    {
      GetVectorRef(LEADING_B)                = &b_jets.at(0);
      leading_b_index_                       = b_jets_indices.at(0);
      GetVectorRef(SCND_LEADING_B)           = &b_jets.at(1);
      scnd_leading_b_index_                  = b_jets_indices.at(1);
    }
  else
    {
      GetVectorRef(LEADING_B)                = &b_jets.at(1);
      leading_b_index_                       = b_jets_indices.at(1);
      GetVectorRef(SCND_LEADING_B)           = &b_jets.at(0);
      scnd_leading_b_index_                  = b_jets_indices.at(0);
    }
    
  //  printf("leading b %p sec lb %p hb %p lb%p\n", GetVector(LEADING_B) , GetVector(SCND_LEADING_B), GetVector(HAD_B), GetVector(LEPT_B) );
}

void CFAT_Core_cmssw::AddVector(VectorCode_t code, const TLorentzVector * other)
{
  
  GetVectorRef(code) =  other;
}

