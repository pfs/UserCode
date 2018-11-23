#define JetAnalysis_cxx
// The class definition in JetAnalysis.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//   Begin():      called every time a loop on the tree starts,
//              a convenient place to create your histograms.
//   SlaveBegin():  called after Begin(), when on PROOF called only on the
//              slave servers.
//   Process():    called for each event, in this function you decide what
//              to read and fill your histograms.
//   SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//              called only on the slave servers.
//   Terminate():   called at the end of the loop on the tree,
//              a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("JetAnalysis.C")
// root> T->Process("JetAnalysis.C","some options")
// root> T->Process("JetAnalysis.C+")
//


#include "JetAnalysis.h"
#include "TRandom.h"
#include "TParameter.h"

void JetAnalysis::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

  // const double jetPtBinning[] = {0, 10, 20, 30, 50, 70, 90, 110, 140, 170, 200, 250, 300, 400};
  // const int nPtBins = sizeof(jetPtBinning)/sizeof(double)-1;
  // const char * ptVar = "Jet p_{T}^{EM} [GeV]";

  const double jetPtBinning[] = {30.0, 36.0, 43.0, 52.0, 63.0, 75.0, 91.0, 109.0, 131.0, 158.0, 190.0, 228.0, 274.0, 397.0, 574.0, 830.0, 1200.0};
  // const double jetPtBinning[] = {30.0, 40.0, 50.0, 65.0, 80.0, 100.0, 125.0, 150.0, 175.0, 200.0, 225.0, 250.0, 275.0, 300.0, 400.0, 500.0, 600.0};
  const int nPtBins = sizeof(jetPtBinning)/sizeof(double)-1;
  const char * ptVar = "Jet p_{T} [GeV]";

  const double jetEtaBinning[] = {0, 1., 1.5, 2., 2.25, 2.5, 2.75, 3., 3.5};
  // const double jetEtaBinning[] = {-3.5, -3., -2.75, -2.5, -2.25, -2., -1.5, -1., 0, 1., 1.5, 2., 2.25, 2.5, 2.75, 3., 3.5};
  const int nEtaBins = sizeof(jetEtaBinning)/sizeof(double)-1;

  hJetPtEtaEGeffDenom_ = newOutput<TH2D>("denom", Form(";Jet |#eta|;%s;Counts", ptVar), nEtaBins, jetEtaBinning, nPtBins, jetPtBinning);
  hJetPtEtaEGeffNum_bxm2_ = newOutput<TH2D>("num_bxm2", Form(";Jet |#eta|;%s;Counts", ptVar), nEtaBins, jetEtaBinning, nPtBins, jetPtBinning);
  hJetPtEtaEGeffNum_bxm1_ = newOutput<TH2D>("num_bxm1", Form(";Jet |#eta|;%s;Counts", ptVar), nEtaBins, jetEtaBinning, nPtBins, jetPtBinning);
  hJetPtEtaEGeffNum_bx0_ =  newOutput<TH2D>("num_bx0",  Form(";Jet |#eta|;%s;Counts", ptVar), nEtaBins, jetEtaBinning, nPtBins, jetPtBinning);
  hJetPtEtaEGeffNum_bx1_ =  newOutput<TH2D>("num_bx1",  Form(";Jet |#eta|;%s;Counts", ptVar), nEtaBins, jetEtaBinning, nPtBins, jetPtBinning);
  hJetPtEtaEGeffNum_bx2_ =  newOutput<TH2D>("num_bx2",  Form(";Jet |#eta|;%s;Counts", ptVar), nEtaBins, jetEtaBinning, nPtBins, jetPtBinning);

  hJet30EGEtaPhi_ = newOutput<TH2D>("jet30EGEtaPhi", ";L1EG #eta;L1EG #phi;Counts", 57, -3, 3, 72, -3.1415, 3.1415);

  hJetPtEtaFinOReffDenom_ = newOutput<TH2D>("denomFinOR", Form(";Jet |#eta|;%s;Counts", ptVar), nEtaBins, jetEtaBinning, nPtBins, jetPtBinning);
  hJetPtEtaFinOReffNum_bxm2_ = newOutput<TH2D>("numFinOR_bxm2", Form(";Jet |#eta|;%s;Counts", ptVar), nEtaBins, jetEtaBinning, nPtBins, jetPtBinning);
  hJetPtEtaFinOReffNum_bxm1_ = newOutput<TH2D>("numFinOR_bxm1", Form(";Jet |#eta|;%s;Counts", ptVar), nEtaBins, jetEtaBinning, nPtBins, jetPtBinning);
  hJetPtEtaFinOReffNum_bx0_ =  newOutput<TH2D>("numFinOR_bx0",  Form(";Jet |#eta|;%s;Counts", ptVar), nEtaBins, jetEtaBinning, nPtBins, jetPtBinning);
  hJetPtEtaFinOReffNum_bx1_ =  newOutput<TH2D>("numFinOR_bx1",  Form(";Jet |#eta|;%s;Counts", ptVar), nEtaBins, jetEtaBinning, nPtBins, jetPtBinning);
  hJetPtEtaFinOReffNum_bx2_ =  newOutput<TH2D>("numFinOR_bx2",  Form(";Jet |#eta|;%s;Counts", ptVar), nEtaBins, jetEtaBinning, nPtBins, jetPtBinning);

  hJetL1ADenom_ = newOutput<TH1D>("denomL1A", Form(";%s;Counts", ptVar), nPtBins, jetPtBinning);
  hJetL1ANum_bxm2_ = newOutput<TH1D>("numL1A_bxm2", Form(";%s;Counts", ptVar), nPtBins, jetPtBinning);
  hJetL1ANum_bxm1_ = newOutput<TH1D>("numL1A_bxm1", Form(";%s;Counts", ptVar), nPtBins, jetPtBinning);
  hJetL1ANum_bx0_  = newOutput<TH1D>("numL1A_bx0",  Form(";%s;Counts", ptVar), nPtBins, jetPtBinning);
  hJetL1ANum_bx1_  = newOutput<TH1D>("numL1A_bx1",  Form(";%s;Counts", ptVar), nPtBins, jetPtBinning);
  hJetL1ANum_bx2_  = newOutput<TH1D>("numL1A_bx2",  Form(";%s;Counts", ptVar), nPtBins, jetPtBinning);

  hJetEGm1thrDenom_ = newOutput<TH1D>("denomJetEGthr", Form(";%s;L1EG in BX -1 Efficiency", ptVar), nPtBins, jetPtBinning);
  hJetEGm1thrNum_EGlow_ = newOutput<TH1D>("numJetEGthr_eglow", Form("L1IsoEGlow;%s;Counts", ptVar), nPtBins, jetPtBinning);
  hJetEGm1thrNum_EGmed_ = newOutput<TH1D>("numJetEGthr_egmed", Form("L1IsoEGmed;%s;Counts", ptVar), nPtBins, jetPtBinning);
  hJetEGm1thrNum_EGhigh_ = newOutput<TH1D>("numJetEGthr_eghigh", Form("L1IsoEGhigh;%s;Counts", ptVar), nPtBins, jetPtBinning);

  hJetEGdeltaR_bxm1_ = newOutput<TH1D>("jetEGdeltaR_bxm1", "L1IsoEG20 in BX -1;#DeltaR(j,EG);Counts", 50, 0, 1);
  hJetEGdeltaR_bx0_ =  newOutput<TH1D>("jetEGdeltaR_bx0", "L1IsoEG20 in BX 0;#DeltaR(j,EG);Counts", 50, 0, 1);
  hJetEGdeltaR_bx1_ =  newOutput<TH1D>("jetEGdeltaR_bx1", "L1IsoEG20 in BX 1;#DeltaR(j,EG);Counts", 50, 0, 1);
}

void JetAnalysis::SlaveBegin(TTree * /*tree*/)
{ 
  hjetKinReweight_=nullptr;
  gLumiReweight_=nullptr;
  useEMfrac_=false;

  /*
  hjetKinReweight_ =(TH2F *) GetInputList()->FindObject("jetKinReweight");
  gLumiReweight_ = (TGraph *) GetInputList()->FindObject("lumiReweight");
  */
  auto getFlag = [this] (const char * name) -> bool {
    auto flag = (TParameter<bool> *) GetInputList()->FindObject(name);
    if ( flag != nullptr ) return flag->GetVal();
    return false;
  };

 // useEMfrac_ = getFlag("useEMfraction");
 
  TString option = GetOption();
  
}

Bool_t JetAnalysis::Process(Long64_t entry)
{
  using namespace ROOT::Math::VectorUtil;
  auto bit = [] (int i) { return 1<<i; };
  fReader.SetEntry(entry);

  // Central jet trigger
  // bool centralJetTrigger{false};
  // for(size_t iJet=0; iJet<jet_p4.GetSize(); ++iJet) {
  //   // HLT_PFJet500
  //   if ( std::abs(jet_p4[iJet].Eta()) < 1.5 and (jet_id[iJet]&8)==8 ) {
  //     centralJetTrigger = true;
  //   }
  // }
  // if ( not centralJetTrigger ) return true;
  
  // Drop 2016H hotspot
  if ( *run >= 280919 and *run <= 284044 ) {
    for(size_t iJet=0; iJet<jet_p4.GetSize(); ++iJet) {
      if ( std::abs(jet_p4[iJet].Eta()+2.81) < 0.2 and std::abs(jet_p4[iJet].Phi()-2.07) < 0.2 ) {
        return true;
      }
    }
  }

  float weight = 1.;
  if ( gLumiReweight_ != nullptr ) {
    if ( runIndexCached_ < 0 or gLumiReweight_->GetX()[runIndexCached_] != *run ) {
      runIndexCached_ = TMath::BinarySearch(gLumiReweight_->GetN(), gLumiReweight_->GetX(), (double) *run);
      if ( gLumiReweight_->GetX()[runIndexCached_] != *run ) {
        std::cout << "Run not in lumi reweight table: " << *run << std::endl;
        return true;
      }
    }
    weight *= gLumiReweight_->GetY()[runIndexCached_];
  }

  auto etaBinCut = [] (LorentzVector jet) -> bool {
    return std::abs(jet.Eta()) >= 2.5 and std::abs(jet.Eta()) < 3.;
  };

  bool vetoEvent{false};
  std::vector<LorentzVector> forwardJets;
  for(size_t iJet=0; iJet<jet_p4.GetSize(); ++iJet) {
    auto& jet = jet_p4[iJet];

    if(jet_pdgid[iJet]!=0) continue;

    // Consider only EM
    if ( useEMfrac_ ) jet *= jet_neutralEmFrac[iJet];
    
    if ( jet_muonFrac[iJet] > .5 ) continue;

    // Select jets for FinOR
    if ( std::abs(jet.Eta()) >= 2. and std::abs(jet.Eta()) < 3.25 ) {
      forwardJets.emplace_back(jet);
    }
    else if ( std::abs(jet.Eta()) >= 2.0 ) {
      // Veto events with any other forward jet not in 2.5-3
      vetoEvent = true;
    }

    float jweight = weight;
    if ( hjetKinReweight_ != nullptr ) {
      jweight *= hjetKinReweight_->GetBinContent(hjetKinReweight_->FindBin(std::abs(jet.Eta()), jet.Pt()));
    }

    int match_bx{0};
    int match_thr{0};
    for(size_t iEG=0; iEG<L1EG_p4.GetSize(); ++iEG) {
      if ( L1EG_p4[iEG].Pt() > 30. and (L1EG_iso[iEG] & 0x1) and DeltaR(L1EG_p4[iEG], jet) < 0.4 ) {
        match_bx |= bit(L1EG_bx[iEG]+2);
        if ( jet.Pt() > 30. and L1EG_bx[iEG] == -1 ) {
          hJet30EGEtaPhi_->Fill(L1EG_p4[iEG].Eta(), L1EG_p4[iEG].Phi(), jweight);
        }
      }
      if ( etaBinCut(jet) and L1EG_bx[iEG] == -1 and (L1EG_iso[iEG] & 0x1) ) {
        if ( DeltaR(L1EG_p4[iEG], jet) < 0.4 ) {
          if ( L1EG_p4[iEG].Pt() > 20. ) match_thr |= bit(0);
          if ( L1EG_p4[iEG].Pt() > 30. ) match_thr |= bit(1);
          if ( L1EG_p4[iEG].Pt() > 40. ) match_thr |= bit(2);
        }
      }
      if ( etaBinCut(jet) and (L1EG_iso[iEG] & 0x1) and L1EG_p4[iEG].Pt() > 20. ) {
        if ( L1EG_bx[iEG] == -1 ) hJetEGdeltaR_bxm1_->Fill(DeltaR(L1EG_p4[iEG], jet), jweight);
        if ( L1EG_bx[iEG] ==  0 ) hJetEGdeltaR_bx0_->Fill(DeltaR(L1EG_p4[iEG], jet), jweight);
        if ( L1EG_bx[iEG] ==  1 ) hJetEGdeltaR_bx1_->Fill(DeltaR(L1EG_p4[iEG], jet), jweight);
      }
    }

    hJetPtEtaEGeffDenom_->Fill(std::abs(jet.Eta()), jet.Pt(), jweight);
    if ( match_bx & bit(-2+2) ) hJetPtEtaEGeffNum_bxm2_->Fill(std::abs(jet.Eta()), jet.Pt(), jweight);
    if ( match_bx & bit(-1+2) ) hJetPtEtaEGeffNum_bxm1_->Fill(std::abs(jet.Eta()), jet.Pt(), jweight);
    if ( match_bx & bit( 0+2) ) hJetPtEtaEGeffNum_bx0_ ->Fill(std::abs(jet.Eta()), jet.Pt(), jweight);
    if ( match_bx & bit( 1+2) ) hJetPtEtaEGeffNum_bx1_ ->Fill(std::abs(jet.Eta()), jet.Pt(), jweight);
    if ( match_bx & bit( 2+2) ) hJetPtEtaEGeffNum_bx2_ ->Fill(std::abs(jet.Eta()), jet.Pt(), jweight);

    if ( etaBinCut(jet) ) {
      hJetEGm1thrDenom_->Fill(jet.Pt(), jweight);
      if ( match_thr & bit(0) ) hJetEGm1thrNum_EGlow_->Fill(jet.Pt(), jweight);
      if ( match_thr & bit(1) ) hJetEGm1thrNum_EGmed_->Fill(jet.Pt(), jweight);
      if ( match_thr & bit(2) ) hJetEGm1thrNum_EGhigh_->Fill(jet.Pt(), jweight);
    }
  }

  if ( forwardJets.size() == 1 and not vetoEvent ) {
    LorentzVector jet = forwardJets[0];
    if ( hjetKinReweight_ != nullptr ) {
      weight *= hjetKinReweight_->GetBinContent(hjetKinReweight_->FindBin(std::abs(jet.Eta()), jet.Pt()));
    }
    hJetPtEtaFinOReffDenom_->Fill(std::abs(jet.Eta()), jet.Pt(), weight);
    if ( etaBinCut(jet) ) hJetL1ADenom_->Fill(jet.Pt(), weight);
    if ( L1GtBx[0] ) {
      hJetPtEtaFinOReffNum_bxm2_->Fill(std::abs(jet.Eta()), jet.Pt(), weight);
      if ( etaBinCut(jet) ) hJetL1ANum_bxm2_->Fill(jet.Pt(), weight);
    }
    if ( L1GtBx[1] ) {
      hJetPtEtaFinOReffNum_bxm1_->Fill(std::abs(jet.Eta()), jet.Pt(), weight);
      if ( etaBinCut(jet) ) hJetL1ANum_bxm1_->Fill(jet.Pt(), weight);
    }
    if ( L1GtBx[2] ) {
      hJetPtEtaFinOReffNum_bx0_->Fill(std::abs(jet.Eta()), jet.Pt(), weight);
      if ( etaBinCut(jet) ) hJetL1ANum_bx0_->Fill(jet.Pt(), weight);
    }
    if ( L1GtBx[3] ) {
      hJetPtEtaFinOReffNum_bx1_->Fill(std::abs(jet.Eta()), jet.Pt(), weight);
      if ( etaBinCut(jet) ) hJetL1ANum_bx1_->Fill(jet.Pt(), weight);
    }
    if ( L1GtBx[4] ) {
      hJetPtEtaFinOReffNum_bx2_->Fill(std::abs(jet.Eta()), jet.Pt(), weight);
      if ( etaBinCut(jet) ) hJetL1ANum_bx2_->Fill(jet.Pt(), weight);
    }
  }

  return kTRUE;
}

void JetAnalysis::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

}

void JetAnalysis::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

}
