
#ifndef JetAnalysis_h
#define JetAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/VectorUtil.h"
#include <TH2D.h>
#include <TGraph.h>
#include <TMath.h>

#include <vector>



class JetAnalysis : public TSelector {
public :
  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> LorentzVector;
  TTreeReader    fReader;  //!the tree reader
  TTree       *fChain = 0;  //!pointer to the analyzed TTree or TChain

  // Readers to access the data (delete the ones you do not need).
  TTreeReaderValue<Long64_t> run = {fReader, "run"};
  TTreeReaderValue<Long64_t> lumi = {fReader, "lumi"};
  TTreeReaderValue<Long64_t> event = {fReader, "event"};
  TTreeReaderValue<Int_t> bunchCrossing = {fReader, "bunchCrossing"};
  TTreeReaderArray<LorentzVector> jet_p4 = {fReader, "jet_p4"};
  TTreeReaderArray<float> jet_neutralEmFrac = {fReader, "jet_neutralEmFrac"};
  TTreeReaderArray<float> jet_neutralHadFrac = {fReader, "jet_neutralHadFrac"};
  TTreeReaderArray<float> jet_muonFrac = {fReader, "jet_muonFrac"};
  TTreeReaderArray<int> jet_pdgid = {fReader, "jet_pdgid"};
  TTreeReaderArray<int> jet_id = {fReader, "jet_id"};
  TTreeReaderArray<int> L1EG_bx = {fReader, "L1EG_bx"};
  TTreeReaderArray<LorentzVector> L1EG_p4 = {fReader, "L1EG_p4"};
  TTreeReaderArray<int> L1EG_iso = {fReader, "L1EG_iso"};
  TTreeReaderArray<int> L1GtBx = {fReader, "L1GtBx"};


  JetAnalysis(TTree * /*tree*/ =0) { }
  virtual ~JetAnalysis() { }
  virtual Int_t  Version() const { return 2; }
  virtual void   Begin(TTree *tree);
  virtual void   SlaveBegin(TTree *tree);
  virtual void   Init(TTree *tree);
  virtual Bool_t  Notify();
  virtual Bool_t  Process(Long64_t entry);
  virtual Int_t  GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
  virtual void   SetOption(const char *option) { fOption = option; }
  virtual void   SetObject(TObject *obj) { fObject = obj; }
  virtual void   SetInputList(TList *input) { fInput = input; }
  virtual TList  *GetOutputList() const { return fOutput; }
  virtual void   SlaveTerminate();
  virtual void   Terminate();

  ClassDef(JetAnalysis,0);

private :
  TH2F * hjetKinReweight_;
  bool useEMfrac_;
  TGraph * gLumiReweight_;
  Long64_t runIndexCached_{-1};


  TH2D * hJetPtEtaEGeffDenom_;
  TH2D * hJetPtEtaEGeffNum_bxm2_;
  TH2D * hJetPtEtaEGeffNum_bxm1_;
  TH2D * hJetPtEtaEGeffNum_bx0_;
  TH2D * hJetPtEtaEGeffNum_bx1_;
  TH2D * hJetPtEtaEGeffNum_bx2_;

  TH2D * hJet30EGEtaPhi_;

  TH2D * hJetPtEtaFinOReffDenom_;
  TH2D * hJetPtEtaFinOReffNum_bxm2_;
  TH2D * hJetPtEtaFinOReffNum_bxm1_;
  TH2D * hJetPtEtaFinOReffNum_bx0_;
  TH2D * hJetPtEtaFinOReffNum_bx1_;
  TH2D * hJetPtEtaFinOReffNum_bx2_;

  TH1D * hJetL1ADenom_;
  TH1D * hJetL1ANum_bxm2_;
  TH1D * hJetL1ANum_bxm1_;
  TH1D * hJetL1ANum_bx0_;
  TH1D * hJetL1ANum_bx1_;
  TH1D * hJetL1ANum_bx2_;

  TH1D * hJetEGm1thrDenom_;
  TH1D * hJetEGm1thrNum_EGlow_;
  TH1D * hJetEGm1thrNum_EGmed_;
  TH1D * hJetEGm1thrNum_EGhigh_;

  TH1D * hJetEGdeltaR_bxm1_;
  TH1D * hJetEGdeltaR_bx0_;
  TH1D * hJetEGdeltaR_bx1_;

  template<typename T, typename... Args>
    T * newOutput(Args... args) {
      T * out = new T(args...);
      fOutput->Add(out);
      return out;
    };
};

#endif

#ifdef JetAnalysis_cxx
void JetAnalysis::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the reader is initialized.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  fReader.SetTree(tree);
}

Bool_t JetAnalysis::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}


#endif // #ifdef JetAnalysis_cxx
