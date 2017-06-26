//
// -*- C++ -*-
//
// Package:    TopLJets2015/TopAnalysis
// Class:      MiniAnalyzer
// 
/**\class MiniAnalyzer MiniAnalyzer.cc Test/MiniAnalyzer/plugins/MiniAnalyzer.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Qamar Ul Hassan
//         Created:  Sun, 13 Jul 2014 06:22:18 GMT
//
//
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TopLJets2015/TopAnalysis/interface/MyIPTools.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Ptr.h"

#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

#include "TopQuarkAnalysis/BFragmentationAnalyzer/interface/BFragmentationAnalyzerUtils.h"

#include "KaMuCa/Calibration/interface/KalmanMuonCalibrator.h"

#include "TLorentzVector.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"

#include <vector>
#include <unordered_map>
#include <memory>
#include <cmath>
#include <iostream>
#include <string>

using namespace edm;
using namespace std;
using namespace reco;
using namespace pat; 

//
// class declaration
//

class MiniAnalyzer : public edm::EDAnalyzer {
public:
  explicit MiniAnalyzer(const edm::ParameterSet&);
  ~MiniAnalyzer();  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  virtual void endRun(const edm::Run&,const edm::EventSetup&);  
private:
  virtual void beginJob() override;
  int genAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  int recAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  float getMiniIsolation(edm::Handle<pat::PackedCandidateCollection> pfcands,
			 const reco::Candidate* ptcl,  
                         float r_iso_min, float r_iso_max, float kt_scale,
                         bool charged_only);

  bool isSoftMuon(const reco::Muon & recoMu,const reco::Vertex &vertex);
  bool isMediumMuon2016ReReco(const reco::Muon & recoMu);

  // member data 
  edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
  edm::EDGetTokenT<GenEventInfoProduct> generatorevtToken_;
  edm::EDGetTokenT<LHEEventProduct> generatorlheToken_;
  edm::EDGetTokenT<LHERunInfoProduct> generatorRunInfoToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > puToken_;
  edm::EDGetTokenT<std::vector<reco::GenJet>  > genLeptonsToken_,   genJetsToken_;
  edm::EDGetTokenT<reco::METCollection> genMetsToken_;
  edm::EDGetTokenT<pat::PackedGenParticleCollection> genParticlesToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> prunedGenParticlesToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> pseudoTopToken_;

  edm::EDGetTokenT<edm::TriggerResults> triggerBits_,metFilterBits_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<double> rhoToken_;
  edm::EDGetTokenT<pat::MuonCollection> muonToken_;
  edm::EDGetTokenT<edm::View<pat::Electron>  >  electronToken_;
  edm::EDGetTokenT<edm::View<pat::Jet> > jetToken_;
  edm::EDGetTokenT<pat::METCollection> metToken_, puppiMetToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
  
  //Electron Decisions
  edm::EDGetTokenT<edm::ValueMap<float> > eleMvaIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleVetoIdMapToken_,eleLooseIdMapToken_,eleMediumIdMapToken_,eleTightIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > eleVetoIdFullInfoMapToken_,eleLooseIdFullInfoMapToken_,eleMediumIdFullInfoMapToken_,eleTightIdFullInfoMapToken_;
  unsigned int evetoIsoBit_, elooseIsoBit_, emediumIsoBit_, etightIsoBit_;
  edm::EDGetTokenT<bool> BadChCandFilterToken_,BadPFMuonFilterToken_;

  //  edm::EDGetTokenT<edm::ValueMap<float> > petersonFragToken_;

  std::unordered_map<std::string,TH1*> histContainer_;

  PFJetIDSelectionFunctor pfjetIDLoose_;

  std::vector<std::string> triggersToUse_,metFiltersToUse_;

  bool saveTree_,savePF_;
  TTree *tree_;
  MiniEvent_t ev_;
  
  KalmanMuonCalibrator *muonCor_;
  bool useRawLeptons_;
  edm::Service<TFileService> fs;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//


//
// constructors and destructor
//
MiniAnalyzer::MiniAnalyzer(const edm::ParameterSet& iConfig) :
  generatorToken_(consumes<GenEventInfoProduct>(edm::InputTag("generator"))),
  generatorevtToken_(consumes<GenEventInfoProduct>(edm::InputTag("generator",""))),
  generatorlheToken_(consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer",""))),
  generatorRunInfoToken_(consumes<LHERunInfoProduct,edm::InRun>({"externalLHEProducer"})),
  puToken_(consumes<std::vector<PileupSummaryInfo>>(edm::InputTag("slimmedAddPileupInfo"))),  
  genLeptonsToken_(consumes<std::vector<reco::GenJet> >(edm::InputTag("pseudoTop:leptons"))),
  genJetsToken_(consumes<std::vector<reco::GenJet> >(edm::InputTag("pseudoTop:jets"))),
  genMetsToken_(consumes<reco::METCollection>(edm::InputTag("pseudoTop:mets"))),
  genParticlesToken_(consumes<pat::PackedGenParticleCollection>(edm::InputTag("packedGenParticles"))),
  prunedGenParticlesToken_(consumes<reco::GenParticleCollection>(edm::InputTag("prunedGenParticles"))),
  pseudoTopToken_(consumes<reco::GenParticleCollection>(edm::InputTag("pseudoTop"))),
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerBits"))),
  metFilterBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("metFilterBits"))),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
  muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  jetToken_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jets"))),  
  metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
  puppiMetToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("puppimets"))),
  pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
  eleMvaIdMapToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("eleMvaIdMap"))),
  eleVetoIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleVetoIdMap"))),
  eleLooseIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleLooseIdMap"))),
  eleMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumIdMap"))),
  eleTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap"))),
  eleVetoIdFullInfoMapToken_(consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleVetoIdFullInfoMap"))),
  eleLooseIdFullInfoMapToken_(consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleLooseIdFullInfoMap"))),
  eleMediumIdFullInfoMapToken_(consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleMediumIdFullInfoMap"))),
  eleTightIdFullInfoMapToken_(consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleTightIdFullInfoMap"))),
  evetoIsoBit_(999), elooseIsoBit_(999), emediumIsoBit_(999), etightIsoBit_(999),
  BadChCandFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("badChCandFilter"))),
  BadPFMuonFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("badPFMuonFilter"))),
  //petersonFragToken_(consumes<edm::ValueMap<float> >(edm::InputTag("bfragWgtProducer:PetersonFrag"))),
  pfjetIDLoose_( PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::LOOSE ),  
  saveTree_( iConfig.getParameter<bool>("saveTree") ),
  savePF_( iConfig.getParameter<bool>("savePF") ),
  muonCor_(0),
  useRawLeptons_( iConfig.getParameter<bool>("useRawLeptons") )
{
  //now do what ever initialization is needed
  electronToken_      = mayConsume<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electrons"));
  triggersToUse_      = iConfig.getParameter<std::vector<std::string> >("triggersToUse");
  metFiltersToUse_  = iConfig.getParameter<std::vector<std::string> >("metFiltersToUse");

  histContainer_["triggerList"] = fs->make<TH1F>("triggerList", ";Trigger bits;",triggersToUse_.size(),0,triggersToUse_.size());
  for(size_t i=0; i<triggersToUse_.size(); i++) histContainer_["triggerList"] ->GetXaxis()->SetBinLabel(i+1,triggersToUse_[i].c_str());
  histContainer_["counter"]    = fs->make<TH1F>("counter", ";Counter;Events",2,0,2);
  histContainer_["fidcounter"] = (TH1 *)fs->make<TH2F>("fidcounter",    ";Variation;Events", 1000, 0., 1000.,11,0,11); 
  histContainer_["pu"]         = fs->make<TH1F>("pu",      ";Pileup observed;Events / 1",100,0,100);
  histContainer_["putrue"]     = fs->make<TH1F>("putrue",  ";Pileup true;Events / 0.1",100,0,100);
  for(std::unordered_map<std::string,TH1*>::iterator it=histContainer_.begin();   it!=histContainer_.end();   it++) it->second->Sumw2();

  //create a tree for the selected events
  if(saveTree_)
    {
      tree_ = fs->make<TTree>("data","data");
      createMiniEventTree(tree_,ev_);
    }
}


//
MiniAnalyzer::~MiniAnalyzer()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//
int MiniAnalyzer::genAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //
  // PILEUP
  //
  edm::Handle<std::vector <PileupSummaryInfo> > PupInfo;
  iEvent.getByToken(puToken_,PupInfo);
  std::vector<PileupSummaryInfo>::const_iterator ipu;
  for (ipu = PupInfo->begin(); ipu != PupInfo->end(); ++ipu) 
    {
      if ( ipu->getBunchCrossing() != 0 ) continue; // storing detailed PU info only for BX=0
      ev_.g_pu=ipu->getPU_NumInteractions();
      ev_.g_putrue=ipu->getTrueNumInteractions();
    }
  histContainer_["pu"]->Fill(ev_.g_pu);
  histContainer_["putrue"]->Fill(ev_.g_putrue);
  
  //
  // GENERATOR WEIGHTS
  //
  ev_.g_nw=0; ev_.g_w[0]=1.0;
  edm::Handle<GenEventInfoProduct> evt;
  iEvent.getByToken( generatorToken_,evt);
  if(evt.isValid())
    {
      ev_.g_w[0] = evt->weight();
      ev_.g_nw++;

      //PDF info
      ev_.g_qscale = evt->pdf()->scalePDF;
      ev_.g_x1     = evt->pdf()->x.first;
      ev_.g_x2     = evt->pdf()->x.second;
      ev_.g_id1    = evt->pdf()->id.first;
      ev_.g_id2    = evt->pdf()->id.second;
    }
  histContainer_["counter"]->Fill(1,ev_.g_w[0]);
  
  //alternative weights for systematics 
  edm::Handle<LHEEventProduct> evet;
  iEvent.getByToken(generatorlheToken_, evet);
  if(evet.isValid())
    {
      double asdd=evet->originalXWGTUP();
      for(unsigned int i=0  ; i<evet->weights().size();i++){
	double asdde=evet->weights()[i].wgt;
	ev_.g_w[ev_.g_nw]=ev_.g_w[0]*asdde/asdd;
	ev_.g_nw++;
      }
    }
     
  //
  // GENERATOR LEVEL EVENT
  //
  ev_.ng=0; 
  edm::Handle<std::vector<reco::GenJet> > genJets;
  iEvent.getByToken(genJetsToken_,genJets);  
  std::map<const reco::Candidate *,int> jetConstsMap;
  //edm::Handle<edm::ValueMap<float> > petersonFrag;
  //iEvent.getByToken(petersonFragToken_,petersonFrag);
  int ngjets(0),ngbjets(0);
  for(auto genJet=genJets->begin(); genJet!=genJets->end(); ++genJet)
    {

      edm::Ref<std::vector<reco::GenJet> > genJetRef(genJets,genJet-genJets->begin());

      //map the gen particles which are clustered in this jet
      JetFragInfo_t jinfo=analyzeJet(*genJet);

      std::vector< const reco::Candidate * > jconst=genJet->getJetConstituentsQuick();
      for(size_t ijc=0; ijc <jconst.size(); ijc++) jetConstsMap[ jconst[ijc] ] = ev_.ng;

      ev_.g_tagCtrs[ev_.ng]       = (jinfo.nbtags&0xf) | ((jinfo.nctags&0xf)<<4) | ((jinfo.ntautags&0xf)<<8);
      ev_.g_xb[ev_.ng]            = jinfo.xb;
      ev_.g_bid[ev_.ng]           = jinfo.leadTagId;
      ev_.g_isSemiLepBhad[ev_.ng] = jinfo.hasSemiLepDecay;
      //cout << "xb=" << jinfo.xb << " petersonFrag=" << (*petersonFrag)[genJetRef] << endl;
      ev_.g_id[ev_.ng]   = genJet->pdgId();
      ev_.g_pt[ev_.ng]   = genJet->pt();
      ev_.g_eta[ev_.ng]  = genJet->eta();
      ev_.g_phi[ev_.ng]  = genJet->phi();
      ev_.g_m[ev_.ng]    = genJet->mass();       
      ev_.ng++;
      
      //gen level selection
      if(genJet->pt()>25 && fabs(genJet->eta())<2.5)
	{
	  ngjets++;	
	  if(abs(genJet->pdgId())==5) ngbjets++;
	}
    }

  //leptons
  int ngleptons(0);
  edm::Handle<std::vector<reco::GenJet> > dressedLeptons;  
  iEvent.getByToken(genLeptonsToken_,dressedLeptons);
  for(auto genLep = dressedLeptons->begin();  genLep != dressedLeptons->end(); ++genLep)
    {
      //map the gen particles which are clustered in this lepton
      std::vector< const reco::Candidate * > jconst=genLep->getJetConstituentsQuick();
      for(size_t ijc=0; ijc <jconst.size(); ijc++) jetConstsMap[ jconst[ijc] ] = ev_.ng;
      
      ev_.g_pt[ev_.ng]   = genLep->pt();
      ev_.g_id[ev_.ng]   = genLep->pdgId();
      ev_.g_eta[ev_.ng]  = genLep->eta();
      ev_.g_phi[ev_.ng]  = genLep->phi();
      ev_.g_m[ev_.ng]    = genLep->mass();       
      ev_.ng++;

      //gen level selection
      if(genLep->pt()>20 && fabs(genLep->eta())<2.5) ngleptons++;
    }
  
  
  //final state particles 
  ev_.ngpf=0;
  edm::Handle<pat::PackedGenParticleCollection> genParticles;
  iEvent.getByToken(genParticlesToken_,genParticles);
  for (size_t i = 0; i < genParticles->size(); ++i)
    {
      const pat::PackedGenParticle & genIt = (*genParticles)[i];

      //this shouldn't be needed according to the workbook
      //if(genIt.status()!=1) continue;
      if(genIt.pt()<0.5 || fabs(genIt.eta())>2.5) continue;
      
      ev_.gpf_id[ev_.ngpf]     = genIt.pdgId();
      ev_.gpf_c[ev_.ngpf]      = genIt.charge();
      ev_.gpf_g[ev_.ngpf]=-1;
      for(std::map<const reco::Candidate *,int>::iterator it=jetConstsMap.begin();
	  it!=jetConstsMap.end();
	  it++)
	{
	  if(it->first->pdgId()!=genIt.pdgId()) continue;
	  if(deltaR( *(it->first), genIt)>0.01) continue; 
	  ev_.gpf_g[ev_.ngpf]=it->second;
	  break;
	}
      ev_.gpf_pt[ev_.ngpf]     = genIt.pt();
      ev_.gpf_eta[ev_.ngpf]    = genIt.eta();
      ev_.gpf_phi[ev_.ngpf]    = genIt.phi();
      ev_.gpf_m[ev_.ngpf]      = genIt.mass();
      ev_.ngpf++;    
    }


 //Bhadrons and top quarks (lastCopy)
  edm::Handle<reco::GenParticleCollection> prunedGenParticles;
  iEvent.getByToken(prunedGenParticlesToken_,prunedGenParticles);
  ev_.ngtop=0; 
  for (size_t i = 0; i < prunedGenParticles->size(); ++i)
    {
      const reco::GenParticle & genIt = (*prunedGenParticles)[i];
      int absid=abs(genIt.pdgId());
  
      //top quarks
      if(absid==6 && genIt.isLastCopy()) 
	{
	  ev_.gtop_id[ ev_.ngtop ]  = genIt.pdgId();
	  ev_.gtop_pt[ ev_.ngtop ]  = genIt.pt();
	  ev_.gtop_eta[ ev_.ngtop ] = genIt.eta();
	  ev_.gtop_phi[ ev_.ngtop ] = genIt.phi();
	  ev_.gtop_m[ ev_.ngtop ]   = genIt.mass();
	  ev_.ngtop++;
	}
    }
  
  //pseudo-tops 
  edm::Handle<reco::GenParticleCollection> pseudoTop;
  iEvent.getByToken(pseudoTopToken_,pseudoTop);
  for (size_t i = 0; i < pseudoTop->size(); ++i)
    {
      const GenParticle & genIt = (*pseudoTop)[i];
      ev_.gtop_id[ ev_.ngtop ]  = genIt.pdgId()*1000;
      ev_.gtop_pt[ ev_.ngtop ]  = genIt.pt();
      ev_.gtop_eta[ ev_.ngtop ] = genIt.eta();
      ev_.gtop_phi[ ev_.ngtop ] = genIt.phi();
      ev_.gtop_m[ ev_.ngtop ]   = genIt.mass();
      ev_.ngtop++;
    }

  //gen met
  edm::Handle<reco::METCollection> genMet;
  iEvent.getByToken(genMetsToken_,genMet);
  ev_.gtop_id[ ev_.ngtop ]  = 0;
  ev_.gtop_pt[ ev_.ngtop ]  = (*genMet)[0].pt();
  ev_.gtop_eta[ ev_.ngtop ] = 0;
  ev_.gtop_phi[ ev_.ngtop ] = (*genMet)[0].phi();
  ev_.gtop_m[ ev_.ngtop ]   = 0;
  ev_.ngtop++;

  //fiducial counters
  for(Int_t iw=0; iw<ev_.g_nw; iw++)
    {
      Double_t x(iw);
      Double_t wgt(ev_.g_w[iw]);
      TH2F *fidCounter=(TH2F *)histContainer_["fidcounter"];
      fidCounter->Fill(x,0.,wgt);
      if(ngleptons>0)               fidCounter->Fill(x, 1., wgt);
      if(ngleptons>1)               fidCounter->Fill(x, 2., wgt);
      if(ngleptons>0 && ngjets>0)   fidCounter->Fill(x, 3., wgt);
      if(ngleptons>1 && ngjets>0)   fidCounter->Fill(x, 4., wgt);
      if(ngleptons>0 && ngjets>1)   fidCounter->Fill(x, 5., wgt);
      if(ngleptons>1 && ngjets>1)   fidCounter->Fill(x, 6., wgt);
      if(ngleptons>0 && ngjets>2)   fidCounter->Fill(x, 7., wgt);
      if(ngleptons>1 && ngjets>2)   fidCounter->Fill(x, 8., wgt);
      if(ngleptons>0 && ngjets>3)   fidCounter->Fill(x, 9., wgt);
      if(ngleptons>1 && ngjets>3)   fidCounter->Fill(x, 10.,wgt);
    }
  
  //return number of leptons in fiducial range
  return ngleptons;
}


//
int MiniAnalyzer::recAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //VERTICES
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  if (vertices->empty()) return 0; // skip the event if no PV found
  const reco::Vertex &primVtx = vertices->front();
  reco::VertexRef primVtxRef(vertices,0);
   ev_.nvtx=vertices->size();
  if(ev_.nvtx==0) return 0;

  //RHO
  edm::Handle< double > rhoH;
  iEvent.getByToken(rhoToken_,rhoH);
  float rho=*rhoH;
  ev_.rho=rho;

  //TRIGGER INFORMATION
  edm::Handle<edm::TriggerResults> h_trigRes;
  iEvent.getByToken(triggerBits_, h_trigRes);
  std::vector<string> triggerList;
  Service<service::TriggerNamesService> tns;
  tns->getTrigPaths(*h_trigRes,triggerList);
  ev_.triggerBits=0;
  for (unsigned int i=0; i< h_trigRes->size(); i++) 
    {	
      if( !(*h_trigRes)[i].accept() ) continue;
      for(size_t itrig=0; itrig<triggersToUse_.size(); itrig++)
	{
	  if (triggerList[i].find(triggersToUse_[itrig])==string::npos) continue;
	  ev_.triggerBits |= (1 << itrig);
	  histContainer_["triggerList"]->Fill(itrig);
	}
    }
  bool passTrigger(ev_.isData ? ev_.triggerBits!=0 : true);
  if(!passTrigger) return 0;

  //PF candidates
  edm::Handle<pat::PackedCandidateCollection> pfcands;
  iEvent.getByToken(pfToken_,pfcands);

  //
  //LEPTON SELECTION
  //
  int nrecleptons(0);
  ev_.nl=0; 

  //MUON SELECTION: cf. https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2
  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);
  float maxMuPtForCor(200.);
  if(!useRawLeptons_ && muonCor_==0) muonCor_=new KalmanMuonCalibrator( ev_.isData? "DATA_80X_13TeV" : "MC_80X_13TeV");
  for (const pat::Muon &mu : *muons) 
    { 

      //correct the 4-momentum
      TLorentzVector p4;

      //apply correction
      float pt=mu.pt();
      if(pt<2) continue; //no need to care about very low pt muons here... (corrections will tend to be meaningless)

      float eta = mu.eta();
      float phi = mu.phi();
      float ptUnc = mu.bestTrack()->ptError();
      if(muonCor_ && mu.muonBestTrackType() == 1 && pt < maxMuPtForCor)
        {
	  bool corrApplied(true);
          if( !ev_.isData )
	    {
	      pt = muonCor_->getCorrectedPt(pt, eta, phi, mu.charge());
	      pt = muonCor_->smear(pt, eta);
	    }
	  else if (pt > 2. && fabs(eta) < 2.4)
	    {
	      pt = muonCor_->getCorrectedPt(pt, eta, phi, mu.charge());
	    }
	  else
	    {
	      corrApplied=false;
	    }
          
          if(pt<0.1 || isnan(pt))
            {
              cout << "After correction this muon pt is " << pt << " !!" << endl
                   << "Before correction this muon pt was " << mu.pt() << " !!" << endl
                   << "...will be skipped" << endl;
              continue;
            }
          
          if(corrApplied)
	    {
	      int n=muonCor_->getN();
	      float uncUp(0),uncDn(0);
	      for(int i=0; i<n; i++)
		{
		  muonCor_->vary(i,+1);
		  uncUp += pow(muonCor_->getCorrectedPt(mu.pt(),eta,phi,mu.charge())-pt,2);
		  muonCor_->vary(i,-1);
		  uncDn += pow(muonCor_->getCorrectedPt(mu.pt(),eta,phi,mu.charge())-pt,2);
		}
	      muonCor_->reset();
	      muonCor_->varyClosure(+1);
	      float vClose=muonCor_->getCorrectedPt(mu.pt(),eta,phi,mu.charge());
	      uncUp += pow(vClose-pt,2);
	      uncDn += pow(vClose-pt,2);
	      ptUnc = 0.5*(TMath::Sqrt(uncUp)+TMath::Sqrt(uncDn));
	    }
        }
      p4.SetPtEtaPhiM(pt,eta,phi,mu.mass());

      //kinematics
      bool passPt(p4.Pt() > 10);
      bool passEta(fabs(p4.Eta()) < 2.5);
      if(!passPt || !passEta) continue;

      //ID
      bool isSoft(isSoftMuon(mu,primVtx));
      bool isLoose(muon::isLooseMuon(mu));
      bool isMedium(muon::isMediumMuon(mu));
      bool isMedium2016ReReco(isMediumMuon2016ReReco(mu));
      bool isTight(muon::isTightMuon(mu,primVtx));
      if(!isSoft && !isLoose) continue;

      //save info
      const reco::GenParticle * gen=mu.genLepton(); 
      ev_.l_isPromptFinalState[ev_.nl] = gen ? gen->isPromptFinalState() : false;
      ev_.l_isDirectPromptTauDecayProductFinalState[ev_.nl] = gen ? gen->isDirectPromptTauDecayProductFinalState() : false;
      ev_.l_id[ev_.nl]=13;
      ev_.l_g[ev_.nl]=-1;
      for(int ig=0; ig<ev_.ng; ig++)
	{
	  if(abs(ev_.g_id[ig])!=ev_.l_id[ev_.nl]) continue;
	  if(deltaR( mu.eta(),mu.phi(), ev_.g_eta[ig],ev_.g_phi[ig])>0.4) continue;
	  ev_.l_g[ev_.nl]=ig;
	  break;
	}	 
      ev_.l_charge[ev_.nl]   = mu.charge();
      ev_.l_pt[ev_.nl]       = p4.Pt();
      ev_.l_eta[ev_.nl]      = p4.Eta();
      ev_.l_phi[ev_.nl]      = p4.Phi();
      ev_.l_mass[ev_.nl]     = p4.M();
      ev_.l_scaleUnc[ev_.nl] = ptUnc;
      ev_.l_mva[ev_.nl]      = 0;
      ev_.l_pid[ev_.nl]      = (isSoft | (isLoose<<1) | (isMedium<<2) | (isMedium2016ReReco<<3) | (isTight<<4));
      ev_.l_chargedHadronIso[ev_.nl] = mu.pfIsolationR04().sumChargedHadronPt;
      ev_.l_miniIso[ev_.nl]  = getMiniIsolation(pfcands,&mu,0.05,0.2, 10., false);
      ev_.l_relIso[ev_.nl]   = (
				mu.pfIsolationR04().sumChargedHadronPt 
				+ max(0., mu.pfIsolationR04().sumNeutralHadronEt + mu.pfIsolationR04().sumPhotonEt - 0.5*mu.pfIsolationR04().sumPUPt)
				) / p4.Pt();
      ev_.l_ip3d[ev_.nl]    = -9999.;
      ev_.l_ip3dsig[ev_.nl] = -9999;
      if(mu.innerTrack().get())
	{
	  std::pair<bool,Measurement1D> ip3dRes = getImpactParameter<reco::TrackRef>(mu.innerTrack(), primVtxRef, iSetup, true);
	  ev_.l_ip3d[ev_.nl]    = ip3dRes.second.value();
	  ev_.l_ip3dsig[ev_.nl] = ip3dRes.second.significance();
	}   
      ev_.nl++;    

      if( p4.Pt()>20 && fabs(p4.Eta())<2.5 && isLoose) nrecleptons++;
    }
  
  // ELECTRON SELECTION: cf. https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
  edm::Handle<edm::View<pat::Electron> > electrons;
  iEvent.getByToken(electronToken_, electrons);    
  edm::Handle<edm::ValueMap<bool> > veto_id, loose_id, medium_id, tight_id;
  iEvent.getByToken(eleVetoIdMapToken_,   veto_id);
  iEvent.getByToken(eleLooseIdMapToken_,  loose_id);
  iEvent.getByToken(eleMediumIdMapToken_, medium_id);
  iEvent.getByToken(eleTightIdMapToken_,  tight_id);
  edm::Handle<edm::ValueMap<vid::CutFlowResult> > veto_cuts, loose_cuts, medium_cuts, tight_cuts;
  iEvent.getByToken(eleVetoIdFullInfoMapToken_,   veto_cuts);
  iEvent.getByToken(eleLooseIdFullInfoMapToken_,  loose_cuts);
  iEvent.getByToken(eleMediumIdFullInfoMapToken_, medium_cuts);
  iEvent.getByToken(eleTightIdFullInfoMapToken_,  tight_cuts);
  edm::Handle<edm::ValueMap<float> > emva_id;
  iEvent.getByToken(eleMvaIdMapToken_, emva_id);
  Int_t nele(0);
  for (const pat::Electron &el : *electrons) 
    {        
      const auto e = electrons->ptrAt(nele); 
      nele++;

      //kinematics cuts
      bool passPt(e->pt() > 15.0);
      bool passEta(fabs(e->eta()) < 2.5 && (fabs(e->superCluster()->eta()) < 1.4442 || fabs(e->superCluster()->eta()) > 1.5660));
      if(!passPt || !passEta) continue;
      
      //take out id information alone
      bool passVetoId(true), passLooseId(true), passMediumId(true), passTightId(true);
      vid::CutFlowResult vetoCutBits   = (*veto_cuts)[e];
      for(size_t icut = 0; icut<vetoCutBits.cutFlowSize(); icut++)  
	{ 
	  if(evetoIsoBit_==999)
	    {
	      const std::string &cutName=vetoCutBits.getNameAtIndex(icut);
	      if(cutName.find("PFIso")!=std::string::npos)
		{
		  cout << "e veto isolation bit set to:" << icut << " from " << cutName << endl;
		  evetoIsoBit_=icut;
		}
	    }
	  if(icut!=evetoIsoBit_ && !vetoCutBits.getCutResultByIndex(icut)) passVetoId=false; 
	}
      vid::CutFlowResult looseCutBits  = (*loose_cuts)[e];
      for(size_t icut = 0; icut<looseCutBits.cutFlowSize(); icut++)  
	{
	  if(elooseIsoBit_==999)
            {
              const std::string &cutName=looseCutBits.getNameAtIndex(icut);
              if(cutName.find("PFIso")!=std::string::npos)
                {
		  cout <<"e loose isolation bit set to:" << icut << " from " << cutName << endl;
                  elooseIsoBit_=icut;
                }
            }
	  if(icut!=elooseIsoBit_ && !looseCutBits.getCutResultByIndex(icut)) passLooseId=false; 
	}
      vid::CutFlowResult mediumCutBits = (*medium_cuts)[e];
      for(size_t icut = 0; icut<mediumCutBits.cutFlowSize(); icut++) 
	{ 
	  if(emediumIsoBit_==999)
            {
              const std::string &cutName=mediumCutBits.getNameAtIndex(icut);
              if(cutName.find("PFIso")!=std::string::npos)
                {
                  cout <<"e medium isolation bit set to:" << icut << " from " << cutName << endl;
                  emediumIsoBit_=icut;
                }
            }
	  if(icut!=emediumIsoBit_ && !mediumCutBits.getCutResultByIndex(icut)) passMediumId=false; 
	}
      vid::CutFlowResult tightCutBits  = (*tight_cuts)[e];
      for(size_t icut = 0; icut<tightCutBits.cutFlowSize(); icut++) 
	{
	  if(etightIsoBit_==999)
            {
              const std::string &cutName=tightCutBits.getNameAtIndex(icut);
              if(cutName.find("PFIso")!=std::string::npos)
                {
                  cout <<"e tight isolation bit set to:" << icut << " from " << cutName << endl;
                  etightIsoBit_=icut;
                }
            }
	  if(icut!=etightIsoBit_ && !tightCutBits.getCutResultByIndex(icut)) passTightId=false;
	}
      if(!passVetoId) continue;

      //for(size_t icut = 0; icut<vetoCutBits.cutFlowSize(); icut++)
      //	  std::cout << vetoCutBits.getNameAtIndex (icut) << std::endl;
	  
      //full id+iso decisions
      bool isVeto( (*veto_id)[e] );
      bool isLoose( (*loose_id)[e] );      
      bool isMedium( (*medium_id)[e] );
      bool isTight( (*tight_id)[e] );

      //impact parameter cuts
      //see details in https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
      bool passIpCuts(true);
      if(el.gsfTrack().isNonnull())
	{
	  float dxy(fabs(el.gsfTrack()->dxy(primVtx.position())));
	  float dz(fabs(el.gsfTrack()->dz(primVtx.position())));
	  if(fabs(e->superCluster()->eta()) < 1.4442)
	    {
	      if(dxy>0.05 || dz>0.10) passIpCuts=false;
	    }
	  else
	    {
	      if(dxy>0.10 || dz>0.20) passIpCuts=false;
	    }	  
	}
      else
	{
	  passIpCuts=false;
	}
	  
      //cout << e->pt() << " " << e->eta() << " " << isVeto << " " << isLoose << " " << isMedium << " " << isTight << " " << passIpCuts << endl;

      //save the electron
      const reco::GenParticle * gen=el.genLepton(); 
      ev_.l_isPromptFinalState[ev_.nl] = gen ? gen->isPromptFinalState() : false;
      ev_.l_isDirectPromptTauDecayProductFinalState[ev_.nl] = gen ? gen->isDirectPromptTauDecayProductFinalState() : false;
      ev_.l_id[ev_.nl]=11;
      ev_.l_g[ev_.nl]=-1;
      for(int ig=0; ig<ev_.ng; ig++)
	{
	  if(abs(ev_.g_id[ig])!=ev_.l_id[ev_.nl]) continue;
	  if(deltaR( el.eta(),el.phi(), ev_.g_eta[ig],ev_.g_phi[ig])>0.4) continue;
	  ev_.l_g[ev_.nl]=ig;
	  break;
	}	      
      ev_.l_mva[ev_.nl]=(*emva_id)[e];

      ev_.l_pid[ev_.nl]=0;
      ev_.l_pid[ev_.nl]= (passVetoId | (isVeto<<1) 
			  | (passLooseId<<2) | (isLoose<<3) 
			  | (passMediumId<<4) | (isMedium<<5) 
			  | (passTightId<<6) | (isTight<<7)
			  | (passIpCuts<<8)
			 );
      ev_.l_charge[ev_.nl]   = el.charge();
      ev_.l_pt[ev_.nl]       = e->pt();
      ev_.l_eta[ev_.nl]      = e->eta();
      ev_.l_phi[ev_.nl]      = e->phi();
      ev_.l_mass[ev_.nl]     = e->mass();
      ev_.l_scaleUnc[ev_.nl] = e->correctedEcalEnergyError();
      ev_.l_miniIso[ev_.nl]  = getMiniIsolation(pfcands,&el,0.05, 0.2, 10., false);
      ev_.l_relIso[ev_.nl]   = (e->chargedHadronIso()+ max(0., e->neutralHadronIso() + e->photonIso()  - 0.5*e->puChargedHadronIso()))/e->pt();     
      ev_.l_chargedHadronIso[ev_.nl] = e->chargedHadronIso();
      ev_.l_ip3d[ev_.nl]     = -9999.;
      ev_.l_ip3dsig[ev_.nl]  = -9999;
      if(el.gsfTrack().get())
	{
	  std::pair<bool,Measurement1D> ip3dRes = getImpactParameter<reco::GsfTrackRef>(el.gsfTrack(), primVtxRef, iSetup, true);
	  ev_.l_ip3d[ev_.nl]    = ip3dRes.second.value();
	  ev_.l_ip3dsig[ev_.nl] = ip3dRes.second.significance();
	}
      ev_.nl++;
      
      if( e->pt()>20 && passEta && passLooseId ) nrecleptons++;
    }

  // JETS
  ev_.nj=0; 
  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByToken(jetToken_,jets);
  std::vector< std::pair<const reco::Candidate *,int> > clustCands;
  for(auto j = jets->begin();  j != jets->end(); ++j)
    {
      //kinematics
      if(j->pt()<20 || fabs(j->eta())>4.7) continue;
      
      // PF jet ID
      pat::strbitset retpf = pfjetIDLoose_.getBitTemplate();
      retpf.set(false);
      bool passLoose=pfjetIDLoose_( *j, retpf );
      if(!passLoose) continue;
     
      //save jet
      const reco::Candidate *genParton = j->genParton();
      ev_.j_area[ev_.nj]  = j->jetArea();
      ev_.j_rawsf[ev_.nj] = j->correctedJet("Uncorrected").pt()/j->pt();
      ev_.j_pt[ev_.nj]    = j->pt();
      ev_.j_mass[ev_.nj]  = j->mass();
      ev_.j_eta[ev_.nj]   = j->eta();
      ev_.j_phi[ev_.nj]   = j->phi();
      ev_.j_g[ev_.nj]     = -1;
      for(int ig=0; ig<ev_.ng; ig++)
	{
	  if(abs(ev_.g_id[ig])==11 || abs(ev_.g_id[ig])==13) continue;
	  if(deltaR( j->eta(),j->phi(), ev_.g_eta[ig],ev_.g_phi[ig])>0.4) continue;
	  ev_.j_g[ev_.nj]=ig;
	  ev_.g_xbp[ig]  = genParton   ? ev_.g_xb[ig]*ev_.g_pt[ig]/genParton->pt() : 0.;
	  break;
	}	 
      ev_.j_csv[ev_.nj]=j->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
      ev_.j_btag[ev_.nj]       = (ev_.j_csv[ev_.nj]>0.8484);
      //ev_.j_deepcsvl[ev_.nj]   = j->bDiscriminator("deepFlavourJetTags:probudsg");
      //ev_.j_deepcsvc[ev_.nj]   = j->bDiscriminator("deepFlavourJetTags:probc")+j->bDiscriminator("deepFlavourJetTags:probcc");
      //ev_.j_deepcsvb[ev_.nj]   = j->bDiscriminator("deepFlavourJetTags:probb")+j->bDiscriminator("deepFlavourJetTags:probbb");

      if( j->hasTagInfo("pfInclusiveSecondaryVertexFinder") )
	{
	  const reco::CandSecondaryVertexTagInfo *candSVTagInfo = j->tagInfoCandSecondaryVertex("pfInclusiveSecondaryVertexFinder");
	  if( candSVTagInfo->nVertices() >= 1 ) 
	    {
	      math::XYZTLorentzVectorD vp4 = candSVTagInfo->secondaryVertex(0).p4();
	      ev_.j_vtxpx[ev_.nj]          = vp4.px();
	      ev_.j_vtxpy[ev_.nj]          = vp4.py();
	      ev_.j_vtxpz[ev_.nj]          = vp4.pz();
	      ev_.j_vtxmass[ev_.nj]        = vp4.mass();
	      ev_.j_vtxNtracks[ev_.nj]     = candSVTagInfo->nVertexTracks(0);
	      ev_.j_vtx3DVal[ev_.nj]       = candSVTagInfo->flightDistance(0).value();
	      ev_.j_vtx3DSig[ev_.nj]       = candSVTagInfo->flightDistance(0).significance();
	    }
	}

      ev_.j_flav[ev_.nj]       = j->partonFlavour();
      ev_.j_hadflav[ev_.nj]    = j->hadronFlavour();
      ev_.j_pid[ev_.nj]        = genParton ? genParton->pdgId() : 0;
      ev_.nj++;

      //save all PF candidates central jet
      if(fabs(j->eta())>2.5) continue;
      for(size_t ipf=0; ipf<j->numberOfDaughters(); ipf++)
	{
	  const reco::Candidate *pf=j->daughter(ipf);
	  clustCands.push_back(std::pair<const reco::Candidate *,int>(pf,ev_.nj-1));
	}
    }
      
  // MET
  ev_.nmet=2;
  for(int i=0; i<2; i++)
    {
      edm::Handle<pat::METCollection> mets;
      if(i==0) iEvent.getByToken(metToken_, mets);
      if(i==1) iEvent.getByToken(puppiMetToken_, mets);
      ev_.met_pt[i]  = mets->at(0).pt();
      ev_.met_phi[i] = mets->at(0).phi();
      ev_.met_sig[i] = mets->at(0).significance();
    }

  //MET filter bits
  ev_.met_filterBits=0;
  edm::Handle<edm::TriggerResults> h_metFilters;
  iEvent.getByToken(metFilterBits_, h_metFilters);
  std::vector<string> metFilterNames;
  Service<service::TriggerNamesService> mfns;
  mfns->getTrigPaths(*h_metFilters,metFilterNames);
  for (unsigned int i=0; i< h_metFilters->size(); i++) 
    {	
      if( !(*h_metFilters)[i].accept() ) continue;
      for(size_t itrig=0; itrig<metFiltersToUse_.size(); itrig++)
	{
	  if (metFilterNames[i].find(metFiltersToUse_[itrig])==string::npos) continue;
	  ev_.met_filterBits |= (1<<itrig);
	}
    }

  try{
    edm::Handle<bool> ifilterbadChCand;
    iEvent.getByToken(BadChCandFilterToken_, ifilterbadChCand);
    bool  filterbadChCandidate = *ifilterbadChCand;
    ev_.met_filterBits |= (filterbadChCandidate<<metFiltersToUse_.size());
  }
  catch(...){
  }
  
  try{
    edm::Handle<bool> ifilterbadPFMuon;
    iEvent.getByToken(BadPFMuonFilterToken_, ifilterbadPFMuon);
    bool filterbadPFMuon = *ifilterbadPFMuon;
    ev_.met_filterBits |= (filterbadPFMuon<<(metFiltersToUse_.size()+1));
  }
  catch(...){
  }

  //PF candidates
  ev_.npf=0;
  for(auto pf = pfcands->begin();  pf != pfcands->end(); ++pf)
    {
      if(ev_.npf>=5000) continue;

      ev_.pf_j[ev_.npf] = -1;
      for(size_t i=0; i<clustCands.size(); i++)
	{
	  if(pf->pdgId()!=clustCands[i].first->pdgId()) continue;
	  if(deltaR(*pf,*(clustCands[i].first))>0.01) continue;
	  ev_.pf_j[ev_.npf]=clustCands[i].second;
	  break;
	}

      //if particle is not associated to jet and is neutral, discard
      if(pf->charge()==0 && ev_.pf_j[ev_.npf]==-1) continue;

      //if particle is charged require association to prim vertex
      if(pf->charge()!=0)
	{
	  //if(pf->pvAssociationQuality()<pat::PackedCandidate::CompatibilityDz) continue;
	  //if(pf->vertexRef().key()!=0) continue;

	  if(pf->pt()<0.9 || fabs(pf->eta())>2.5) continue;
	  const pat::PackedCandidate::PVAssoc pvassoc=pf->fromPV();
	  if(pvassoc< pat::PackedCandidate::PVTight) continue;
	  
	  ev_.pf_dxy[ev_.npf]      = pf->dxy();
	  ev_.pf_dz[ev_.npf]       = pf->dz();
	  //ev_.pf_dxyUnc[ev_.npf]   = pf->dxyError();
	  //ev_.pf_dzUnc[ev_.npf]    = pf->dzError();
	  //ev_.pf_vtxRef[ev_.npf]   = pf->vertexRef().key();
	  //ev_.pf_pvAssoc[ev_.npf]  = pf->fromPV() + 10*(pf->pvAssociationQuality());
	}
      
      ev_.pf_id[ev_.npf]       = pf->pdgId();
      ev_.pf_c[ev_.npf]        = pf->charge();
      ev_.pf_pt[ev_.npf]       = pf->pt();
      ev_.pf_eta[ev_.npf]      = pf->eta();
      ev_.pf_phi[ev_.npf]      = pf->phi();
      ev_.pf_m[ev_.npf]        = pf->mass();
      ev_.pf_puppiWgt[ev_.npf] = pf->puppiWeight();      
      ev_.npf++;
    }

  return nrecleptons;
}

//cf. https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2#Soft_Muon
bool MiniAnalyzer::isSoftMuon(const reco::Muon & recoMu,const reco::Vertex &vertex)
{

  bool isGood(muon::isGoodMuon(recoMu, muon::TMOneStationTight));
  bool passLayersWithMeas(recoMu.innerTrack().isNonnull()
			  && recoMu.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5
			  && recoMu.innerTrack()->hitPattern().pixelLayersWithMeasurement() > 0 );
  bool matchesVertex(recoMu.innerTrack().isNonnull()
		     && fabs(recoMu.innerTrack()->dxy(vertex.position())) < 0.3 
		     && fabs(recoMu.innerTrack()->dz(vertex.position())) < 20. );
  return (isGood && passLayersWithMeas && matchesVertex);
}

//cf. https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2#Standard_MediumID_to_be_used_wit
bool MiniAnalyzer::isMediumMuon2016ReReco(const reco::Muon & recoMu) 
{
  bool goodGlob = recoMu.isGlobalMuon() && 
    recoMu.globalTrack()->normalizedChi2() < 3 && 
    recoMu.combinedQuality().chi2LocalPosition < 12 && 
    recoMu.combinedQuality().trkKink < 20; 
  bool isMedium = muon::isLooseMuon(recoMu) && 
    recoMu.innerTrack()->validFraction() > 0.8 && 
    muon::segmentCompatibility(recoMu) > (goodGlob ? 0.303 : 0.451); 
  return isMedium; 
}



// ------------ method called for each event  ------------
void MiniAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  histContainer_["counter"]->Fill(0);

  //analyze the event
  int ngleptons(0),nrecleptons(0);
  if(!iEvent.isRealData()) ngleptons=genAnalysis(iEvent,iSetup);
  nrecleptons=recAnalysis(iEvent,iSetup);
  
  //save event if at least one object at gen or reco level
  if((ngleptons==0 && nrecleptons==0) || !saveTree_) return;  
  ev_.run     = iEvent.id().run();
  ev_.lumi    = iEvent.luminosityBlock();
  ev_.event   = iEvent.id().event(); 
  ev_.isData  = iEvent.isRealData();
  if(!savePF_) { ev_.ngpf=0; ev_.npf=0; }
  tree_->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
MiniAnalyzer::beginJob(){
}

//
void 
MiniAnalyzer::endRun(const edm::Run& iRun,
		     const EventSetup& iSetup) 
{
  try{

    cout << "[MiniAnalyzer::endRun]" << endl;

    edm::Handle<LHERunInfoProduct> lheruninfo;
    typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
    iRun.getByToken(generatorRunInfoToken_, lheruninfo );
    //iRun.getByLabel( "externalLHEProducer", lheruninfo );

    LHERunInfoProduct myLHERunInfoProduct = *(lheruninfo.product());
    for (headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); 
	 iter!=myLHERunInfoProduct.headers_end(); 
	 iter++)
      {
	std::string tag("generator");
	if(iter->tag()!="") tag+="_"+iter->tag();
	
	std::vector<std::string> lines = iter->lines();
	std::vector<std::string> prunedLines;
	for (unsigned int iLine = 0; iLine<lines.size(); iLine++) 
	  {
	    if(lines.at(iLine)=="") continue;
	    if(lines.at(iLine).find("weightgroup")!=std::string::npos) continue;
	    prunedLines.push_back( lines.at(iLine) );
	  }
	
	if(histContainer_.find(tag)==histContainer_.end()) 
	  {
	    std::cout << "Starting histo for " << tag << std::endl;
	    histContainer_[tag]=fs->make<TH1F>(tag.c_str(),tag.c_str(),prunedLines.size(),0,prunedLines.size());
	  }
	for (unsigned int iLine = 0; iLine<prunedLines.size(); iLine++) 
	  histContainer_[tag]->GetXaxis()->SetBinLabel(iLine+1,prunedLines.at(iLine).c_str());  
      }
  }
  catch(std::exception &e){
    std::cout << e.what() << endl
	      << "Failed to retrieve LHERunInfoProduct" << std::endl;
  }
}

//-------------
//cf. https://twiki.cern.ch/twiki/bin/view/CMS/MiniIsolationSUSY
float MiniAnalyzer::getMiniIsolation(edm::Handle<pat::PackedCandidateCollection> pfcands,
				     const reco::Candidate* ptcl,  
				     float r_iso_min, float r_iso_max, float kt_scale,
				     bool charged_only) 
{

    if (ptcl->pt()<5.) return 99999.;

    float deadcone_nh(0.), deadcone_ch(0.), deadcone_ph(0.), deadcone_pu(0.);
    if(ptcl->isElectron()) {
      if (fabs(ptcl->eta())>1.479) {deadcone_ch = 0.015; deadcone_pu = 0.015; deadcone_ph = 0.08;}
    } else if(ptcl->isMuon()) {
      deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01;  
    } else {
      //deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01; // maybe use muon cones??
    }

    float iso_nh(0.), iso_ch(0.), iso_ph(0.), iso_pu(0.);
    float ptThresh(0.5);
    if(ptcl->isElectron()) ptThresh = 0;
    float r_iso = (float)TMath::Max((float)r_iso_min,
				    (float)TMath::Min((float)r_iso_max, (float)(kt_scale/ptcl->pt())));
    for (const pat::PackedCandidate &pfc : *pfcands) {
      if (abs(pfc.pdgId())<7) continue;
      
      float dr = deltaR(pfc, *ptcl);
      if (dr > r_iso) continue;
      
      //////////////////  NEUTRALS  /////////////////////////
      if (pfc.charge()==0){
        if (pfc.pt()>ptThresh) {
          /////////// PHOTONS ////////////
          if (abs(pfc.pdgId())==22) {
            if(dr < deadcone_ph) continue;
            iso_ph += pfc.pt();
	    /////////// NEUTRAL HADRONS ////////////
          } else if (abs(pfc.pdgId())==130) {
            if(dr < deadcone_nh) continue;
            iso_nh += pfc.pt();
          }
        }
        //////////////////  CHARGED from PV  /////////////////////////
      } else if (pfc.fromPV()>1){
        if (abs(pfc.pdgId())==211) {
          if(dr < deadcone_ch) continue;
          iso_ch += pfc.pt();
        }
        //////////////////  CHARGED from PU  /////////////////////////
      } else {
        if (pfc.pt()>ptThresh){
          if(dr < deadcone_pu) continue;
          iso_pu += pfc.pt();
        }
      }
    }
    float iso(0.);
    if (charged_only){
      iso = iso_ch;
    } else {
      iso = iso_ph + iso_nh;
      iso -= 0.5*iso_pu;
      if (iso>0) iso += iso_ch;
      else iso = iso_ch;
    }
    iso = iso/ptcl->pt();

    return iso;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MiniAnalyzer::endJob() 
{
  std::cout << "[MiniAnalyzer::endJob]" << endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MiniAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniAnalyzer);
