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
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
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
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSDetId.h"
#include "DataFormats/ProtonReco/interface/ProtonTrack.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetResolutionObject.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TopLJets2015/TopAnalysis/interface/MyIPTools.h"
#include "TopLJets2015/TopAnalysis/interface/JetShapes.h"
#include "TopLJets2015/TopAnalysis/interface/RoccoR.h"
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
#include "FWCore/Framework/interface/Run.h"
#include "TopQuarkAnalysis/BFragmentationAnalyzer/interface/BFragmentationAnalyzerUtils.h"

#include "TLorentzVector.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TRandom2.h"
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

typedef math::XYZTLorentzVector LorentzVector;

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
  void beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const& iSetup);
  void genAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void recAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup);
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
  edm::EDGetTokenT<GenLumiInfoHeader> generatorLumiHeaderToken_;
  edm::EDGetTokenT<GenEventInfoProduct> generatorEventInfoToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > puToken_;
  edm::EDGetTokenT<std::vector<reco::GenParticle>  > genPhotonsToken_;
  edm::EDGetTokenT<std::vector<reco::GenJet>  > genLeptonsToken_, genJetsToken_;
  edm::EDGetTokenT<reco::METCollection> genMetsToken_;
  edm::EDGetTokenT<pat::PackedGenParticleCollection> genParticlesToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> prunedGenParticlesToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> particleLevelToken_;

  edm::EDGetTokenT<edm::TriggerResults> triggerBits_,metFilterBits_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_,l1triggerPrescales_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<double> rhoToken_;
  edm::EDGetTokenT<pat::MuonCollection> muonToken_;
  edm::EDGetTokenT<edm::View<pat::Electron>  >  electronToken_;
  edm::EDGetTokenT<edm::View<pat::Photon>  >  photonToken_;
  edm::EDGetTokenT<edm::View<pat::Jet> > jetToken_;
  edm::EDGetTokenT<pat::METCollection> metToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
  edm::EDGetTokenT<std::vector<CTPPSLocalTrackLite> > ctppsToken_;
  edm::EDGetTokenT<std::vector<reco::ProtonTrack>> tokenRecoProtons_;

  //
  edm::EDGetTokenT<bool> BadChCandFilterToken_,BadPFMuonFilterToken_;

  //  edm::EDGetTokenT<edm::ValueMap<float> > petersonFragToken_;

  std::unordered_map<std::string,TH1*> histContainer_;

  std::string jetIdToUse_;
  std::vector<JetCorrectionUncertainty *> jecCorrectionUncs_;

  std::vector<std::string> triggersToUse_,metFiltersToUse_;

  bool saveTree_,savePF_;
  TTree *tree_;
  MiniEvent_t ev_;
  
  RoccoR *muonRC_;

  edm::Service<TFileService> fs;

  //counters
  int nrecleptons_,nrecphotons_,ngleptons_,ngphotons_;

  //apply filter to save tree
  bool applyFilt_;

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
  generatorLumiHeaderToken_(consumes<GenLumiInfoHeader,edm::InLumi>(edm::InputTag("generator"))),
  generatorEventInfoToken_(consumes<GenEventInfoProduct>(edm::InputTag("generator"))),
  puToken_(consumes<std::vector<PileupSummaryInfo>>(edm::InputTag("slimmedAddPileupInfo"))),  
  genPhotonsToken_(consumes<std::vector<reco::GenParticle> >(edm::InputTag("particleLevel:photons"))),
  genLeptonsToken_(consumes<std::vector<reco::GenJet> >(edm::InputTag("particleLevel:leptons"))),
  genJetsToken_(consumes<std::vector<reco::GenJet> >(edm::InputTag("particleLevel:jets"))),
  genMetsToken_(consumes<reco::METCollection>(edm::InputTag("particleLevel:mets"))),
  genParticlesToken_(consumes<pat::PackedGenParticleCollection>(edm::InputTag("packedGenParticles"))),
  prunedGenParticlesToken_(consumes<reco::GenParticleCollection>(edm::InputTag("prunedGenParticles"))),
  particleLevelToken_(consumes<reco::GenParticleCollection>(edm::InputTag("particleLevel"))),
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerBits"))),
  metFilterBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("metFilterBits"))),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
  l1triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("l1prescales"))),
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
  muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  jetToken_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jets"))),  
  metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),  
  pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
  ctppsToken_(consumes<std::vector<CTPPSLocalTrackLite> >(iConfig.getParameter<edm::InputTag>("ctppsLocalTracks"))),
  tokenRecoProtons_(consumes<std::vector<reco::ProtonTrack>>(iConfig.getParameter<InputTag>("tagRecoProtons"))),
  BadChCandFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("badChCandFilter"))),
  BadPFMuonFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("badPFMuonFilter"))),
  saveTree_( iConfig.getParameter<bool>("saveTree") ),
  savePF_( iConfig.getParameter<bool>("savePF") ),
  applyFilt_( iConfig.getParameter<bool>("applyFilt") )
{
  //now do what ever initialization is needed
  electronToken_      = mayConsume<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electrons"));
  photonToken_        = mayConsume<edm::View<pat::Photon> >(iConfig.getParameter<edm::InputTag>("photons"));
  triggersToUse_      = iConfig.getParameter<std::vector<std::string> >("triggersToUse");
  metFiltersToUse_    = iConfig.getParameter<std::vector<std::string> >("metFiltersToUse");
  jetIdToUse_         = iConfig.getParameter<std::string>("jetIdToUse");  
  std::string jecUncFile(iConfig.getParameter<std::string>("jecUncFile"));
  //std::string jecUncFile(edm::FileInPath(iConfig.getParameter<std::string>("jecUncFile")).fullPath());
  for(auto name : iConfig.getParameter<std::vector<std::string> >("jecUncSources") ) {
    JetCorrectorParameters *p = new JetCorrectorParameters(jecUncFile,name.c_str());
    jecCorrectionUncs_.push_back(new JetCorrectionUncertainty(*p));
  }

  muonRC_ = new RoccoR();
  muonRC_->init(iConfig.getParameter<std::string>("RoccoR"));
  //muonRC_->init(edm::FileInPath(iConfig.getParameter<std::string>("RoccoR")).fullPath());

  histContainer_["triggerList"] = fs->make<TH1F>("triggerList", ";Trigger bits;",triggersToUse_.size(),0,triggersToUse_.size());
  histContainer_["triggerPrescale"] = fs->make<TH1D>("triggerPrescale", ";Trigger prescale sum;",triggersToUse_.size(),0,triggersToUse_.size());
  for(size_t i=0; i<triggersToUse_.size(); i++) histContainer_["triggerList"] ->GetXaxis()->SetBinLabel(i+1,triggersToUse_[i].c_str());
  histContainer_["counter"]    = fs->make<TH1F>("counter", ";Counter;Events",2,0,2);
  histContainer_["fidcounter"] = (TH1 *)fs->make<TH2F>("fidcounter",    ";Variation;Events", 1500, 0., 1500.,11,0,11); 
  histContainer_["pu"]         = fs->make<TH1F>("pu",      ";Pileup observed;Events / 1",100,0,100);
  histContainer_["putrue"]     = fs->make<TH1F>("putrue",  ";Pileup true;Events / 0.1",100,0,100);
  for(std::unordered_map<std::string,TH1*>::iterator it=histContainer_.begin();   it!=histContainer_.end();   it++) it->second->Sumw2();

  //create a tree for the selected events
  if(saveTree_)
    {
      tree_ = fs->make<TTree>("data","data");
      createMiniEventTree(tree_,ev_,jecCorrectionUncs_.size());
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
void MiniAnalyzer::genAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
  ev_.g_nw=0; ev_.g_npsw=0; ev_.g_w[0]=1.0;
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
      const auto & hepeup = evet->hepeup();
      ev_.g_nup=hepeup.NUP;
    }

  //parton shower weights
  edm::Handle<GenEventInfoProduct> genEventInfoProduct;
  iEvent.getByToken(generatorEventInfoToken_, genEventInfoProduct);      
  if(genEventInfoProduct.isValid()){    
    for(unsigned int i=0; i< genEventInfoProduct->weights().size(); i++){
      ev_.g_psw[ev_.g_npsw]=genEventInfoProduct->weights().at(i);
      ev_.g_npsw++;
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
  if(genJets.isValid()){
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
  }

  //leptons
  edm::Handle<std::vector<reco::GenJet> > dressedLeptons;  
  iEvent.getByToken(genLeptonsToken_,dressedLeptons);
  if(dressedLeptons.isValid()) {
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
        if(genLep->pt()>20 && fabs(genLep->eta())<2.5) ngleptons_++;
      }
  }
  
  edm::Handle<std::vector<reco::GenParticle> > genPhotons;
  iEvent.getByToken(genPhotonsToken_,genPhotons);
  if(genPhotons.isValid()){
    for(auto genPhoton = genPhotons->begin();  genPhoton != genPhotons->end(); ++genPhoton)
      {
        if(genPhoton->pt()<15) continue;
        if(fabs(genPhoton->eta())>2.5) continue;
        
        ev_.g_pt[ev_.ng]   = genPhoton->pt();
        ev_.g_id[ev_.ng]   = genPhoton->pdgId();
        ev_.g_eta[ev_.ng]  = genPhoton->eta();
        ev_.g_phi[ev_.ng]  = genPhoton->phi();
        ev_.g_m[ev_.ng]    = genPhoton->mass();       
        ev_.ng++;
        
        //gen level selection
        if(genPhoton->pt()>20 && fabs(genPhoton->eta())<2.5) ngphotons_++;
      }
  }
  
  //final state particles 
  ev_.g_nchPV=0;
  ev_.g_sumPVChPt=0; 
  ev_.g_sumPVChPz=0; 
  ev_.g_sumPVChHt=0; 
  edm::Handle<pat::PackedGenParticleCollection> genParticles;
  iEvent.getByToken(genParticlesToken_,genParticles);
  LorentzVector pvP4(0,0,0,0);
  if(genParticles.isValid()){
    for (size_t i = 0; i < genParticles->size(); ++i)
      {
        const pat::PackedGenParticle & genIt = (*genParticles)[i];
        if(genIt.pt()<0.5) continue;
        if(genIt.charge()==0) continue;
        ev_.g_nchPV++;
        pvP4+=genIt.p4();
        ev_.g_sumPVChPz+=fabs(genIt.pz()); 
        ev_.g_sumPVChHt+=genIt.pt(); 
      }
  }
  ev_.g_sumPVChPt=pvP4.Pt();

  //Bhadrons and top quarks (lastCopy)
  edm::Handle<reco::GenParticleCollection> prunedGenParticles;
  iEvent.getByToken(prunedGenParticlesToken_,prunedGenParticles);
  ev_.ngtop=0; 
  if(prunedGenParticles.isValid()){
    for (size_t i = 0; i < prunedGenParticles->size(); ++i)
      {
        const reco::GenParticle & genIt = (*prunedGenParticles)[i];
        int absid=abs(genIt.pdgId());
        bool outGoingProton( absid==2212 && genIt.status()==1 && fabs(genIt.eta())>4.7 );
        bool topLastCopy(absid==6 && genIt.isLastCopy());
        bool wLastCopy(absid==24 && genIt.isLastCopy());
        if(outGoingProton || topLastCopy)
          {
            ev_.gtop_id[ ev_.ngtop ]  = genIt.pdgId();
            ev_.gtop_pt[ ev_.ngtop ]  = genIt.pt();
            ev_.gtop_eta[ ev_.ngtop ] = genIt.eta();
            ev_.gtop_phi[ ev_.ngtop ] = genIt.phi();
            ev_.gtop_m[ ev_.ngtop ]   = genIt.mass();
            ev_.ngtop++;
          }

        //save top daughters
        if(topLastCopy || wLastCopy) {
          size_t ndau=genIt.numberOfDaughters();
          for(size_t idau=0; idau<ndau; idau++){
            const reco::Candidate *d=genIt.daughter(idau);
            ev_.gtop_id[ ev_.ngtop ]  = d->pdgId();
            ev_.gtop_pt[ ev_.ngtop ]  = d->pt();
            ev_.gtop_eta[ ev_.ngtop ] = d->eta();
            ev_.gtop_phi[ ev_.ngtop ] = d->phi();
            ev_.gtop_m[ ev_.ngtop ]   = d->mass();
            ev_.ngtop++;
          }                
        }
      }
  }
  
  //pseudo-tops 
/*  edm::Handle<reco::GenParticleCollection> particleLevel;
  iEvent.getByToken(particleLevelToken_,particleLevel);
  for (size_t i = 0; i < particleLevel->size(); ++i)
    {
      const GenParticle & genIt = (*particleLevel)[i];
      ev_.gtop_id[ ev_.ngtop ]  = genIt.pdgId()*1000;
      ev_.gtop_pt[ ev_.ngtop ]  = genIt.pt();
      ev_.gtop_eta[ ev_.ngtop ] = genIt.eta();
      ev_.gtop_phi[ ev_.ngtop ] = genIt.phi();
      ev_.gtop_m[ ev_.ngtop ]   = genIt.mass();
      ev_.ngtop++;
    }
*/
  //gen met
  edm::Handle<reco::METCollection> genMet;
  iEvent.getByToken(genMetsToken_,genMet);
  if(genMet.isValid()){
    ev_.gtop_id[ ev_.ngtop ]  = 0;
    ev_.gtop_pt[ ev_.ngtop ]  = (*genMet)[0].pt();
    ev_.gtop_eta[ ev_.ngtop ] = 0;
    ev_.gtop_phi[ ev_.ngtop ] = (*genMet)[0].phi();
    ev_.gtop_m[ ev_.ngtop ]   = 0;
    ev_.ngtop++;
  }

  //fiducial counters
  for(Int_t iw=0; iw<ev_.g_nw; iw++)
    {
      Double_t x(iw);
      Double_t wgt(ev_.g_w[iw]);
      TH2F *fidCounter=(TH2F *)histContainer_["fidcounter"];
      fidCounter->Fill(x,0.,wgt);
      if(ngleptons_>0)               fidCounter->Fill(x, 1., wgt);
      if(ngleptons_>1)               fidCounter->Fill(x, 2., wgt);
      if(ngleptons_>0 && ngjets>0)   fidCounter->Fill(x, 3., wgt);
      if(ngleptons_>1 && ngjets>0)   fidCounter->Fill(x, 4., wgt);
      if(ngleptons_>0 && ngjets>1)   fidCounter->Fill(x, 5., wgt);
      if(ngleptons_>1 && ngjets>1)   fidCounter->Fill(x, 6., wgt);
      if(ngleptons_>0 && ngjets>2)   fidCounter->Fill(x, 7., wgt);
      if(ngleptons_>1 && ngjets>2)   fidCounter->Fill(x, 8., wgt);
      if(ngleptons_>0 && ngjets>3)   fidCounter->Fill(x, 9., wgt);
      if(ngleptons_>1 && ngjets>3)   fidCounter->Fill(x, 10.,wgt);
    }
  
}


//
void MiniAnalyzer::recAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //VERTICES
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  if (vertices->empty()) return; // skip the event if no PV found
  const reco::Vertex &primVtx = vertices->front();
  reco::VertexRef primVtxRef(vertices,0);
   ev_.nvtx=vertices->size();
  if(ev_.nvtx==0) return;

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
  edm::Handle<pat::PackedTriggerPrescales> h_trigPrescale;
  iEvent.getByToken(triggerPrescales_, h_trigPrescale);
  edm::Handle<pat::PackedTriggerPrescales> h_l1trigPrescale;
  iEvent.getByToken(l1triggerPrescales_, h_l1trigPrescale);
  ev_.triggerBits=0;
  ev_.addTriggerBits=0;
  for (unsigned int i=0; i< h_trigRes->size(); i++) 
    {	
      if( !(*h_trigRes)[i].accept() ) continue;
      for(size_t itrig=0; itrig<triggersToUse_.size(); itrig++)
	{
	  if (triggerList[i].find(triggersToUse_[itrig])==string::npos) continue;          
          int prescale=h_trigPrescale->getPrescaleForIndex(i); 
          int l1prescale=h_l1trigPrescale->getPrescaleForIndex(i); 
          if(itrig<32)
            ev_.triggerBits |= (1 << itrig);
          else
            ev_.addTriggerBits |= (1 << (itrig-32));
	  histContainer_["triggerList"]->Fill(itrig);
          histContainer_["triggerPrescale"]->Fill(itrig,prescale);
          bool isZeroBias(triggerList[i].find("ZeroBias")!=string::npos);
          if(!isZeroBias) continue;
          ev_.zeroBiasPS=prescale*l1prescale;
	}
    }
  bool passTrigger(ev_.isData ? ev_.triggerBits!=0 : true);
  if(!passTrigger) return;

  //PF candidates
  edm::Handle<pat::PackedCandidateCollection> pfcands;
  iEvent.getByToken(pfToken_,pfcands);

  //
  //CTPPS local tracks (only present in data)
  //
  ev_.nfwdtrk=0;
  edm::Handle<vector<reco::ProtonTrack>> recoProtons;
  iEvent.getByToken(tokenRecoProtons_, recoProtons);
  if(recoProtons.isValid()){
    try{
      for (const auto & proton : *recoProtons)
        {
          if(!proton.valid()) continue;

          CTPPSDetId detid(* proton.contributingRPIds.begin());
          ev_.fwdtrk_pot[ev_.nfwdtrk]       = 100*detid.arm()+10*detid.station()+detid.rp();
          ev_.fwdtrk_chisqnorm[ev_.nfwdtrk] = proton.fitChiSq;
          ev_.fwdtrk_method[ev_.nfwdtrk]    = proton.method;
          ev_.fwdtrk_ex[ev_.nfwdtrk]        = proton.direction().x();
          ev_.fwdtrk_ey[ev_.nfwdtrk]        = proton.direction().y();
          ev_.fwdtrk_ez[ev_.nfwdtrk]        = proton.direction().z();
          ev_.fwdtrk_y[ev_.nfwdtrk]         = proton.vertex().y();

          float xi=proton.xi();
          ev_.fwdtrk_xi[ev_.nfwdtrk]        = xi;

          float th_x = proton.direction().x() / proton.direction().mag();
          float th_y = proton.direction().y() / proton.direction().mag();
          float mp = 0.938; // GeV
          float Eb = 6500.; // GeV
          float t0 = 2.*pow(mp,2) + 2.*pow(Eb,2)*(1.-xi) - 2.*sqrt( (pow(mp,2) + pow(Eb,2)) * (pow(mp,2) + pow(Eb,2)*pow(1.-xi,2)) );
          float th = sqrt(th_x * th_x + th_y * th_y);
          float S = sin(th/2.);
          ev_.fwdtrk_t[ev_.nfwdtrk] = t0 - 4. * pow(Eb,2)* (1.-xi) * S*S;
                  
          ev_.nfwdtrk++;
        }
    }
    catch(...){
    }
  }

  //
  //LEPTON SELECTION 
  ev_.nl=0; 
  
  //MUON SELECTION: cf. https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2
  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);
  ev_.nrawmu=0;
  std::vector<Double_t> muStatUncReplicas(100,0);
  for (const pat::Muon &mu : *muons) 
    { 

      //raw muon info
      ev_.rawmu_pt[ev_.nrawmu]=(Short_t)mu.pt();
      ev_.rawmu_eta[ev_.nrawmu]=(Short_t)10*mu.eta();
      ev_.rawmu_phi[ev_.nrawmu]=(Short_t)10*mu.phi();
      ev_.rawmu_pid[ev_.nrawmu]= mu.selectors();
      ev_.nrawmu++;


      //apply correction
      float pt  = mu.pt();
      if(pt<2) continue; //no need to care about very low pt muons here... (corrections will tend to be meaningless)
      float eta = mu.eta();
      float phi = mu.phi();
      float q   = mu.charge();
      const reco::GenParticle * gen=mu.genLepton(); 

      //rochester corrections
      float sf(1.0),smearSeed(-1);
      float statUnc(0),zptUnc(0),ewkUnc(0),deltamUnc(0);
      if(iEvent.isRealData())
        {
          sf = muonRC_->kScaleDT(q, pt, eta, phi);
          for(int i=0; i<100; i++)
            muStatUncReplicas[i] = muonRC_->kScaleDT(q, pt, eta, phi, 1,i);
          statUnc = (1.0+TMath::StdDev(muStatUncReplicas.begin(),muStatUncReplicas.end()))/sf;
          zptUnc = muonRC_->kScaleDT(q, pt, eta, phi, 2)/sf;
          ewkUnc = muonRC_->kScaleDT(q, pt, eta, phi, 3)/sf;
          deltamUnc = muonRC_->kScaleDT(q, pt, eta, phi, 4)/sf;
        }
      else
        {
          smearSeed=gRandom->Rndm();
          int tlwm=(mu.innerTrack().isNonnull() ? mu.innerTrack()->hitPattern().trackerLayersWithMeasurement() : 0);
          sf = gen ?
            muonRC_->kSpreadMC(q, pt, eta, phi, gen->pt()) :
            muonRC_->kSmearMC(q, pt, eta, phi, tlwm , smearSeed);
          for(int i=0; i<100; i++)
            muStatUncReplicas[i]=gen ?
              muonRC_->kSpreadMC(q, pt, eta, phi, gen->pt(),1,i) :
              muonRC_->kSmearMC(q, pt, eta, phi, tlwm, smearSeed,1,i);
          statUnc=(1.0+TMath::StdDev(muStatUncReplicas.begin(),muStatUncReplicas.end()))/sf;
          zptUnc=(gen ?
                  muonRC_->kSpreadMC(q, pt, eta, phi, gen->pt(),2) :
                  muonRC_->kSmearMC(q, pt, eta, phi, tlwm, smearSeed,2))/sf;
          ewkUnc=(gen ?
                  muonRC_->kSpreadMC(q, pt, eta, phi, gen->pt(),3) :
                  muonRC_->kSmearMC(q, pt, eta, phi, tlwm, smearSeed,3))/sf;
          deltamUnc=(gen ?
                     muonRC_->kSpreadMC(q, pt, eta, phi, gen->pt(),4) :
                     muonRC_->kSmearMC(q, pt, eta, phi, tlwm, smearSeed,4))/sf;
        }      
      
      auto p4  = mu.p4() * sf;

      //kinematics
      bool passPt(p4.Pt() > 10);
      bool passEta(fabs(p4.Eta()) < 2.5);
      if(!passPt || !passEta) continue;

      //ID
      bool isLoose(muon::isLooseMuon(mu));
      if(!isLoose) continue;

      //save info
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
      
      ev_.l_charge[ev_.nl]   = q;
      ev_.l_pt[ev_.nl]       = p4.Pt();
      ev_.l_eta[ev_.nl]      = p4.Eta();
      ev_.l_phi[ev_.nl]      = p4.Phi();
      ev_.l_mass[ev_.nl]     = p4.M();
      ev_.l_scaleUnc1[ev_.nl]= statUnc;
      ev_.l_scaleUnc2[ev_.nl]= zptUnc;
      ev_.l_scaleUnc3[ev_.nl]= ewkUnc;
      ev_.l_scaleUnc4[ev_.nl]= deltamUnc;
      ev_.l_highpt[ev_.nl]   = mu.tunePMuonBestTrack()->pt();
      ev_.l_scaleUnc5[ev_.nl]= mu.tunePMuonBestTrack()->ptError();
      ev_.l_mva[ev_.nl]      = 0;
      ev_.l_pid[ev_.nl]      = mu.selectors();
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

      if( p4.Pt()>20 && fabs(p4.Eta())<2.5 && isLoose) nrecleptons_++;
    }
  
  // ELECTRON SELECTION: cf. https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
  edm::Handle<edm::View<pat::Electron> > electrons;
  iEvent.getByToken(electronToken_, electrons);    
  for (const pat::Electron &e : *electrons) 
    {        

      auto corrP4  = e.p4() * e.userFloat("ecalTrkEnergyPostCorr") / e.energy();

      //kinematics cuts
      bool passPt(corrP4.pt() > 15.0);
      bool passEta(fabs(corrP4.eta()) < 2.5 && (fabs(e.superCluster()->eta()) < 1.4442 || fabs(e.superCluster()->eta()) > 1.5660));
      if(!passPt || !passEta) continue;
      
      //full id+iso decisions
      bool isVeto( e.electronID("cutBasedElectronID-Fall17-94X-V1-veto") );
      int vetoBits( e.userInt("cutBasedElectronID-Fall17-94X-V1-veto")  );
      bool passVetoId( (vetoBits | 0xc0)== 0x3ff);  //mask isolation cuts and require all bits active      
      bool isLoose( e.electronID("cutBasedElectronID-Fall17-94X-V1-loose") );
      int looseBits( e.userInt("cutBasedElectronID-Fall17-94X-V1-loose")  );
      bool passLooseId( (looseBits | 0xc0)== 0x3ff);  //mask isolation cuts and require all bits active
      bool isMedium( e.electronID("cutBasedElectronID-Fall17-94X-V1-medium") );
      int mediumBits( e.userInt("cutBasedElectronID-Fall17-94X-V1-medium")  );
      bool passMediumId( (mediumBits | 0xc0)== 0x3ff);  //mask isolation cuts and require all bits active
      bool isTight( e.electronID("cutBasedElectronID-Fall17-94X-V1-tight") );
      int tightBits( e.userInt("cutBasedElectronID-Fall17-94X-V1-tight") );
      bool passTightId( (tightBits | 0xc0)== 0x3ff);  //mask isolation cuts and require all bits active

      bool mvawp80(e.electronID("mvaEleID-Fall17-iso-V2-wp80"));
      bool mvawp90(e.electronID("mvaEleID-Fall17-iso-V2-wp90"));
      bool mvawploose(e.electronID("mvaEleID-Fall17-iso-V2-wpLoose"));
      bool mvanonisowp80(e.electronID("mvaEleID-Fall17-noIso-V2-wp80"));
      bool mvanonisowp90(e.electronID("mvaEleID-Fall17-noIso-V2-wp90"));
      bool mvanonisowploose(e.electronID("mvaEleID-Fall17-noIso-V2-wpLoose"));
      bool passHEEP(e.electronID("heepElectronID-HEEPV70"));

      //impact parameter cuts
      //see details in https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
      bool passIpCuts(true);
      if(e.gsfTrack().isNonnull())
	{
	  float dxy(fabs(e.gsfTrack()->dxy(primVtx.position())));
	  float dz(fabs(e.gsfTrack()->dz(primVtx.position())));
	  if(fabs(e.superCluster()->eta()) < 1.4442)
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

      //save the electron
      const reco::GenParticle * gen=e.genLepton(); 
      ev_.l_isPromptFinalState[ev_.nl] = gen ? gen->isPromptFinalState() : false;
      ev_.l_isDirectPromptTauDecayProductFinalState[ev_.nl] = gen ? gen->isDirectPromptTauDecayProductFinalState() : false;
      ev_.l_id[ev_.nl]=11;
      ev_.l_g[ev_.nl]=-1;
      for(int ig=0; ig<ev_.ng; ig++)
	{
	  if(abs(ev_.g_id[ig])!=ev_.l_id[ev_.nl]) continue;
	  if(deltaR( corrP4.eta(),corrP4.phi(), ev_.g_eta[ig],ev_.g_phi[ig])>0.4) continue;
	  ev_.l_g[ev_.nl]=ig;
	  break;
	}	      
      ev_.l_mva[ev_.nl]=e.userFloat("ElectronMVAEstimatorRun2Fall17IsoV2Values");
      ev_.l_mvaCats[ev_.nl]=e.userInt("ElectronMVAEstimatorRun2Fall17IsoV2Categories");

      ev_.l_pid[ev_.nl]=0;
      ev_.l_pid[ev_.nl]= (passVetoId | (isVeto<<1) 
			  | (passLooseId<<2) | (isLoose<<3) 
			  | (passMediumId<<4) | (isMedium<<5) 
			  | (passTightId<<6) | (isTight<<7)
			  | (passIpCuts<<8) 
                          | (mvawp80<<9) | (mvawp90<<10) | (mvawploose<<11)
                          | (mvanonisowp80<<12) | (mvanonisowp90<<13) | (mvanonisowploose<<14)
                          | (passHEEP<<15)
			 );

      ev_.l_charge[ev_.nl]   = e.charge();
      ev_.l_pt[ev_.nl]       = corrP4.pt();
      ev_.l_highpt[ev_.nl]   = corrP4.pt();
      ev_.l_eta[ev_.nl]      = corrP4.eta();
      ev_.l_phi[ev_.nl]      = corrP4.phi();
      ev_.l_mass[ev_.nl]     = corrP4.mass();
      ev_.l_scaleUnc1[ev_.nl] = 0.5*(e.userFloat("energyScaleStatUp")-e.userFloat("energyScaleStatDown"));
      ev_.l_scaleUnc2[ev_.nl] = 0.5*(e.userFloat("energyScaleGainUp")-e.userFloat("energyScaleGainDown"));
      ev_.l_scaleUnc3[ev_.nl] = 0.5*(e.userFloat("energyScaleSystUp")-e.userFloat("energyScaleSystDown"));
      ev_.l_scaleUnc4[ev_.nl] = 0.5*(e.userFloat("energySigmaUp")-e.userFloat("energySigmaDown"));
      ev_.l_scaleUnc5[ev_.nl] = 0.5*(e.userFloat("energySigmaPhiUp")-e.userFloat("energySigmaPhiDown"));
      ev_.l_scaleUnc6[ev_.nl] = 0.5*(e.userFloat("energySigmaRhoUp")-e.userFloat("energySigmaRhoDown"));
      ev_.l_scaleUnc7[ev_.nl] = 0.5*(e.userFloat("energyScaleEtUp")-e.userFloat("energyScaleEtDown"));
      ev_.l_miniIso[ev_.nl]  = getMiniIsolation(pfcands,&e,0.05, 0.2, 10., false);
      ev_.l_relIso[ev_.nl]   = (e.chargedHadronIso()+ max(0., e.neutralHadronIso() + e.photonIso()  - 0.5*e.puChargedHadronIso()))/corrP4.pt();     
      ev_.l_chargedHadronIso[ev_.nl] = e.chargedHadronIso();
      ev_.l_ip3d[ev_.nl]     = -9999.;
      ev_.l_ip3dsig[ev_.nl]  = -9999;
      if(e.gsfTrack().get())
	{
	  std::pair<bool,Measurement1D> ip3dRes = getImpactParameter<reco::GsfTrackRef>(e.gsfTrack(), primVtxRef, iSetup, true);
	  ev_.l_ip3d[ev_.nl]    = ip3dRes.second.value();
	  ev_.l_ip3dsig[ev_.nl] = ip3dRes.second.significance();
	}
      ev_.nl++;
      
      if( corrP4.pt()>20 && passEta && passLooseId ) nrecleptons_++;
    }

  // PHOTON SELECTION: cf. https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2
  ev_.ngamma=0;
  edm::Handle<edm::View<pat::Photon> > photons;
  iEvent.getByToken(photonToken_, photons);    
  for (const pat::Photon &g : *photons)
    {        
      auto corrP4  = g.p4() * g.userFloat("ecalEnergyPostCorr") / g.energy();

      //kinematics cuts
      bool passPt(corrP4.pt() > 30.0);
      float eta=corrP4.eta();
      bool passEta(fabs(eta) < 2.5 && (fabs(eta) < 1.4442 || fabs(eta) > 1.5660));
      if(!passPt || !passEta) continue;

      //full id+iso decisions
      //bool isLoose( g.electronID("cutBasedPhotonID-Fall17-94X-V1-loose") );
      int looseBits( g.userInt("cutBasedPhotonID-Fall17-94X-V1-loose") );
      //bool passLooseId( (looseBits & 0x3) == 0x3 ); //require first two bits (h/e + sihih)
      //bool isMedium( g.electronID("cutBasedPhotonID-Fall17-94X-V1-medium") );
      int mediumBits( g.userInt("cutBasedPhotonID-Fall17-94X-V1-medium") );
      //bool passMediumId( (mediumBits & 0x3)== 0x3); //require first two bits (h/e + sihih)
      //bool isTight( g.photonID("acutBasedPhotonID-Fall17-94X-V1-tight") );
      int tightBits( g.userInt("cutBasedPhotonID-Fall17-94X-V1-tight") );
      //bool passTightId( (tightBits & 0x3)== 0x3);  //require first two bits (h/e + sihih)
      bool ismvawp80( g.photonID("mvaPhoID-RunIIFall17-v2-wp80"));
      bool ismvawp90( g.photonID("mvaPhoID-RunIIFall17-v2-wp90"));

      //save the photon
      const reco::GenParticle * gen=(const reco::GenParticle *)g.genPhoton(); 
      ev_.gamma_isPromptFinalState[ev_.ngamma] = gen ? gen->isPromptFinalState() : false;
      ev_.gamma_g[ev_.ngamma]=-1;
      for(int ig=0; ig<ev_.ng; ig++)
	{
	  if(abs(ev_.g_id[ig])!=22) continue;
	  if(deltaR( corrP4.eta(),corrP4.phi(), ev_.g_eta[ig],ev_.g_phi[ig])>0.4) continue;
	  ev_.gamma_g[ev_.ngamma]=ig;
	  break;
	}	      
      
      ev_.gamma_mva[ev_.ngamma]=g.userFloat("PhotonMVAEstimatorRunIIFall17v1Values");
      ev_.gamma_mvaCats[ev_.ngamma]=g.userInt("PhotonMVAEstimatorRunIIFall17v1Categories");
      ev_.gamma_idFlags[ev_.ngamma]= g.passElectronVeto() | (g.hasPixelSeed()<<1) | (ismvawp80<<2) | (ismvawp90<<3);
      ev_.gamma_pid[ev_.ngamma]= ( (looseBits & 0x3ff)
                                   | ((mediumBits & 0x3ff)<<10)
                                   | ((tightBits & 0x3ff)<<20));
      ev_.gamma_pt[ev_.ngamma]  = corrP4.pt();
      ev_.gamma_eta[ev_.ngamma] = corrP4.eta();
      ev_.gamma_phi[ev_.ngamma] = corrP4.phi();   
      ev_.gamma_scaleUnc1[ev_.ngamma] = 0.5*(g.userFloat("energyScaleStatUp")-g.userFloat("energyScaleStatDown"));
      ev_.gamma_scaleUnc2[ev_.ngamma] = 0.5*(g.userFloat("energyScaleGainUp")-g.userFloat("energyScaleGainDown"));
      ev_.gamma_scaleUnc3[ev_.ngamma] = 0.5*(g.userFloat("energyScaleSystUp")-g.userFloat("energyScaleSystDown"));
      ev_.gamma_scaleUnc4[ev_.ngamma] = 0.5*(g.userFloat("energySigmaUp")-g.userFloat("energySigmaDown"));
      ev_.gamma_scaleUnc5[ev_.ngamma] = 0.5*(g.userFloat("energySigmaPhiUp")-g.userFloat("energySigmaPhiDown"));
      ev_.gamma_scaleUnc6[ev_.ngamma] = 0.5*(g.userFloat("energySigmaRhoUp")-g.userFloat("energySigmaRhoDown"));
      ev_.gamma_scaleUnc7[ev_.ngamma] = 0.5*(g.userFloat("energyScaleEtUp")-g.userFloat("energyScaleEtDown"));
      ev_.gamma_chargedHadronIso[ev_.ngamma] = g.chargedHadronIso();
      ev_.gamma_neutralHadronIso[ev_.ngamma] = g.neutralHadronIso();
      ev_.gamma_photonIso[ev_.ngamma]        = g.photonIso();
      ev_.gamma_hoe[ev_.ngamma]              = g.hadTowOverEm();
      ev_.gamma_sieie[ev_.ngamma]            = g.full5x5_sigmaIetaIeta();
      ev_.gamma_r9[ev_.ngamma]               = g.full5x5_r9();
      ev_.ngamma++;
      if(ev_.ngamma>50) break;
      if( corrP4.pt()>30 && passEta) nrecphotons_++;
    }

  // JETS
  ev_.nj=0; 
  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByToken(jetToken_,jets);
  JME::JetResolution jerResolution  = JME::JetResolution::get(iSetup, "AK4PFchs_pt");
  JME::JetResolutionScaleFactor jerResolutionSF = JME::JetResolutionScaleFactor::get(iSetup, "AK4PFchs");
  std::vector< std::pair<const reco::Candidate *,int> > clustCands;
  for(auto j = jets->begin();  j != jets->end(); ++j)
    {
      //base kinematics
      if(j->pt()<15 || fabs(j->eta())>4.7) continue;
      
      //resolution corrections
      float jerSF[]={1.0,1.0,1.0};
      ev_.j_g[ev_.nj] = -1;
      if(!iEvent.isRealData())
        {
          //match to gen level jet/parton
          float genj_pt(0);
          const reco::Candidate *genParton = j->genParton();
          ev_.j_flav[ev_.nj]       = j->partonFlavour();
          ev_.j_hadflav[ev_.nj]    = j->hadronFlavour();
          ev_.j_pid[ev_.nj]        = genParton ? genParton->pdgId() : 0;
          for(int ig=0; ig<ev_.ng; ig++)
            {
              if(abs(ev_.g_id[ig])==11 || abs(ev_.g_id[ig])==13) continue;
              if(deltaR( j->eta(),j->phi(), ev_.g_eta[ig],ev_.g_phi[ig])>0.4) continue;
              genj_pt=ev_.g_pt[ig];
              ev_.j_g[ev_.nj]=ig;
              ev_.g_xbp[ig]  = genParton   ? ev_.g_xb[ig]*genj_pt/genParton->pt() : 0.;
              break;
            }	 
          
          //jet energy resolution
          JME::JetParameters jerParams = {{JME::Binning::JetPt, j->pt()}, 
                                          {JME::Binning::JetEta, j->eta()},
                                          {JME::Binning::Rho, rho}};
          float r = jerResolution.getResolution(jerParams);
          jerSF[0] = jerResolutionSF.getScaleFactor(jerParams);
          jerSF[1] = jerResolutionSF.getScaleFactor(jerParams, Variation::UP);
          jerSF[2] = jerResolutionSF.getScaleFactor(jerParams, Variation::DOWN);
          for(int i=0; i<3; i++) {
            //use stochasting smearing for unmatched jets
            if(genj_pt<=0)
              {
                float sigma = r * std::sqrt(std::max(float( pow(jerSF[i],2)- 1.0),float(0.)));
                jerSF[i] = std::max(float(1.0 + gRandom->Gaus(0, sigma)),float(0.));
              }
            else {           
              float dPt = j->pt()-genj_pt;
              jerSF[i] = std::max(float(1.0 + (jerSF[i] - 1.) * dPt / j->pt()),float(0.));
            }
           
            //make up/down variations relative
            if(jerSF[0]>0) { jerSF[1]/=jerSF[0]; jerSF[2]/=jerSF[0]; }
          }
        }

      auto corrP4  = j->p4() * jerSF[0];

      //jet id cf.
      //2016 https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016
      //2017 https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2017
      float NHF  = j->neutralHadronEnergyFraction();
      float NEMF = j->neutralEmEnergyFraction();
      float CHF  = j->chargedHadronEnergyFraction();
      float MUF  = j->muonEnergyFraction();
      float CEMF = j->chargedEmEnergyFraction();
      float NumChargedParticles = j->chargedMultiplicity();
      float NumNeutralParticles = j->neutralMultiplicity();
      float NumConst = NumChargedParticles+NumNeutralParticles;
      float CHM = j->chargedMultiplicity();

      bool tightLepVeto(true),looseJetID(true);//,tightJetId(true);
      if(abs(j->eta())<2.4) {
        looseJetID = (NHF<0.99 && NEMF<0.99 && NumConst>1 && CHF>0 && CHM>0 && CEMF<0.99);
        //tightJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1 && CHF>0 && CHM>0 && CEMF<0.99);
        tightLepVeto = (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8 && CHF>0 && CHM>0 && CEMF<0.80);
      }
      else if(abs(j->eta())<2.7) {
        looseJetID = (NHF<0.99 && NEMF<0.99 && NumConst>1);
        //tightJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1);
        tightLepVeto = (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8);
      }
      else if(abs(j->eta())<3.0) {
        looseJetID = (NHF<0.98 && NEMF>0.01 && NumNeutralParticles>2);
        //tightJetID = (NHF<0.98 && NEMF>0.01 && NumNeutralParticles>2);
        tightLepVeto = (NHF<0.98 && NEMF>0.01 && NumNeutralParticles>2);
      }
      else {
        looseJetID = (NEMF<0.90 && NumNeutralParticles>10);
        //tightJetID = (NEMF<0.90 && NumNeutralParticles>10);
        tightLepVeto = (NEMF<0.90 && NumNeutralParticles>10);
      }

      if(jetIdToUse_=="tightLepVeto") { if(!tightLepVeto) continue; }
      else { if(!looseJetID) continue; }

      //save jet
      ev_.j_area[ev_.nj]    = j->jetArea();
      ev_.j_jerUp[ev_.nj]   = jerSF[1];
      ev_.j_jerDn[ev_.nj]   = jerSF[2];
      for(size_t iunc=0; iunc<jecCorrectionUncs_.size(); iunc++){
        jecCorrectionUncs_[iunc]->setJetPt(j->pt());
        jecCorrectionUncs_[iunc]->setJetEta(j->eta());
        ev_.j_jecUp[iunc][ev_.nj]=1.+jecCorrectionUncs_[iunc]->getUncertainty(true);
        jecCorrectionUncs_[iunc]->setJetPt(j->pt());
        jecCorrectionUncs_[iunc]->setJetEta(j->eta());
        ev_.j_jecDn[iunc][ev_.nj]=1.+jecCorrectionUncs_[iunc]->getUncertainty(false);
      }
      ev_.j_rawsf[ev_.nj]   = j->correctedJet("Uncorrected").pt()/j->pt();
      ev_.j_pt[ev_.nj]      = corrP4.pt();
      ev_.j_mass[ev_.nj]    = corrP4.mass();
      ev_.j_eta[ev_.nj]     = corrP4.eta();
      ev_.j_phi[ev_.nj]     = corrP4.phi();
      ev_.j_qg[ev_.nj]      = j->userFloat("QGTagger:qgLikelihood");
      ev_.j_pumva[ev_.nj]   = j->userFloat("pileupJetId:fullDiscriminant");
      ev_.j_id[ev_.nj]      = j->userInt("pileupJetId:fullId");
      ev_.j_csv[ev_.nj]     = j->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
      ev_.j_deepcsv[ev_.nj] = j->bDiscriminator("pfDeepCSVJetTags:probb") + j->bDiscriminator("pfDeepCSVJetTags:probbb");
      ev_.j_btag[ev_.nj]    = (ev_.j_deepcsv[ev_.nj]>0.4941);
      ev_.j_emf[ev_.nj]     = CEMF+NEMF;

      //jet shape variables
      ev_.j_c2_00[ev_.nj]    = getC(2, 0.0, &(*j), true, 0.9);
      ev_.j_c2_02[ev_.nj]    = getC(2, 0.2, &(*j), true, 0.9);
      ev_.j_c2_05[ev_.nj]    = getC(2, 0.5, &(*j), true, 0.9);
      ev_.j_c2_10[ev_.nj]    = getC(2, 1.0, &(*j), true, 0.9);
      ev_.j_c2_20[ev_.nj]    = getC(2, 2.0, &(*j), true, 0.9);
      ev_.j_zg[ev_.nj]       = getZg(&(*j),true,0.9)[0];
      ev_.j_mult[ev_.nj]     = calcGA(0,0,&(*j),true,0.9);
      ev_.j_gaptd[ev_.nj]    = calcGA(0,2,&(*j),true,0.9);
      ev_.j_gawidth[ev_.nj]  = calcGA(1,1,&(*j),true,0.9);
      ev_.j_gathrust[ev_.nj] = calcGA(2,1,&(*j),true,0.9);
      //this function throws an exception in rare events
      try{
        ev_.j_tau32[ev_.nj]    = getTau(3,2,&(*j),true,0.9);
        ev_.j_tau21[ev_.nj]    = getTau(2,1,&(*j),true,0.9);
      }
      catch(...){
        ev_.j_tau32[ev_.nj]=-99;
        ev_.j_tau21[ev_.nj]=-99;
      }
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
  edm::Handle<pat::METCollection> mets;
  iEvent.getByToken(metToken_, mets);
  ev_.met_pt  = mets->at(0).pt();
  ev_.met_phi = mets->at(0).phi();
  ev_.met_sig = mets->at(0).significance();
  for(size_t i=0; i<14; i++){
    ev_.met_ptShifted[i]  = mets->at(0).shiftedPt(pat::MET::METUncertainty(i));
    ev_.met_phiShifted[i] = mets->at(0).shiftedPhi(pat::MET::METUncertainty(i));
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
  LorentzVector vtxPt(0,0,0,0);
  ev_.nchPV=0; ev_.sumPVChPt=0; ev_.sumPVChPz=0; ev_.sumPVChHt=0;  
  for(int i=0; i<8; i++){
    ev_.nPFCands[i]=0;
    ev_.sumPFHt[i]=0;
    ev_.sumPFEn[i]=0;
    ev_.sumPFPz[i]=0;
    ev_.nPFChCands[i]=0;
    ev_.sumPFChHt[i]=0;
    ev_.sumPFChEn[i]=0;
    ev_.sumPFChPz[i]=0;
  }
  for(auto pf = pfcands->begin();  pf != pfcands->end(); ++pf)
    {
      int ieta(-1);
      if(pf->eta()>-4.7) ieta=0;
      if(pf->eta()>-3)   ieta=1;
      if(pf->eta()>-2.5) ieta=2;
      if(pf->eta()>-1.5) ieta=3;
      if(pf->eta()>0)    ieta=4;
      if(pf->eta()>1.5)  ieta=5;
      if(pf->eta()>2.5)  ieta=6;
      if(pf->eta()>3.0)  ieta=7;
      if(pf->eta()>4.7)  ieta=-1;
      if(ieta<0) continue;
      ev_.nPFCands[ieta]++;
      ev_.sumPFHt[ieta] += pf->pt();
      ev_.sumPFEn[ieta] += pf->energy();
      ev_.sumPFPz[ieta] += fabs(pf->pz());
      if(pf->charge()!=0){       
        ev_.nPFChCands[ieta]++;
        ev_.sumPFChHt[ieta] += pf->pt();
        ev_.sumPFChEn[ieta] += pf->energy();
        ev_.sumPFChPz[ieta] += fabs(pf->pz());
        bool passChargeSel(pf->pt()>0.9 && fabs(pf->eta())<2.5);
        const pat::PackedCandidate::PVAssoc pvassoc=pf->fromPV();
        if(passChargeSel && pvassoc>=pat::PackedCandidate::PVTight){
          vtxPt+=pf->p4();
          ev_.nchPV++;
          ev_.sumPVChPz+=fabs(pf->pz()); 
          ev_.sumPVChHt+=pf->pt();
        }
      }
    }

  ev_.sumPVChPt=vtxPt.pt(); 
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

  ngleptons_=0;   ngphotons_=0;
  nrecleptons_=0; nrecphotons_=0;
  ev_.g_nw=0; ev_.g_npsw=0; ev_.ng=0; ev_.ngtop=0;
  ev_.nl=0; ev_.ngamma=0; ev_.nj=0; ev_.nfwdtrk=0; ev_.nrawmu=0;
  
  //analyze the event
  if(!iEvent.isRealData()) genAnalysis(iEvent,iSetup);
  recAnalysis(iEvent,iSetup);
  
  //save event if at least one object at gen or reco level
  if(applyFilt_)
    if((ngleptons_==0 && ngphotons_==0 && nrecleptons_==0 && nrecphotons_==0) || !saveTree_) return;  
  ev_.run     = iEvent.id().run();
  ev_.lumi    = iEvent.luminosityBlock();
  ev_.event   = iEvent.id().event(); 
  ev_.isData  = iEvent.isRealData();
  tree_->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
MiniAnalyzer::beginJob(){
}

void MiniAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const& iSetup)
{

  if(histContainer_.find("genlumiwgts")!=histContainer_.end()) return;
  try{
    cout << "[MiniAnalyzer::beginLuminosityBlock]" << endl;
    edm::Handle<GenLumiInfoHeader> gen_header;
    iLumi.getByToken(generatorLumiHeaderToken_, gen_header);
    
    unsigned int nwgts(gen_header->weightNames().size());
    cout << nwgts << " weights are found" << endl;
    histContainer_["genlumiwgts"]=fs->make<TH1F>("genlumiwgts",";Weight name;",nwgts,0,nwgts);
    int i=1;
    for(auto it : gen_header->weightNames() ) {
      histContainer_["genlumiwgts"]->GetXaxis()->SetBinLabel(i,it.c_str());
      i++;
    }
  }
  catch(std::exception &e){
    std::cout << e.what() << endl
	      << "Failed to retrieve GenLumiInfoHeader" << std::endl;
  }

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
//    iRun.getByLabel( "externalLHEProducer", lheruninfo );

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
