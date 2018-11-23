/** 
    @class PrefiringVBFAna
    dumps a small tree with jets, photons and L1 objects to analyse the pre-fire probability
    largely based on https://github.com/nsmith-/PrefireAnalysis/
*/

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/transform.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"

#include "TTree.h"

#include "Math/LorentzVector.h"
#include "Math/PxPyPzE4D.h"
// == reco::LorentzVector
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> LorentzVector;

namespace {
  struct EventStruct {    
    int triggerRule;
    std::vector<LorentzVector> jet_p4;
    std::vector<float> jet_neutralEmFrac,jet_emFrac,jet_muFrac;
    // pdgid=0 if jet 22 if gamma
    // bits   jet: 0=loose ID,  1=tight ID, 2-4=match to HLT_PFJet 450,500,550
    //      gamma: 0=medium ID, 1=tight ID
    std::vector<int> jet_pdgid,jet_id;

    std::vector<int> L1EG_bx;
    std::vector<LorentzVector> L1EG_p4;
    std::vector<int> L1EG_iso;

    std::vector<int> L1Jet_bx;
    std::vector<LorentzVector> L1Jet_p4;

    std::vector<int> L1GtBx;
  };

  std::vector<const pat::TriggerObjectStandAlone*> getMatchedObjs(const float eta, const float phi, const std::vector<pat::TriggerObjectStandAlone>& trigObjs, const float maxDeltaR=0.1)
  {
    std::vector<const pat::TriggerObjectStandAlone*> matchedObjs;
    const float maxDR2 = maxDeltaR*maxDeltaR;
    for(auto& trigObj : trigObjs){
      const float dR2 = reco::deltaR2(eta,phi,trigObj.eta(),trigObj.phi());
      if(dR2<maxDR2) matchedObjs.push_back(&trigObj);
    }
    return matchedObjs;
  }
}

class PrefiringVBFAna : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:
    explicit PrefiringVBFAna(const edm::ParameterSet&);
    ~PrefiringVBFAna();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    // virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    // virtual void endJob() override;

    edm::EDGetTokenT<int> triggerRuleToken_;
    edm::EDGetTokenT<pat::JetCollection> jetToken_;
    edm::EDGetTokenT<pat::ElectronCollection>  electronToken_;
    edm::EDGetTokenT<pat::PhotonCollection> photonToken_;
    edm::EDGetTokenT<pat::METCollection> metToken_;
    edm::EDGetTokenT<BXVector<l1t::EGamma>> l1egToken_;
    edm::EDGetTokenT<BXVector<l1t::Jet>> l1jetToken_;
    edm::EDGetTokenT<BXVector<GlobalAlgBlk>> l1GtToken_;
    edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjectsToken_;
    edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescalesToken_;

    TTree * tree_;
    EventStruct event_;
};

PrefiringVBFAna::PrefiringVBFAna(const edm::ParameterSet& iConfig):
  triggerRuleToken_(consumes<int>(iConfig.getParameter<edm::InputTag>("triggerRule"))),  
  jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jetSrc"))),
  electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electronSrc"))),
  photonToken_(consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photonSrc"))),
  l1egToken_(consumes<BXVector<l1t::EGamma>>(iConfig.getParameter<edm::InputTag>("l1egSrc"))),
  l1jetToken_(consumes<BXVector<l1t::Jet>>(iConfig.getParameter<edm::InputTag>("l1jetSrc"))),
  l1GtToken_(consumes<BXVector<GlobalAlgBlk>>(iConfig.getParameter<edm::InputTag>("l1GtSrc"))),
  triggerObjectsToken_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
  triggerPrescalesToken_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("triggerPrescales")))
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  tree_ = fs->make<TTree>("l1prefire","Event Summary");
  tree_->Branch("triggerRule", &event_.triggerRule);
  tree_->Branch("jet_p4", &event_.jet_p4);
  tree_->Branch("jet_neutralEmFrac", &event_.jet_neutralEmFrac);
  tree_->Branch("jet_emFrac", &event_.jet_emFrac);
  tree_->Branch("jet_muFrac", &event_.jet_muFrac);
  tree_->Branch("jet_pdgid", &event_.jet_pdgid);  
  tree_->Branch("jet_id", &event_.jet_id);
  tree_->Branch("L1EG_bx", &event_.L1EG_bx);
  tree_->Branch("L1EG_p4", &event_.L1EG_p4);
  tree_->Branch("L1EG_iso", &event_.L1EG_iso);
  tree_->Branch("L1Jet_bx", &event_.L1Jet_bx);
  tree_->Branch("L1Jet_p4", &event_.L1Jet_p4);
  tree_->Branch("L1GtBx", &event_.L1GtBx);
}


//
PrefiringVBFAna::~PrefiringVBFAna()
{
}

//
void PrefiringVBFAna::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  Handle<int> triggerRuleHandle;
  iEvent.getByToken(triggerRuleToken_, triggerRuleHandle);
  event_.triggerRule = *triggerRuleHandle;

  Handle<pat::JetCollection> jetHandle;
  iEvent.getByToken(jetToken_, jetHandle);

  Handle<pat::ElectronCollection> elecHandle;
  iEvent.getByToken(electronToken_, elecHandle);

  Handle<pat::PhotonCollection> phoHandle;
  iEvent.getByToken(photonToken_, phoHandle);

  Handle<BXVector<l1t::EGamma>> l1egHandle;
  iEvent.getByToken(l1egToken_, l1egHandle);

  Handle<BXVector<l1t::Jet>> l1jetHandle;
  iEvent.getByToken(l1jetToken_, l1jetHandle);

  Handle<BXVector<GlobalAlgBlk>> l1GtHandle;
  iEvent.getByToken(l1GtToken_, l1GtHandle);

  Handle<pat::TriggerObjectStandAloneCollection> triggerObjectsHandle;
  iEvent.getByToken(triggerObjectsToken_, triggerObjectsHandle);

  Handle<pat::PackedTriggerPrescales> triggerPrescalesHandle;
  iEvent.getByToken(triggerPrescalesToken_, triggerPrescalesHandle);
  const edm::TriggerResults& triggerResults = triggerPrescalesHandle->triggerResults();

  event_.jet_p4.clear();
  event_.jet_neutralEmFrac.clear();
  event_.jet_emFrac.clear();
  event_.jet_muFrac.clear();
  event_.jet_pdgid.clear();
  event_.jet_id.clear();

  //jets
  std::string jetfilters[]={"hltDiCaloJet30MJJ300AllJetsDEta3Filter", 
                            "hltDiPFJet30MJJ300AllJetsDEta3Filter",
                            "hltDiCaloJet30MJJ600AllJetsDEta3Filter", 
                            "hltDiPFJet30MJJ600AllJetsDEta3Filter",
                            "hltSinglePFJet450",
                            "hltSinglePFJet500",
                            "hltSinglePFJet550"};
  for(const auto& jet : *jetHandle) {
    pat::Jet jetNew(jet);

    float NHF  = jetNew.neutralHadronEnergyFraction();
    float NEMF = jetNew.neutralEmEnergyFraction();
    float CHF  = jetNew.chargedHadronEnergyFraction();
    float MUF  = jetNew.muonEnergyFraction();
    float CEMF = jetNew.chargedEmEnergyFraction();
    float NumChargedParticles = jetNew.chargedMultiplicity();
    float NumNeutralParticles = jetNew.neutralMultiplicity();
    float NumConst = NumChargedParticles+NumNeutralParticles;
    float CHM = jetNew.chargedMultiplicity();
    
    bool tightLepVeto(true),looseJetID(true),tightJetID(true);
    if(abs(jetNew.eta())<2.4) {
      looseJetID = (NHF<0.99 && NEMF<0.99 && NumConst>1 && CHF>0 && CHM>0 && CEMF<0.99);
      tightJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1 && CHF>0 && CHM>0 && CEMF<0.99);
      tightLepVeto = (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8 && CHF>0 && CHM>0 && CEMF<0.80);
    }
    else if(abs(jetNew.eta())<2.7) {
      looseJetID = (NHF<0.99 && NEMF<0.99 && NumConst>1);
      tightJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1);
      tightLepVeto = (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8);
    }
    else if(abs(jetNew.eta())<3.0) {
      looseJetID = (NHF<0.98 && NEMF>0.01 && NumNeutralParticles>2);
      tightJetID = (NHF<0.98 && NEMF>0.01 && NumNeutralParticles>2);
      tightLepVeto = (NHF<0.98 && NEMF>0.01 && NumNeutralParticles>2);
    }
    else {
      looseJetID = (NEMF<0.90 && NumNeutralParticles>10);
      tightJetID = (NEMF<0.90 && NumNeutralParticles>10);
      tightLepVeto = (NEMF<0.90 && NumNeutralParticles>10);
    }

    int packedIds = (looseJetID <<0) | (tightJetID <<1) | (tightLepVeto<<2);

    if ( !looseJetID && !tightLepVeto) continue;
    std::vector<const pat::TriggerObjectStandAlone*> matchedTrigObjs = getMatchedObjs(jet.eta(), jet.phi(), *triggerObjectsHandle, 0.3);
    for (const auto trigObjConst : matchedTrigObjs) {
      pat::TriggerObjectStandAlone trigObj(*trigObjConst);
      trigObj.unpackNamesAndLabels(iEvent, triggerResults);
      for(size_t fidx=0; fidx<sizeof(jetfilters)/sizeof(std::string); fidx++)
        if ( trigObj.hasFilterLabel(jetfilters[fidx]) ) packedIds |= 1<<(3+fidx);
    }
    event_.jet_p4.push_back( jetNew.p4() );
    event_.jet_neutralEmFrac.push_back( NEMF );
    event_.jet_emFrac.push_back( NEMF+CEMF );
    event_.jet_muFrac.push_back( MUF );
    event_.jet_pdgid.push_back( 0 );
    event_.jet_id.push_back( packedIds );
  }

  //electrons
  std::string egfilters[]={"hltEG35L1SingleEGOrEtFilter",                           
                           "hltPreEle23Ele12CaloIdLTrackIdLIsoVLDZ",
                           "hltPreEle23Ele12CaloIdLTrackIdLIsoVL"};
  for (const pat::Electron &g : *elecHandle) {
    auto corrP4  = g.p4() * g.userFloat("ecalEnergyPostCorr") / g.energy();
    bool looseID(g.electronID("cutBasedElectronID-Fall17-94X-V1-loose"));
    bool mediumID(g.electronID("cutBasedElectronID-Fall17-94X-V1-medium"));
    bool tightID(g.electronID("cutBasedElectronID-Fall17-94X-V1-tight"));
    int packedIds = (looseID <<0)  | (mediumID<<1) | (tightID <<2);
    if ( !looseID || corrP4.pt()<20) continue;
    std::vector<const pat::TriggerObjectStandAlone*> matchedTrigObjs = getMatchedObjs(corrP4.eta(), corrP4.phi(), *triggerObjectsHandle, 0.3);
    for (const auto trigObjConst : matchedTrigObjs) {
      pat::TriggerObjectStandAlone trigObj(*trigObjConst);
      trigObj.unpackNamesAndLabels(iEvent, triggerResults);
      for(size_t fidx=0; fidx<sizeof(egfilters)/sizeof(std::string); fidx++)
        if ( trigObj.hasFilterLabel(egfilters[fidx]) ) packedIds |= 1<<(3+fidx);
    }    
    
    event_.jet_p4.push_back( corrP4 );
    event_.jet_neutralEmFrac.push_back( 1.0 );
    event_.jet_emFrac.push_back( 1.0 );
    event_.jet_muFrac.push_back( 0.0 );
    event_.jet_pdgid.push_back( 11 );
    event_.jet_id.push_back( packedIds );
  }


  
  //photons
  std::string pho_egfilters[]={"hltEG200EtFilter",
                           "hltEG200HEFilter",
                           "hltL1sSingleEG40to50"};
  for (const pat::Photon &g : *phoHandle) {
    auto corrP4  = g.p4() * g.userFloat("ecalEnergyPostCorr") / g.energy();
    bool looseID(g.photonID("cutBasedPhotonID-Fall17-94X-V1-loose"));
    bool mediumID(g.photonID("cutBasedPhotonID-Fall17-94X-V1-medium"));
    bool tightID(g.photonID("cutBasedPhotonID-Fall17-94X-V1-tight"));
    int packedIds = (looseID <<0)  | (mediumID<<1) | (tightID <<2);
    if ( !looseID || corrP4.pt()<20) continue;
    std::vector<const pat::TriggerObjectStandAlone*> matchedTrigObjs = getMatchedObjs(corrP4.eta(), corrP4.phi(), *triggerObjectsHandle, 0.3);
    for (const auto trigObjConst : matchedTrigObjs) {
      pat::TriggerObjectStandAlone trigObj(*trigObjConst);
      trigObj.unpackNamesAndLabels(iEvent, triggerResults);
      for(size_t fidx=0; fidx<sizeof(pho_egfilters)/sizeof(std::string); fidx++)
        if ( trigObj.hasFilterLabel(pho_egfilters[fidx]) ) packedIds |= 1<<(3+fidx);
    }    
      
    event_.jet_p4.push_back( corrP4 );
    event_.jet_neutralEmFrac.push_back( 1.0 );
    event_.jet_emFrac.push_back( 1.0 );
    event_.jet_muFrac.push_back( 0.0 );
    event_.jet_pdgid.push_back( 22 );
    event_.jet_id.push_back( packedIds );
  }
  
  event_.L1EG_bx.clear();
  event_.L1EG_p4.clear();
  event_.L1EG_iso.clear();
  for (auto bx=l1egHandle->getFirstBX(); bx<l1egHandle->getLastBX()+1; ++bx) {
    for (auto itL1=l1egHandle->begin(bx); itL1!=l1egHandle->end(bx); ++itL1) {
      event_.L1EG_bx.push_back(bx);
      event_.L1EG_p4.push_back(itL1->p4());
      event_.L1EG_iso.push_back(itL1->hwIso());
    }
  }

  event_.L1Jet_bx.clear();
  event_.L1Jet_p4.clear();
  for (auto bx=l1jetHandle->getFirstBX(); bx<l1jetHandle->getLastBX()+1; ++bx) {
    for (auto itL1=l1jetHandle->begin(bx); itL1!=l1jetHandle->end(bx); ++itL1) {
      event_.L1Jet_bx.push_back(bx);
      event_.L1Jet_p4.push_back(itL1->p4());
    }
  }

  event_.L1GtBx.clear();
  for (auto bx=l1GtHandle->getFirstBX(); bx<l1GtHandle->getLastBX()+1; ++bx) {
    event_.L1GtBx.push_back(l1GtHandle->begin(bx)->getFinalOR());
  }

  tree_->Fill();
}


void PrefiringVBFAna::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


DEFINE_FWK_MODULE(PrefiringVBFAna);
