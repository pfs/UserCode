import FWCore.ParameterSet.Config as cms

analysis = cms.EDAnalyzer("MiniAnalyzer",
                          saveTree        = cms.bool(True),
                          triggerBits     = cms.InputTag("TriggerResults","","HLT"),
                          prescales       = cms.InputTag("patTrigger"),
                          rho             = cms.InputTag("fixedGridRhoFastjetAll"),
                          vertices        = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          muons           = cms.InputTag("slimmedMuons"),
                          electrons       = cms.InputTag("slimmedElectrons"),
                          jets            = cms.InputTag("slimmedJetsReapplyJEC"),
                          mets            = cms.InputTag("slimmedMETs"),                          
                          puppimets       = cms.InputTag('slimmedMETsPuppi'),
                          pfCands         = cms.InputTag('packedPFCandidates'),                          
                          eleVetoIdMap    = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
                          eleTightIdMap   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
                          eleTightIdFullInfoMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
                          muTriggersToUse = cms.vstring('IsoMu20_v', 'IsoTkMu20_v'),
                          elTriggersToUse = cms.vstring('Ele22_eta2p1_WPLoose_Gsf_v')
                          )
