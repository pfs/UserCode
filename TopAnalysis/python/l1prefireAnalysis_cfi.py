import FWCore.ParameterSet.Config as cms

def defineL1PrefireAnalysis(process,era):

    
    process.prefireVetoFilter = cms.EDFilter("TriggerRulePrefireVetoFilter",
                                             l1AcceptRecordLabel = cms.InputTag("scalersRawToDigi"),
                                             )

    process.prefiringVBFAna = cms.EDAnalyzer("PrefiringVBFAna",
                                             triggerRule = cms.InputTag("prefireVetoFilter:ruleIndex"),
                                             jetSrc = cms.InputTag("slimmedJets"),
                                             electronSrc = cms.InputTag("slimmedElectrons"),
                                             photonSrc = cms.InputTag("slimmedPhotons"),
                                             l1egSrc = cms.InputTag("caloStage2Digis:EGamma"),
                                             l1jetSrc = cms.InputTag("caloStage2Digis:Jet"),
                                             l1GtSrc = cms.InputTag("gtStage2Digis"),
                                             triggerObjects = cms.InputTag("selectedPatTrigger" if '2016' in era else "slimmedPatTrigger"),
                                             triggerPrescales = cms.InputTag("patTrigger"),
                                             )

    process.bxm1_pass = cms.EDAnalyzer("L1EGCheck",
                                       l1egSrc = cms.InputTag("caloStage2Digis:EGamma"),
                                       bx = cms.int32(-1),
                                       )    
    process.l1prefirePath = cms.Path(process.prefireVetoFilter+process.bxm1_pass+process.prefiringVBFAna)
    
    #process.bxm1_fail = cms.EDAnalyzer("L1EGCheck",
    #                                   l1egSrc = cms.InputTag("caloStage2Digis:EGamma"),
    #                                   bx = cms.int32(-1),
    #                                   )
    #process.checkPath = cms.Path(~process.prefireVetoFilter+process.bxm1_fail)
