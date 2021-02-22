import FWCore.ParameterSet.Config as cms
process = cms.Process('GEN')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)

process.load("GeneratorInterface.RivetInterface.mergedGenParticles_cfi")
process.load("GeneratorInterface.RivetInterface.genParticles2HepMC_cff")
process.load("GeneratorInterface.RivetInterface.rivetAnalyzer_cfi")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(200000)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'/store/relval/CMSSW_10_6_0/RelValTTbar_13/MINIAODSIM/PUpmx25ns_106X_upgrade2018_realistic_v4-v1/10000/DB4F12C2-3B53-4247-A1C3-BEB74F177362.root',
        #'/store/mc/RunIIAutumn18MiniAOD/TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/00000/0126270E-DB20-FB42-A2B6-CE59183BEB12.root',
#        '/store/mc/RunIIAutumn18MiniAOD/TT_noSC_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/60000/26D8AE17-B828-C948-A961-A785C6CB71CE.root',
        #'/store/mc/RunIIAutumn18MiniAOD/W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/1A3A8B55-3655-114D-B933-DA41260548F4.root',
        '/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/0693E0E7-97BE-E611-B32F-0CC47A78A3D8.root',
        '/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/0806AB92-99BE-E611-9ECD-0025905A6138.root',
         '/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/165F54A0-A3BE-E611-B3F7-0025905A606A.root',
        '/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/18E31463-B3BE-E611-B6A3-0CC47A4D7678.root',
        '/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/22B7C5F1-A4BE-E611-AC0A-0025905B85AE.root',
        '/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/26ABF488-A0BE-E611-BEEB-0CC47A4D7640.root'
    ),
)

process.options = cms.untracked.PSet()

process.genParticles2HepMC.genParticles = cms.InputTag("mergedGenParticles")
process.rivetAnalyzer.HepMCCollection = cms.InputTag("genParticles2HepMC:unsmeared")

#process.rivetAnalyzer.AnalysisNames = cms.vstring('CMS_2013_I1224539_DIJET', 'MC_TTBAR', 'CMS_2018_I1662081')
process.rivetAnalyzer.AnalysisNames = cms.vstring('CMS_2016_Viesturs_new')
process.rivetAnalyzer.OutputFile = cms.string('test.yoda')
process.rivetAnalyzer.CrossSection    = 831.76 # NNLO (arXiv:1303.6254)

process.rivetAnalyzer.useLHEweights = False

process.p = cms.Path(process.mergedGenParticles * process.genParticles2HepMC * process.rivetAnalyzer)
