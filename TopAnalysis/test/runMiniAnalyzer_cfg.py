import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('runOnData', False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Run this on real data"
                 )
options.register('runL1PrefireAna', False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Run L1 pre-firing analysis"
                 )
options.register('noParticleLevel', False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Do not run the particleLevel sequence"
                 )
options.register('era', 'era2017',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "era to use (configurations defined in python/EraConfig.py)"
                 )
options.register('outFilename', 'MiniEvents.root',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Output file name"
                 )
options.register('baseJetCollection','slimmedJets',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Base jet collection"
                 )
options.register('inputFile', None,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "input file to process"
                 )
options.register('lumiJson', None,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 'apply this lumi json file'
                 )
options.register('saveTree', True,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "save summary tree"
                 )
options.register('savePF', False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 'save PF candidates'
                 )
options.parseArguments()

#start process
process = cms.Process("MiniAnalysis")

#get the configuration to apply
from TopLJets2015.TopAnalysis.EraConfig import getEraConfiguration
globalTag, jecTag, jecDB, jerTag, jerDB = getEraConfiguration(era=options.era,isData=options.runOnData)

# Load the standard set of configuration modules
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

#EGM customization
from TopLJets2015.TopAnalysis.customizeEGM_cff import customizeEGM
customizeEGM(process=process,era=options.era)

# global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globalTag)

process.load("GeneratorInterface.RivetInterface.mergedGenParticles_cfi")
process.load("GeneratorInterface.RivetInterface.genParticles2HepMC_cfi")
process.genParticles2HepMC.genParticles = cms.InputTag("mergedGenParticles")
process.load("GeneratorInterface.RivetInterface.particleLevel_cfi") 

# particle level definitions
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.mergedGenParticles = cms.EDProducer("MergedGenParticleProducer",
                                            inputPruned = cms.InputTag("prunedGenParticles"),
                                            inputPacked = cms.InputTag("packedGenParticles"),
                                            )
process.load('GeneratorInterface.RivetInterface.genParticles2HepMC_cfi')
process.genParticles2HepMC.genParticles = cms.InputTag("mergedGenParticles")

#jet energy corrections
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
from JetMETCorrections.Configuration.DefaultJEC_cff import *
from JetMETCorrections.Configuration.JetCorrectionServices_cff import *
from TopLJets2015.TopAnalysis.customizeJetTools_cff import *
customizeJetTools(process=process,
                  jecDB=jecDB,
                  jecTag=jecTag,
                  jerDB=jerDB,
                  jerTag=jerTag,
                  baseJetCollection=options.baseJetCollection,
                  runOnData=options.runOnData)

#message logger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# set input to process
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('/store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/90000/FE7ABAEB-4A42-E811-87A3-0CC47AD98D26.root'),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck') 
                            )
if '2016' in options.era:
      process.source.fileNames = cms.untracked.vstring('/store/mc/RunIISummer16MiniAODv3/DYToLL_0J_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/00000/004ABB41-17C1-E811-8194-0CC47AB35C32.root')

if options.runOnData:
      if '2017' in options.era:
            process.source.fileNames = cms.untracked.vstring('/store/data/Run2017F/SinglePhoton/MINIAOD/31Mar2018-v1/00000/0033B9DF-4338-E811-AB11-1CB72C1B6CC6.root')
            if options.runL1PrefireAna:
                  print 'Adding secondary filenames to run L1 prefire analysis'
                  process.source.secondaryFileNames = cms.untracked.vstring(['/store/data/Run2017F/SinglePhoton/AOD/17Nov2017-v1/70000/BC8110F6-3CE0-E711-9210-02163E014564.root',
                                                                             '/store/data/Run2017F/SinglePhoton/AOD/17Nov2017-v1/70000/F68021F2-3CE0-E711-85E6-02163E014702.root',
                                                                             '/store/data/Run2017F/SinglePhoton/AOD/17Nov2017-v1/60000/90EE7B12-BCDE-E711-8C47-A4BF0112BC06.root',
                                                                             '/store/data/Run2017F/SinglePhoton/AOD/17Nov2017-v1/60000/F8CD978A-17DF-E711-A05F-A0369F836334.root',
                                                                             '/store/data/Run2017F/SinglePhoton/AOD/17Nov2017-v1/60000/3E1B045C-28DF-E711-8B81-002590A37114.root',
                                                                             '/store/data/Run2017F/SinglePhoton/AOD/17Nov2017-v1/60000/506AF751-29DF-E711-9E1F-001E67397F2B.root',
                                                                             '/store/data/Run2017F/SinglePhoton/AOD/17Nov2017-v1/60000/1CDC7A53-20DF-E711-A272-1CC1DE1CDD20.root',
                                                                             '/store/data/Run2017F/SinglePhoton/AOD/17Nov2017-v1/60000/28249C05-0DDF-E711-BA56-008CFA0527CC.root',
                                                                             '/store/data/Run2017F/SinglePhoton/AOD/17Nov2017-v1/60000/142F6B6F-28DF-E711-9D48-A4BF0112BDF8.root'])
      else:
            process.source.fileNames = cms.untracked.vstring('/store/data/Run2016B/SingleMuon/MINIAOD/17Jul2018_ver2-v1/00000/0219DD4A-7D8C-E811-B7EB-001E67248A43.root')
            if options.runL1PrefireAna:
                  print 'Adding secondary filenames to run L1 prefire analysis'
                  process.source.secondaryFileNames = cms.untracked.vstring(['/store/data/Run2016B/SingleMuon/AOD/07Aug17_ver2-v1/110000/B6479036-9486-E711-A6CF-003048FFD734.root',
                                                                             '/store/data/Run2016B/SingleMuon/AOD/07Aug17_ver2-v1/110001/0A8A7032-9086-E711-BD1B-0025905A60C6.root',
                                                                             '/store/data/Run2016B/SingleMuon/AOD/07Aug17_ver2-v1/110001/22FFA42A-9086-E711-9B49-0CC47A7C3472.root',
                                                                             '/store/data/Run2016B/SingleMuon/AOD/07Aug17_ver2-v1/70000/767AE9A5-BA81-E711-9A8E-0CC47A745282.root',
                                                                             '/store/data/Run2016B/SingleMuon/AOD/07Aug17_ver2-v1/70000/D67D18A5-BA81-E711-B6C2-0CC47A7C361E.root',
                                                                             '/store/data/Run2016B/SingleMuon/AOD/07Aug17_ver2-v1/70000/E06164A7-BA81-E711-B612-0025905B858A.root',
                                                                             '/store/data/Run2016B/SingleMuon/AOD/07Aug17_ver2-v1/70001/302E2840-C481-E711-8913-002618FDA259.root',
                                                                             '/store/data/Run2016B/SingleMuon/AOD/07Aug17_ver2-v1/70001/687FA66D-BC81-E711-B41D-0CC47A745282.root',
                                                                             '/store/data/Run2016B/SingleMuon/AOD/07Aug17_ver2-v1/70001/6890466C-BC81-E711-9C62-0CC47A4C8E1C.root',
                                                                             '/store/data/Run2016B/SingleMuon/AOD/07Aug17_ver2-v1/70001/E64E346F-BC81-E711-93BC-0CC47A4C8F18.root',
                                                                             '/store/data/Run2016B/SingleMuon/AOD/07Aug17_ver2-v1/70001/EC3632F7-C081-E711-B219-0CC47A7452DA.root',
                                                                             '/store/data/Run2016B/SingleMuon/AOD/07Aug17_ver2-v1/70001/F6499DE0-C281-E711-9835-0CC47A78A33E.root'])

if options.inputFile:
      fileList=[]
      if '.root' in options.inputFile :
            fileList=[options.inputFile]
      else:
            import os
            inDir=options.inputFile
            fileList = ['file:'+os.path.join(inDir,f) for f in os.listdir(inDir)]
      process.source.fileNames = cms.untracked.vstring(fileList)
print  "Processing",process.source.fileNames

#apply lumi json, if passed in command line
if options.lumiJson:
    print 'Lumi sections will be selected with',options.lumiJson
    from FWCore.PythonUtilities.LumiList import LumiList
    myLumis = LumiList(filename = options.lumiJson).getCMSSWString().split(',')
    process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange()
    process.source.lumisToProcess.extend(myLumis)


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(options.outFilename)
                                  )

#analysis
from TopLJets2015.TopAnalysis.miniAnalyzer_cfi import  ANALYSISJETIDS,ANALYSISTRIGGERLISTS
process.load('TopLJets2015.TopAnalysis.miniAnalyzer_cfi')
print 'MiniAnalyzer configuration is as follows:'
process.analysis.saveTree=cms.bool(options.saveTree)
process.analysis.savePF=cms.bool(options.savePF)
print '\t save tree=',options.saveTree,' save PF=',options.savePF
if 'era2017' in options.era:
      process.analysis.jetIdToUse=ANALYSISJETIDS[2017]
      process.analysis.triggersToUse=ANALYSISTRIGGERLISTS[2017]
      print '\t Using 2017 triggers/jet ids'
else:
      process.analysis.jetIdToUse=ANALYSISJETIDS[2016]
      process.analysis.triggersToUse=ANALYSISTRIGGERLISTS[2016]
      print '\t Using 2016 triggers/jet ids'
if options.runOnData:
      process.analysis.metFilterBits = cms.InputTag("TriggerResults","","RECO")
      print '\t will save met filter bits'

#schedule execution
toSchedule=[]
if process.egammaPostReco:
      toSchedule.append( process.egammaPostReco )
if process.updatedPatJetsUpdatedJECBTag:
      process.custom_jec_seq=cms.Sequence(process.patJetCorrFactorsUpdatedJECBTag * process.updatedPatJetsUpdatedJECBTag)
      process.custom_jec=cms.Path(process.custom_jec_seq)
      toSchedule.append( process.custom_jec)
if process.fullPatMetSequenceModifiedMET:
      process.custom_met=cms.Path(process.fullPatMetSequenceModifiedMET)
      toSchedule.append(process.custom_met)
if not (options.runOnData or options.noParticleLevel):
      process.mctruth=cms.Path(process.mergedGenParticles*process.genParticles2HepMC*process.particleLevel)
      toSchedule.append( process.mctruth )
process.ana=cms.Path(process.analysis)
toSchedule.append( process.ana )
if options.runOnData and options.runL1PrefireAna:
      print 'Prefire analysis is scheduled to be executed'
      from TopLJets2015.TopAnalysis.l1prefireAnalysis_cfi import *
      defineL1PrefireAnalysis(process,options.era)
      toSchedule.append(process.l1prefirePath)

process.schedule=cms.Schedule( (p for p in toSchedule) )
