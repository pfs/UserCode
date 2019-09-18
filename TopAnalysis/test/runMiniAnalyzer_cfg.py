import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('runOnData', False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Run this on real data"
                 )
options.register('runProtonFastSim', None,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 "Run proton fastsim for this angle"
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
options.register('secInputFile', None,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "secondary input file to process"
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
options.register('applyFilt', True,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 'save PF candidates'
                 )
options.parseArguments()

#start process
if options.runProtonFastSim:
      from Configuration.StandardSequences.Eras import eras
      process = cms.Process("MiniAnalysis", eras.ctpps_2016)      
else:
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
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# set input to process
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('/store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/90000/FE7ABAEB-4A42-E811-87A3-0CC47AD98D26.root'),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck') 
                            )
if '2016' in options.era:
      process.source.fileNames = cms.untracked.vstring('/store/mc/RunIISummer16MiniAODv3/TT_TuneCUETP8M2T4_PSweights_13TeV-powheg-pythia8/MINIAODSIM/94X_mcRun2_asymptotic_v3-v1/100000/B28E96ED-08A6-E911-8D2C-FA163EEC7BB7.root')

if options.runOnData:
      if '2017' in options.era:
            #process.source.fileNames = cms.untracked.vstring('/store/data/Run2017F/ZeroBias/MINIAOD/31Mar2018-v1/30000/8604984F-DF37-E811-ACFE-008CFA197B74.root')
            process.source.fileNames = cms.untracked.vstring('/store/data/Run2017B/SingleMuon/MINIAOD/31Mar2018-v1/90000/46CE7E24-EA37-E811-95CB-0025905A6132.root')
            if options.runL1PrefireAna:
                  print 'Adding secondary filenames to run L1 prefire analysis'
                  process.source.secondaryFileNames = cms.untracked.vstring(['/store/data/Run2017B/SingleMuon/AOD/17Nov2017-v1/70003/C4FDB602-77D8-E711-B975-02163E011E63.root',
                                                                             '/store/data/Run2017B/SingleMuon/AOD/17Nov2017-v1/70003/089281B2-74D8-E711-96D9-02163E019CF1.root',
                                                                             '/store/data/Run2017B/SingleMuon/AOD/17Nov2017-v1/40001/74593CEA-11D9-E711-A49E-FA163E62ECD8.root',
                                                                             '/store/data/Run2017B/SingleMuon/AOD/17Nov2017-v1/40001/4CC8D127-EED8-E711-91CC-02163E016CCF.root',
                                                                             '/store/data/Run2017B/SingleMuon/AOD/17Nov2017-v1/40001/3AFF8583-F0D8-E711-ADB8-FA163EA84AFA.root',
                                                                             '/store/data/Run2017B/SingleMuon/AOD/17Nov2017-v1/40001/24A6AA06-21D9-E711-B5AD-FA163E3BC51E.root',
                                                                             '/store/data/Run2017B/SingleMuon/AOD/17Nov2017-v1/40001/1A02BCE9-FDD8-E711-BB87-FA163E408AED.root',
                                                                             '/store/data/Run2017B/SingleMuon/AOD/17Nov2017-v1/40001/0A90781C-EBD8-E711-9198-FA163E9875C8.root',
                                                                             '/store/data/Run2017B/SingleMuon/AOD/17Nov2017-v1/40001/0481E927-F2D8-E711-A11A-FA163EFB3309.root',
                                                                             '/store/data/Run2017B/SingleMuon/AOD/17Nov2017-v1/40000/FC54E9BB-D8D8-E711-8076-FA163EC00925.root',
                                                                             '/store/data/Run2017B/SingleMuon/AOD/17Nov2017-v1/40000/DA1366B2-D2D8-E711-ACEE-FA163E24EC20.root',
                                                                             '/store/data/Run2017B/SingleMuon/AOD/17Nov2017-v1/40000/8CCE5537-D3D8-E711-B316-FA163E653C95.root',
                                                                             '/store/data/Run2017B/SingleMuon/AOD/17Nov2017-v1/40000/8492E4C5-C2D8-E711-9631-FA163EA53599.root',
                                                                             '/store/data/Run2017B/SingleMuon/AOD/17Nov2017-v1/40000/82906F7E-C4D8-E711-8A90-FA163EDC0DF0.root',
                                                                             '/store/data/Run2017B/SingleMuon/AOD/17Nov2017-v1/40000/38AFFEDF-C7D8-E711-85A5-FA163EE5AC62.root',
                                                                             '/store/data/Run2017B/SingleMuon/AOD/17Nov2017-v1/40000/16C2CF2F-CAD8-E711-8630-FA163E9C5649.root'
                                                                             ])
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
      secFileList=[]
      if '.root' in options.inputFile :
            fileList=[options.inputFile]
            if options.secInputFile and '.root' in options.secInputFile:
                  secFileList=[options.secInputFile]
                  print 'Will run also on',secFileList
      else:
            import os
            inDir=options.inputFile
            fileList = ['file:'+os.path.join(inDir,f) for f in os.listdir(inDir)]
      process.source.fileNames = cms.untracked.vstring(fileList)     
      process.source.secondaryFileNames = cms.untracked.vstring(secFileList)
      
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
process.analysis.saveTree  = cms.bool(options.saveTree)
process.analysis.savePF    = cms.bool(options.savePF)
process.analysis.applyFilt = cms.bool(options.applyFilt)
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
      process.custom_jec_seq=cms.Sequence(process.QGTagger * process.patJetCorrFactorsUpdatedJECBTag * process.updatedPatJetsUpdatedJECBTag)
      process.custom_jec=cms.Path(process.custom_jec_seq)
      toSchedule.append( process.custom_jec)
if process.fullPatMetSequenceModifiedMET:
      process.custom_met=cms.Path(process.fullPatMetSequenceModifiedMET)
      toSchedule.append(process.custom_met)
if not (options.runOnData or options.noParticleLevel):
      process.mctruth=cms.Path(process.mergedGenParticles*process.genParticles2HepMC*process.particleLevel)
      toSchedule.append( process.mctruth )
if options.runOnData:
      if 'era2017' in options.era:
            process.load("RecoCTPPS.ProtonReconstruction.year_2017_OF.ctppsProtonReconstructionOF_cfi")
            process.ppsReco=cms.Path(process.ctppsProtonReconstructionOFDB)
      else:
            process.load("RecoCTPPS.ProtonReconstruction.year_2016.ctppsProtonReconstruction_cfi")
            process.ppsReco=cms.Path(process.ctppsProtonReconstruction)
      toSchedule.append(process.ppsReco)

if options.runProtonFastSim:
      from TopLJets2015.TopAnalysis.protonReco_cfg import setupProtonReco
      setupProtonReco(process,options.runProtonFastSim)
      toSchedule.append(process.pps_fastsim)
      toSchedule.append(process.pps_simulation_step)
      toSchedule.append(process.pps_reco_step)

      process.analysis.tagRecoProtons = cms.InputTag('ctppsProtonReconstructionOFDB')


process.ana=cms.Path(process.analysis)
toSchedule.append( process.ana )
if options.runOnData:

      if 'era2017' in options.era:
            process.analysis.tagRecoProtons = cms.InputTag('ctppsProtonReconstructionOFDB')
      else:
            process.analysis.tagRecoProtons = cms.InputTag('ctppsProtonReconstruction')

      #if options.runL1PrefireAna:
      #      print 'Prefire analysis is scheduled to be executed'
      #      from TopLJets2015.TopAnalysis.l1prefireAnalysis_cfi import *
      #      defineL1PrefireAnalysis(process,options.era)
      #      toSchedule.append(process.l1prefirePath)
print process.analysis.tagRecoProtons
                      
process.schedule=cms.Schedule( (p for p in toSchedule) )
print process.schedule
