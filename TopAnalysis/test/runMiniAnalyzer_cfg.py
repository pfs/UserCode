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
options.register('runWithAOD', False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "run with AOD"
                 )
options.register('redoProtonRecoFromRAW', False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "run proton reco from scratch"
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
options.register('globalTag', None,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Override with this global tag, if none use default"
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
from Configuration.StandardSequences.Eras import eras
process = cms.Process("MiniAnalysis", eras.Run2_2017)      

#get the configuration to apply
from TopLJets2015.TopAnalysis.EraConfig import getEraConfiguration
globalTag, jecTag, jecDB, jerTag, jerDB = getEraConfiguration(era=options.era,isData=options.runOnData)
if options.globalTag:
      globalTag=options.globalTag
      print 'Forcing global tag to',globalTag


# Load the standard set of configuration modules
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

#EGM customization
from TopLJets2015.TopAnalysis.customizeEGM_cff import customizeEGM
customizeEGM(process=process,era=options.era,runWithAOD=options.runWithAOD)

# global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globalTag)
print 'Global tag is',globalTag

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
                            fileNames = cms.untracked.vstring(),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck') 
                            )
      
if options.inputFile:
      import os

      fileList=[]
      if '.root' in options.inputFile :
            fileList=[options.inputFile]
      else:
            inDir=options.inputFile
            fileList = ['file:'+os.path.join(inDir,f) for f in os.listdir(inDir)]
      print 'Will run on',fileList
      process.source.fileNames = cms.untracked.vstring(fileList)     

      if options.secInputFile:
            secFileList=[]
            if '.root' in options.secInputFile:
                  secFileList=options.secInputFile.split(',')
            else:
                  inDir=options.secInputFile
                  secFileList = ['file:'+os.path.join(inDir,f) for f in os.listdir(inDir)]
            print 'Will run also on',secFileList
            process.source.secondaryFileNames = cms.untracked.vstring(secFileList)
else:
      #use standard test files
      from TopLJets2015.TopAnalysis.customizeInputFiles import *
      customTestInputFiles(process,options.era,options.runOnData,True if options.runWithAOD or options.redoProtonRecoFromRAW else False)

print  "Processing",process.source.fileNames
if hasattr(process.source,'secondaryFileNames'):
      print  "+",process.source.secondaryFileNames

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
elif 'era2018' in options.era:
      process.analysis.jetIdToUse=ANALYSISJETIDS[2018]
      process.analysis.triggersToUse=ANALYSISTRIGGERLISTS[2018]
      print '\t Using 2018 triggers/jet ids'
else:
      process.analysis.jetIdToUse=ANALYSISJETIDS[2016]
      process.analysis.triggersToUse=ANALYSISTRIGGERLISTS[2016]
      print '\t Using 2016 triggers/jet ids'
if options.runOnData:
      process.analysis.metFilterBits = cms.InputTag("TriggerResults","","RECO")
      print '\t will save met filter bits'
      process.analysis.tagRecoProtons = cms.InputTag('ctppsProtons','singleRP')      

#schedule execution
toSchedule=[]
try:
      toSchedule.append( process.egammaPostReco )
except:
      print '\te/g post reco not found, will take e/g as it is'
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

#if options.runOnData:
#      from TopLJets2015.TopAnalysis.protonReco_cfg import ctppsCustom
#      ctppsCustom(process,options.era)
#
#      if options.redoProtonRecoFromRAW:
#            process.load('EventFilter.CTPPSRawToDigi.ctppsRawToDigi_cff')
#            process.load('RecoCTPPS.Configuration.recoCTPPS_cff')
#            process.ppsReco=cms.Path(process.ctppsRawToDigi*process.recoCTPPS*process.ppsSeq)
#      else:
#            process.ppsReco=cms.Path(process.ppsSeq)
#      toSchedule.append(process.ppsReco)

if options.runProtonFastSim:
      from TopLJets2015.TopAnalysis.protonReco_cfg import setupProtonReco
      setupProtonReco(process,options.runProtonFastSim)
      toSchedule.append(process.pps_fastsim)
      toSchedule.append(process.pps_simulation_step)
      toSchedule.append(process.pps_reco_step)
      process.analysis.tagRecoProtons = cms.InputTag('ctppsProtonReconstructionOFDB')

process.ana=cms.Path(process.analysis)
toSchedule.append( process.ana )
                           
process.schedule=cms.Schedule( (p for p in toSchedule) )
print process.schedule
