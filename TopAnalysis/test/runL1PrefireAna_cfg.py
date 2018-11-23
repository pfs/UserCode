import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('runOnData', True,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Run this on real data"
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

# global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globalTag)


#EGM customization
from TopLJets2015.TopAnalysis.customizeEGM_cff import customizeEGM
process.egammaPostRecoSeq=customizeEGM(process=process,era=options.era)

#JetMET customization
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
from JetMETCorrections.Configuration.DefaultJEC_cff import *
from JetMETCorrections.Configuration.JetCorrectionServices_cff import *
from TopLJets2015.TopAnalysis.customizeJetTools_cff import *
customizeJetTools(process=process,
                  jecDB=jecDB,
                  jecTag=jecTag,
                  jerDB=jerDB,
                  jerTag=jerTag,
                  baseJetCollection='slimmedJets',
                  runOnData=options.runOnData)

#message logger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# set input to process
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('/store/data/Run2017F/SinglePhoton/MINIAOD/31Mar2018-v1/00000/0033B9DF-4338-E811-AB11-1CB72C1B6CC6.root'),
                            secondaryFileNames = cms.untracked.vstring(['/store/data/Run2017F/SinglePhoton/AOD/17Nov2017-v1/70000/BC8110F6-3CE0-E711-9210-02163E014564.root',
                                                                        '/store/data/Run2017F/SinglePhoton/AOD/17Nov2017-v1/70000/F68021F2-3CE0-E711-85E6-02163E014702.root',
                                                                        '/store/data/Run2017F/SinglePhoton/AOD/17Nov2017-v1/60000/90EE7B12-BCDE-E711-8C47-A4BF0112BC06.root',
                                                                        '/store/data/Run2017F/SinglePhoton/AOD/17Nov2017-v1/60000/F8CD978A-17DF-E711-A05F-A0369F836334.root',
                                                                        '/store/data/Run2017F/SinglePhoton/AOD/17Nov2017-v1/60000/3E1B045C-28DF-E711-8B81-002590A37114.root',
                                                                        '/store/data/Run2017F/SinglePhoton/AOD/17Nov2017-v1/60000/506AF751-29DF-E711-9E1F-001E67397F2B.root',
                                                                        '/store/data/Run2017F/SinglePhoton/AOD/17Nov2017-v1/60000/1CDC7A53-20DF-E711-A272-1CC1DE1CDD20.root',
                                                                        '/store/data/Run2017F/SinglePhoton/AOD/17Nov2017-v1/60000/28249C05-0DDF-E711-BA56-008CFA0527CC.root',
                                                                        '/store/data/Run2017F/SinglePhoton/AOD/17Nov2017-v1/60000/142F6B6F-28DF-E711-9D48-A4BF0112BDF8.root']),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck') 
                            )
print  "Processing",process.source.fileNames

#L1 pre-fire analysis
from TopLJets2015.TopAnalysis.l1prefireAnalysis_cfi import *
defineL1PrefireAnalysis(process,options.era)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(options.outFilename)
                                  )
#schedule execution
toSchedule=[]
if process.egammaPostRecoSeq:
      process.egm=cms.Path(process.egammaPostRecoSeq) 
      toSchedule.append( process.egm )
toSchedule.append(process.l1prefirePath)
process.schedule=cms.Schedule( (p for p in toSchedule) )
