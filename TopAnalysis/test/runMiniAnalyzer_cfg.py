import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('runOnData', False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Run this on real data"
                 )
options.register('outFilename', 'MiniEvents.root',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Output file name"
                 )
options.register('inputDir', '',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "input directory with files to process"
                 )
options.register('saveTree', True,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "save summary tree"
                 )
options.register('savePF', True,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 'save PF candidates'
                 )
options.parseArguments()

process = cms.Process("MiniAnalysis")


# Load the standard set of configuration modules
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')


# global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_Prompt_ICHEP16JEC_v0' if options.runOnData else '80X_mcRun2_asymptotic_2016_miniAODv2_v1')

#message logger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


# set input to process
from TopLJets2015.TopAnalysis.Compressed_T2tt_cfi import *
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",
                            #fileNames = cms.untracked.vstring('/store/mc/RunIISpring16MiniAODv1/TTToSemiLeptonic_13TeV_ScaleDown-powheg/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/067717F2-220A-E611-A553-0090FAA58D84.root')
                            fileNames = cms.untracked.vstring(Compressed_T2tt_2),
                            
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck') 
                            )
if options.runOnData:
    process.source.fileNames = cms.untracked.vstring('/store/data/Run2016B/SingleElectron/MINIAOD/PromptReco-v2/000/273/158/00000/0A7BD549-131A-E611-8287-02163E0134FC.root')
    #process.source.fileNames = cms.untracked.vstring('root://cmseos.fnal.gov//store/data/Run2015D/SingleMuon/MINIAOD/16Dec2015-v1/10000/00006301-CAA8-E511-AD39-549F35AD8BC9.root')
    #/store/data/Run2016B/SingleElectron/MINIAOD/PromptReco-v2/000/273/158/00000/06277EC1-181A-E611-870F-02163E0145E5.root')

#this make the process crash ?!
#if options.inputDir!='': 
#    from TopLJets2015.TopAnalysis.storeTools import getEOSlslist 
#    print 'Will process files from',options.inputDir
#    process.source.fileNames=cms.untracked.vstring(getEOSlslist(directory=options.inputDir))


#analysis
process.load('TopLJets2015.TopAnalysis.miniAnalyzer_cfi')
if not options.saveTree:
    print 'Summary tree won\'t be saved'
    process.analysis.saveTree=cms.bool(False)
if not options.savePF:
    print 'Summary PF info won\'t be saved'
    process.analysis.savePF=cms.bool(False)

#pseudo-top
if not options.runOnData:
    process.load('TopQuarkAnalysis.TopEventProducers.producers.pseudoTop_cfi')
    process.pseudoTop.leptonMinPt=cms.double(20)
    process.pseudoTop.leptonMaxEta=cms.double(2.5)
    process.pseudoTop.jetMaxEta=cms.double(5.0)

#EGM smearer https://twiki.cern.ch/twiki/bin/view/CMS/EGMSmearer
process.load('Configuration.StandardSequences.Services_cff')
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                   calibratedPatElectrons  = cms.PSet( initialSeed = cms.untracked.uint32(81),
                                                                                       engineName = cms.untracked.string('TRandom3'),
                                                                                       )
                                                   )
process.load('EgammaAnalysis.ElectronTools.calibratedElectronsRun2_cfi')
process.calibratedPatElectrons.correctionFile = 'EgammaAnalysis/ElectronTools/data/ScalesSmearings/80X_ichepV1_2016_ele'
print 'Using calibrated electrons with corrections from',process.calibratedPatElectrons.correctionFile
process.calibratedPatElectrons.isMC=False if options.runOnData else True

# Set up electron ID (VID framework)
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_Trig_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff'] 
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)


#jet energy corrections
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
from JetMETCorrections.Configuration.DefaultJEC_cff import *
from JetMETCorrections.Configuration.JetCorrectionServices_cff import *
from TopLJets2015.TopAnalysis.customizeJetTools_cff import *
jecLevels=['L1FastJet','L2Relative','L3Absolute']
jecFile='Spring16_25nsV6_MC.db'
jecTag='Spring16_25nsV6_MC_AK4PFchs'
if options.runOnData : 
    print 'Warning we\'re still using Spring16 MC corrections for data - to be updated'
    jecLevels.append( 'L2L3Residual' )
    jecFile='Spring16_25nsV6_DATA.db'
    jecTag='Spring16_25nsV6_DATA_AK4PFchs'
customizeJetTools(process=process,jecLevels=jecLevels,jecFile=jecFile,jecTag=jecTag)

#tfile service
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('histo.root')
                                   )

if options.runOnData:
    process.p = cms.Path(process.calibratedPatElectrons*process.egmGsfElectronIDSequence*process.customizeJetToolsSequence*process.analysis)
else:
    process.p = cms.Path(process.calibratedPatElectrons*process.egmGsfElectronIDSequence*process.customizeJetToolsSequence*process.pseudoTop*process.analysis)
