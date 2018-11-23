from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

#from WMCore.Configuration import Configuration
#config = Configuration()

config.section_("General")
config.General.requestName = 'SingleMuonv2'

#config.General.requestName = 'DYJets_ToLL_M50_SEP10'

config.General.workArea = 'crabprojects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
#config.JobType.inputFiles = ["UserCode/TopAnalysis/jec/Summer15_50nsV4_DATA.db"]
#config.JobType.inputFiles = ["UserCode/TopAnalysis/jec/Summer15_50nsV4_MC.db"]
config.JobType.psetName = 'runMiniAnalyzer_cfg.py'
config.JobType.pyCfgParams = ["runOnData=True"]

config.section_("Data")
config.Data.inputDataset = '/SingleMuon/Run2017B-PromptReco-v1/MINIAOD'

#config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3/MINIAODSIM'

config.Data.inputDBS = 'global'
#config.Data.splitting = 'FileBased'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 20
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-300575_13TeV_PromptReco_Collisions17_JSON.txt'
config.Data.runRange =  '294927-300575' # '193093-194075'
config.Data.publication = False
#config.Data.publishDbsUrl = 'phys03'
#config.Data.publishDataName = 'CSA14_wjets'

config.section_("Site")
config.Site.blacklist = ['T2_CH_CERN']
config.Site.storageSite = 'T2_IT_Legnaro'
