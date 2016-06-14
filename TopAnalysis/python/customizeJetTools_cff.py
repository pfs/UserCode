import FWCore.ParameterSet.Config as cms

def customizeJetTools(process,jecLevels,jecFile,jecTag):

	process.load("CondCore.DBCommon.CondDBCommon_cfi")
	from CondCore.DBCommon.CondDBSetup_cfi import CondDBSetup
	process.jec = cms.ESSource("PoolDBESSource",
				   DBParameters = cms.PSet(
			messageLevel = cms.untracked.int32(0)
			),
				   timetype = cms.string('runnumber'),
				   toGet = cms.VPSet(
			cms.PSet(
				record = cms.string('JetCorrectionsRecord'),
				tag    = cms.string('JetCorrectorParametersCollection_%s'%jecTag),
				label  = cms.untracked.string('AK4PFchs')
				),


      	## here you add as many jet types as you need
      	## note that the tag name is specific for the particular sqlite file 
     	 ), 
      	connect = cms.string('sqlite_file:%s'%jecFile)
	)
	## add an es_prefer statement to resolve a possible conflict from simultaneous connection to a global tag
	process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

	from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetCorrFactorsUpdated
	process.patJetCorrFactorsReapplyJEC = patJetCorrFactorsUpdated.clone(
		src = cms.InputTag("slimmedJets"),
		levels = jecLevels,
		payload = 'AK4PFchs' ) # Make sure to choose the appropriate levels and payload here!
	
	from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetsUpdated	
	process.slimmedJetsReapplyJEC = patJetsUpdated.clone(
		jetSource = cms.InputTag("slimmedJets"),
		jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
  	)

	#sequence to include
	process.customizeJetToolsSequence = cms.Sequence(process.patJetCorrFactorsReapplyJEC 
							 + process. slimmedJetsReapplyJEC)
