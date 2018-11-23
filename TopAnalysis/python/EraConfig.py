import os

def getEraConfiguration(era,isData):
    
    """ defines global tags, JEC/R corrections, etc. depending on the era """

    globalTags = {
        'era2016':('94X_mcRun2_asymptotic_v3', '94X_dataRun2_v10'),
        'era2017':('94X_mc2017_realistic_v14', '94X_dataRun2_v6')
        }
    jecFiles    = {
        'era2016':('Summer16_07Aug2017_V11_MC',  'Summer16_07Aug2017All_V11_DATA', 'Summer16_07Aug2017_V11_MC_UncertaintySources_AK4PFchs'),
        'era2017':('Fall17_17Nov2017_V6_MC',     'Fall17_17Nov2017BCDEF_V6_DATA',  'Fall17_17Nov2017_V8_MC_UncertaintySources_AK4PFchs')
        }
    jerFiles    = {
        'era2016':('Summer16_25nsV1_MC',         'Summer16_25nsV1_DATA'),
        'era2017':('Summer16_25nsV1_MC',         'Summer16_25nsV1_DATA'),
        }
    muonFiles   = {
        'era2016':'RoccoR2016.txt',
        'era2017':'RoccoR2017.txt'

        }
    globalTag = globalTags[era][isData]
    jecFile   = jecFiles[era][isData]
    jecTag    = '_'.join( jecFile.split('_')[0:-1] )
    jecDB     = 'jec_DATA.db'  if isData else 'jec_MC.db'
    jerFile   = jerFiles[era][isData]
    jerTag    = '_'.join( jerFile.split('_')[0:-1] )
    jerDB     = 'jer_DATA.db'  if isData else 'jer_MC.db'

    #copy correction files to a common CMSSW search path directory
    os.system('cp ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/%s/%s.db %s'%(era,jecFile,jecDB))
    os.system('cp ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/%s/%s.db %s'%(era,jerFile,jerDB))
    os.system('cp ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/%s/%s.txt jecUncSources.txt'%(era,jecFiles[era][2]))
    os.system('cp ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/%s/%s muoncorr_db.txt'%(era,muonFiles[era]))

    print 'JEC tag: ',jecTag,'to be read from',jecDB
    print 'JER tag: ',jerTag,'to be read from',jerDB
    print 'Muon corrections to be read from muoncorr_db.txt'

    return globalTag, jecTag, jecDB, jerTag, jerDB
    
