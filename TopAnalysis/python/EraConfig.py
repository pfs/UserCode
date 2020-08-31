import os

def getEraConfiguration(era,isData,dataset=None):
    
    """ defines global tags, JEC/R corrections, etc. depending on the era """

    globalTags = {
        'era2016':('102X_mc2016_asymptotic_v8',       '102X_dataRun2_v13'),
        'era2017':('102X_mc2017_realistic_v8',       '102X_dataRun2_v13'),
        'era2018':('102X_upgrade2018_realistic_v21', '102X_dataRun2_Prompt_v16')
        }
    jecFiles    = {
        'era2016':('Summer16_07Aug2017_V11_MC', 'Summer16_07Aug2017All_V11_DATA', 'Summer16_07Aug2017_V11_MC_UncertaintySources_AK4PFchs'),
        'era2017':('Summer19UL17_V5_MC',        'Summer19UL17_RunBCDEF_V5_DATA',  'Summer19UL17_V5_MC_UncertaintySources_AK4PFchs'),
        'era2018':('Autumn18_V19_MC',           'Autumn18_RunABCD_V19_DATA',      'Autumn18_V19_MC_UncertaintySources_AK4PFchs')
        }
    jerFiles    = {
        'era2016':('Summer16_25nsV1_MC',   'Summer16_25nsV1_DATA'),
        'era2017':('Summer19UL17_JRV2_MC', 'Summer19UL17_JRV2_DATA'),
        'era2018':('Autumn18_V7b_MC',      'Autumn18_V7b_DATA'),
        }
    muonFiles   = {
        'era2016':'RoccoR2016.txt',
        'era2017':'RoccoR2017.txt',
        'era2018':'RoccoR2018.txt'
        }
    globalTag = globalTags[era][isData]
    jecFile   = jecFiles[era][isData]
    jecTag    = '_'.join( jecFile.split('_')[0:-1] )
    jecDB     = 'jec_DATA.db'  if isData else 'jec_MC.db'
    jerFile   = jerFiles[era][isData]
    jerTag    = '_'.join( jerFile.split('_')[0:-1] )
    jerDB     = 'jer_DATA.db'  if isData else 'jer_MC.db'
    qgDBFile  = 'QGL_AK4chs_94X.db'
    #ctppsDBFile= 'CTPPSRPAlignment_real_offline_v7.db'
    muonDBFile = muonFiles[era]

    #customize based on dataset
    if dataset:
        if isData and '2018D' in dataset:
            globalTag='102X_dataRun2_Prompt_v14'
            print '[Warn] Global tag has been overriden to',globalTag,'for',dataset,'in EraConfig.py'

    #copy correction files to a common CMSSW search path directory
    os.system('cp ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/%s/%s.db %s'%(era,jecFile,jecDB))
    os.system('cp ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/%s/%s.db %s'%(era,jecFile,jecDB))
    os.system('cp ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/%s/%s.db %s'%(era,jerFile,jerDB))
    os.system('cp ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/%s/%s.txt jecUncSources.txt'%(era,jecFiles[era][2]))
    os.system('cp ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/%s qg_db.db'%(qgDBFile))
    #os.system('cp ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/%s ctpps_db.db'%(ctppsDBFile))
    os.system('cp ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/%s/%s muoncorr_db.txt'%(era,muonDBFile))


    print 'JEC tag: ',jecTag,'to be read from',jecDB
    print 'JER tag: ',jerTag,'to be read from',jerDB
    print 'Muon corrections to be read from muoncorr_db.txt'
    print 'q/g discriminator to be read from qg_db.db'
    #print 'CTPPS config to be read from ctpps_db.db'

    return globalTag, jecTag, jecDB, jerTag, jerDB
    
