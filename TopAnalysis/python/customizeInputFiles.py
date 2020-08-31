import FWCore.ParameterSet.Config as cms

def customTestInputFiles(process,era,runOnData,runWithAOD):
    if '2016' in era:
        if runOnData:
            process.source.fileNames = cms.untracked.vstring('/store/data/Run2016B/SingleMuon/MINIAOD/17Jul2018_ver2-v1/00000/0219DD4A-7D8C-E811-B7EB-001E67248A43.root')
            if runWithAOD:
                print 'Adding secondary filenames'
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
            else:
                process.source.fileNames = cms.untracked.vstring('/store/mc/RunIISummer16MiniAODv3/ST_t-channel_antitop_4f_mtop1715_inclusiveDecays_13TeV-powhegV2-madspin-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/120000/16CEB785-3FE6-E811-AAE8-FA163E9D74F8.root')
    elif '2017' in era:
        if runOnData:
            process.source.fileNames = cms.untracked.vstring('/store/data/Run2017F/DoubleMuon/MINIAOD/09Aug2019_UL2017-v1/270000/E6579C3C-A837-6045-AEAF-8957079DF1AE.root')
            if runWithAOD:
                print 'Adding secondary filenames'
                process.source.secondaryFileNames = cms.untracked.vstring([
                    '/store/data/Run2017F/DoubleMuon/RAW/v1/000/305/282/00000/E8E6CF6E-F9B4-E711-AC84-02163E01360E.root',
                    '/store/data/Run2017F/DoubleMuon/RAW/v1/000/305/282/00000/A6EB9360-F9B4-E711-BB38-02163E011C66.root',
                    '/store/data/Run2017F/DoubleMuon/RAW/v1/000/305/282/00000/74671608-F9B4-E711-B718-02163E0144E3.root',
                    '/store/data/Run2017F/DoubleMuon/RAW/v1/000/305/282/00000/288C227D-F9B4-E711-A614-02163E014105.root',
                    '/store/data/Run2017F/DoubleMuon/RAW/v1/000/305/282/00000/18AD8180-F9B4-E711-A43A-02163E019DB5.root',
                    '/store/data/Run2017F/DoubleMuon/RAW/v1/000/305/282/00000/14C48D98-F9B4-E711-A9D9-02163E014732.root',
                    '/store/data/Run2017F/DoubleMuon/RAW/v1/000/305/112/00000/FAC0BFBE-7EB3-E711-ACE2-02163E011FB0.root',
                    '/store/data/Run2017F/DoubleMuon/RAW/v1/000/305/112/00000/F820E455-7EB3-E711-8169-02163E0141FF.root',
                    '/store/data/Run2017F/DoubleMuon/RAW/v1/000/305/112/00000/E655E2E3-7EB3-E711-A106-02163E019C3D.root',
                    '/store/data/Run2017F/DoubleMuon/RAW/v1/000/305/112/00000/BA6AC59A-7EB3-E711-B9F8-02163E01A450.root',
                    '/store/data/Run2017F/DoubleMuon/RAW/v1/000/305/112/00000/ACE01C8B-7EB3-E711-9F27-02163E0138A4.root',
                    '/store/data/Run2017F/DoubleMuon/RAW/v1/000/305/112/00000/AC631519-7FB3-E711-9B94-02163E01A749.root',
                    '/store/data/Run2017F/DoubleMuon/RAW/v1/000/305/112/00000/AA00D009-7FB3-E711-9E0B-02163E019D2B.root',
                    '/store/data/Run2017F/DoubleMuon/RAW/v1/000/305/112/00000/A40CBDC1-7EB3-E711-9D3E-02163E01341F.root',
                    '/store/data/Run2017F/DoubleMuon/RAW/v1/000/305/112/00000/96FFC79F-7EB3-E711-94CB-02163E019BF3.root',
                    '/store/data/Run2017F/DoubleMuon/RAW/v1/000/305/112/00000/6A48BD33-7EB3-E711-AD33-02163E01A27C.root',
                    '/store/data/Run2017F/DoubleMuon/RAW/v1/000/305/112/00000/50180243-7EB3-E711-A516-02163E01A4DD.root'
                ])
        else:
            process.source.fileNames = cms.untracked.vstring('/store/mc/RunIISummer19UL17MiniAOD/TTTo2L2Nu_TuneCP5CR1_13TeV-powheg-pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v1/10000/04E43816-D0FA-F346-9308-1481C9C2CA99.root')
    elif '2018' in era:
        if runOnData:
            process.source.fileNames = cms.untracked.vstring('/store/data/Run2018C/SingleMuon/MINIAOD/17Sep2018-v1/00000/AB61FB4F-3A42-4B4F-93E2-78CD7E7CF0A4.root')
            if runWithAOD:
                print 'Adding secondary filenames'
                process.source.secondaryFileNames = cms.untracked.vstring([
                    '/store/data/Run2018C/SingleMuon/RAW/v1/000/320/065/00000/C45A09EB-328E-E811-BFEC-FA163EDD3A25.root',
                    '/store/data/Run2018C/SingleMuon/RAW/v1/000/320/065/00000/0A3EDD77-128E-E811-8885-FA163E08D4C1.root',
                    '/store/data/Run2018C/SingleMuon/RAW/v1/000/320/065/00000/1E0CFF73-128E-E811-9B2B-FA163E6DCA76.root',
                    '/store/data/Run2018C/SingleMuon/RAW/v1/000/320/065/00000/9659EBAF-108E-E811-93CB-FA163E82172A.root',
                    '/store/data/Run2018C/SingleMuon/RAW/v1/000/320/025/00000/F05577B5-568D-E811-A9FA-FA163E513903.root',
                    '/store/data/Run2018C/SingleMuon/RAW/v1/000/320/026/00000/26D29DBA-5F8D-E811-9D7B-FA163E600F07.root',
                    '/store/data/Run2018C/SingleMuon/RAW/v1/000/320/026/00000/9239C9D4-5F8D-E811-BD2C-FA163ECED2A8.root',
                    '/store/data/Run2018C/SingleMuon/RAW/v1/000/320/065/00000/BEED33AC-328E-E811-8305-FA163EDE7DC6.root',
                    '/store/data/Run2018C/SingleMuon/RAW/v1/000/320/025/00000/C662E202-578D-E811-8E1E-FA163E7FB452.root',
                    '/store/data/Run2018C/SingleMuon/RAW/v1/000/320/024/00000/F4B08B5F-528D-E811-965A-FA163ED40F48.root',
                    '/store/data/Run2018C/SingleMuon/RAW/v1/000/320/024/00000/0863DA3E-548D-E811-B8E1-A4BF01277823.root',
                    '/store/data/Run2018C/SingleMuon/RAW/v1/000/320/009/00000/0218B251-F28C-E811-81CD-FA163ED486E3.root',
                    '/store/data/Run2018C/SingleMuon/RAW/v1/000/320/010/00000/768D18BE-F98C-E811-A332-FA163E60ED2F.root',
                    '/store/data/Run2018C/SingleMuon/RAW/v1/000/320/007/00000/D208B88E-E98C-E811-9B7F-FA163E61ADEA.root',
                    '/store/data/Run2018C/SingleMuon/RAW/v1/000/320/007/00000/DAA30463-E98C-E811-918B-FA163EF70A72.root',
                    '/store/data/Run2018C/SingleMuon/RAW/v1/000/320/010/00000/7A33E7D8-F78C-E811-AF32-FA163E3AFBF6.root',
                    '/store/data/Run2018C/SingleMuon/RAW/v1/000/319/854/00000/A4D6882F-768A-E811-9978-FA163E11A4DA.root',
                    '/store/data/Run2018C/SingleMuon/RAW/v1/000/320/007/00000/B4B0ADA6-E88C-E811-B65F-FA163E105698.root'
                ])
        else:
            process.source.fileNames = cms.untracked.vstring('/store/mc/RunIIAutumn18MiniAOD/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15_ext2-v1/50001/FD88E288-BFD2-F444-AD4C-844482719AA2.root')
            
