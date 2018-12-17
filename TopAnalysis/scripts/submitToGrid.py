import optparse
import os,sys
import json
import commands

"""
creates the crab cfg and submits the job
"""
def submitProduction(tag,lfnDirBase,dataset,isData,cfg,workDir,lumiMask,era='era2017',submit=False,addParents=False):
    
    from TopLJets2015.TopAnalysis.EraConfig import getEraConfiguration
    globalTag, jecTag, jecDB, jerTag, jerDB = getEraConfiguration(era=era,isData=bool(isData))

    #cmssw=os.environ['CMSSW_BASE']+'/src'

    os.system("rm -rvf %s/*%s* "%(workDir,tag))
    crabConfigFile=workDir+'/'+tag+'_cfg.py'
    config_file=open(crabConfigFile,'w')
    config_file.write('from WMCore.Configuration import Configuration\n')
    config_file.write('import os\n')
    config_file.write('config = Configuration()\n')
    config_file.write('\n')
    config_file.write('config.section_("General")\n')
    config_file.write('config.General.requestName = "%s"\n' % tag)
    config_file.write('config.General.workArea = "%s"\n' % workDir)
    config_file.write('config.General.transferOutputs=True\n')
    #config_file.write('config.General.transferLogs=True\n')
    config_file.write('\n')
    config_file.write('config.section_("JobType")\n')
    config_file.write('config.JobType.pluginName = "Analysis"\n')
    config_file.write('config.JobType.psetName = "'+cfg+'"\n')
    config_file.write('config.JobType.disableAutomaticOutputCollection = False\n')
    config_file.write('config.JobType.pyCfgParams = [\'runOnData=%s\',\'era=%s\']\n' % (bool(isData),era))
    #config_file.write('config.JobType.inputFiles = [\'{0}/{1}\',\'{0}/{2}\',\'{0}/muoncorr_db.txt\',\'{0}/jecUncSources.txt\']\n'.format(cmssw,jecDB,jerDB))
    config_file.write('config.JobType.inputFiles = [\'{0}\',\'{1}\',\'muoncorr_db.txt\',\'jecUncSources.txt\',\'qg_db.db\']\n'.format(jecDB,jerDB))
    config_file.write('\n')
    config_file.write('config.section_("Data")\n')
    config_file.write('config.Data.inputDataset = "%s"\n' % dataset)
    config_file.write('config.Data.inputDBS = "global"\n')
    config_file.write('config.Data.useParent = %s\n'% bool(addParents) )
    if isData :         
        config_file.write('config.Data.splitting = "Automatic"\n')
        #config_file.write('config.Data.splitting = "LumiBased"\n')
        #config_file.write('config.Data.unitsPerJob = 100\n')
        config_file.write('config.Data.lumiMask = \'%s\'\n' %lumiMask)
    else : 
        config_file.write('config.Data.splitting = "FileBased"\n')
        config_file.write('config.Data.unitsPerJob = 4\n')
     
    config_file.write('config.Data.ignoreLocality = False\n')    
    config_file.write('config.Data.publication = False\n')
    config_file.write('config.Data.outLFNDirBase = \"%s\"\n' % lfnDirBase)
    config_file.write('\n')
    config_file.write('config.section_("Site")\n')
    config_file.write('config.Site.storageSite = "T2_CH_CERN"\n')
    config_file.close()
    
    if submit : os.system('alias crab=\'/cvmfs/cms.cern.ch/crab3/crab-env-bootstrap.sh\' && crab submit -c %s' % crabConfigFile )

"""
steer the script
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-c', '--cfg',         dest='cfg'   ,      help='cfg to be sent to grid',       default=None,    type='string')
    parser.add_option(      '--addParents',  dest='addParents',  help='analyze paents as well',       default=False,   action='store_true')
    parser.add_option('-j', '--json',        dest='json'  ,      help='json with list of files',      default=None,    type='string')
    parser.add_option('-o', '--only',        dest='only'  ,      help='submit only these (csv)',      default=None,    type='string')
    parser.add_option('-l', '--lumi',        dest='lumiMask',    help='json with list of good lumis', default='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt')
    parser.add_option(      '--era',         dest='era',         help='era to use (sub-dir in data/)', default='era2017',  type='string')
    parser.add_option('-w', '--workDir',     dest='workDir',     help='working directory',             default='grid',     type='string')
    parser.add_option(      '--lfn',         dest='lfn',         help='base lfn to store outputs',     default='/store/group/cmst3/group/top/psilva', type='string')
    parser.add_option('-s', '--submit',      dest='submit',      help='submit jobs',                   default=False,   action='store_true')
    (opt, args) = parser.parse_args()

    #read list of samples
    jsonFile = open(opt.json,'r')
    samplesList=json.load(jsonFile,encoding='utf-8').items()
    jsonFile.close()

    githash=commands.getstatusoutput('git log --pretty=format:\'%h\' -n 1')[1]
    lfnDirBase='%s/%s/' % (opt.lfn,githash)
    
    onlyList=[]
    try:
        onlyList=opt.only.split(',')
    except:
        pass

    #submit jobs
    os.system("mkdir -p %s" % opt.workDir)
    for tag,sample in samplesList: 
        submit=False if len(onlyList)>0 else True
        for filtTag in onlyList:
            if filtTag in tag :
                submit=True
        if not submit : continue
        if len(sample[2])==0 : 
            print 'Ignoring',tag,'no dataset is set in the json file'
            continue 
        submitProduction(tag=tag,
                         lfnDirBase=lfnDirBase,
                         dataset=sample[2],
                         isData=sample[1],
                         lumiMask=opt.lumiMask,
                         cfg=opt.cfg,
                         era=opt.era,
                         workDir=opt.workDir,
                         submit=opt.submit,
                         addParents=opt.addParents)

    print 'crab cfg files have been created under %s' % opt.workDir
    if opt.submit:
        print 'Jobs have been submitted'
        print 'Monitor jobs with http://dashb-cms-job.cern.ch/dashboard/templates/task-analysis/ or locally with crab'
        print 'Output will be found in %s' % lfnDirBase

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
