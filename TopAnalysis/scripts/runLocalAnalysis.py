import os
import sys
import optparse
import ROOT
import pickle
from collections import OrderedDict
import json
import re
import commands
from TopLJets2015.TopAnalysis.storeTools import *
from TopLJets2015.TopAnalysis.batchTools import *

"""
Wrapper to be used when run in parallel
"""
def RunMethodPacked(args):

    method,inF,outF,channel,charge,flag,runSysts,systVar,era,tag,debug,CR,QCDTemp,SRfake,mvatree,genWeights,xsec=args
    print args
    print 'Running ',method,' on ',inF
    print 'Output file',outF
    print 'Selection ch=',channel,' charge=',charge,' flag=',flag,' systs=',runSysts
    print 'Normalization applied from tag=',tag,'from',genWeights,'xsec=',xsec
    print 'Corrections will be retrieved for era=',era
    print 'Make the mva tree? ', mvatree
    print 'Make CR for photon fake rate? ', CR
    print 'Is QCD?', QCDTemp
    print 'Prepare the region to apply fake ratio?', SRfake

    try:
        cmd='analysisWrapper --era %s --normTag %s --in %s --out %s --method %s --charge %d --channel %d --flag %d --systVar %s --genWeights %s --xsec %f'\
            %(era, tag, inF, outF, method, charge, channel, flag, systVar,genWeights,xsec)
        if runSysts : cmd += ' --runSysts'
        if debug : cmd += ' --debug'
        if mvatree : cmd += ' --mvatree'
        if CR : cmd += ' --CR'
        if QCDTemp : cmd += ' --QCDTemp'
        if SRfake : cmd += ' --SRfake'
        print(cmd)
        os.system(cmd)

    except :
        print 50*'<'
        print "  Problem  (%s) with %s continuing without"%(sys.exc_info()[1],inF)
        print 50*'<'
        return False
    return True

"""
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-m', '--method',      dest='method',      help='method to run [%default]',                   default='TOP-16-006::RunTop16006',  type='string')
    parser.add_option('-i', '--in',          dest='input',       help='input directory with files or single file [%default]',  default=None,       type='string')
    parser.add_option('-o', '--out',         dest='output',      help='output directory (or file if single file to process)  [%default]',  default='analysis', type='string')
    parser.add_option(      '--only',        dest='only',        help='csv list of samples to process  [%default]',             default=None,       type='string')
    parser.add_option(      '--skip',        dest='skip',        help='csv list of samples to skip  [%default]',             default=None,       type='string')
    parser.add_option(      '--runSysts',    dest='runSysts',    help='run systematics  [%default]',                            default=False,      action='store_true')
    parser.add_option(      '--systVar',     dest='systVar',     help='single systematic variation  [%default]',   default='nominal',       type='string')
    parser.add_option(      '--debug',       dest='debug',       help='debug mode  [%default]',                            default=False,      action='store_true')
    parser.add_option(      '--flag',        dest='flag',        help='job specific flag  [%default]',   default=0,          type=int)
    parser.add_option(      '--xsec',        dest='xsec',        help='use this xsec value instead of the json one  [%default]',   default=None,          type=float)
    parser.add_option(      '--ch',          dest='channel',     help='channel  [%default]',                                    default=13,         type=int)
    parser.add_option(      '--charge',      dest='charge',      help='charge  [%default]',                                     default=0,          type=int)
    parser.add_option(      '--era',         dest='era',         help='era to use for corrections/uncertainties  [%default]',   default='era2016',       type='string')
    parser.add_option(      '--tag',         dest='tag',         help='normalize from this tag  [%default]',                    default=None,       type='string')
    parser.add_option('-q', '--queue',       dest='queue',       help='if not local send to batch with condor. queues are now called flavours, see http://batchdocs.web.cern.ch/batchdocs/local/submit.html#job-flavours   [%default]',     default='local',    type='string')    
    parser.add_option('-n', '--njobs',       dest='njobs',       help='# jobs to run in parallel  [%default]',                  default=0,    type='int')
    parser.add_option(      '--skipexisting',dest='skipexisting',help='skip jobs with existing output files  [%default]',       default=False,      action='store_true')
    parser.add_option(      '--exactonly',   dest='exactonly',   help='match only exact sample tags to process  [%default]',    default=False,      action='store_true')
    parser.add_option(      '--outputonly',  dest='outputonly',  help='filter job submission for a csv list of output files  [%default]',             default=None,       type='string')
    parser.add_option(      '--farmappendix',dest='farmappendix',help='Appendix to condor FARM directory [%default]',             default='',       type='string')
    parser.add_option(      '--genWeights',  dest='genWeights',  help='genWeights to get the normalization from (found within data/era directory) [%default]',             default='genweights.root',       type='string')
    parser.add_option(      '--mvatree',     dest='mvatree',     help='make mva tree  [%default]',                            default=False,      action='store_true'),
    parser.add_option(      '--CR',          dest='CR',          help='provide control region for photon FR  [%default]',                            default=False,      action='store_true')
    parser.add_option(      '--QCDTemp',     dest='QCDTemp',     help='provide template for fake photons [%default]',                            default=False,      action='store_true')
    parser.add_option(      '--SRfake',      dest='SRfake',     help='provide the region to apply the fake ratio [%default]',                            default=False,      action='store_true')
    (opt, args) = parser.parse_args()

    #parse selection lists
    onlyList=[]
    onlyListXsec={}
    try:
        for t in opt.only.split(','):
            if '.json' in t:
                print "The json file is %s" %t
                try:
                    myfile = open(t, 'r') # or "a+", whatever you need
                except IOError:
                    print "Could not open json file",t
                jsonFile = open(t,'r')
                samplesList = json.load(jsonFile, encoding='utf-8', object_pairs_hook=OrderedDict).items()
                for s,desc in samplesList: 
                    onlyList.append(s)
                    onlyListXsec[s]=desc[0]
                jsonFile.close()
            else:
                onlyList.append(t)
        print onlyListXsec
    except:
        pass
    skipList=[]
    try:
        skipList=opt.skip.split(',')
    except:
        pass
    outputOnlyList=[]
    try:
        outputOnlyList=opt.outputonly.split(',')
    except:
        pass
    #parse list of systematic variations
    varList=[]
    if opt.systVar == 'all':
        allSystVars = ['jec_CorrelationGroupMPFInSitu', 'jec_RelativeFSR',
                       'jec_CorrelationGroupUncorrelated', 'jec_FlavorPureGluon', 'jec_FlavorPureQuark',
                       'jec_FlavorPureCharm', 'jec_FlavorPureBottom', 'jer',
                       'btag_heavy', 'btag_light', 'csv_heavy', 'csv_light', 'tracking']
        for var in allSystVars:
            varList.append(var+'_up')
            varList.append(var+'_down')
    else:
        try:
            varList=opt.systVar.split(',')
        except:
            pass
    print 'Running following variations: ', varList

    #prepare output if a directory
    if not '.root' in opt.output :
        print opt.output
        if '/store/' in opt.output:
            os.system('eos mkdir %s'%opt.output)
            os.system('eos mkdir %s/Chunks'%opt.output)            
        else:
            os.system('mkdir -p %s/Chunks'%opt.output)

    #correct location of corrections to be used using cmsswBase, if needed
    cmsswBase=os.environ['CMSSW_BASE']
    if not cmsswBase in opt.era : opt.era=cmsswBase+'/src/TopLJets2015/TopAnalysis/data/'+opt.era

    #process tasks
    task_list = []
    processedTags=[]
    if '.root' in opt.input:
        inF=opt.input
        if '/store/' in inF and not 'root:' in inF : inF='root://eoscms//eos/cms'+opt.input              
        for systVar in varList:
            outF=opt.output
            xsec=1.
            if opt.tag in onlyListXsec: xsec=onlyListXsec[opt.tag]
            if opt.xsec: xsec=opt.xsec
            if systVar != 'nominal' and not systVar in opt.output: outF=opt.output[:-5]+'_'+systVar+'.root'
            task_list.append( (opt.method,inF,outF,opt.channel,opt.charge,opt.flag,opt.runSysts,systVar,opt.era,opt.tag,opt.debug, opt.CR, opt.QCDTemp, opt.SRfake, opt.mvatree,opt.genWeights,xsec) )
    else:

        inputTags=getEOSlslist(directory=opt.input,prepend='')
        for baseDir in inputTags:

            tag=os.path.basename(baseDir)
            if tag=='backup' : continue
            xsec=1.
            if tag in onlyListXsec: xsec=onlyListXsec[tag]
            if opt.xsec: xsec=opt.xsec

            #filter tags            
            if len(onlyList)>0:
                processThisTag=False
                for itag in onlyList:
                    if ((opt.exactonly and itag == tag) or 
                        (not opt.exactonly and itag in tag)):
                        processThisTag=True
                        break
                if not processThisTag : continue
            if len(skipList)>0:
                processThisTag=True
                for itag in skipList:
                    if itag in tag:
                        processThisTag=False
                        break
                if not processThisTag : continue

            input_list=getEOSlslist(directory='%s/%s' % (opt.input,tag) )
            
            for systVar in varList:
                nexisting = 0
                for ifile in xrange(0,len(input_list)):
                    inF=input_list[ifile]
                
                    outF=os.path.join(opt.output,'Chunks','%s_%d.root' %(tag,ifile))
                    if systVar != 'nominal' and not systVar in tag: outF=os.path.join(opt.output,'Chunks','%s_%s_%d.root' %(tag,systVar,ifile))
                    if (opt.skipexisting and os.path.isfile(outF)):
                        nexisting += 1
                        continue
                    if (len(outputOnlyList) > 1 and not outF in outputOnlyList):
                        continue
                    task_list.append( (opt.method,inF,outF,opt.channel,opt.charge,opt.flag,opt.runSysts,systVar,opt.era,tag,opt.debug, opt.CR, opt.QCDTemp, opt.SRfake, opt.mvatree,opt.genWeights,xsec) )
                if (opt.skipexisting and nexisting): print '--skipexisting: %s - skipping %d of %d tasks as files already exist'%(systVar,nexisting,len(input_list))

    #run the analysis jobs
    if opt.queue=='local':
        print 'launching %d tasks in %d parallel jobs'%(len(task_list),opt.njobs)
        if opt.njobs == 0:
            for args in task_list: RunMethodPacked(args)
        else:
            from multiprocessing import Pool
            pool = Pool(opt.njobs)
            pool.map(RunMethodPacked, task_list)
    else:
        
        FarmDirectory = '%s/FARM%s%s'%(cmsswBase,os.path.basename(opt.output),opt.farmappendix)
        os.system('mkdir -p %s'%FarmDirectory)
        
        print 'Preparing %d tasks to submit to the batch'%len(task_list)
        print 'Executables and condor wrapper are stored in %s'%FarmDirectory

        with open ('%s/condor.sub'%FarmDirectory,'w') as condor:

            condor.write('executable = {0}/$(cfgFile).sh\n'.format(FarmDirectory))
            condor.write('output     = {0}/output_$(cfgFile).out\n'.format(FarmDirectory))
            condor.write('error      = {0}/output_$(cfgFile).err\n'.format(FarmDirectory))
            condor.write('log        = {0}/output_$(cfgFile).log\n'.format(FarmDirectory))
            condor.write('+JobFlavour = "{0}"\n'.format(opt.queue))

            jobNb=0
            for method,inF,outF,channel,charge,flag,runSysts,systVar,era,tag,debug,CR,QCDTemp,SRfake,mvatree,genWeights,xsec in task_list:

                jobNb+=1
                cfgFile='%s'%(os.path.splitext(os.path.basename(outF))[0])

                condor.write('cfgFile=%s\n'%cfgFile)
                condor.write('queue 1\n')
                
                with open('%s/%s.sh'%(FarmDirectory,cfgFile),'w') as cfg:

                    cfg.write('#!/bin/bash\n')
                    cfg.write('WORKDIR=`pwd`\n')
                    cfg.write('echo "Working directory is ${WORKDIR}"\n')
                    cfg.write('cd %s\n'%cmsswBase)
                    cfg.write('eval `scram r -sh`\n')
                    cfg.write('cd ${WORKDIR}\n')
                    localOutF=os.path.basename(outF)
                    runOpts='-i %s -o ${WORKDIR}/%s --charge %d --ch %d --era %s --tag %s --flag %d --method %s --systVar %s --genWeights %s --xsec %f'\
                        %(inF, localOutF, charge, channel, era, tag, flag, method, systVar,genWeights,xsec)
                    if runSysts : runOpts += ' --runSysts'
                    if debug :    runOpts += ' --debug'
                    if mvatree :  runOpts += ' --mvatree'                    
                    if CR :       runOpts += ' --CR'
                    if QCDTemp :  runOpts += ' --QCDTemp'
                    if SRfake :  runOpts += ' --SRfake'
                    cfg.write('python %s/src/TopLJets2015/TopAnalysis/scripts/runLocalAnalysis.py %s\n'%(cmsswBase,runOpts))
                    if '/store' in outF:
                        cfg.write('xrdcp --force ${WORKDIR}/%s root://eoscms//%s\n'%(localOutF,outF))
                        cfg.write('rm ${WORKDIR}/%s'%localOutF)
                    elif outF!=localOutF:
                        cfg.write('  mv -v ${WORKDIR}/%s %s\n'%(localOutF,outF))

                os.system('chmod u+x %s/%s.sh'%(FarmDirectory,cfgFile))

        print 'Submitting jobs to condor, flavour "%s"'%(opt.queue)
        os.system('condor_submit %s/condor.sub'%FarmDirectory)
        

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
