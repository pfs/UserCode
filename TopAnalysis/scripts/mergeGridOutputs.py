import os
import sys
import optparse
from TopLJets2015.TopAnalysis.storeTools import *


def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--inDir',      dest='inDir',       help='input directory with files',  default=None,       type='string')
    parser.add_option('-o', '--outDir',     dest='outDir',      help='output directory with files', default=None,       type='string')
    parser.add_option('-s', '--chunkSize',  dest='chunkSize',   help='size of the output chunk (Gb) [%default]', default=2, type=float)
    parser.add_option(      '--only',       dest='only',        help='only this tag',               default=None    ,   type='string')
    parser.add_option(      '--farm',       dest='farm',        help='farm tag',                    default=None    ,   type='string')
    parser.add_option(      '--localProd',  dest='localProd',   help='local production',            default=False, action='store_true')
    parser.add_option('-q', '--queue',      dest='queue',       help='batch queue [%default]',      default='workday',  type='string')
    parser.add_option(      '--dry',        dest='dry',         help='do not submit/create output directory [%default]',
                            default=False, action='store_true')
    (opt, args) = parser.parse_args()

    #prepare output directory
    if opt.outDir is None: opt.outDir=opt.inDir

    cmsswBase=os.environ['CMSSW_BASE']
    FARMDIR='%s/INTEGRITYFARM'%cmsswBase
    if opt.farm: FARMDIR += opt.farm
    os.system('mkdir -p %s'%FARMDIR)
    os.system('rm  %s/*'%FARMDIR)
    print 'Submissions scripts @ ',FARMDIR

    onlyList=opt.only.split(',') if opt.only else []

    with open('%s/custom_merge.sh'%FARMDIR,'w') as shell:
        shell.write('#!/bin/bash\n')
        shell.write('WORKDIR=`pwd`\n')
        shell.write('outfile=${1}\n')
        shell.write('local_outfile=`basename ${outfile}`\n')
        shell.write('shift\n')
        shell.write('infiles=$*\n')
        shell.write('echo "Working directory is ${WORKDIR}"\n')
        shell.write('cd %s\n'%cmsswBase)
        shell.write('eval `scram r -sh`\n')
        shell.write('cd ${WORKDIR}\n')
        shell.write('echo "CMSSW_BASE=${CMSSW_BASE}"\n')
        shell.write('echo "Hadding files"\n')
        shell.write('echo "hadd -f -k ${local_outfile} ${infiles}"\n')
        shell.write('hadd -f -k ${local_outfile} ${infiles}\n')
        shell.write('echo "Copying to final destination and removing local output"\n')
        shell.write('xrdcp -f ${local_outfile} ${outfile}\n')
        shell.write('rm ${local_outfile}\n')
    os.system('chmod u+x %s/custom_merge.sh'%FARMDIR)

    condor=open('%s/condor.sub'%FARMDIR,'w')
    condor.write('executable = %s/custom_merge.sh\n'%FARMDIR)
    condor.write('output = %s/job_$(ProcId).out\n'%FARMDIR)
    condor.write('error  = %s/job_$(ProcId).err\n'%FARMDIR)
    condor.write('log    = %s/job_$(ProcId).log\n'%FARMDIR)
    condor.write('+JobFlavour ="%s"\n'%opt.queue)

    if not opt.localProd:
        dset_list=getEOSlslist(directory=opt.inDir,prepend='')
        for dset in dset_list:
            dsetname=dset.split('/')[-1]

            pub_list=getEOSlslist(directory=dset,prepend='')
            for pubDir in pub_list:

                if not 'crab' in pubDir:
                    print 'Ambiguity found @ <publication-name> for <primary-dataset>=%s , bailing out'%dsetname
                    continue
                pub=pubDir.split('/crab_')[-1]

                if len(onlyList)>0:
                    found=False
                    for tag in onlyList:
                        if not tag in pub: continue
                        found=True
                        break
                    if not found : 
                        print 'Skipping %s, not in process only list'%pub
                        continue

                #check if it's an extension
                pubExt=None
                try:
                    extSplit=pub.split('_ext')
                    pubExt='ext%d'%(len(extSplit)-1)
                    pub=extSplit[0]
                    print 'Extension will be postfixed with ',pubExt
                except:
                    print 'Core sample (no extension)'
                
                time_list=getEOSlslist(directory=pubDir,prepend='')
                if len(time_list)!=1:
                    print 'Ambiguity found @ <time-stamp> for <primary-dataset>=%s , bailing out'%dsetname
                    continue
                time_stamp=time_list[0].split('/')[-1]

                out_list=[]
                count_list=getEOSlslist(directory=time_list[0],prepend='')

                chunkList=getChunksInSizeOf(chunkSize=opt.chunkSize,directoryList=count_list,prepend='/eos/cms/')
                print pub,'will be hadded in',len(chunkList),'chunks of approx %fGb'%opt.chunkSize
                for ichunk in xrange(0,len(chunkList)):
                    outFile='/eos/cms/{0}/{1}/Chunk_{2}_{3}.root'.format(opt.outDir,pub,ichunk,pubExt)
                    condor.write('arguments = %s %s\n'%(outFile,' '.join(chunkList[ichunk])))
                    condor.write('queue 1\n')

                #prepare output directory
                if not opt.dry: os.system('mkdir -p /eos/cms/{0}/{1}'.format(opt.outDir,pub))

    else:
        print 'Local production'
        dset_list=getEOSlslist(directory=opt.inDir,prepend='')
        for dset in dset_list:
            pub=os.path.basename(dset)
            if len(onlyList)>0:
                found=False
                for tag in onlyList:
                    if not tag in pub: continue
                    found=True
                    break
                if not found : 
                    print 'Skipping %s, not in process only list'%pub
                    continue

            chunkList=getChunksInSizeOf(chunkSize=opt.chunkSize,directoryList=[dset],prepend='/eos/cms/')
            print pub,'will be hadded in',len(chunkList),'chunks of approx %fGb'%opt.chunkSize
            for ichunk in xrange(0,len(chunkList)):
                outFile='/eos/cms/{0}/{1}/Chunk_{2}_ext0.root'.format(opt.outDir,pub,ichunk)
                condor.write('arguments = %s %s\n'%(outFile,' '.join(chunkList[ichunk])))
                condor.write('queue 1\n')                

            if not opt.dry: os.system('mkdir -p /eos/cms/{0}/{1}'.format(opt.outDir,pub))

    condor.close()
    if not opt.dry:
        os.system('condor_submit %s/condor.sub'%FARMDIR)


if __name__ == "__main__":
    sys.exit(main())

