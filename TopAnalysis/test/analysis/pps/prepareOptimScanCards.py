import ROOT
import os
import sys
import argparse
import itertools
from generateBinnedWorkspace import VALIDLHCXANGLES
import numpy as np

KINEMATICS = '(((cat==169 || cat==121) && l1pt>30 && l2pt>20 && bosonpt>40) || (cat==22 && bosonpt>95))'
RPSEL      = 'csi1>0.035 && csi2>0.035'
CATEGS     = [list(itertools.product( ['protonCat==%d'%protonCat],
                                      ['']+['xangle==%d'%x for x in VALIDLHCXANGLES],
                                      ['']+['nvtx<20','nvtx>=20','nch<15','njets==0']))
              for protonCat in [1,2,3,4] ]
CATEGS= [filter(None,list(y)) for x in CATEGS for y in x]
OPTIMLIST=[ ' && '.join( [KINEMATICS]+[RPSEL]+x ) for x in CATEGS]
MASSPOINTS=[ [600,660,720,780],
             [800,840,900,960],
             [1000,1020,1080,1140],
             [1200,1260,1320,1380],
             [1400,1440,1500,1560,1600],
             [1620],
             [1060,1100,1120,1160],
             [1180,1220,1240,1280],
             [1300,1340,1360,1420],
             [1460,1480,1520,1540,1580],
             [620,640,680,700],
             [740,760,820,860],
             [880,920,940,980],
]

def main(args):

    parser = argparse.ArgumentParser(description='usage: %prog [options]')
    parser.add_argument('-i', '--input',
                        dest='input',   
                        default='/eos/cms/store/cmst3/user/psilva/ExclusiveAna/final/ab05162/analysis_0p05',
                        help='input directory with the files [default: %default]')
    parser.add_argument('--signed',
                        dest='signed',
                        default=False,
                        help='used rapidity-signed missing mass [%default]',
                        action='store_true')
    parser.add_argument('-o', '--output',
                        dest='output', 
                        default='ppvx_analysis',
                        help='Output directory [default: %default]')
    parser.add_argument('--just',
                        dest='just',
                        default=None,
                        help='Run only this optimization points (CSV list) [default: %default]')
    opt=parser.parse_args(args)

    #build a list of the points to run
    if opt.just: opt.just=[int(x) for x in opt.just.split(',')]

    finalOptimList=[]
    for ipt,ana in enumerate(OPTIMLIST):
        for finalState in ['169','121','22']:
            finalOptimList.append( (ipt,ana,finalState) )

    optimJobs=[]
    combJobs=[]
    for ipt,ana,finalState in finalOptimList:

        if opt.just and not ipt in opt.just: continue

        workDir='%s/optim_%d'%(opt.output,ipt)
        os.system('mkdir -p ' + workDir)

        optimJobs.append((ipt,finalState))
        with open('%s/optimJob_%s.sh'%(workDir,finalState),'w') as script:
            script.write('#!/bin/bash\n')

            #initialization
            script.write('\n')
            script.write('doBackground=${1}\n')
            script.write('massList=${2}\n')
            script.write('input=%s\n'%opt.input)
            script.write('cuts="%s"\n'%ana)
            script.write('finalState="%s"\n'%finalState)

            if opt.signed:
                script.write('extraOpt=--signed\n')
            else:
                script.write('extraOpt=""\n')
        
            script.write('if [ "${doBackground}" == "1" ]; then\n')
            script.write('\t extraOpt="${extraOpt} --doBackground"\n')
            script.write('fi\n')

            script.write('if [ "${massList}" != "-1" ]; then\n')
            script.write('\t extraOpt="${extraOpt} --massList ${massList}"\n')
            script.write('fi\n')

            script.write('output=%s\n'%os.path.abspath(workDir))     
            script.write('\n')

            #environment
            script.write('echo "Setting up environment"\n')
            script.write('cd %s/src\n'%os.environ['CMSSW_BASE'])
            script.write('eval `scram r -sh`\n')
            script.write('cd -\n')
            script.write('\n')

            #create datacard
            script.write('echo "Running datacard creation"\n')
            deltaM=50
            if 'protonCat==1' in ana : deltaM=25
            if 'protonCat==4' in ana : deltaM=75
            script.write('python ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/pps/generateBinnedWorkspace.py -i ${input} -o ${output} --cuts "${cuts}" ${extraOpt} --finalState ${finalState} --mMin 0 --mMax 2000 --mBin %s\n'%(deltaM))
            script.write('\n')

    #submit optimization points to crab
    print 'Will submit %d optimization scan points'%len(OPTIMLIST)
    with open('%s/zxstatana_scan.sub'%opt.output,'w') as condor:
        condor.write("executable  = %s/optim_$(optimId)/optimJob_$(finalState).sh\n"%os.path.abspath(opt.output))
        condor.write("arguments   = $(doBackground) $(massList)\n")
        condor.write("output      = zxstatana_scan.out\n")
        condor.write("error       = zxstatana_scan.err\n")
        condor.write("log         = zxstatana_scan.log\n")
        condor.write("+JobFlavour = \"nextweek\"\n")
        condor.write("request_cpus = 4\n")

        for ipt,finalState in optimJobs:
            condor.write("\n")
            condor.write("optimId  = %d\n"%ipt)
            condor.write("finalState  = %s\n"%finalState)
            
            #do background
            condor.write('doBackground = 1\n')
            condor.write('massList     = -1\n')
            condor.write("queue 1\n\n")

            #do mass points
            condor.write('doBackground = 0\n')
            condor.write('queue massList from (\n') 
            for x in MASSPOINTS:
                xstr=','.join([ str(m) for m in x ])
                condor.write('\t%s\n'%xstr)
            condor.write(")\n")

    os.system('condor_submit %s/zxstatana_scan.sub'%opt.output)
    

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
