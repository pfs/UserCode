import ROOT
import os
import sys
import argparse
import itertools
"""
KINEMATICS = [('bosonpt>20', 'bosonpt>95'),
              ('bosonpt>30', 'bosonpt>100'),
              ('bosonpt>40', 'bosonpt>110'),
              ('bosonpt>50', 'bosonpt>120'),
              ('bosonpt>60', 'bosonpt>130')]
RPSEL      = ['csi1>0.05 && csi2>0.05',
              'csi1>0.05 && csi2>0.055',
              'csi1>0.055 && csi2>0.05',
              'csi1>0.06 && csi2>0.06']
CATEGS     = ['nvtx>15',
              'nvtx<15,nvtx>=15',
              'nvtx<20,nvtx>=20',
              'nvtx<25,nvtx>=25',
              'nvtx<30,nvtx>=30']
"""


KINEMATICS = [('bosonpt>30',                       'bosonpt>95'),
              ('bosonpt>40',                       'bosonpt>100'),
              ('bosonpt>40 && l1pt>40 && l2pt>30', 'bosonpt>110'),             
              ('bosonpt>50',                       'bosonpt>120'),
              ('bosonpt>50 && l1pt>40 && l2pt>30', 'bosonpt>130'),             
              ('bosonpt>60',                       'bosonpt>140')]
RPSEL      = ['csi1>0.04 && csi2>0.04',
              'csi1>0.05 && csi2>0.05',              
              'csi1>0.06 && csi2>0.06']
CATEGS     = ['nvtx<20,nvtx>=20',
              'nvtx<25,nvtx>=25',
              'nvtx<30,nvtx>=30']

OPTIMLIST=list(itertools.product(KINEMATICS, RPSEL,CATEGS))

def main(args):

    parser = argparse.ArgumentParser(description='usage: %prog [options]')
    parser.add_argument('-i', '--input',
                        dest='input',   
                        default='/eos/cms/store/cmst3/user/psilva/ExclusiveAna/final/ab05162/analysis_0p05',
                        help='input directory with the files [default: %default]')
    parser.add_argument('--injectMass',
                        dest='injectMass',
                        default=None,
                        help='mass to inject in pseudo-data [%default]')
    parser.add_argument('-o', '--output',
                        dest='output', 
                        default='ppvx_analysis',
                        help='Output directory [default: %default]')
    parser.add_argument('--just',
                        dest='just',
                        default=None,
                        help='Run only this optimization points (CSV list) [default: %default]')
    parser.add_argument('--unblind',
                        dest='unblind', 
                        default=False,
                        action='store_true',
                        help='Use non-mixed data in the final fit [default: %default]')
    opt=parser.parse_args(args)

    #build a list of the points to run
    if opt.just: opt.just=[int(x) for x in opt.just.split(',')]

    ipt=0
    n2sub=[]
    for ana in OPTIMLIST:
        ipt+=1

        if opt.just and not ipt in opt.just: continue
        n2sub.append(ipt)

        kin,rpsel,cats=ana

        workDir='%s/optim_%d'%(opt.output,ipt)
        print workDir

        os.system('mkdir -p ' + workDir)
        with open('%s/optimJob.sh'%workDir,'w') as script:
            script.write('#!/bin/bash\n')

            #initialization
            script.write('\n')
            script.write('input=%s\n'%opt.input)
            script.write('preselZ="%s && %s"\n'%(kin[0],rpsel))
            script.write('preselGamma="%s && %s"\n'%(kin[1],rpsel))        
            script.write('categs="%s"\n'%(cats))
            if opt.injectMass:
                script.write('injectMass="--injectMass %s"\n'%opt.injectMass)
            else:
                script.write('injectMass=""\n')
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
            script.write('python ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/pps/generateBinnedWorkspace.py -i ${input} -o ${output} --preselZ "${preselZ}" --preselGamma "${preselGamma}" --categs "${categs}" ${injectMass}\n')
            script.write('\n')

            #combine cards
            script.write('echo "Combining datacards"\n')
            script.write('cd $output\n')        
            script.write('for v in z g; do\n')
            script.write('\ta=(`ls shapes-parametric.datacard_${v}*.dat`)\n')
            script.write('\tcombstr=""\n')
            script.write('\tncards=${#a[@]}\n')
            script.write('\tfor (( i=0; i<${ncards}; i++ ));\n')
            script.write('\tdo\n')
            script.write('\t\tcombstr="${combstr} ${v}cat${i}=${a[$i]}"\n')
            script.write('\tdone\n') 
            script.write('\tcombineCards.py $combstr > ${v}_datacard.dat\n')
            script.write('done\n')
            script.write('\n')
        
            #run statistical analysis
            script.write('echo "Running combine"\n')
            script.write('for b in z g; do\n')
            script.write('\t for m in 780 800 840 900 960 1000 1020 1080 1140 1200 1260 1320 1380 1400 1440 1500 1560 1600; do \n')
            script.write('\t\t pfix=${b}_m${m}\n')
            script.write('\t\t text2workspace.py ${b}_datacard.dat -m ${m} -o ${pfix}_workspace.root\n')
            script.write('\t\t combine ${pfix}_workspace.root -n PP${b}X.obs -M AsymptoticLimits -m ${m}\n')
            script.write('\t\t combine ${pfix}_workspace.root -n PP${b}X     -M AsymptoticLimits -m ${m} -t -1 --expectSignal=1 --setParameters mu_outfidsig=1\n')
            script.write('\t\t combine ${pfix}_workspace.root -n PP${b}X.obs -M Significance -m ${m}\n')
            script.write('\t\t combine ${pfix}_workspace.root -n PP${b}X     -M Significance -m ${m} -t -1 --expectSignal=1 --setParameters mu_outfidsig=1\n')
            script.write('\t\t #combine ${pfix}_workspace.root -n PP${b}X -M FitDiagnostics -m ${m} -t -1 --expectSignal=1 --setParameters mu_outfidsig=1\n')
            script.write('\t done\n')
            script.write('done\n')
            script.write('cd -\n')        

    #submit optimization points to crab
    print 'Will submit %d optimization scan points'%len(n2sub)
    with open('zxstatana_scan.sub','w') as condor:
        condor.write("executable  = %s/optim_$(point)/optimJob.sh\n"%os.path.abspath(opt.output))
        condor.write("output       = zxstatana_scan.out\n")
        condor.write("error        = zxstatana_scan.err\n")
        condor.write("log          = zxstatana_scan.log\n")
        condor.write("+JobFlavour = \"tomorrow\"\n")
        condor.write("request_cpus = 4\n")
        for i in n2sub:
            condor.write("point=%d\n"%i)
            condor.write("queue 1\n")
    os.system('condor_submit zxstatana_scan.sub')


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
