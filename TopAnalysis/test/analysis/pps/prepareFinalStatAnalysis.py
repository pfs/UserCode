import ROOT
import os
import sys
import argparse
import itertools
import re
from prepareOptimScanCards import OPTIMLIST, MASSPOINTS
from collections import defaultdict


def main(args):

    parser = argparse.ArgumentParser(description='usage: %prog [options]')
    parser.add_argument('-i', '--input',
                        dest='input',   
                        default=None,
                        help='input directory with the files [default: %default]')
    parser.add_argument('-t', '--tag',
                        dest='tag',   
                        default='exp',
                        help='cards tag [default: %default]')
    parser.add_argument('--cmssw',
                        dest='cmssw',   
                        default='/afs/cern.ch/work/p/psilva/HIN-19-001/CMSSW_10_2_13/src/',
                        help='input directory with the files [default: %default]')
    opt=parser.parse_args(args)

    isSingleRP=True if 'excsingle' in opt.input else False


    statAna=defaultdict(list)
    for i in range(len(OPTIMLIST)):
        ana=OPTIMLIST[i]
        protonCat      = int(re.search('protonCat==(\d+)',ana).group(1))
        pcatRec        = 'pcat%d'%protonCat
        requiresXangle      = True if 'xangle' in ana else False
        requiresNvtx        = True if 'nvtx'   in ana else False
        requiresNjets       = True if 'njets'  in ana else False
        requiresNch         = True if 'nch'    in ana else False
        
        if isSingleRP and protonCat!=4 : continue

        if not requiresXangle and not requiresNvtx and not requiresNch and not requiresNjets:
            statAna['inc'].append(i)
            statAna[pcatRec].append(i)

        if     requiresXangle and not requiresNvtx and not requiresNch and not requiresNjets:
            statAna['inc_xangle'].append(i)
            statAna[pcatRec+'_xangle'].append(i)

        if not requiresXangle and     requiresNvtx and not requiresNch and not requiresNjets:
            statAna['inc_nvtx'].append(i)
            statAna[pcatRec+'_nvtx'].append(i)

        if     requiresXangle and     requiresNvtx and not requiresNch and not requiresNjets:
            statAna['inc_xangle_nvtx'].append(i)
            statAna[pcatRec+'_xangle_nvtx'].append(i)

        if not requiresXangle and not requiresNvtx and     requiresNch and not requiresNjets:
            pfix='loose' if '15' in ana else '10'
            statAna['inc_nch'+pfix].append(i)
            statAna[pcatRec+'_nch'+pfix].append(i)

        if not requiresXangle and not requiresNvtx and not requiresNch and requiresNjets:
            statAna['inc_jveto'].append(i)
            statAna[pcatRec+'_jveto'].append(i)

        if requiresXangle and not requiresNch and requiresNjets:
            statAna['inc_xangle_jveto'].append(i)
            statAna[pcatRec+'_xangle_jveto'].append(i)

    for key,optimPts in statAna.items():

        odir='%s/%s/%s'%(opt.input,opt.tag,key)
        os.system('mkdir -p '+odir)
        
        #a human-readable summary
        baseDataCards=[]
        with open(os.path.join(odir,'info.dat'),'w') as f:
            f.write(key+' analysis\n')
            f.write('tag={}\n'.format(opt.tag))
            f.write('-'*50+'\n')
            for ioptim in optimPts:
                icat=len(baseDataCards)
                f.write('%3d == %3d %s\n'%(icat,ioptim,OPTIMLIST[ioptim]))
                baseDataCards.append('../../optim_%d/shapes-parametric.datacard_{ch}_%s.dat'%(ioptim,opt.tag))
            f.write('-'*50+'\n')
            f.write('Final statistical analysis will be composed for the different categories from\n')
            f.write(' '.join(baseDataCards))
        
        #the steering script for combine
        with open(os.path.join('%s/statAnaJob.sh'%odir),'w') as script:

            script.write('#!/bin/bash\n')

            #environment
            script.write('echo "Setting up environment"\n')
            script.write('cd %s\n'%opt.cmssw)
            script.write('eval `scram r -sh`\n')
            script.write('output=%s\n'%os.path.abspath(odir))     
            script.write('cd ${output}\n')
            script.write('\n')

            script.write('m=${1}\n')
            script.write('b=${2}\n')
            script.write('pfix=${b}_m${m}\n')
            script.write('echo "Running combine for b=${2} m=${m} with CMSSW_BASE=${CMSSW_BASE}"\n')
            script.write('\n')

            cardsToCombine   = ' '.join(baseDataCards).format(ch='${b}')
            zCardsToCombine  = ' '.join(baseDataCards).format(ch='zee')
            zCardsToCombine += ' '+' '.join(baseDataCards).format(ch='zmm')
            script.write('if [ "${b}" = "z" ]; then\n')
            script.write('    cardsList=(%s)\n'%zCardsToCombine)
            script.write('else\n')
            script.write('    cardsList=(%s)\n'%cardsToCombine)
            script.write('fi\n')
            script.write('\n')

            script.write('mkdir -p ${b}_cards\n')
            script.write('combStr=""\n')
            script.write('for icat in "${!cardsList[@]}"; do\n')
            script.write('     dc=${cardsList[$icat]}\n')
            script.write('     full_dc=`readlink -f ${dc}`\n')
            script.write('     full_dc=`dirname ${full_dc}`\n')
            script.write('     full_dc=${full_dc//\//\\\\/}\n')
            script.write('     regex="_(.+)_(.+)\.dat"\n')
            script.write('     if [[ "`basename $dc`" =~ $regex ]]; then\n')
            script.write('         ch="${BASH_REMATCH[1]}"\n')
            script.write('     fi\n')
            script.write('     mod_dc=${b}_cards/datacard_cat${icat}.dat\n')
            script.write('     sed -e "s/mu_bkg/mu_bkgCat${icat}/" -e "s/shapes_/${full_dc}\/shapes_/" ${dc} > ${mod_dc}\n')
            script.write('     echo "nuisance edit rename bkg * ${ch}_bkgShapeEM         ${ch}_bkgShapeEMCat${icat}" >> ${mod_dc}\n')
            script.write('     echo "nuisance edit rename bkg * ${ch}_bkgShapeSingleDiff ${ch}_bkgShapeSingleDiffCat${icat}" >> ${mod_dc}\n')
            script.write('     combStr="${combStr} cat${icat}=${mod_dc}"\n')
            script.write('done\n')
            script.write('\n')

            script.write('combineCards.py ${combStr} > ${b}_datacard.dat\n')
            script.write('sed -i -e "s/${b}_cards\///" ${b}_datacard.dat\n')
            script.write('\n')

            script.write('text2workspace.py ${b}_datacard.dat -m ${m} -o ${pfix}_workspace.root --channel-masks\n')
            script.write('\n')

            script.write('baseCmd=\"combine ${pfix}_workspace.root -m ${m} --X-rtd MINIMIZER_analytic\"\n')
            script.write('${baseCmd} -n PP${b}X.obs   -M AsymptoticLimits\n')
            script.write('${baseCmd} -n PP${b}X       -M AsymptoticLimits --run blind\n')
            script.write('${baseCmd} -n PP${b}X.obs   -M Significance\n')
            script.write('${baseCmd} -n PP${b}X       -M Significance     -t -1 --expectSignal=0.5 --setParameters mu_outfidsig=0.5\n')
            script.write('if [ "${m}" = "1000" ]; then\n')
            script.write('  ${baseCmd} --cminDefaultMinimizerStrategy 0 -n PP${b}X.m${m} -M FitDiagnostics\n')
            script.write('fi\n')
            script.write('cd -\n')

        #run with condor
        with open('%s/zxstatana_run.sub'%odir,'w') as condor:
            condor.write("executable  = %s/statAnaJob.sh\n"%os.path.abspath(odir))
            condor.write("output       = zxstatana_run.out\n")
            condor.write("error        = zxstatana_run.err\n")
            condor.write("log          = zxstatana_run.log\n")
            condor.write("+JobFlavour = \"longlunch\"\n")
            condor.write("request_cpus = 4\n")        
            for mass in  [m for sublist in MASSPOINTS for m in sublist]:
                for boson in ['z','g','zmm','zee']:
                    condor.write("arguments=%d %s\n"%(mass,boson))
                    condor.write("queue 1\n")
        print 'Submitting jobs for',odir
        os.system('condor_submit %s/zxstatana_run.sub'%odir)
        

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
