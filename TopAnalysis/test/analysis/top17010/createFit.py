#!/usr/bin/env python

import os,sys
import optparse

def main():

    #configuration
    usage = 'usage: %prog [options] tag1=datacard1 tag2=datacard2 ... (for single datacards no tag is needed)'
    parser = optparse.OptionParser(usage)
    parser.add_option('-c', '--combine',
                      dest='combine',       
                      help='higgs combination tool location [%default]',
                      default='/afs/cern.ch/user/p/psilva/work/CMSSW_10_3_0_pre4/src/HiggsAnalysis/CombinedLimit',
                      type='string')
    parser.add_option('-o', '--outdir',          
                      dest='outdir',       
                      help='output directory [%default]',
                      default='fit_results',
                      type='string')
    parser.add_option('-t', '--toys',
                      dest='toys',       
                      help='# toys to run [%default]',
                      default=0,
                      type=int)
    parser.add_option('-a', '--asimov',
                      dest='asimov',       
                      help='run the asimov',
                      default=False,
                      action='store_true')
    parser.add_option('--tag',
                      dest='tag',
                      help='use this to tag the fit script',
                      default=None,
                      type='string')

    (opt, args) = parser.parse_args()

    if len(args)==0: 
        parser.print_help()
        sys.exit(-1)

    script=open('runFit.sh','w')
    script.write('#!/bin/bash\n')

    script.write('\n')
    script.write('#setup environment\n')
    script.write('cd %s\n'%opt.combine)
    script.write('eval `scram r -sh`\n')
    script.write('source /afs/cern.ch/user/b/bendavid/work/cmspublic/pythonvenv/tensorflowfit_h5py/bin/activate\n')
    script.write('cd -\n')
    
    script.write('\n')
    if len(args)>1:
        script.write('#combine cards\n')
        for a in args:
            tag,f=a.split('=')
            if not os.path.isfile(f):
                print 'No file found for',tag,':',f


        absPathArgs=[ '{0}={1}'.format(a.split('=')[0],os.path.abspath(a.split('=')[1])) for a in args ]
        script.write('combineCards.py --noDirPrefix %s > datacard.dat\n'%' '.join(absPathArgs))
    else:
        script.write('#make local copy of the datacard\n')
        script.write('cat %s > datacard.dat\n'%os.path.abspath(args[0]))

    script.write('\n')
    script.write('#convert to HDF5 and run TF-based fits\n')
    script.write('text2hdf5.py datacard.dat\n')
    fitresultsName='fitresults' if not opt.tag else 'fitresults_%s'%opt.tag
    script.write('combinetf.py datacard.dat.hdf5 -o %s.root\n'%fitresultsName)
    if opt.asimov:
        script.write('combinetf.py datacard.dat.hdf5 -t -1 -o %s_asimov.root\n'%fitresultsName)
    if opt.toys>0:
        script.write('combinetf.py datacard.dat.hdf5 -t %d -o %s_toys.root\n'%(opt.toys,fitresultsName))
    script.write('cp -v fitresults*root %s\n'%opt.outdir)

    script.close()

    #mv script to output
    os.system('mkdir -p %s'%opt.outdir)
    os.system('mv -v runFit.sh %s'%os.path.join(opt.outdir,'runFit.sh' if not opt.tag else 'runFit_%s.sh'%opt.tag))
    

if __name__ == "__main__":
    sys.exit(main())
