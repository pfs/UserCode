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


    fGammaData, fJetData, fJetQCD, fGammaMC=args
    print args

    try:
        cmd='runFRcalculation --fGdata %s --fJdata %s --fJQCD %s --fGMC %s' %(fGammaData,fJetData,fJetQCD,fGammaMC)
        print(cmd)
        os.system(cmd)
    except :
        print 50*'<'
        print " Something is wrong! "
        print 50*'<'
        return False
    return True

"""
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('', '--fGdata',      dest='fGammaData',        help='Photon data from CR [%default]',                                         default='Data13TeV_SinglePhoton_2017.root',                    type='string')
    parser.add_option('', '--fJdata',      dest='fJetData',     help='Jet data from CR [%default]',                                    default='Data13TeV_JetHT_2017.root',                           type='string')
    parser.add_option('', '--fJQCD',       dest='fJetQCD',       help='Data for QCD photon template  [%default]',                                         default='Data13TeV_JetHTQCD_2017.root',                            type='string')
    parser.add_option('', '--fGMC',        dest='fGammaMC',         help='MC for prompt photon template  [%default]',             default='MC13TeV_GJets.root',       type='string')



    (opt, args) = parser.parse_args()
    print opt
    print args

    args = [opt.fGammaData,opt.fJetData,opt.fJetQCD,opt.fGammaMC]
    RunMethodPacked(args)
    
"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
