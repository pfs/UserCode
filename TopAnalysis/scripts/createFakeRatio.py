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

    gData,jData,jQCD,gMC,cats,binvar,oer=args
    print args
    print 'Running with photon data ',gData,' and photon MC ',gMC 
    print 'The JetData is ', jData, ' and ', jQCD, ' has teh QCD template'
    print 'Fake rates will be provided for these categories', cats
    print 'FR will be calculated in bins of ', binvar
    print 'oer is ', oer
    try:
        cmd='runFRcalculation --fGdata %s --fJdata %s --fJQCD %s --fGMC %s --cats %s --binvar %s --oer %s' %(gData,jData,jQCD,gMC,cats,binvar,oer)        
        print '>>>>>>'
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
    parser.add_option(''  , '--fGdata',     dest='gData',        help='Photon data in signal region [%default]',      default='Data13TeV_SinglePhoton_2016.root',          type='string')
    parser.add_option(''  , '--fJdata',     dest='jData',        help='Jet data from control region [%default]',      default='Data13TeV_JetHT_2016.root',                 type='string')
    parser.add_option(''  , '--fJQCD',      dest='jQCD',         help='Jet data for QCD templates [%default]',        default='Data13TeV_JetHTQCD_2016.root',              type='string')
    parser.add_option(''  , '--fGMC',       dest='gMC',          help='Photon MC in signal region [%default]',        default='MC13TeV_GJets.root',                        type='string')
    parser.add_option(''  , '--cats',       dest='cats',         help='selection categories [%default]',              default='HighVPtA',                                  type='string')
    parser.add_option(''  , '--binvar',     dest='binvar',       help='FR in bins of [%default]',                     default='Mjj',                                       type='string')
    parser.add_option(''  , '--oer',        dest='oer',          help='fake underestimate regulator [%default]',      default='1',                                         type='string')



    (opt, args) = parser.parse_args()
    print opt
    print args

    args = [opt.gData,opt.jData,opt.jQCD,opt.gMC,opt.cats,opt.binvar,opt.oer]
    RunMethodPacked(args)
    
"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
