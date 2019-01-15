
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

    hist,channel,boson,YEs,nBin,sigPH,year=args
    print args
    print 'Running for ',hist,' variable in the ',channel,boson, ' channel '
    print 'Systematics are ',YEs
    print 'Number of bins are ', nBin
    print 'Running on %s data sample' %(year)
    if sigPH:
        print 'Connecting MM and A signals with paremetricHists'
    else:
        print 'Signals are kept unrelated'
    try:
        cmd='runWorkspaceProvider --Hist %s --Chan %s --V %s --YieldErr %s --nBin %d --doSignalPH --year %s' %(hist,channel,boson,YEs,nBin,year)
        if not sigPH:
            cmd='runWorkspaceProvider --Hist %s --Chan %s --V %s --YieldErr %s --nBin %d --year %s' %(hist,channel,boson,YEs,nBin,year)
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
    parser.add_option(''  , '--Hist',        dest='hist',        help='variable name [%default]',                                         default='vbfmva',                       type='string')
    parser.add_option('-C', '--Chan',        dest='channel',     help='selection category [%default]',                                    default=None,                           type='string')
    parser.add_option('-v', '--V',           dest='boson',       help='vector boson  [%default]',                                         default='A',                            type='string')
    parser.add_option('-y', '--YieldErr',    dest='YEs',         help='nuisance for signal yields: name,down,up  [%default]',             default='lumi_13TeV,1.024,1.024',       type='string')
    parser.add_option('-b', '--nBin',        dest='nBin',        help='number of final bins  [%default]',                                 default=4,                              type='int')
    parser.add_option('-p', '--doSignalPH',  dest='sigPH',       help='Make a parametric signal relation between MM and A',               default=False,                          action='store_true')
    parser.add_option(''  , '--year',        dest='year',        help='dataset year [%default]',                                          default='2017',                         type='string')

    (opt, args) = parser.parse_args()
    print opt
    print args

    args = [opt.hist,opt.channel,opt.boson,opt.YEs,opt.nBin,opt.sigPH,opt.year]
    RunMethodPacked(args)
    
"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
