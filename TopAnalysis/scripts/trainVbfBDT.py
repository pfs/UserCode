
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

    options,extention,signal,background,category,card=args
    print args
    print 'Running with options ',options,' in ',category, ' category'
    print 'The extention is ', extention, ' and the card is ', card
    print 'Signal input is ',signal, ' and background input is ', background
    try:
        cmd='runTMVAClassification --vbf %s --ext %s --sig %s --bkg %s --cat %s --card %s' %(options,extention,signal,background,category,card)        
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
    parser.add_option(''  , '--vbf',        dest='options',      help='options for BDT training [%default]',                default='nt=50:mns=5:md=3:abb=0.6:nc=1',                       type='string')
    parser.add_option(''  , '--ext',        dest='extention',    help='extention for the weight file [%default]',           default='LowVPt',                                       type='string')
    parser.add_option(''  , '--sig',        dest='signal',       help='signal file name  [%default]',                       default='signal.root',                                         type='string')
    parser.add_option(''  , '--bkg',        dest='background',  help='background file name [%default]',                    default='backgrounds.root',                                    type='string')
    parser.add_option(''  , '--cat',        dest='category',     help='selection category [%default]',                      default='A:VBF',                                               type='string')
    parser.add_option(''  , '--card',       dest='card',         help='card name including input variables [%default]',     default='LowVPtCard',                                   type='string')


    (opt, args) = parser.parse_args()
    print opt
    print args

    args = [opt.options,opt.extention,opt.signal,opt.background,opt.category,opt.card]
    RunMethodPacked(args)
    
"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
