import ROOT
ROOT.gROOT.SetBatch(True)
import optparse
import os,sys
import json
import re
from collections import OrderedDict
from math import sqrt

debug = True

"""
steer the script
"""
def main():
    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--inDir',       dest='inDir' ,      help='input directory',                default='unfolding/fill',              type='string')
    (opt, args) = parser.parse_args()
    
    observables = ["mult", "width", "ptd", "ptds", "ecc", "tau21", "tau32", "tau43", "zg", "zgxdr", "zgdr", "ga_width", "ga_lha", "ga_thrust", "c1_02", "c1_05", "c1_10", "c1_20", "c2_02", "c2_05", "c2_10", "c2_20", "c3_02", "c3_05", "c3_10", "c3_20"]
    
    flavors = ['all', 'light', 'bottom', 'gluon']

    for obs in observables:
        for flavor in flavors:
            fIn = ROOT.TFile.Open('%s/MC13TeV_TTJets.root'%(opt.inDir))
            h2dForTest = fIn.Get('%s_charged_%s_responsematrix'%(obs, flavor))
            h2dForTest.RebinY(2)
            
            nx,ny=h2dForTest.GetNbinsX(),h2dForTest.GetNbinsY()
            matrix=ROOT.TMatrixD(ny+1,nx)
            for xbin in xrange(1,nx+1):
                for ybin in xrange(0,ny+1):
                    matrix[ybin][xbin-1]=h2dForTest.GetBinContent(xbin,ybin)
            svd=ROOT.TDecompSVD(matrix)
            sig=svd.GetSig()
            maxSig,minSig=sig.Max(),sig.Min()
            condK=-1 if minSig==0 else maxSig/max(0,minSig)
            print(obs,flavor,maxSig,minSig,condK) 
        
"""
for execution from another script
"""
if __name__ == "__main__":
    main()
    #sys.exit(main())

