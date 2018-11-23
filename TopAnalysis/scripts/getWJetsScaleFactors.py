import os
import sys
import optparse
import ROOT
import pickle
import json
from TopLJets2015.TopAnalysis.storeTools import *
from runQCDEstimation import getTemplates

"""
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('--shape',  
                      dest='shape',  
                      help='name of the plotter with the sample used for shape',
                      default='syst_plotter.root',
                      type='string')
    parser.add_option('--norm',  
                      dest='norm',
                      help='name of the plotter with the sample used for normalization',
                      default='plotter.root',
                      type='string')
    parser.add_option('--out',   
                      dest='outdir',      
                      help='output directory',                   
                      default='./',       
                      type='string')
    (opt, args) = parser.parse_args()

    #prepare output
    os.system('mkdir -p %s'%opt.outdir)

    #prepare fitter
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    #ROOT.gROOT.SetBatch(True)
    #ROOT.AutoLibraryLoader.enable()
    ROOT.FWLiteEnabler.enable()

    #open inputs
    fNorm=ROOT.TFile.Open(opt.norm)
    fShape=ROOT.TFile.Open(opt.shape)

    #perform a fit to a variable of interest
    wjetsSF={}
    for sel in ['1j','2j','3j','4j']: 

        #data in the sideband
        #_,_,normH = getTemplates(fIn=fNorm,   dist='nbtags_%s'%sel, tag='norm',  refName='_W+')
        _,_,normH = getTemplates(fIn=fNorm,   dist='nbtags_%s'%sel, tag='norm',  refName='_W')

        #data in the signal region
        #_,_,shapeH = getTemplates(fIn=fShape, dist='nbtags_%s'%sel, tag='shape', refName='_W+')
        _,_,shapeH = getTemplates(fIn=fShape, dist='nbtags_%s'%sel, tag='shape', refName='_W')

        expUnc=ROOT.Double()
        expNorm=shapeH.IntegralAndError(1,shapeH.GetNbinsX(),expUnc)
        targetUnc=ROOT.Double()
        targetNorm=normH.IntegralAndError(1,normH.GetNbinsX(),targetUnc)
        wjetsSF[sel]=( 
            (targetNorm/expNorm,
             ROOT.TMath.Sqrt((targetUnc*expNorm)**2+(targetNorm*expUnc)**2)/(expNorm**2) )
            )

        #free mem
        normH.Delete()
        shapeH.Delete()

    fShape.Close()
    fNorm.Close()
        
    #dump to file    
    print 'W+jets scaling factors'
    print wjetsSF
    cachefile=open('%s/.wjetsscalefactors.pck'%opt.outdir,'w')
    pickle.dump(wjetsSF, cachefile, pickle.HIGHEST_PROTOCOL)
    cachefile.close()

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
