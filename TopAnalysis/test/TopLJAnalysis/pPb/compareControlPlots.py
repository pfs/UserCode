import ROOT
from TopLJets2015.TopAnalysis.Plot import Plot
import sys
import os
import optparse

"""
"""
def getPlots(dist,argList):
    plots={}

    colors=[ROOT.kAzure+4, ROOT.kRed+1, ROOT.kGray, ROOT.kGreen+3]
    ci=0
    for arg in argList:
        era,url=arg.split(':')
        
        #get distribution from file
        fIn=ROOT.TFile.Open(url)
        h=fIn.Get(dist)
        try:
            plots[era]=h.Clone(era)
            plots[era].SetTitle(era)
            plots[era].SetDirectory(0)
            plots[era].SetLineColor(colors[ci])
        except:
            pass
        fIn.Close()
        ci+=1

    return plots


"""                                                                                                                                                                                                 
"""
def main():

    #configuration                                                                                                                                                                                   
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option(      '--outDir', dest='outDir', default='./',    type='string',       help='data flag [%default]')
    parser.add_option(      '--ref',    dest='ref',    default=None,    type='string',       help='plotter to use as reference [%default]')
    parser.add_option(      '--data',   dest='data',   default=False,   action="store_true", help='data flag [%default]')
    parser.add_option(      '--shape',  dest='shape',   default=False,   action="store_true", help='shape only flag [%default]')
    (opt, args) = parser.parse_args()

    ROOT.gROOT.SetBatch(True) #False)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)

    os.system('mkdir -p %s'%opt.outDir)

    for ch in ['e','mu']:
        for cat in ['1l4j2q','1l4j1b1q','1l4j2b']:
            for p in ['mjj','mthad','mtlep','mtw','met','pt_l','pt_b','y_b','y_l']:
                dist='%s_%s%s'%(p,ch,cat)
            
                plots=getPlots(dist,args)
                if len(plots)==0 : continue
                p=Plot(dist)
                for era in plots:
                    if opt.shape: plots[era].Scale(1./plots[era].Integral())
                    isRef = True if era == opt.ref else False
                    ci=plots[era].GetLineColor()
                    if isRef : ci=1
                    p.add(plots[era],era,ci,isRef,False,False)
                p.ratiorange=(0.7,1.3)
                p.doPoissonErrorBars=False
                p.show(outDir=opt.outDir,lumi=174e-3,noStack=True,noRatio=False) 


"""
"""
if __name__ == "__main__":
    main()
