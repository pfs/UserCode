import ROOT
import sys
import pickle

from runDataFit import lumi
from runQCDestimation import getHistos,QCDNORM,PROCNORM
from TopLJets2015.TopAnalysis.Plot import Plot

def showPlotFor(args,var):
    
    histos=getHistos(args,var)

    for ch in ['e','mu']:

        qcdKey=(ch,'1l4j2q')        
        for cat in ['1l4j2q','1l4j1b1q','1l4j2b']:
            key=(ch,cat)

            #scale up processes
            histos['ttbar'][key].Scale(lumi[0]*PROCNORM[('ttbar',ch)])
            histos['wjets'][key].Scale(lumi[0]*PROCNORM[('wjets',ch)])
            histos['dy'][key].Scale(lumi[0]*PROCNORM[('dy',ch)])
            qcdHisto=histos['qcd'][qcdKey].Clone('qcd_%s%s'%(ch,cat))
            qcdHisto.Scale(QCDNORM[key][0]/qcdHisto.Integral(0,qcdHisto.GetNbinsX()+1))
            for xbin in xrange(1,qcdHisto.GetNbinsX()+1):
                totalUnc=qcdHisto.GetBinError(xbin)**2+(qcdHisto.GetBinContent(xbin)*QCDNORM[key][1]/QCDNORM[key][0])**2
                qcdHisto.SetBinError(xbin,ROOT.TMath.Sqrt(totalUnc))


            p=Plot('%s%s_%s_control'%(ch,cat,var))
            p.add(histos['ttbar'][key],'t#bar{t}',0,False,False,False)
            p.add(histos['wjets'][key],'W',ROOT.TColor.GetColor('#fee090'),False,False,False)
            p.add(histos['dy'][key],'DY',ROOT.TColor.GetColor('#fc8d59'),False,False,False)
            p.add(qcdHisto,'QCD',ROOT.TColor.GetColor('#e0f3f8'),False,False,False)
            p.add(histos['data'][key],'Data',1,True,False,False)
        
            p.ratiorange=(0.4,1.7)
            p.com='8.16 TeV'
            p.show(outDir='./',lumi=174*1e-3)
            p.reset()

"""
"""
def main():

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(True)

    #load qcd normalization
    global QCDNORM
    with open('qcdnorm.pck','r') as fIn:
        QCDNORM=pickle.load(fIn)

    args=sys.argv[1:]
    for var in ['met','mtw','minmlb','maxmlb','mjj','drjj','ntracks','ntrackshp','mthad','mtlep']: 
        showPlotFor(args,var)


if __name__ == "__main__":
    main()
