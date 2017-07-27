import ROOT
import sys
from roofitTools import *

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gROOT.SetBatch(True)

fIn=ROOT.TFile.Open(sys.argv[1])
w=fIn.Get('w')

for ch in ['e','mu']:
    redData = w.data('data').reduce(ROOT.RooFit.Cut("sample==sample::%s1f4j2q"%ch))
    for varName in ['mjj','mthad','mtlep']:
        showFitResult(fitVar=varName,
                      data=redData,
                      pdf=w.pdf('QCD_%s_%s'%(varName,ch)),
                      categs=[''],
                      w=w,
                      showComponents=[],
                      rangeX=(0,400),
                      tagTitle='%s_QCD'%ch,
                      outDir='./')
        
