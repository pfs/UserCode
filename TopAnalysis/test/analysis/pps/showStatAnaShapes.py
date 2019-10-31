from compareOptimResults import showShapes
import sys
import ROOT

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gROOT.SetBatch(True)

showShapes(resultsDir=sys.argv[1],
           name=sys.argv[2],
           title=sys.argv[3],
           mass=int(sys.argv[4]),
           boson=sys.argv[5],
           r95=0,
           sig=0,
           lumi=37.5)
