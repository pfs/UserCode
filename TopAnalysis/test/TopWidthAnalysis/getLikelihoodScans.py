#!/usr/bin/env python
import re
from sys import argv
import os.path
import ROOT
from pprint import pprint
from optparse import OptionParser
from tdrStyle import setTDRStyle
import CMS_lumi
parser = OptionParser(
    usage="%prog [options] [label=datacard.txt | datacard.txt]",
    epilog="Collects quantiles information from signal statistics output and turns it into a nice TGraph and LaTeX table. Format of .txt files is stats__<wid>_<lfs>.txt"
    )
parser.add_option("-i",    type="string", dest="indir"  , default="./", help="directory to look for stats files in")
parser.add_option("--wid", type="string", dest="widList", default="0p5w,1p5w,2p0w,2p5w,3p0w,3p5w,4p0w,4p5w,5p0w", help="a list of widths to look for in stats filenames")
parser.add_option("--lfs", type="string", dest="lfsList", default="step1,step2,step3,step4,step5", help="a list of lepton final states to look for in stats filenames")
parser.add_option("-o",    type="string", dest="outdir" , default="./", help="the base filename for the quantiles plot")
parser.add_option("--scanType",    type="string", dest="scanType" , default="scan_Asimov", help="the base filename for the quantiles plot")

(options, args) = parser.parse_args()

# get lists to loop over
rawWidList=options.widList.split(',')
rawLfsList=options.lfsList.split(',')

blueShift = [1,2,3,0,6]
orngShift = [1,3,0,2,9]

for wid in rawWidList :

    # create canvas
    c=ROOT.TCanvas()
    setTDRStyle()
    c.SetGrid()
    c.cd()

    # set the bin and axis labels

    # produce legend
    #leg=ROOT.TLegend(0.15,0.67,0.24,0.85)

    lfsPlotList=[]
    for lfs in rawLfsList :
        f0=ROOT.TFile("%s/higgsCombinex0_%s_%s_%s.MultiDimFit.mH172.5.root"%(options.indir,options.scanType,wid,lfs))
        f1=ROOT.TFile("%s/higgsCombinex1_%s_%s_%s.MultiDimFit.mH172.5.root"%(options.indir,options.scanType,wid,lfs))

        tf0=ROOT.TF1("x0_%s"%lfs, "pol2",0,1)
        tf1=ROOT.TF1("x1_%s"%lfs, "pol2",0,1)

        h0=ROOT.TH2F("h0","h0",100,0,1,500,0,5)
        h1=ROOT.TH2F("h1","h1",100,0,1,500,0,5)

        f0.Get('limit').Draw("2*deltaNLL:x>>h0")
        f1.Get('limit').Draw("2*deltaNLL:x>>h1")

        h0.Fit(tf0);
        h1.Fit(tf1);

        lfsPlotList += [tf0] + [tf1]

    i=0
    for lfsPlot in lfsPlotList :
        if i==0 :
            xax=lfsPlot.GetXaxis()
            xax.SetTitle("Hypothesis sample fraction (x)")
            yax=lfsPlot.GetYaxis()
            yax.SetTitle("-2 #times ln(L_{alt}/L_{null})")
            lfsPlot.SetTitle("")

        lfsPlot.SetLineColor((ROOT.kOrange - orngShift[i/2] if i%2==0 else ROOT.kBlue - blueShift[(i-1)/2]))
        c.cd()
        lfsPlot.Draw(("SAME" if i > 0 else ""))
        i+=1


    # add legend
    # leg.Draw()

    # CMS text
    #CMSLine="CMS"
    #CP=ROOT.TLatex(0.12,0.92, CMSLine)
    #CP.SetNDC(ROOT.kTRUE)
    #CP.SetTextSize(0.05)
    #CP.Draw()

    ## Lumi
    #CMSLineLumi="#sqrt{s}=13 TeV, 2.1 fb^{-1}"
    #CP1=ROOT.TLatex(0.67,0.92, CMSLineLumi)
    #CP1.SetNDC(ROOT.kTRUE)
    #CP1.SetTextSize(0.04)
    #CP1.Draw()

    ## ExtraText
    #CMSLineExtra="#bf{#it{Preliminary}}"
    #CP2=ROOT.TLatex(0.195,0.92, CMSLineExtra)
    #CP2.SetNDC(ROOT.kTRUE)
    #CP2.SetTextSize(0.04)
    #CP2.Draw()

    CMS_lumi.relPosX = 0.180
    CMS_lumi.lumiTextSize = 0.55
    CMS_lumi.CMS_lumi(canvas,4,0)

    # save plots
    c.Modified()
    c.Update()
    c.SaveAs(options.outdir+"LikelihoodScan_%s.pdf"%wid)
    c.SaveAs(options.outdir+"LikelihoodScan_%s.png"%wid)
