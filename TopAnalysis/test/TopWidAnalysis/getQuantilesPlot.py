#!/usr/bin/env python
import re
from sys import argv
import os.path
import ROOT
from pprint import pprint
from optparse import OptionParser
parser = OptionParser(
    usage="%prog [options] [label=datacard.txt | datacard.txt]",
    epilog="Collects quantiles information from signal statistics output and turns it into a nice TGraph and LaTeX table. Format of .txt files is stats__<wid>_<lfs>.txt"
    )
parser.add_option("-i",    type="string", dest="indir"  , default="./"            ,   help="directory to look for stats files in")
parser.add_option("--wid", type="string", dest="widList", default="0p5w,2p0w,4p0w",   help="a list of widths to look for in stats filenames")
parser.add_option("--lfs", type="string", dest="lfsList", default=""              ,   help="a list of lepton final states to look for in stats filenames")
parser.add_option("-o",    type="string", dest="outdir" , default="./"   ,   help="the base filename for the quantiles plot")

(options, args) = parser.parse_args()

# get lists to loop over
rawWidList=options.widList.split(',')
rawLfsList=options.lfsList.split(',')

# create base arrays for eventual tgraph
nPoints = 2*len(rawLfsList)*len(rawWidList)

x    =ROOT.TVector(nPoints)
y    =ROOT.TVector(nPoints)
ex   =ROOT.TVector(nPoints)
eyl1N=ROOT.TVector(nPoints) # 1 sigma deviations
eyu1N=ROOT.TVector(nPoints)
eyl1A=ROOT.TVector(nPoints) # 1 sigma deviations
eyu1A=ROOT.TVector(nPoints)
eyl2N=ROOT.TVector(nPoints) # 2 sigma
eyu2N=ROOT.TVector(nPoints)
eyl2A=ROOT.TVector(nPoints) # 2 sigma
eyu2A=ROOT.TVector(nPoints)
eyl3N=ROOT.TVector(nPoints) # 3 sigma
eyu3N=ROOT.TVector(nPoints)
eyl3A=ROOT.TVector(nPoints) # 3 sigma
eyu3A=ROOT.TVector(nPoints)

eyexp=ROOT.TVector(nPoints)

# initialize standard arrays
for i in xrange(0,nPoints) :
    x[i]     = 0.25 + 0.5*i
    ex[i]    = 0.2
    eyexp[i] = 0.25

# loop over widths, lfs, parse array info
i=0
for wid,lfs in [(wid,lfs) for wid in rawWidList for lfs in rawLfsList]:
    statsFileName="%s/stats__%s_%s.txt"%(options.indir,wid,lfs)
    for line in open(statsFileName,"r"):
        if "nulquant" in line :
            tline = map(float,line.split(";")[1:8]);
            eyl3N[i] = ROOT.TMath.Abs(tline[3]-tline[0])
            eyl2N[i] = ROOT.TMath.Abs(tline[3]-tline[1])
            eyl1N[i] = ROOT.TMath.Abs(tline[3]-tline[2])
            y[i]     = tline[3]
            eyu1N[i] = ROOT.TMath.Abs(tline[4]-tline[3])
            eyu2N[i] = ROOT.TMath.Abs(tline[5]-tline[3])
            eyu3N[i] = ROOT.TMath.Abs(tline[6]-tline[3])
        elif "altquant" in line :
            tline = map(float,line.split(";")[1:8]);
            eyl3A[i+1] = ROOT.TMath.Abs(tline[3]-tline[0])
            eyl2A[i+1] = ROOT.TMath.Abs(tline[3]-tline[1])
            eyl1A[i+1] = ROOT.TMath.Abs(tline[3]-tline[2])
            y[i+1]     = tline[3]
            eyu1A[i+1] = ROOT.TMath.Abs(tline[4]-tline[3])
            eyu2A[i+1] = ROOT.TMath.Abs(tline[5]-tline[3])
            eyu3A[i+1] = ROOT.TMath.Abs(tline[6]-tline[3])
        else : continue
    i+=2

for i in xrange(0,nPoints) :
    print y[i]

# create graphs
quantGraph1sigN = ROOT.TGraphAsymmErrors(x,y,ex,ex,eyl1N,eyu1N);
quantGraph1sigA = ROOT.TGraphAsymmErrors(x,y,ex,ex,eyl1A,eyu1A);
quantGraph2sigN = ROOT.TGraphAsymmErrors(x,y,ex,ex,eyl2N,eyu2N);
quantGraph2sigA = ROOT.TGraphAsymmErrors(x,y,ex,ex,eyl2A,eyu2A);
quantGraph3sigN = ROOT.TGraphAsymmErrors(x,y,ex,ex,eyl3N,eyu3N);
quantGraph3sigA = ROOT.TGraphAsymmErrors(x,y,ex,ex,eyl3A,eyu3A);
quantGraphExp   = ROOT.TGraphAsymmErrors(x,y,ex,ex,eyexp,eyexp);
#quantGraphData

# create canvas
c=ROOT.TCanvas()
c.SetGrid()
c.cd()

# format all graphs: color
quantGraph1sigN.SetFillColor(ROOT.kBlue)
quantGraph1sigA.SetFillColor(ROOT.kOrange+7)
quantGraph2sigN.SetFillColor(ROOT.kBlue-7)
quantGraph2sigA.SetFillColor(ROOT.kOrange+1)
quantGraph3sigN.SetFillColor(ROOT.kBlue-9)
quantGraph3sigA.SetFillColor(ROOT.kOrange)

quantGraphExp.SetFillColor(ROOT.kBlack)

# draw as a multigraph
totalGraph=ROOT.TMultiGraph()
totalGraph.Add(quantGraph3sigN)
totalGraph.Add(quantGraph3sigA)
totalGraph.Add(quantGraph2sigN)
totalGraph.Add(quantGraph2sigA)
totalGraph.Add(quantGraph1sigN)
totalGraph.Add(quantGraph1sigA)
totalGraph.Add(quantGraphExp)
totalGraph.Draw("a2")

# CMS text
CMSLine="CMS"
CP=ROOT.TLatex(0.12,0.92, CMSLine)
CP.SetNDC(ROOT.kTRUE)
CP.SetTextSize(0.05)
CP.Draw()

# Lumi
CMSLineLumi="#sqrt{s}=13 TeV, 2.1 fb^{-1}"
CP1=ROOT.TLatex(0.67,0.92, CMSLineLumi)
CP1.SetNDC(ROOT.kTRUE)
CP1.SetTextSize(0.04)
CP1.Draw()

# ExtraText
CMSLineExtra="#bf{#it{Preliminary}}"
CP2=ROOT.TLatex(0.195,0.92, CMSLineExtra)
CP2.SetNDC(ROOT.kTRUE)
CP2.SetTextSize(0.04)
CP2.Draw()

# set the bin and axis labels
xax=totalGraph.GetXaxis()
xax.SetTitle("")
i=0
for wid,lfs in [(wid,lfs) for wid in rawWidList for lfs in rawLfsList]:
    bin_index = xax.FindBin(0.5+i)
    label = "%s %s"%(wid,lfs)
    xax.SetBinLabel(bin_index,label)
    i+=1

yax=totalGraph.GetYaxis()
yax.SetTitle("-2 #times ln(L_{alt}/L_{null})")

# add legend
leg=ROOT.TLegend(0.15,0.67,0.24,0.85)
leg.AddEntry(quantGraph1sigN,"Null, 1#sigma","f")
leg.AddEntry(quantGraph2sigN,"Null, 2#sigma","f")
leg.AddEntry(quantGraph3sigN,"Null, 3#sigma","f")
leg.AddEntry(quantGraph1sigA,"Alt,  1#sigma","f")
leg.AddEntry(quantGraph2sigA,"Alt,  2#sigma","f")
leg.AddEntry(quantGraph3sigA,"Alt,  3#sigma","f")
leg.AddEntry(quantGraphExp  ,"Median"       ,"l")
leg.Draw()

# save plots
c.Modified()
c.Update()
c.SaveAs(options.outdir+options.lfsList.replace(',','')+"quantiles.pdf")
c.SaveAs(options.outdir+options.lfsList.replace(',','')+"quantiles.png")
