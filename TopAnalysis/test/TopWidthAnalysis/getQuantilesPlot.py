#!/usr/bin/env python
import re
from sys import argv
import os.path
import ROOT

ROOT.gROOT.SetBatch(True)

from pprint import pprint
from optparse import OptionParser
from tdrStyle import setTDRStyle
import CMS_lumi
parser = OptionParser(
    usage="%prog [options]",
    epilog="Collects quantiles information from signal statistics output and turns it into a nice TGraph and LaTeX table. Format of .txt files is stats__<wid>_<lfs>_<dist>.txt"
    )
parser.add_option("-i",    type="string", dest="indir"  , default="./",     help="directory to look for stats files in")
parser.add_option("-o",    type="string", dest="outdir" , default="./",     help="the base filename for the quantiles plot")
parser.add_option("--lfs", type="string", dest="lfsList", default="",       help="a list of lepton final states to look for in stats filenames")
parser.add_option("--dist",type="string", dest="dist"   , default="incmlb", help="the observable distribution to look at")
parser.add_option("--wid", type="string", dest="widList", default="0p5w,1p0w,1p5w,2p0w,2p5w,3p0w,3p5w,4p0w,4p5w,5p0w", help="a list of widths to look for in stats filenames")
parser.add_option("--axisOverwrite", type="string", dest="aoverList", default="", help="Axis labels to use if desired")

(options, args) = parser.parse_args()

# get axis labels from axis overwrite
axisLabels=options.aoverList.split(',')

# get lists to loop over
rawWidList=options.widList.split(',')
rawLfsList=options.lfsList.split(',')

# create base arrays for eventual tgraph
nPoints = 2*len(rawLfsList)*len(rawWidList)
if options.aoverList != "" and len(axisLabels)*2 != nPoints :
    print "ERROR: axisOverwrite does not write the correct number of labels! Exiting..."
    quit()

distToTitle={
        "mlb": ("Minimum M_{lb}",0.185,0.23),
        "minmlb": ("Minimum M_{lb}",0.185,0.23),
        "mdrmlb": ("#DeltaR-Filtered M_{lb}",0.185,0.23),
        "incmlb": ("Inclusive M_{lb}",0.185,0.23),
        "sncmlb": ("Semi-Inclusive M_{lb}",0.185,0.23),
        "mt2mlb": ("M_{T2}^{lb} Strategy",0.185,0.23)
       }

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
    if wid == "1p0w" :
        eyl3N[i] = 0
        eyl2N[i] = 0
        eyl1N[i] = 0
        y[i]     = 0
        eyu1N[i] = 0
        eyu2N[i] = 0
        eyu3N[i] = 0
        eyl3N[i+1] = 0
        eyl2N[i+1] = 0
        eyl1N[i+1] = 0
        y[i+1]     = 0
        eyu1N[i+1] = 0
        eyu2N[i+1] = 0
        eyu3N[i+1] = 0
        i+=2
        continue

    statsFileName="%s/stats__%s_%s_%s.txt"%(options.indir,wid,lfs,options.dist)
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
setTDRStyle()
c.SetRightMargin(c.GetRightMargin()/4)
c.SetLeftMargin(c.GetLeftMargin()*1.5)
c.SetBottomMargin(c.GetBottomMargin()*1.5)
c.SetGrid(1,1)
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

# draw dist information if available
if options.dist in [key for key in distToTitle] :
    DistInfo,xpos,ypos=distToTitle[options.dist]
    DistLaTeX=ROOT.TLatex(xpos,ypos, DistInfo)
    DistLaTeX.SetNDC(ROOT.kTRUE)
    DistLaTeX.SetTextSize(0.04)
    DistLaTeX.Draw()

    TMassInfo=ROOT.TLatex(.185,.185,"m_{t} = 172.5 GeV")
    TMassInfo.SetNDC(ROOT.kTRUE)
    TMassInfo.SetTextSize(0.03)
    TMassInfo.Draw()

# CMS text
#CMSLine="CMS"
#CP=ROOT.TLatex(0.12,0.92, CMSLine)
#CP.SetNDC(ROOT.kTRUE)
#CP.SetTextSize(0.05)
#CP.Draw()
#
## Lumi
#CMSLineLumi="#sqrt{s}=13 TeV, 2.3 fb^{-1}"
#CP1=ROOT.TLatex(0.67,0.92, CMSLineLumi)
#CP1.SetNDC(ROOT.kTRUE)
#CP1.SetTextSize(0.04)
#CP1.Draw()
#
## ExtraText
#CMSLineExtra="#bf{#it{Preliminary}}"
#CP2=ROOT.TLatex(0.195,0.92, CMSLineExtra)
#CP2.SetNDC(ROOT.kTRUE)
#CP2.SetTextSize(0.04)
#CP2.Draw()

# set the bin and axis labels
xax=totalGraph.GetXaxis()
xax.SetTitle("")
i=0
for wid,lfs in [(wid,lfs) for wid in rawWidList for lfs in rawLfsList]:
    bin_index = xax.FindBin(0.5+i)
    label = "%s#times#Gamma_{SM} %s"%(wid.replace('p','.').replace('w',''),lfs)
    if options.aoverList != "" and len(axisLabels) == nPoints/2 :
        label = axisLabels[i].replace('_',' ');
    xax.SetBinLabel(bin_index,label)
    i+=1

yax=totalGraph.GetYaxis()
yax.SetTitle("-2 #times ln(L_{alt}/L_{null})")
yax.SetTitleOffset(1.1);

# add legend
leg=ROOT.TLegend(0.17,0.63,0.32,0.88)
leg.AddEntry(quantGraph1sigN,"Null, 1#sigma","f")
leg.AddEntry(quantGraph2sigN,"Null, 2#sigma","f")
leg.AddEntry(quantGraph3sigN,"Null, 3#sigma","f")
leg.AddEntry(quantGraph1sigA,"Alt,  1#sigma","f")
leg.AddEntry(quantGraph2sigA,"Alt,  2#sigma","f")
leg.AddEntry(quantGraph3sigA,"Alt,  3#sigma","f")
leg.AddEntry(quantGraphExp  ,"Median"       ,"l")
leg.Draw()

#massAxis=ROOT.TGaxis(-.5,-2.5,10.5,-2.5,.25*1.324,5.25*1.324,505,"")
#massAxis.SetTitle("Generator-level #Gamma_{t}")
#massAxis.Draw()

CMS_lumi.relPosX = 0.180
CMS_lumi.lumiTextSize = 0.55
CMS_lumi.CMS_lumi(c,4,0)

# save plots
c.Modified()
c.Update()
c.SaveAs(options.outdir+"quantiles.pdf")
c.SaveAs(options.outdir+"quantiles.png")
