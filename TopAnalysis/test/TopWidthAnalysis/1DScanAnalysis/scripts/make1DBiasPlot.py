#!/usr/bin/env python

import optparse
import ROOT
import os,sys
import glob

ROOT.gROOT.SetBatch(True)

parser=optparse.OptionParser("Get plots of 1D bias for TOP-16-019")
parser.add_option('-i', '--input', dest='input', help='input datacards', type='string',
        default='./biasplots/')
parser.add_option('-o', '--outdir', dest='outdir', help='output directory', type='string',
        default='./biasplots/')
opt,args=parser.parse_args()

wids=[x*0.2 for x in range(1,16)]
wids=[("%.1fw"%x).replace('.','p') for x in wids]

nPoints = len(glob.glob('%s/higgsCombineTopWid1D_w*_m*.root'%opt.input))
x   = ROOT.TVector(nPoints)
y   = ROOT.TVector(nPoints)
xerr= ROOT.TVector(nPoints)
yerr= ROOT.TVector(nPoints)

iFile=0
for fName in glob.glob('%s/higgsCombineTopWid1D_w*_m*.root'%opt.input) :
    fIn = ROOT.TFile(fName);
    limit=fIn.Get('limit')
    limit.Draw('limit:limitErr','quantileExpected==-1')
    x[iFile]=float(fName.split('_')[1].replace('p','.').replace('w',''))
    y[iFile]=(limit.GetV1()[0]-x[iFile])*1.324
    yerr[iFile]=limit.GetV2()[0]

    fIn.Close()
    iFile+=1


canvas=ROOT.TCanvas()
canvas.cd()
canvas.SetRightMargin(0.05)
canvas.SetTopMargin(0.05);
canvas.SetGrid()

gr=ROOT.TGraphErrors(x,y,xerr,yerr)
gr.Draw("AP")
gr.SetMarkerStyle(20)
gr.GetXaxis().SetTitleOffset(1.25)
gr.GetXaxis().SetTitle("Injected width (#Gamma/#Gamma_{SM})")
gr.GetYaxis().SetTitle("Bias [GeV]")
gr.SetTitle('')
canvas.SaveAs("bla.pdf")
