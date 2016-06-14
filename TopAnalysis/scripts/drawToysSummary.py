import ROOT
import numpy as np
import sys

inDir=sys.argv[1]

gr1s=ROOT.TGraphAsymmErrors()
gr1s.SetFillStyle(1001)
gr1s.SetFillColor(19)
gr1s.SetLineColor(19)
gr1s.SetMarkerColor(1)
gr1s.SetMarkerStyle(24)
gr1s.SetName('gr1s')
gr1s.SetTitle('68%CL')
gr2s=ROOT.TGraphAsymmErrors()
gr2s.SetFillStyle(1001)
gr2s.SetFillColor(18)
gr2s.SetLineColor(18)
gr2s.SetMarkerColor(18)
gr2s.SetName('gr2s')
gr2s.SetTitle('95%CL')

for rstr in ['0.9', '1', '1.1']:
    fIn=ROOT.TFile.Open('%s/higgsCombineTest.MaxLikelihoodFit.mH%s.123456.root'%(inDir,rstr))
    limit=fIn.Get('limit')

    toys=[]
    for i in xrange(0,limit.GetEntriesFast()):
        limit.GetEntry(i)
        toys.append(limit.limit)
    fIn.Close()

    print len(toys)
    median=np.percentile(toys, 50)
    lo68,up68=np.percentile(toys,16), np.percentile(toys,84)
    lo95,up95=np.percentile(toys,2.5),np.percentile(toys,97.5)
    print lo95,lo68,median,up68,up95

    r=float(rstr)
    npt=gr1s.GetN()
    gr1s.SetPoint(npt,r,median)
    gr1s.SetPointError(npt,0.0,0.0,median-lo68,up68-median)
    gr2s.SetPoint(npt,r,median)
    gr2s.SetPointError(npt,0.0,0.0,median-lo95,up95-median)
    
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gROOT.SetBatch(True)
   
c=ROOT.TCanvas('c','c',500,500)
c.SetBottomMargin(0.1)
c.SetTopMargin(0.05)
c.SetLeftMargin(0.12)
c.SetRightMargin(0.05)

gr2s.Draw('a3')
gr2s.GetYaxis().SetRangeUser(0.7,1.3)
gr2s.GetXaxis().SetRangeUser(0.8,1.2)
gr2s.GetYaxis().SetTitle('Fitted #mu')
gr2s.GetXaxis().SetTitle('Generated #mu')
gr1s.Draw('3')
gr1s.Draw('p')
gr1s.Fit('pol1')

pol1=gr1s.GetFunction('pol1')

leg=ROOT.TLegend(0.15,0.9,0.4,0.8)
leg.SetTextFont(42)
leg.SetTextSize(0.035)
leg.SetBorderSize(0)
leg.SetFillStyle(1001)
leg.SetFillColor(0)
leg.AddEntry(gr1s,gr1s.GetTitle(),'fp')
leg.AddEntry(gr2s,gr2s.GetTitle(),'f')
leg.Draw()

cmsLabel=ROOT.TLatex()
cmsLabel.SetTextFont(42)
cmsLabel.SetTextSize(0.035)
cmsLabel.SetNDC()
cmsLabel.DrawLatex(0.15,0.92,'#bf{CMS} #it{simulation}')
cmsLabel.DrawLatex(0.15,0.2,'#scale[0.8]{#it{offset=%3.2f#pm%3.2f}}'%(pol1.GetParameter(0),pol1.GetParError(0)))
cmsLabel.DrawLatex(0.15,0.15,'#scale[0.8]{#it{slope =%3.2f#pm%3.2f}}'%(pol1.GetParameter(1),pol1.GetParError(1)))

c.SaveAs('toyssummary.pdf')

