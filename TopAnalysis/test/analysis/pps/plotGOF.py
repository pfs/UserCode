import os
import sys
import ROOT
import numpy as np

def getTestStatFrom(url):
    print url
    fIn=ROOT.TFile.Open(url)
    limit=fIn.Get('limit')
    vals=[]
    for i in range(limit.GetEntriesFast()):
        limit.GetEntry(i)
        vals.append(limit.limit)
    return vals

url=sys.argv[1]
obs=getTestStatFrom(url)[0]
exp=getTestStatFrom(url.replace('.gof.','.gof.toy.').replace('.root','.123456.root'))
p=np.arange(0,110,10)
q=np.percentile(exp, p )


pfix=os.path.basename(url).split('.gof')[0]
lumi=37193 if 'z' in pfix else 2288

gr=ROOT.TGraph()
for i in range(len(p)):
    gr.SetPoint(i,q[i],p[i])
pval=gr.Eval(obs)
gr.Print('v')
print(pval)
if pval>50.: pval=100-pval
print(pval)
pval/=100.
print(pval)

title='pp#rightarrow ppZX'
if 'PPzmmX' in url: title='pp#rightarrow ppZ#mu#muX'
if 'PPzeeX' in url: title='pp#rightarrow ppZeeX'
if 'PPgX' in url  : title='pp#rightarrow pp#gammaX'


htstat=ROOT.TH1F('tstat',';Test statistics;Toys',25,q[0]*0.95,q[-1]*1.05)
for t in exp: htstat.Fill(t)

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gROOT.SetBatch(True)

c=ROOT.TCanvas('c','c',500,500)
c.SetLeftMargin(0.12)
c.SetRightMargin(0.03)
c.SetTopMargin(0.05)
c.SetBottomMargin(0.1)

htstat.Draw('hist')
htstat.GetYaxis().SetRangeUser(0,htstat.GetMaximum()*1.2)
l=ROOT.TLine()
l.SetLineWidth(2)
l.SetLineColor(ROOT.kBlue)
l.DrawLine(obs,0,obs,htstat.GetMaximum())

tex=ROOT.TLatex()
tex.SetTextFont(42)
tex.SetTextSize(0.04)
tex.SetNDC()
tex.DrawLatex(0.12,0.96,'#bf{CMS} #it{preliminary}')
tex.SetTextAlign(31)
tex.DrawLatex(0.95,0.96,'#scale[0.7]{ %3.1f fb^{-1} (13 TeV)}'%(lumi/1000.))
tex.SetTextAlign(ROOT.kHAlignLeft+ROOT.kVAlignCenter)
tex.DrawLatex(0.15,0.9,title)
tex.DrawLatex(0.15,0.84,'GOF p-val:{:3.3f}'.format(pval))
c.Modified()
c.Update()   
c.SaveAs('gof_{}.pdf'.format(pfix))
c.SaveAs('gof_{}.png'.format(pfix))

