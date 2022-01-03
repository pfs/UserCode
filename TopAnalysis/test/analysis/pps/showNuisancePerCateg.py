#!/usr/bin/env python
import ROOT
import sys

categ_dict={'multi-multi':[34,32,36,33,37,38,35,39,4,6,3,2,0,5,7,1],
            'multi-single':[40,45,46,44,47,43,42,41,10,8,9,11,13,15,14,12],
            'single-multi':[49,52,55,48,54,50,51,53,20,19,18,23,17,16,21,22],
            'single-single':[60,58,63,57,62,59,56,61,24,31,27,30,25,26,28,29]}

#gamma
categ_dict={'multi-multi':[3,6,0,4,7,1,5,2],
            'multi-single':[15,11,9,10,14,13,8,12],
            'single-multi':[18,17,21,16,23,19,22,20],
            'single-single':[30,31,28,29,27,26,24,25]}

url=sys.argv[1]
inF=ROOT.TFile.Open(url)
fres=inF.Get('fit_s')
pars=fres.floatParsFinal()
iter = pars.createIterator()
var = iter.Next()

gr_dict={}
while var :
    name=var.GetName()
    try:

        categ=int(name.split('Cat')[-1])

        reco=None
        for reco_categ,clist in categ_dict.items():
            if not categ in clist: continue
            reco=reco_categ
            break

        syst=None
        if 'bkgShapeEM' in name :
            syst='Shape HP'
        if 'bkgShapeSingleDiff' in name : 
            syst='Shape SD'

        val=var.getVal()
        ehi=var.getErrorHi()
        elo=var.getErrorLo()

        key=(reco,syst)
        if not key in gr_dict:
            gr_dict[key]=ROOT.TGraphAsymmErrors()
            gr_dict[key].SetName('nuisances_{}'.format(key))
            gr_dict[key].SetTitle('{} ({})'.format(syst,reco))
            gr_dict[key].SetMarkerStyle(20)
        n=gr_dict[key].GetN()
        gr_dict[key].SetPoint(n,val,n)
        gr_dict[key].SetPointError(n,abs(elo),abs(ehi),0,0)

    except Exception as e:
        pass

    var=iter.Next()

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gROOT.SetBatch(True)
c=ROOT.TCanvas('c','c',600,600)
c.SetBottomMargin(0.1)
c.SetRightMargin(0.05)
c.SetTopMargin(0.06)
c.SetLeftMargin(0.08)
npars=len(categ_dict['multi-multi'])
frame=ROOT.TH2F('frame',';#hat{#theta};Nuisance parameter',1,-3,3,npars,-0.5,npars-0.5)
frame.GetYaxis().SetLabelSize(0)

#1-sigma
gr1s=ROOT.TGraph()
gr1s.SetName('gr1s')
gr1s.SetMarkerStyle(1)
gr1s.SetMarkerColor(19)
gr1s.SetLineColor(19)
gr1s.SetFillStyle(1001)
gr1s.SetFillColor(19)
gr1s.SetPoint(0,-1,-0.5)
gr1s.SetPoint(1,-1,npars)
gr1s.SetPoint(2,1,npars)
gr1s.SetPoint(3,1,-0.5)
gr1s.SetPoint(4,-1,-0.5)

#2-sigma
gr2s=gr1s.Clone('gr2s')
gr2s.SetMarkerColor(18)
gr2s.SetLineColor(18)
gr2s.SetFillStyle(1001)
gr2s.SetFillColor(18)
gr2s.SetPoint(0,-2,-0.5)
gr2s.SetPoint(1,-2,npars)
gr2s.SetPoint(2,2,npars)
gr2s.SetPoint(3,2,-0.5)
gr2s.SetPoint(4,-2,-0.5)

for k,gr in gr_dict.items():
    c.Clear()
    frame.Draw()
    gr2s.Draw('f')
    gr1s.Draw('f')
    gr.Draw('p')
    name='{}_{}_grouped_nuisances.png'.format(*k)
    name=name.replace(' ','')
    name=name.replace('-','_')

    txt=ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(42)
    txt.SetTextSize(0.05)
    txt.SetTextAlign(12)
    txt.DrawLatex(0.08,0.955,'#bf{CMS} #it{Preliminary}')

    txt.SetTextAlign(ROOT.kHAlignRight+ROOT.kVAlignCenter)
    txt.SetTextSize(0.04)
    txt.DrawLatex(0.97,0.955,gr.GetTitle())

    l=ROOT.TLine(0,-0.5,0,npars-0.5)
    l.SetLineColor(ROOT.kAzure+3)
    l.Draw()

    c.RedrawAxis()
    c.Modified()
    c.Update()
    c.SaveAs(name)
    

