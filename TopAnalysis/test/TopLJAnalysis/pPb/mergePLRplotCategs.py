#!/usr/bin/env python

import ROOT
import sys

lumi=(174.5,8.725)
#from runDataFit import lumi

plots={'pdfcomp1':None,
       'pdfcomp0':None,
       'data':None,
       'totalpdf':None}

var=sys.argv[1]
varTitle='M(jj) [GeV]'
varRange=(25,300)
if 'mthad' in var:
    varTitle='M_{top} [GeV]'
    varRange=(100,400)
if 'mtlep' in var:
    varTitle='M_{t_{lep}} [GeV]'
    varRange=(100,400)

cat=sys.argv[2]
catTitle = 'e^{#pm}/mu^{#pm} + #geq4j'
if '2q' in cat   : catTitle += ' (=0b)'
if '1b1q' in cat : catTitle += ' (=1b)'
if '2b' in cat  : catTitle += ' (#geq2b)'

maxY=0
for ch in ['e','mu']:
    fIn=ROOT.TFile.Open('%s_%s%s_final.root'%(var,ch,cat))
    for p in fIn.Get('c').GetPrimitive('p1').GetListOfPrimitives():
        pname=p.GetName()
        if not pname in plots : continue

        print ch,pname,p.Integral()

        if plots[pname] is None :
            plots[pname]=p.Clone()
        else:

            x,y=ROOT.Double(0),ROOT.Double(0)
            if pname!='data':
                
                for n in xrange(0,plots[pname].GetN()):
                    plots[pname].GetPoint(n,x,y)
                    plots[pname].SetPoint(n,x,y+p.Eval(x))
            else:

                xcur,ycur=ROOT.Double(0),ROOT.Double(0)
                for n in xrange(0,p.GetN()):
                    p.GetPoint(n,x,y)
                    
                    for k in xrange(0,plots[pname].GetN()):
                        
                        plots[pname].GetPoint(k,xcur,ycur)
                        if xcur!=x : continue 
                    
                        exhi=p.GetErrorXhigh(n)
                        exlo=p.GetErrorXlow(n)
                        eyhi=p.GetErrorYhigh(n)
                        eylo=p.GetErrorYlow(n)

                        exhicur=plots[pname].GetErrorXhigh(k)
                        exlocur=plots[pname].GetErrorXlow(k)
                        eyhicur=plots[pname].GetErrorYhigh(k)
                        eylocur=plots[pname].GetErrorYlow(k)

                        maxY=max(y+ycur,maxY)
                        plots[pname].SetPoint(k,x,y+ycur)
                        plots[pname].SetPointError(k,
                                                   ROOT.TMath.Sqrt(exlo**2+exlocur**2),
                                                   ROOT.TMath.Sqrt(exhi**2+exhicur**2),
                                                   ROOT.TMath.Sqrt(eylo**2+eylocur**2),
                                                   ROOT.TMath.Sqrt(eyhi**2+eyhicur**2))

                        

    fIn.Close()


ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gROOT.SetBatch(True)

c=ROOT.TCanvas('c','c',800,800)
c.SetTopMargin(0)
c.SetLeftMargin(0)
c.SetRightMargin(0)
c.SetBottomMargin(0)
p1 = ROOT.TPad('p1','p1',0.0,0.85,1.0,0.0)
p1.SetRightMargin(0.05)
p1.SetLeftMargin(0.12)
p1.SetTopMargin(0.008)
p1.SetBottomMargin(0.15)
p1.Draw()
p1.cd()
plots['totalpdf'].Draw('al')
plots['pdfcomp1'].Draw('lf')
plots['pdfcomp0'].Draw('lf')
plots['data'].Draw('ep')
plots['totalpdf'].GetYaxis().SetRangeUser(0,maxY*1.23)
plots['totalpdf'].GetYaxis().SetTitle("Events")
plots['totalpdf'].GetYaxis().SetTitleOffset(1.0)
plots['totalpdf'].GetYaxis().SetTitleSize(0.06)
plots['totalpdf'].GetYaxis().SetLabelSize(0.05)
plots['totalpdf'].GetXaxis().SetTitleSize(0.06)
plots['totalpdf'].GetXaxis().SetLabelSize(0.05)
plots['totalpdf'].GetXaxis().SetTitleOffset(0.9)
plots['totalpdf'].GetXaxis().SetTitle(varTitle)
plots['totalpdf'].GetXaxis().SetRangeUser(varRange[0],varRange[1])

label = ROOT.TLatex()
label.SetNDC()
label.SetTextFont(42)
label.SetTextSize(0.05)
label.DrawLatex(0.5,0.92,'#bf{CMS} #it{preliminary}')
label.DrawLatex(0.5,0.86,'%3.0f nb^{-1} (#sqrt{s_{NN}}=8.16 TeV)'%lumi[0])
label.DrawLatex(0.5,0.8,'#it{%s}'%catTitle)
label.DrawLatex(0.5,0.74,'#chi^{2}/ndof=%3.2f'%plots['totalpdf'].chiSquare(plots['data'],1))

leg=ROOT.TLegend(0.57,0.7,0.95,0.5)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.SetTextFont(42)
leg.SetTextSize(0.035)
leg.AddEntry(plots['data'],'data','ep')
leg.AddEntry(plots['totalpdf'],'background','f')
leg.AddEntry(plots['pdfcomp1'],'t#bar{t} correct permutations','f')
leg.AddEntry(plots['pdfcomp0'],'t#bar{t} wrong permutations','f')
leg.Draw()


c.cd()
p2 = ROOT.TPad('p2','p2',0.0,0.85,1.0,1.0)
p2.SetBottomMargin(0.005)
p2.SetRightMargin(0.05)
p2.SetLeftMargin(0.12)
p2.SetTopMargin(0.05)
p2.SetGridy(True)
p2.Draw()
p2.cd()
p2.Clear()
pullGr=plots['data'].makeResidHist(plots['totalpdf'],True,True)
pullGr.Draw('ap')
pullGr.GetYaxis().SetRangeUser(-2.4,2.4)
pullGr.GetYaxis().SetTitle("Pull")
pullGr.GetYaxis().SetTitleSize(0.35)
pullGr.GetYaxis().SetLabelSize(0.3)
pullGr.GetXaxis().SetTitleSize(0)
pullGr.GetXaxis().SetLabelSize(0)
pullGr.GetYaxis().SetTitleOffset(0.17)
pullGr.GetYaxis().SetNdivisions(4)
pullGr.GetYaxis().SetRangeUser(-3.1,3.1)
pullGr.GetXaxis().SetTitleOffset(0.8)
pullGr.GetXaxis().SetRangeUser(varRange[0],varRange[1])
c.cd()
c.Modified()
c.Update()
for ext in ['png','pdf']:
    c.SaveAs('%s_comb_%s.%s'%(var,cat,ext))
             
