#!/usr/bin/env python
import os,sys
import ROOT
from optparse import OptionParser
from collections import defaultdict
import numpy as np
from rounding import *
import optparse
import pandas as pd

_skipstat=True

def doNuisanceReport(url_list,tag='',npergroup=20,poi='r'):

    """compare postfit nuisances"""

    fitResults=defaultdict(dict)
    poiResults={}
    for title,fit,url in url_list:
        inF=ROOT.TFile.Open(url)
        fres=inF.Get('fit_%s'%fit)
        pars=fres.floatParsFinal()
        iter = pars.createIterator()
        var = iter.Next()
        fitResults[title]=defaultdict(dict)
        while var :
            name=var.GetName()

            if _skipstat and name.find('prop_bincat')==0:
                var = iter.Next()
                continue

            if name.find('mu_bkg')==0:
                var = iter.Next()
                continue
           
            print name,'->',
            for t,rt in [
                    ('_',' '),
                    ('prop bincat','stat. cat'),
                    ('zmm','#mu#mu'),
                    ('zee','ee'),
                    #('mm','#mu#mu'),
                    ('mu outfidsig','#nu out.fid.sig.'),
                    ('bkgShapeEM','shape HP'),
                    ('bkgShapeSingleDiff','shape SD'),
                    ('mu bkg',' #lambda(bkg)'),
                    ('Cat',' cat.'),
                    ('sigPPSEff','PPS eff.'),
                    ('sigCalib','Time dependency'),
                    ('sigShapeEM','Signal shape HP'),
                    ('g shape','#gamma shape'),
            ]:            
                name=name.replace(t,rt)
            print name
            
            val=var.getVal()
            ehi=var.getErrorHi()
            elo=var.getErrorLo()
            if name!=poi:
                corr=fres.correlation(var.GetName(),poi) if fit=='s' else 0.
                print(fit,var.GetName(),poi,corr)
                fitResults[title][name]=(val,ehi,elo,corr)
            else:
                poiResults[title]=(val,ehi,elo)

            var = iter.Next()

    varNames=[]
    for title,_,url in url_list:
        varNames += [x for x in fitResults[title].keys()]
    varNames=sorted(list(set(varNames)))
    print('{} nuisances to be plotted',len(varNames))

    #save in a pandas dataframe
    summary=[]
    for t,fr in fitResults.items():
        for name,vals in fr.items():
            summary.append([t,name,vals[0],vals[1],vals[2]])
    pd.DataFrame(summary,columns=['title','name','val','ehi','elo']).to_hdf('nuisances%s.h5'%tag,key='nuisances')

    #show nuisances
    c=ROOT.TCanvas('c','c',500,500)
    c.SetLeftMargin(0)
    c.SetTopMargin(0)
    c.SetRightMargin(0)
    c.SetBottomMargin(0)

    c.cd()
    p1=ROOT.TPad('p1','p1',0,0,0.7,1)
    p1.SetLeftMargin(0.3)
    p1.SetTopMargin(0.15)
    p1.SetRightMargin(0.03)
    p1.SetBottomMargin(1.0)
    p1.SetGridy(True)    
    p1.Draw()

    c.cd()
    p2=ROOT.TPad('p2','p2',0.7,0,1.0,1)
    p2.SetLeftMargin(0.02)
    p2.SetTopMargin(0.15)
    p2.SetRightMargin(0.03)
    p2.SetBottomMargin(1.0)
    p2.SetGridy(True)    
    p2.Draw()

    ngroups=len(varNames)/npergroup
    if ngroups*npergroup<len(varNames) : 
        ngroups+=1

    for ig in range(ngroups):

        first=npergroup*ig
        last=min(npergroup*(ig+1),len(varNames))
        varList=varNames[first:last]
        print(varList)
    
        #prepare a frame
        npars=len(varList)

        p1.cd()
        frame=ROOT.TH2F('frame',';#hat{#theta};Nuisance parameter',1,-3,3,npars,0,npars)
        frame.SetDirectory(0)
        for ybin in range(npars):
            label=varList[ybin]
            frame.GetYaxis().SetBinLabel(ybin+1,'#color[%d]{%s}'%((ybin%2)*10+1,label))
        frame.Draw()
        frame.GetXaxis().SetTitleSize(0.05)
        frame.GetYaxis().SetTitleSize(0.05)
        frame.GetXaxis().SetLabelSize(0.04)
        frame.GetYaxis().SetLabelSize(0.04)
        frame.GetXaxis().SetTitleOffset(0.8)
        frame.GetXaxis().CenterTitle()
        frame.GetYaxis().SetTitleOffset(3.2)

        #1-sigma
        gr1s=ROOT.TGraph()
        gr1s.SetName('gr1s')
        gr1s.SetMarkerStyle(1)
        gr1s.SetMarkerColor(19)
        gr1s.SetLineColor(19)
        gr1s.SetFillStyle(1001)
        gr1s.SetFillColor(19)
        gr1s.SetPoint(0,-1,0)
        gr1s.SetPoint(1,-1,npars)
        gr1s.SetPoint(2,1,npars)
        gr1s.SetPoint(3,1,0)
        gr1s.SetPoint(4,-1,0)

        #2-sigma
        gr2s=gr1s.Clone('gr2s')
        gr2s.SetMarkerColor(18)
        gr2s.SetLineColor(18)
        gr2s.SetFillStyle(1001)
        gr2s.SetFillColor(18)
        gr2s.SetPoint(0,-2,0)
        gr2s.SetPoint(1,-2,npars)
        gr2s.SetPoint(2,2,npars)
        gr2s.SetPoint(3,2,0)
        gr2s.SetPoint(4,-2,0)

        gr2s.Draw('f')
        gr1s.Draw('f')

        txt=ROOT.TLatex()
        txt.SetTextFont(42)
        txt.SetTextSize(0.03)
        txt.SetTextColor(ROOT.kGray+3)
        txt.SetNDC(True)
        for delta,title in [(0.72,'-1#sigma'),(0.82,'+2#sigma'),(0.5,'-1#sigma'),(0.4,'-2#sigma')]:
            txt.DrawLatex(delta,0.87,title)

        txt=ROOT.TLatex()
        txt.SetNDC(True)
        txt.SetTextFont(42)
        txt.SetTextSize(0.05)
        txt.SetTextAlign(12)
        txt.DrawLatex(0.05,0.955,'#bf{CMS} #it{Preliminary}')
        
        nuisGrs={}
        corrGrs={}
        markers=[20,24]
        colors=[1,ROOT.kGreen+2]
        for iv,v in enumerate(varList):

            for it,(title,_,url) in enumerate(url_list):

                #start graph if needed
                if not title in nuisGrs:
                    nuisGrs[title]=ROOT.TGraphAsymmErrors()
                    nuisGrs[title].SetTitle(title)
                    marker=markers[it]
                    ci=colors[it]
                    nuisGrs[title].SetMarkerStyle(marker)
                    nuisGrs[title].SetMarkerColor(ci)
                    nuisGrs[title].SetLineColor(ci)
                    nuisGrs[title].SetLineWidth(2)
                    nuisGrs[title].SetFillStyle(0)
                    corrGrs[title]=nuisGrs[title].Clone('corr_{}'.format(it))
                    corrGrs[title].SetFillStyle(1001)
                    corrGrs[title].SetFillColor(ci) #ROOT.kGray)

                npts=nuisGrs[title].GetN()

                val,uncLo,uncHi,rho=-99,0,0,0
                if v in fitResults[title]:
                    val,uncLo,uncHi,rho =fitResults[title][v]
                else:
                    print v,'is missing for',title,'!!!'
                y0=frame.GetYaxis().GetBinCenter(iv+1)
                dy=frame.GetYaxis().GetBinWidth(iv)

                corrGrs[title].SetPoint(npts,rho,y0+it*0.5)
                corrGrs[title].SetPointError(npts,rho,0,0.25*dy,0.25*dy)

                y0=y0+0.2*it-0.1*dy
                nuisGrs[title].SetPoint(npts,val,y0)
                nuisGrs[title].SetPointError(npts,abs(uncLo),abs(uncHi),0.,0.)

        leg=ROOT.TLegend(0.6,0.9,0.95,0.97)
        leg.SetTextFont(42)
        leg.SetTextSize(0.035)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        #leg.SetHeader('#bf{B-only fit}' if fit=='b' else '#bf{S+B fit}')
        for title in nuisGrs:
            nuisGrs[title].Draw('p')
            leg.AddEntry(nuisGrs[title],title,'p')
        leg.Draw()
        p1.RedrawAxis()

        p2.cd()
        frame2=ROOT.TH2F('frame2',';#rho;',1,-0.5,0.5,npars,0,npars)
        frame2.Draw()
        frame2.GetYaxis().SetTitleSize(0)
        frame2.GetYaxis().SetLabelSize(0)
        frame2.GetXaxis().SetTitleSize(0.12)
        frame2.GetXaxis().SetLabelSize(0.09)
        frame2.GetXaxis().SetTitleOffset(0.4)
        frame2.GetXaxis().SetLabelOffset(-0.04)
        frame2.GetXaxis().SetNdivisions(5)
        frame2.GetXaxis().CenterTitle()
        for title in corrGrs:
            corrGrs[title].Draw('2')

        txt2=ROOT.TLatex()
        txt2.SetNDC(True)
        txt2.SetTextFont(42)
        txt2.SetTextSize(0.09)
        txt2.SetTextAlign(22)
        for it,title in enumerate(poiResults):
            poiFit,uncUp,uncDown=poiResults[title]
            txt2.DrawLatex(0.5,0.95-0.05*it,'#hat{%s}_{%s} = %s'%(poi,title,toROOTRounded(poiFit,((uncUp,uncDown),))))

        p2.RedrawAxis()

        c.cd()
        c.Modified()
        c.Update()

        for ext in ['png','pdf']:
            c.SaveAs('nuisances%s_%d.%s'%(tag,ig,ext))

    return

def main():

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    parser = optparse.OptionParser(usage='usage: %prog [opts] ', version='%prog 1.0')
    parser.add_option('-t', '--tag',          
                      dest='tag',       
                      help='tag [%default]',  
                      default='',
                      type='string')

    (opt, args) = parser.parse_args()

    url_list=[]
    for a in args:
        url_list.append(a.split(':'))

    doNuisanceReport(url_list,opt.tag)


if __name__ == "__main__":
    sys.exit(main())
