import ROOT
import os
import sys
import optparse
import numpy as np
from scipy.ndimage import gaussian_filter1d

GTLIST=(0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.85, 1.9, 1.95, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 4.0)
MTLIST=(169.5, 170, 170.5, 171, 171.5, 172, 172.5, 173, 173.5, 174, 174.5, 175, 175.5)

def getSummedDist(url,dist,ch=['ee','em','mm']):

    """sums up all channels for a given distribution"""

    h=None    
    try:
        fIn=ROOT.TFile.Open(url)
        for c in ch:
            if h:
                h.Add(fIn.Get(c+dist))
            else:
                h=fIn.Get(c+dist)
                h.SetDirectory(0)
        fIn.Close()
    except:
        pass

    return h

def estimateLocalSensitivity(dist,opt):

    """ steers the estimation of the local sensitivity for a given distribution """

    #get main histogram
    url=os.path.join(opt.input,'MC13TeV_2016_TTJets.root')
    h=getSummedDist(url,dist)
    
    #local  variation graphs will be approximated by pol1
    gtEvol=[]
    mtEvol=[]
    for xbin in range(h.GetNbinsX()):
        gtEvol.append(ROOT.TGraph())
        mtEvol.append(ROOT.TGraph())
    gfunc=ROOT.TF1('grad','[0]*x+[1]',-100,100)

    #scan in width
    m=172.5
    for g in GTLIST:
        gidx=int((g-0.7)/0.01)
        midx=int((m-169)/0.25)
        flag=((midx<<16)|(gidx))
        url=os.path.join(opt.input,'scenario%d/MC13TeV_2016_TTJets.root'%flag)
        hvar=getSummedDist(url,dist)
        if not hvar: continue
        #hvar.Scale(h.Integral()/hvar.Integral())
        for xbin in range(h.GetNbinsX()):
            rel_diff=100.*(hvar.GetBinContent(xbin+1)/h.GetBinContent(xbin+1)-1.)
            gtEvol[xbin].SetPoint(gtEvol[xbin].GetN(),100.*(g/1.31-1.),rel_diff)
    
    #local sensitivity to the width
    
    #scan in mass
    g=1.2
    for m in MTLIST:
        gidx=int((g-0.7)/0.01)
        midx=int((m-169)/0.25)
        flag=((midx<<16)|(gidx))
        url=os.path.join(opt.input,'scenario%d/MC13TeV_2016_TTJets.root'%flag)
        hvar=getSummedDist(url,dist)
        if not hvar : continue
        #hvar.Scale(h.Integral()/hvar.Integral())
        for xbin in range(h.GetNbinsX()):
            rel_diff=100.*(hvar.GetBinContent(xbin+1)/h.GetBinContent(xbin+1)-1.)
            mtEvol[xbin].SetPoint(mtEvol[xbin].GetN(),100.*(m/172.5-1.),rel_diff)
    
    #local sensitivities
    gls=h.Clone('gls')
    gls.Reset('ICE')
    gls.SetFillStyle(0)
    gls.SetLineWidth(2)
    gls.SetLineColor(ROOT.kAzure-2)
    mls=h.Clone('mls')
    mls.Reset('ICE')
    mls.SetFillStyle(0)
    mls.SetLineWidth(2)
    mls.SetLineColor(ROOT.kRed+1)
    for xbin in range(h.GetNbinsX()):
        gtEvol[xbin].Sort()
        gtEvol[xbin].Fit(gfunc,'MQ+')
        grad=gfunc.GetParameter(0)
        gls.SetBinContent(xbin+1,grad*50)

        mtEvol[xbin].Sort()
        mtEvol[xbin].Fit(gfunc,'MQ+')
        grad=gfunc.GetParameter(0)
        mls.SetBinContent(xbin+1,grad)

    #display results
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetBottomMargin(0.1)
    c.SetLeftMargin(0.15)
    c.SetRightMargin(0.04)
    c.SetGridy()
    gls.Draw('hist')
    gls.SetTitle('#Gamma_{top} #scale[0.8]{(x50)}')
    gls.GetYaxis().SetTitle('d(N/N_{0}) / d(x/x_{0}) (x=m_{t}, #Gamma_{t}) [%]')
    gls.GetYaxis().SetRangeUser(-20,20)
    gls.GetYaxis().SetTitleOffset(1.2)
    gls.GetYaxis().SetTitleSize(0.045)
    mls.Draw('histsame')
    mls.SetTitle('m_{top}')
    leg=c.BuildLegend(0.2,0.95,0.4,0.8)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.05)
    leg.SetFillStyle(0)
    tex=ROOT.TLatex()
    tex.SetTextFont(42)
    tex.SetTextSize(0.04)
    tex.SetNDC()
    tex.DrawLatex(0.15,0.96,'#bf{CMS} #it{simulation preliminary}')
    tex.SetTextAlign(31)
    tex.DrawLatex(0.97,0.96,'13 TeV')  

    c.Modified()
    c.Update()
    c.RedrawAxis()
    outDir=opt.output if opt.output else './'
    for ext in ['png','pdf']:
        c.SaveAs('%slocalsens_%s.%s'%(outDir,dist,ext))

def main():

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--in',          
                      dest='input',       
                      help='input directory [%default]',  
                      default='/eos/cms/store/cmst3/group/top/TOP17010/0c522df',
                      type='string')
    parser.add_option('-o', '--out',          
                      dest='output',
                      help='output directory [%default]',  
                      default=None,
                      type='string')
    parser.add_option('-d', '--dist',          
                      dest='dist',       
                      help='distribution',
                      default='_mlb,highpt2b_mlb,highpt1b_mlb,highpt_mlb,lowpt2b_mlb,lowpt1b_mlb,lowpt_mlb',
                      type='string')
    (opt, args) = parser.parse_args()

    if opt.output: os.system('mkdir -p %s'%opt.output)
    for d in opt.dist.split(','): estimateLocalSensitivity(d,opt)

if __name__ == "__main__":
    sys.exit(main())
