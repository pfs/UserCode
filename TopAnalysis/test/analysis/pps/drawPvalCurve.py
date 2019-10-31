import ROOT
import sys
import os
import re
from compareOptimResults import getSignificanceFrom,getLimitsFrom,showLimits

def fillPvalGraphFrom(files):

    gr=ROOT.TGraph()
    for f in files:    
        mass = int(re.search('.mH(\d+)', f).group(1))
        sig=getSignificanceFrom(f)
        gr.SetPoint(gr.GetN(),mass,sig[1])
    gr.Sort()
    return gr


def fillLimitsFrom(files,getObs=False):
    results=[]
    for f in files:
        mass = int(re.search('.mH(\d+)', f).group(1))
        results.append( (mass,getLimitsFrom(f,getObs)) )
    return results

def showPvalSummary(grpval):

    c=ROOT.TCanvas('c','c',600,600)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.03)
    c.SetTopMargin(0.05)
    c.SetBottomMargin(0.1)
    c.SetGridx()
    c.SetGridy()
    c.SetLogy()

    leg=ROOT.TLegend(0.65,0.15,0.95,0.3)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.03)

    mg=ROOT.TMultiGraph()
    for gr in grpval:
        mg.Add(gr)
        leg.AddEntry(gr,gr.GetTitle(),'lep')

    mg.Draw('apc')
    mg.GetYaxis().SetRangeUser(1e-8,1.5)
    mg.GetYaxis().SetTitle('Local p-value')
    mg.GetYaxis().SetTitleOffset(1.3)
    mg.GetXaxis().SetTitle('m_{X} [GeV]')

    leg.Draw()

    tex=ROOT.TLatex()
    tex.SetTextFont(42)
    tex.SetTextSize(0.04)
    tex.SetNDC()
    tex.DrawLatex(0.12,0.96,'#bf{CMS} #it{preliminary}')
    tex.SetTextAlign(ROOT.kHAlignRight+ROOT.kVAlignBottom)
    tex.DrawLatex(0.95,0.96,'#scale[0.8]{2.64-37.5 fb^{-1} (13 TeV)}')

    lines=[]
    l=ROOT.TLine()
    l.SetLineColor(ROOT.kGray+2)
    l.SetLineStyle(9)
    tex.SetTextSize(0.03)
    tex.SetNDC(False)
    tex.SetTextColor(ROOT.kGray+2)
    for i in range(6):
        pval=0.5*(1-ROOT.TMath.Erf(float(i)/ROOT.TMath.Sqrt(2.)))
        l.DrawLine(mg.GetXaxis().GetXmin(),pval,mg.GetXaxis().GetXmax(),pval)
        if i==0: continue
        tex.DrawLatex(mg.GetXaxis().GetXmax()-25,0.6*pval,'%d#sigma'%i)
    
    c.Modified()
    c.Update()
    for ext in ['png','pdf']:
        c.SaveAs('pval_summary.%s'%ext)


def main():

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(True)

    grpval=[]
    colors=[ROOT.kAzure+3,ROOT.kRed+1]
    lines=[1,9]
    markers=[20,24,21]

    for i in range(1,len(sys.argv)):
        key,resultsDir=sys.argv[i].split(':')
        mass=int(re.search('analysis_(\d+)', resultsDir).group(1))
        print key,resultsDir,mass
        files=[os.path.join(resultsDir,f) for f in os.listdir(resultsDir) if '.obs.Significance' in f and key in f]
        grpval.append( fillPvalGraphFrom(files) )
        grpval[-1].SetName(key+'_%d'%i)
        title=key
        title=title.replace('g','#gamma')
        title=title.replace('z','Z')
        title=title.replace('PP','pp#rightarrowpp')
        title +='(%d)'%mass
        grpval[-1].SetTitle(title)
        grpval[-1].SetLineWidth(2)

        ci=colors[1] if '#gamma' in title else colors[0]
        grpval[-1].SetLineColor(ci)
        grpval[-1].SetMarkerColor(ci)
        mk=markers[i%3]
        grpval[-1].SetMarkerStyle(mk)
        ls=lines[i%2]
        grpval[-1].SetLineStyle(ls)
        grpval[-1].SetFillStyle(0)

        obslimfiles=[os.path.join(resultsDir,f) for f in os.listdir(resultsDir) if '.obs.AsymptoticLimits' in f and key in f]
        obs_limits=fillLimitsFrom(obslimfiles,True)
        explimfiles=[os.path.join(resultsDir,f) for f in os.listdir(resultsDir) if 'AsymptoticLimits' in f and key in f and not '.obs.' in f]
        exp_limits=fillLimitsFrom(explimfiles,False)
        lumi=2.64 if '#gamma' in title else 37.5
        showLimits(exp_limits,'limits_%s_%d'%(key,mass),title,lumi,obs_limits)

    showPvalSummary(grpval)


if __name__ == "__main__":
    main()

