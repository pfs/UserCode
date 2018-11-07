import ROOT
import sys

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gROOT.SetBatch(True)

c=ROOT.TCanvas('c','c',500,500)
c.SetLeftMargin(0)
c.SetRightMargin(0)
c.SetTopMargin(0)
c.SetBottomMargin(0)
c.cd()
p1=ROOT.TPad('p1','p1',0.,0.4,1.0,1.0)
p1.SetLeftMargin(0.12)
p1.SetRightMargin(0.02)
p1.SetTopMargin(0.08)
p1.SetBottomMargin(0.01)
p1.Draw()
c.cd()
p2=ROOT.TPad('p2','p2',0.,0.,1.0,0.4)
p2.SetLeftMargin(0.12)
p2.SetRightMargin(0.02)
p2.SetTopMargin(0.01)
p2.SetBottomMargin(0.12)
p2.Draw()

def getPlotFrom(pname,fIn):
    return fIn.Get(pname+'/'+pname)


def showComparison(a,b,outName,outDir):

    p1.cd()
    a.SetMarkerStyle(20)
    a.SetTitle('Data')
    b.SetFillStyle(1001)
    b.SetFillColor(ROOT.kGray)
    b.SetLineColor(ROOT.kGray)
    b.SetTitle('Event mixing')
    b.GetYaxis().SetRangeUser(1e-1,b.GetMaximum()*1.2)
    b.GetXaxis().SetTitleSize(0)
    b.GetXaxis().SetLabelSize(0)
    b.GetYaxis().SetTitleSize(0.04)
    b.GetYaxis().SetLabelSize(0.04)

    b.Draw('hist')
    a.Draw('e1same')
    leg=ROOT.TLegend(0.6,0.9,0.95,0.75)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.045)
    leg.AddEntry(a,a.GetTitle(),'ep')
    leg.AddEntry(b,b.GetTitle(),'f')
    leg.Draw()

    txt=ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(42)
    txt.SetTextSize(0.06)
    txt.SetTextAlign(12)
    txt.DrawLatex(0.15,0.87,'#bf{CMS} #it{preliminary}')
    txt.SetTextAlign(32)
    txt.DrawLatex(0.96,0.96,'#scale[0.9]{37.5 fb^{-1} (13 TeV)}')
    p1.RedrawAxis()

    p2.cd()
    p2.SetGridy()
    ratio=a.Clone('ratio')
    ratio.Divide(b)
    ratio.Draw('e1')
    ratio.GetYaxis().SetTitle('Ratio')
    ratio.GetYaxis().SetRangeUser(0.62,1.38)
    ratio.GetXaxis().SetTitleSize(0.06)
    ratio.GetXaxis().SetLabelSize(0.06)
    ratio.GetYaxis().SetTitleSize(0.06)
    ratio.GetYaxis().SetLabelSize(0.06)

    p1.RedrawAxis()
    p2.RedrawAxis()
    c.Modified()
    c.Update()

    for ext in ['png','pdf']:
        c.SaveAs('%s/%s.%s'%(outDir,outName,ext))

    ratio.Delete()


fIn=ROOT.TFile.Open(sys.argv[1])
outDir=sys.argv[2]

for dilcat in ['ee','mm','em']:
    for subcat in ['','Z','hpt','hptZ','Zhighpur','hpthighpur','hptZhighpur']:
        if dilcat=='em' and 'Z' in subcat: continue
        ch=dilcat+subcat
        for side in ['neg','pos']:
            for dist in ['csi','ntk']:
                a=getPlotFrom('%s_%s_%s'%(dist,ch,side),fIn)
                b=getPlotFrom('%s_%s_mix%s'%(dist,ch,side),fIn)
                showComparison(a,b,'bkg_%s_%s_%s'%(dist,ch,side),outDir)
                                   
        for dist in ['mpp','mmass']:
            a=getPlotFrom('%s_%s'%(dist,ch),fIn)
            b=getPlotFrom('%s_%s_mix'%(dist,ch),fIn)
            showComparison(a,b,'bkg_%s_%s'%(dist,ch),outDir)
        
