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
p2.SetBottomMargin(0.2)
p2.Draw()

def getPlotFrom(pname,fIn):
    return fIn.Get(pname+'/'+pname)


def showComparison(a,b,data,outName,outDir):

    p1.cd()
    data.SetMarkerStyle(20)
    a.SetLineColor(ROOT.kMagenta-3)
    a.SetLineWidth(2)
    b.SetLineColor(ROOT.kAzure-3)
    b.SetLineWidth(2)
    b.GetYaxis().SetRangeUser(1e-1,b.GetMaximum()*1.2)
    b.GetXaxis().SetTitleSize(0)
    b.GetXaxis().SetLabelSize(0)
    b.GetYaxis().SetTitleSize(0.05)
    b.GetYaxis().SetLabelSize(0.05)

    b.Draw('hist')
    a.Draw('histsame')
    data.Draw('e1same')
    leg=ROOT.TLegend(0.6,0.9,0.95,0.75)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.045)
    leg.AddEntry(data,data.GetTitle(),'ep')
    leg.AddEntry(a,a.GetTitle(),'l')
    leg.AddEntry(b,b.GetTitle(),'l')
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
    ratio=data.Clone('ratio')
    ratio.Divide(b)
    ratio.SetLineColor(b.GetLineColor())
    ratio.SetMarkerColor(b.GetLineColor())
    ratio.Draw('e1')
    ratio2=data.Clone('ratio2')
    ratio2.Divide(a)
    ratio2.Draw('e1same')
    ratio2.SetLineColor(a.GetLineColor())
    ratio2.SetMarkerColor(a.GetLineColor())
    ratio.GetYaxis().SetTitle('Ratio')
    ratio.GetYaxis().SetRangeUser(0.42,1.52)
    ratio.GetYaxis().SetNdivisions(5)
    ratio.GetXaxis().SetTitleSize(0.08)
    ratio.GetXaxis().SetLabelSize(0.08)
    ratio.GetYaxis().SetTitleSize(0.08)
    ratio.GetYaxis().SetLabelSize(0.08)
    ratio.GetYaxis().SetTitleOffset(0.95)

    p1.RedrawAxis()
    p2.RedrawAxis()
    c.Modified()
    c.Update()

    for ext in ['png','pdf']:
        c.SaveAs('%s/%s.%s'%(outDir,outName,ext))

    ratio.Delete()


fIn=ROOT.TFile.Open(sys.argv[1])
outDir=sys.argv[2]

for dilcat in ['mm']:
    for subcat in ['Z','hptZ','hptZhighpur']:
        for xangle in [120,130,140]:
            for side in ['neg','pos']:
                
                for dist in ['csi','ntk']:
                    a=getPlotFrom('{0}_{1}mix2{2}{3}_{4}'.format(dist,dilcat,subcat,xangle,side),fIn)
                    b=getPlotFrom('{0}_{1}mix1{2}{3}_{4}'.format(dist,dilcat,subcat,xangle,side),fIn)
                    data=getPlotFrom('{0}_{1}{2}{3}_{4}'.format(dist,dilcat,subcat,xangle,side),fIn)
                    a.Scale(data.Integral()/a.Integral())
                    b.Scale(data.Integral()/b.Integral())
                    a.SetTitle('2p mix')
                    b.SetTitle('1p mix')
                    data.SetTitle('Data')
                    showComparison(a,b,data,'bkg_'+data.GetName(),outDir)
                                   
            for dist in ['mpp']:
                a=getPlotFrom('{0}_{1}mix2{2}{3}'.format(dist,dilcat,subcat,xangle),fIn)
                b=getPlotFrom('{0}_{1}mix1{2}{3}'.format(dist,dilcat,subcat,xangle),fIn)
                data=getPlotFrom('{0}_{1}{2}{3}'.format(dist,dilcat,subcat,xangle),fIn)
                a.Scale(data.Integral()/a.Integral())
                b.Scale(data.Integral()/b.Integral())
                a.SetTitle('2p mix')
                b.SetTitle('1p mix')
                data.SetTitle('Data')
                showComparison(a,b,data,'bkg_'+data.GetName(),outDir)
        
