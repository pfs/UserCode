import ROOT
import sys

url=sys.argv[1]
tag=url.split('_')[-1]
fIn=ROOT.TFile.Open(url)
nopu=fIn.Get('mmass_eeZelppnopu')
nopu_highpur=fIn.Get('mmass_eeZelpphighPurnopu')
pu=fIn.Get('mmass_eeZelpp')
pu_highpur=fIn.Get('mmass_eeZelpphighPur')

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
c=ROOT.TCanvas('c','c',500,500)
c.SetLeftMargin(0.12)
c.SetRightMargin(0.02)
c.SetTopMargin(0.05)
c.SetBottomMargin(0.1)
leg=ROOT.TLegend(0.2,0.9,0.4,0.7)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.SetTextFont(42)

nopu.SetLineWidth(2)
nopu.Draw('hist')
leg.AddEntry(nopu,'no PU','l')
pu.SetLineWidth(2)
pu.SetLineColor(2)
pu.Draw('histsame')
leg.AddEntry(pu,'#geq1p (PU)','l')
pu_highpur.SetLineWidth(2)
pu_highpur.SetLineColor(2)
pu_highpur.SetLineStyle(6)
pu_highpur.Draw('histsame')
leg.AddEntry(pu_highpur,'=1p (PU)','l')
leg.Draw()
tex=ROOT.TLatex()
tex.SetTextFont(42)
tex.SetTextSize(0.04)
tex.SetNDC()
tex.DrawLatex(0.12,0.96,'#bf{CMS} #it{simulation preliminary}')
c.Modified()
c.Update()   
c.RedrawAxis()
for ext in ['png','pdf']:
    c.SaveAs('mmass_{0}_sig.{1}'.format(tag,ext))

print pu_highpur.Integral()/nopu_highpur.Integral()
