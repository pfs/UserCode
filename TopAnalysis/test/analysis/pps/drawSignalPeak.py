import ROOT
import os
import sys
import re

url=sys.argv[1]
tkns=re.findall(r'\d+',os.path.basename(url))
mass,xangle=tkns[0],tkns[1]
print mass,xangle
fIn=ROOT.TFile.Open(url)
nopu=fIn.Get('mmass_eerpinnopu')
nopu_highpur=fIn.Get('mmass_eerpinhurnopu')
pu=fIn.Get('mmass_eerpin')
pu_highpur=fIn.Get('mmass_eerpinhpur')
data=fIn.Get('data')



ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
c=ROOT.TCanvas('c','c',500,500)
c.SetLeftMargin(0.12)
c.SetRightMargin(0.04)
c.SetTopMargin(0.05)
c.SetBottomMargin(0.1)
leg=ROOT.TLegend(0.65,0.9,0.9,0.7)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.SetTextFont(42)

nopu.SetLineWidth(2)
nopu.GetYaxis().SetTitleOffset(1.3)
nopu.Draw('hist')
nopu.GetYaxis().SetRangeUser(0,1.5*pu.GetMaximum())
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
tex.SetTextAlign(31)
tex.DrawLatex(0.95,0.96,'#scale[0.8]{{m={0} GeV #alpha={1}#murad}}'.format(mass,xangle))
c.Modified()
c.Update()   
c.RedrawAxis()
for ext in ['png','pdf']:
    c.SaveAs('mmass_{0}_xangle{1}_sig.{2}'.format(mass,xangle,ext))


data.Draw("mixmmiss>>hout(50,0,2500)","wgt*((gen_pzpp<-300 || gen_pzpp>500) && mixmmiss>0)")
hout=c.GetPrimitive('hout')
hout.SetDirectory(0)
data.Draw("mixmmiss>>hin(50,0,2500)","wgt*(gen_pzpp>-300 && gen_pzpp<500 && mixmmiss>0)","histsame")
hin=c.GetPrimitive('hin')
hin.SetDirectory(0)
hout.SetLineColor(1)
hin.SetLineColor(2)
hout.SetFillStyle(1001)
hout.SetFillColor(ROOT.kGray)

hout.Draw("hist")
hin.Draw("histsame")

hout.GetYaxis().SetTitleOffset(1.3)
hout.GetYaxis().SetRangeUser(0,1.5*max(hout.GetMaximum(),hin.GetMaximum()))
hout.GetXaxis().SetTitle('Missing mass [GeV]')
hout.GetYaxis().SetTitle('PDF')
hin.SetLineWidth(2)
leg=ROOT.TLegend(0.65,0.9,0.9,0.7)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.SetTextFont(42)

leg.AddEntry(hin,'in acceptance','l')
leg.AddEntry(hout,'out acceptance','f')
leg.Draw()


tex=ROOT.TLatex()
tex.SetTextFont(42)
tex.SetTextSize(0.04)
tex.SetNDC()
tex.DrawLatex(0.12,0.96,'#bf{CMS} #it{simulation preliminary}')
tex.SetTextAlign(31)
tex.DrawLatex(0.95,0.96,'#scale[0.8]{{m={0} GeV #alpha={1}#murad}}'.format(mass,xangle))
c.Modified()
c.Update()   
c.RedrawAxis()
raw_input()
for ext in ['png','pdf']:
    c.SaveAs('mmass_{0}_xangle{1}_sigacc.{2}'.format(mass,xangle,ext))


