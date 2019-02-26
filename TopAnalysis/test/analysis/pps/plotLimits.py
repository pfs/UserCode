import ROOT
import pickle
import sys

massPts=[x.split('=') for x in sys.argv[1].split(',')]

r68=ROOT.TGraphAsymmErrors()
r68.SetFillColor(8)
r68.SetFillStyle(1001)
r68.SetMarkerColor(8)
r68.SetLineColor(1)
r68.SetLineWidth(2)
r68.SetLineStyle(7)
r68.SetName('r68')
r68.SetTitle('Expected (68%)')
r95=ROOT.TGraphAsymmErrors()
r95.SetFillColor(5)
r95.SetFillStyle(1001)
r95.SetMarkerColor(5)
r95.SetLineColor(1)
r95.SetLineWidth(2)
r95.SetLineStyle(7)
r95.SetName('r95')
r95.SetTitle('Expected (95%)')
rmed=ROOT.TGraph()
rmed.SetFillColor(0)
rmed.SetFillStyle(0)
rmed.SetLineColor(1)
rmed.SetLineWidth(2)
rmed.SetLineStyle(7)
sig1=rmed.Clone('sig1')
sig10=rmed.Clone('sig10')
sig10.SetLineColor(ROOT.kGray)
sig10.SetMarkerColor(ROOT.kGray)

bySignificance=False #True

for m,url in massPts:

    m=float(m)
    with open(url,'r') as cache:
        results=pickle.load(cache)

    np=r68.GetN()
    maxSig=0
    minR95=9999.
    bestKey=None
    for key in results:
        
        vals=results[key]
        if bySignificance:
            if vals[6]<maxSig: continue
        else:
            if vals[2]>minR95: continue

        bestKey=key
        minR95=vals[2]
        maxSig=vals[6]
        rmed.SetPoint(np,m,vals[2])
        r68.SetPoint(np,m,vals[2])
        r68.SetPointError(np,10,10,vals[2]-vals[1],vals[3]-vals[2])
        r95.SetPoint(np,m,vals[2])
        r95.SetPointError(np,10,100,vals[2]-vals[0],vals[4]-vals[2])
        sig1.SetPoint(np,m,vals[5])
        sig10.SetPoint(np,m,vals[6])
    print m,bestKey

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
c=ROOT.TCanvas('c','c',500,500)
c.SetLeftMargin(0.12)
c.SetRightMargin(0.03)
c.SetTopMargin(0.05)
c.SetBottomMargin(0.1)
mg=ROOT.TMultiGraph()
mg.Add(r95,'l3')
mg.Add(r68,'l3')
mg.Add(rmed,'l')
frame=ROOT.TH1F('frame',';m_{X} [GeV];95% CL limits on #sigma/#sigma_{fid}',1,800,2000)
frame.GetYaxis().SetRangeUser(0,10)
frame.SetBinContent(1,1)
frame.SetLineWidth(2)
frame.SetLineColor(ROOT.kRed)
frame.Draw('hist')
mg.Draw('l3')
leg=ROOT.TLegend(0.15,0.92,0.4,0.72)
leg.SetTextFont(42)
leg.SetTextSize(0.035)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.AddEntry(r68,r68.GetTitle(),'lf')
leg.AddEntry(r95,r95.GetTitle(),'lf')
leg.Draw()
tex=ROOT.TLatex()
tex.SetTextFont(42)
tex.SetTextSize(0.04)
tex.SetNDC()
tex.DrawLatex(0.12,0.96,'#bf{CMS} #it{preliminary}')
tex.SetTextAlign(ROOT.kHAlignRight+ROOT.kVAlignCenter)
tex.DrawLatex(0.97,0.975,'#scale[0.9]{37.5 fb^{-1} (13 TeV)}')
c.Modified()
c.Update()   
c.RedrawAxis()
for ext in ['png','pdf']:
    c.SaveAs('limits.%s'%ext)

c.Clear()
frame.Draw()
frame.GetYaxis().SetTitle('Signal asymptotic significance (#sigma)')
frame.SetBinContent(1,5)
sig1.Draw('l')
sig10.Draw('l')
leg=ROOT.TLegend(0.15,0.92,0.4,0.72)
leg.SetTextFont(42)
leg.SetTextSize(0.035)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.AddEntry(sig1,'#mu=1','lf')
leg.AddEntry(sig10,'#mu=10','lf')
leg.Draw()
tex=ROOT.TLatex()
tex.SetTextFont(42)
tex.SetTextSize(0.04)
tex.SetNDC()
tex.DrawLatex(0.12,0.96,'#bf{CMS} #it{preliminary}')
tex.SetTextAlign(ROOT.kHAlignRight+ROOT.kVAlignCenter)
tex.DrawLatex(0.97,0.975,'#scale[0.9]{37.5 fb^{-1} (13 TeV)}')
c.Modified()
c.Update()   
c.RedrawAxis()
for ext in ['png','pdf']:
    c.SaveAs('significance.%s'%ext)
