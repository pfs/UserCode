import ROOT
import os
from runExclusiveAnalysis import getTracksPerRomanPot

def cmsHeader():
    txt=ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(42)
    txt.SetTextSize(0.04)
    txt.SetTextAlign(12)
    txt.DrawLatex(0.15,0.97,'#bf{CMS} #it{preliminary}')
    txt.SetTextAlign(32)
    txt.DrawLatex(0.96,0.97,'#scale[0.9]{37.5 fb^{-1} (13 TeV)}')


baseDir='/eos/cms/store/cmst3/user/psilva/ExclusiveAna/Chunks/'
tree=ROOT.TChain('tree')
for f in os.listdir(baseDir):
    if not 'Data13TeV_2017F_DoubleMuon' in f : continue
    tree.AddFile(os.path.join(baseDir,f))
    break

histos={'csi23':ROOT.TH1F('csi23','RP23;#xi;PDF',50,0,0.3),
        'csi123':ROOT.TH1F('csi123','RP123;#xi;PDF',50,0,0.3),
        'csi2d':ROOT.TH2F('csi2d',';#xi_{RP23};#xi_{RP123};PDF',50,0,0.3,50,0,0.3),
        'n23':ROOT.TH1F('n23','RP23;Track multiplicity;PDF',5,0,5),
        'n123':ROOT.TH1F('n123','RP123;Track multiplicity;PDF',5,0,5),
        'n2d':ROOT.TH2F('n2d',';Track multiplicity (RP23); Track multiplicity (RP123); PDF',5,0,5,5,0,5)
        }

for i in range(tree.GetEntries()):
    tree.GetEntry(i)
    rp23,rp123=getTracksPerRomanPot(tree)
    histos['n2d'].Fill(len(rp23),len(rp123))
    histos['n23'].Fill(len(rp23))
    histos['n123'].Fill(len(rp123))
    for x in rp23: histos['csi23'].Fill(x)
    for x in rp123: histos['csi123'].Fill(x)
    if len(rp23)==1 and len(rp123)==1:
        histos['csi2d'].Fill(rp23[0],rp123[0])

for key in histos:
    histos[key].SetDirectory(0)
    histos[key].Sumw2()
    histos[key].Scale(1./histos[key].Integral())

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
c=ROOT.TCanvas('c','c',500,500)
c.SetTopMargin(0.05)
c.SetBottomMargin(0.1)
c.SetLeftMargin(0.15)
for d in ['csi','n']:
    c.SetRightMargin(0.03)
    histos['%s23'%d].Draw('e1')
    histos['%s23'%d].SetLineColor(1)
    histos['%s23'%d].SetMarkerColor(1)
    histos['%s23'%d].SetMarkerStyle(20)

    histos['%s123'%d].Draw('e1same')
    histos['%s123'%d].SetLineColor(2)
    histos['%s123'%d].SetMarkerColor(2)
    histos['%s123'%d].SetMarkerStyle(24)

    leg=c.BuildLegend(0.75,0.94,0.95,0.8)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)

    cmsHeader()

    c.Modified()
    c.Update()
    c.SaveAs('%s.pdf'%d)


    c.SetRightMargin(0.15)
    ROOT.gStyle.SetPaintTextFormat("4.2f");
    histos['%s2d'%d].Draw('colztext' if d =='n' else 'colz')
    cmsHeader()
    c.Modified()
    c.Update()
    c.SaveAs('%s2d.pdf'%d)

    


