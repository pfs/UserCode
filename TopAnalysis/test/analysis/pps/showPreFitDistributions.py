import ROOT
import sys

fIn=ROOT.TFile.Open(sys.argv[1])
tag=sys.argv[2]
title='Z#rightarrow#mu#mu'
if '22' in tag : title='#gamma'
if '121' in tag: title='Z#rightarrowee'
if '120' in tag: title+= ', 120 #murad'
if '130' in tag: title+= ', 130 #murad'
if '140' in tag: title+= ', 140 #murad'
if '150' in tag: title+= ', 150 #murad'

frame=None
histos={}
for key in fIn.GetListOfKeys():
    name=key.GetName()
    tkns=name.split('_')
    cat=None    
    if tkns[-1].isdigit(): cat=int(tkns[-1])
    if tkns[-2].isdigit(): cat=int(tkns[-2])
    if cat is None : continue

    if not cat in histos:
        histos[cat]={'sig':[],'bkg':[]}
    proc='sig' if 'sig' in name else 'bkg'

    ci=ROOT.kRed if proc=='sig' else 1
    lw=1 if 'Shape' in name or 'Calib' in name else 2
    if 'Calib' in name : ci=9
    tit='Signal' if proc=='sig' else 'Background'
    
    histos[cat][proc].append( key.ReadObj() )
    histos[cat][proc][-1].SetTitle(tit)
    histos[cat][proc][-1].SetLineColor(ci)
    histos[cat][proc][-1].SetLineWidth(lw)
    histos[cat][proc][-1].SetMarkerColor(ci)

    if frame: continue
    frame=histos[cat][proc][-1].Clone('frame')
    frame.Reset('ICE')

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gROOT.SetBatch(True)
c=ROOT.TCanvas('c','c',500,500)
c.SetTopMargin(0.05)
c.SetBottomMargin(0.1)
c.SetLeftMargin(0.15)
c.SetRightMargin(0.04)
for cat in histos:

    frame.Draw()
    leg=ROOT.TLegend(0.7,0.9,0.9,0.7)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.035)
    leg.SetFillStyle(0)
    maxY=1
    for proc in histos[cat]:
        histos[cat][proc].sort(key=lambda x: len(x.GetName()),reverse=True)
        nhistos=len(histos[cat][proc])
        for i in range(nhistos):
            if i==nhistos-1 : 
                leg.AddEntry(histos[cat][proc][i],histos[cat][proc][i].GetTitle(),'l')
            histos[cat][proc][i].Draw('histsame')
            maxY=max(maxY,histos[cat][proc][i].GetMaximum())
    frame.GetYaxis().SetRangeUser(0.01,1.2*maxY)
    frame.GetYaxis().SetTitle('Events')
    frame.GetXaxis().SetTitle('Missing mass [GeV]')
    leg.Draw()
    tex=ROOT.TLatex()
    tex.SetTextFont(42)
    tex.SetTextSize(0.04)
    tex.SetNDC()
    tex.DrawLatex(0.15,0.96,'#bf{CMS} #it{preliminary}')
    tex.SetTextAlign(31)
    tex.DrawLatex(0.97,0.96,'13 TeV')  
    tex.DrawLatex(0.95,0.9,title)
    #c.SetLogy()
    c.Modified()
    c.Update()
    c.RedrawAxis()
    for ext in ['png','pdf']:
        c.SaveAs('shapescat%d_%s.%s'%(cat,tag,ext))


fIn.Close()
