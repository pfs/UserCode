import ROOT

def compareResults(hname):
    c=ROOT.TCanvas('c','c',500,500)
    c.SetLeftMargin(0.2)
    c.SetRightMargin(0.02)
    c.SetTopMargin(0.05)
    c.SetBottomMargin(0.1)
    c.SetGridx()
    c.SetGridy()
    leg=ROOT.TLegend(0.2,0.05,0.7,0.03)
    leg.SetNColumns(3)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)

    drawOpt='hbar'
    histos=[]
    for m,ci,offs in [(1000,ROOT.kGray,0.05),
                      (1200,ROOT.kOrange,0.35),
                      (1400,ROOT.kGreen+1,0.65),
                      ] :
        fIn=ROOT.TFile.Open('optimresults_%d.root'%m)
        histos.append(fIn.Get(hname).Clone('%s_%d'%(hname,m)))
        histos[-1].SetDirectory(0)
        fIn.Close()
        histos[-1].SetBarWidth(0.3)
        histos[-1].SetBarOffset(offs)
        histos[-1].SetFillColor(ci)    
        histos[-1].Draw(drawOpt)
        histos[-1].GetXaxis().SetTitleOffset(2.4)    
        histos[-1].GetYaxis().SetTitleOffset(1.3)
        drawOpt='hbarsame'
        leg.AddEntry(histos[-1],'M_{X}=%d GeV'%m,'f')

    leg.Draw()
    tex=ROOT.TLatex()
    tex.SetTextFont(42)
    tex.SetTextSize(0.04)
    tex.SetNDC()
    tex.DrawLatex(0.2,0.96,'#bf{CMS} #it{simulation preliminary}')
    c.Modified()
    c.Update()
    for ext in ['png','pdf']:
        c.SaveAs('%s.%s'%(hname,ext))


ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gROOT.SetBatch(True)
for hname in ['optim','optim_sig1','optim_sig10']: compareResults(hname)
