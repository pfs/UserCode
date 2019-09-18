import ROOT
import os

leg_dict={
'MC13TeV_2016_TTJetsTo2L2Nu_psweights':(1,'CP5'),
'MC13TeV_2016_TTJets_fsrdn':(ROOT.kMagenta,'FSR up'),
'MC13TeV_2016_TTJets_fsrup':(ROOT.kMagenta+2,'FSR dn'),
'MC13TeV_2016_TTJets_hdampup':(ROOT.kAzure+7,'hdamp up'),
'MC13TeV_2016_TTJets_hdampdn':(ROOT.kBlue-7,'hdamp dn'),
'MC13TeV_2016_TTJets_uedn':(ROOT.kViolet+2,'UE dn'),
'MC13TeV_2016_TTJets_ueup':(ROOT.kRed+1,'UE up'),
'MC13TeV_2016_TTJets_erdon':(ROOT.kGray,'ERD on'),
'MC13TeV_2016_TTJets_qcdbased':(ROOT.kGray+1,'QCD based'),
'MC13TeV_2016_TTJets_gluonmove':(ROOT.kGray+2,'Gluon move')
}



#read all corrections
gr={}
fIn=ROOT.TFile.Open('test/analysis/top17010/mc2mc_corrections.root')
for k in fIn.GetListOfKeys():
    d=k.ReadObj()
    for kk in d.GetListOfKeys():
        kkname=kk.GetName()
        if not kkname in gr:
            gr[kkname]={}
        gr[kkname][k.GetName()]=kk.ReadObj()
fIn.Close()


ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
c=ROOT.TCanvas('c','c',500,500)
c.SetLeftMargin(0.12)
c.SetTopMargin(0.05)
c.SetRightMargin(0.03)
c.SetBottomMargin(1.0)
c.SetLogx(True)

#display corrections
for k in gr:

    c.Clear()

    leg=ROOT.TLegend(0.15,0.93,0.95,0.93-0.06*len(gr[k])/3)
    leg.SetTextFont(42)
    leg.SetNColumns(3)
    leg.SetTextSize(0.035)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)

    mg=ROOT.TMultiGraph()
    for kk in gr[k]:
        ci,title=leg_dict[kk]
        gr[k][kk].SetLineColor(ci)
        gr[k][kk].SetLineWidth(2)
        gr[k][kk].SetTitle(title)
        mg.Add(gr[k][kk],'l')
        leg.AddEntry(gr[k][kk],title,'l')
    mg.Draw('ac')
    mg.GetYaxis().SetTitle('Ratio to nominal')
    mg.GetYaxis().SetRangeUser(0.92,1.08)
    mg.GetXaxis().SetRangeUser(30,500)
    mg.GetXaxis().SetTitle('Transverse momentum [GeV]')
    mg.GetXaxis().SetMoreLogLabels()

    leg.Draw()

    txt=ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(42)
    txt.SetTextSize(0.045)
    txt.SetTextAlign(12)
    txt.DrawLatex(0.12,0.97,'#bf{CMS} #it{Simulation Preliminary}')
    txt.SetTextAlign(31)
    txt.DrawLatex(0.98,0.965,'#scale[0.8]{35.6 fb^{-1} (13 TeV)}')
    
    c.RedrawAxis()
    c.Modified()
    c.Update()

    for ext in ['png','pdf']:
        c.SaveAs('%s.%s'%(k,ext))
