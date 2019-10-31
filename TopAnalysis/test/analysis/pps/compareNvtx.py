import ROOT
import os

BASEDIR='/eos/cms/store/cmst3/user/psilva/ExclusiveAna/final/ab05162/Chunks/'
colors=[ROOT.kBlack,ROOT.kGreen+1,ROOT.kAzure+1,ROOT.kRed+1,ROOT.kGray]

def cmsHeader():
    tex=ROOT.TLatex()
    tex.SetTextFont(42)
    tex.SetTextSize(0.04)
    tex.SetNDC()
    tex.DrawLatex(0.15,0.96,'#bf{CMS} #it{preliminary}')
    tex.SetTextAlign(ROOT.kHAlignRight+ROOT.kVAlignBottom)
    tex.DrawLatex(0.97,0.96,'37.5 fb^{-1} (13 TeV)')  


def main():
    
    tree=ROOT.TChain('tree')
    treePerEra={}
    for e in 'BCDEF': treePerEra[e]=ROOT.TChain('tree')
    fList=[os.path.join(BASEDIR,f) for f in os.listdir(BASEDIR) if 'DoubleMuon' in f]
    for f in fList:
        tree.AddFile(f)
        era='B'
        if '2017C' in f : era='C'
        if '2017D' in f : era='D'
        if '2017E' in f : era='E'
        if '2017F' in f : era='F'
        treePerEra[era].AddFile(f)


    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    #ROOT.gROOT.SetBatch(True)
    c=ROOT.TCanvas('c','c',500,500)
    c.SetRightMargin(0.02)
    c.SetTopMargin(0.05)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.1)
    
    def showHistos(histos,pname):
        for i in range(len(histos)):
            histos[i].SetLineWidth(2)
            histos[i].SetLineColor(colors[i])
            histos[i].Draw('hist' if i==0 else 'histsame')
            histos[i].GetXaxis().SetTitle('Vertex multiplicity')
            histos[i].GetYaxis().SetTitle('PDF')
            histos[i].GetYaxis().SetRangeUser(0,0.08)
        
        leg=c.BuildLegend(0.7,0.94,0.95,0.8)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.035)

        cmsHeader()
        c.Modified()
        c.Update()
        raw_input(pname)
        for ext in ['png','pdf']:
            c.SaveAs(pname+'.'+ext)

    #per era
    histos=[]
    for era in 'BCDEF':
       treePerEra[era].Draw('nvtx >> h(60,0,60)','','goff')
       h=treePerEra[era].GetHistogram()
       histos.append(h.Clone('nvtx'+era))
       histos[-1].SetTitle('2017'+era)
       histos[-1].Scale(1./histos[-1].Integral())
       h.Reset('ICE')
    showHistos(histos,'nvtxvsera')


    #per crossing angle
    histos=[]
    for xangle in [120,130,140,150]:
        tree.Draw('nvtx >> h(60,0,60)','beamXangle==%d'%xangle,'goff')
        h=tree.GetHistogram()
        histos.append(h.Clone('nvtx%d'%xangle) )
        histos[-1].SetTitle('%d #murad'%xangle)
        histos[-1].Scale(1./histos[-1].Integral())
        h.Reset('ICE')
    showHistos(histos,'nvtxvsangle')

    histos=[]
    for cat,cut in [('=0 protons','nRPtk==0'),
                    ('=1 proton', 'nRPtk==2 && RPid[0]!=RPid[1]'),
                    ('#geq 1 proton', 'nRPtk>2')]:
        tree.Draw('nvtx >> h(60,0,60)',cut,'goff')
        h=tree.GetHistogram()
        histos.append(h.Clone('nvtxvsprotons%d'%xangle) )
        histos[-1].SetTitle(cat)
        histos[-1].Scale(1./histos[-1].Integral())
        h.Reset('ICE')
    showHistos(histos,'nvtxvsproton')

if __name__ == "__main__":
    main()
