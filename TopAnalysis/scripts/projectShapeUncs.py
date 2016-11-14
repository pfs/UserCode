from plotter import *
import sys
import ROOT

def projectShapeUncs(url,proc,systList):

    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    fIn=ROOT.TFile.Open(url)

    c=ROOT.TCanvas('c','c',500,500)
    c.SetRightMargin(0.05)
    c.SetTopMargin(0.05)
    c.SetBottomMargin(0.1)
    c.SetLeftMargin(0.12)

    nomH=fIn.Get('nom/%s'%proc)
    nomH.SetLineWidth(2)
    nomH.SetLineColor(1)
    nomH.SetFillStyle(0)
    nomH.Draw('hist')
    nomH.GetYaxis().SetTitleOffset(1.2)
    nomH.GetYaxis().SetRangeUser(1e-3,nomH.GetMaximum()*2.0)

    leg= ROOT.TLegend(0.65,0.5,0.85,0.95)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.035)
    leg.SetHeader('#bf{CMS} #it{simulation}')
    leg.AddEntry(nomH,'Nominal','l')

    colors=[2,4,8,9,36,3]
    allGr=[]
    nSysts=len(systList)
    for isyst in xrange(0,nSysts):
        key=systList[isyst]
        color=colors[isyst]
        hup=fIn.Get('%sUp/%s'%(key,proc))
        hdn=fIn.Get('%sDown/%s'%(key,proc))
        allGr.append( ROOT.TGraphAsymmErrors() )
        for xbin in xrange(1,hup.GetNbinsX()+1):
            valUp=hup.GetBinContent(xbin)
            valDn=hdn.GetBinContent(xbin)
            valCen=nomH.GetBinContent(xbin)
            binw=hdn.GetXaxis().GetBinWidth(xbin)
            binc=hdn.GetXaxis().GetBinCenter(xbin)

            dx=(binw*0.9/nSysts)
            if isyst>0:
                if isyst%2:
                    binc -= dx*(isyst/2+1)
                else:
                    binc += dx*(isyst/2)

            allGr[-1].SetPoint(xbin,binc,valCen)
            ylo=ROOT.TMath.Min(valUp-valCen,valDn-valCen)
            yhi=ROOT.TMath.Max(valUp-valCen,valDn-valCen)
            allGr[-1].SetPointError(xbin,dx*0.5,dx*0.5,ROOT.TMath.Abs(ylo),ROOT.TMath.Abs(yhi))

        allGr[-1].SetTitle(key)
        allGr[-1].SetName(key)
        allGr[-1].SetMarkerStyle(0)
        allGr[-1].SetMarkerColor(color)
        allGr[-1].SetFillStyle(1001)
        allGr[-1].SetFillColor(color)
        allGr[-1].Draw('2')
        leg.AddEntry(allGr[-1],key,'f')

    leg.Draw()
    c.Modified()
    c.Update()
    outName=os.path.splitext( os.path.basename(url) )[0]
    dirName=os.path.dirname(url)
    for ext in ['pdf','png']:
        c.SaveAs(os.path.join(dirName,'%s_%s_unc.%s'%(outName,key,ext)))
    fIn.Close()

def main():

    ROOT.gROOT.SetBatch(True)
    url = sys.argv[1]
    systList=sys.argv[2]
    proc='tbart' if len(sys.argv)<4 else sys.argv[3]
    projectShapeUncs(url,proc,systList.split(','))
    
"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
