import ROOT
import numpy as np
import os

def divide2DHistograms(hnum,hden):

    hratio=hnum.Clone('%s_%s'%(hnum.GetName(),hden.GetName()))

    for xbin in range(1,hnum.GetNbinsX()+1):
        xcen=hnum.GetXaxis().GetBinCenter(xbin)
        xbin_den = hden.GetXaxis().FindBin( max(min(xcen,hden.GetXaxis().GetXmax()),hden.GetXaxis().GetXmin()) )

        for ybin in range(1,hnum.GetNbinsY()+1):
            ycen=hnum.GetYaxis().GetBinCenter(ybin)
            ybin_den = hden.GetYaxis().FindBin( max(min(ycen,hden.GetYaxis().GetXmax()),hden.GetYaxis().GetXmin()) )
            
            val_den=hden.GetBinContent(xbin_den,ybin_den)
            unc_den=hden.GetBinError(xbin_den,ybin_den)
            val_num=hnum.GetBinContent(xbin,ybin)
            unc_num=hnum.GetBinError(xbin,ybin)
                    
            ratio=val_num/val_den
            ratio_unc=ratio*np.sqrt((unc_den/val_den)**2+(unc_num/val_num)**2)
            hratio.SetBinContent(xbin,ybin,ratio)
            hratio.SetBinError(xbin,ybin,ratio_unc)

    return hratio


baseDir='/eos/user/p/psilva/www/EXO-19-009/scaleFactors'
os.system('mkdir -p %s'%baseDir)

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetPaintTextFormat(".2f");
c=ROOT.TCanvas('c','c',500,500)
c.SetLeftMargin(0.18)
c.SetTopMargin(0.05)
c.SetBottomMargin(0.12)
c.SetRightMargin(0.15)

for tag,rereco,ul,out in [ ('gamma',
                            'data/era2017/2017_PhotonsMVAwp90.root', 
                            'data/era2017/egammaEffi.txt_EGM2D_PHO_MVA90_UL17.root',
                            'data/era2017/egammaEffi.txt_EGM2D_PHO_MVA90_UL17_MCReReco.root'),
                           ('electron',
                            'data/era2017/2017_ElectronMVA90.root',
                            'data/era2017/egammaEffi.txt_EGM2D_MVA90iso_UL17.root',
                            'data/era2017/egammaEffi.txt_EGM2D_ELE_MVA90_UL17_MCReReco.root'),
                           ]:


    frereco=ROOT.TFile.Open(rereco)
    ful=ROOT.TFile.Open(ul)

    #compute the ratios UL/RE-RECO
    hratios={}
    for h in ['EGamma_EffMC2D','EGamma_SF2D']:
        hrereco=frereco.Get(h)
        hul=ful.Get(h)
        if 'SF2D' in h:
            hratios[h]=hul.Clone()
            hratios[h].SetDirectory(0)
            hratios[h].Multiply(hratios['EGamma_EffMC2D_ULtoReReco'])
            hratios[h].GetZaxis().SetTitle('Data UL / MC Re-Reco')
        else:
            h=h+'_ULtoReReco'
            hratios[h]=divide2DHistograms(hul,hrereco)
            hratios[h].GetZaxis().SetTitle('MC UL / MC Re-Reco')
            
        hratios[h].Draw('colztext')
        hratios[h].GetZaxis().SetTitleOffset(0.)
        c.SaveAs(baseDir+'/'+tag+'_'+h+'.png')

        if 'SF2D' in h:
            hsfratio=divide2DHistograms(hratios[h],hrereco)
            hsfratio.Draw('colztext')
            hsfratio.GetZaxis().SetTitle('SF / SF re-reco')
            hsfratio.GetZaxis().SetTitleOffset(0.)
            c.SaveAs(baseDir+'/'+tag+'_'+h+'_ratio.png')


    fout=ROOT.TFile.Open(out,'RECREATE')
    for h in hratios:
        hratios[h].SetDirectory(fout)
        hratios[h].Write(h)
    fout.Close()


fIn=ROOT.TFile.Open('data/era2017/TriggerSF_2017_UL.root')
for ch in ['ee','mumu']:
    h=fIn.Get('h2D_SF_%s_lepABpt_FullError'%ch)
    h.Draw('colztext')
    c.SetLogx()
    c.SetLogy()
    h.GetZaxis().SetRangeUser(0.9,1)
    c.SaveAs('%s/trig_sf_%s.pdf'%(baseDir,ch))
fIn.Close()

hcombSF=None
for i,hname in [('ID','NUM_TightID_DEN_genTracks_pt_abseta'),
                ('ISO','NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta')]:

    fIn=ROOT.TFile.Open('data/era2017/RunBCDEF_SF_Mu%s.root'%i)
    h=fIn.Get(hname)
    h.Print("all")
    if hcombSF is None:
        hcombSF=h.Clone('muidiso_SF')
        hcombSF.SetDirectory(0)
    else:
        for xbin in range(1,h.GetNbinsX()+1):
            for ybin in range(1,h.GetNbinsY()+1):
                sf=hcombSF.GetBinContent(xbin,ybin)
                sf*=h.GetBinContent(xbin,ybin)
                hcombSF.SetBinContent(xbin,ybin,sf)
    fIn.Close()

hcombSF.Draw('colztext')
c.SetLogx(False)
c.SetLogy(False)
hcombSF.GetZaxis().SetRangeUser(0.9,1)
c.SaveAs('%s/%s.pdf'%(baseDir,hcombSF.GetName()))



