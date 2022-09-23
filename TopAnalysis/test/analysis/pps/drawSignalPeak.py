import ROOT
import os
import sys
import re
from TopLJets2015.TopAnalysis.Plot import *
from generateBinnedWorkspace import SIGNALXSECS,VALIDLHCXANGLES

LUMI=37200.
FPRETS2=14586.4464/41529.3
FPOSTTS2=1.-FPRETS2
FIDUCIALCUTS='gencsi1>0.02 & gencsi1<0.16 && gencsi2>0.03 && gencsi2<0.18 && !isOffZ && (cat==121 || cat==169)'

def main():

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(True)

    
    isSingleRP=True if 'single' in sys.argv[1] else False
    isExc=True if '1exc' in sys.argv[1] else False
    pfix='pxonly' if isSingleRP else ''
    pfix += 'exc' if isExc else ''

    hnopu={}
    hrpinhpur={}
    hout={}
    hin={}
    for xangle in VALIDLHCXANGLES:
        for era in 'preTS2','postTS2':

            url=sys.argv[1].format(xangle,era)
            mass=re.findall(r'\d+',os.path.basename(url))[0]
            fEra=FPRETS2 if era=='preTS2' else FPOSTTS2
            xsecWgt=SIGNALXSECS[xangle]*LUMI*fEra
            print url
            #determine analysis type
            cats=[4] if isSingleRP else [1,2,3,4]

            #open file and build te plots
            fIn=ROOT.TFile.Open(url)
            data=fIn.Get('data')
            for icat in cats:
                
                if not icat in hnopu:
                    hnopu[icat]=fIn.Get('mmass_mm%dnopu'%icat).Clone('hnopu%d'%icat)
                    hnopu[icat].Reset('ICE')
                    hnopu[icat].SetDirectory(0)
                    hrpinhpur[icat]=fIn.Get('mmass_mmrpinhpur%d'%icat).Clone('hrpinhpur%d'%icat)
                    hrpinhpur[icat].Reset('ICE')
                    hrpinhpur[icat].SetDirectory(0)
                for ch in ['ee','mm']:
                    hnopu[icat].Add(fIn.Get('mmass_%s%dnopu'%(ch,icat)),xsecWgt)
                    hrpinhpur[icat].Add(fIn.Get('mmass_%srpinhpur%d'%(ch,icat)),xsecWgt)

                #
                # project in-fiducial vs out-fiducial components
                #
                data.Draw("mmiss>>hout(50,0,2500)",
                          "ppsEff*wgt*(!(%s) && protonCat==%d && mmiss>0 && mixType==1)"%(FIDUCIALCUTS,icat),
                          'goff')
                h=ROOT.gDirectory.Get('hout')
                if not icat in hout:
                    hout[icat]=h.Clone('hout%d'%icat)
                    hout[icat].SetDirectory(0)
                    hout[icat].Reset('ICE')
                hout[icat].Add(h,xsecWgt)

                data.Draw("mmiss>>hin(50,0,2500)",
                          "ppsEff*wgt*(%s && protonCat==%d && mmiss>0 && mixType==1)"%(FIDUCIALCUTS,icat),
                          "goff")
                h=ROOT.gDirectory.Get('hin')
                if not icat in hin:
                    hin[icat]=h.Clone('hin%d'%icat)
                    hin[icat].SetDirectory(0)
                    hin[icat].Reset('ICE')
                hin[icat].Add(h,xsecWgt)

    #final plots
    for icat in hin:

        pcatTitle='single(+z)-single(-z)'
        if not isSingleRP:
            pcatTitle='multi(+z)-multi(-z)'
            if icat==2: pcatTitle='multi(+z)-single(-z)'
            if icat==3: pcatTitle='single(+z)-multi(-z)'
            if icat==4: pcatTitle='single(+z)-single(-z)'
        finalExtraText='m_{{X}}={0} GeV\\{1}'.format(mass,pcatTitle)

        p=Plot('mmass_{0}_sigacc_{1}{2}'.format(mass,icat,pfix),com='13 TeV')
        #p.cmsLabel='#bf{CMS-Totem} #it{Preliminary}'
        p.cmsLabel='#bf{CMS-TOTEM}'
        p.plotformats=['pdf','png','C','root']
        p.saveLog=True
        p.xtit='Missing mass [GeV]'
        p.ytit='Events'

        hin[icat].SetLineWidth(2)
        p.add(hin[icat],  title='In fiducial',  color=ROOT.kRed,   isData=False, spImpose=False, isSyst=False)
        p.add(hout[icat], title='Out fiducial', color=ROOT.kBlack, isData=False, spImpose=True,  isSyst=False)
        p.show(outDir='./', lumi=LUMI, extraText=finalExtraText)
        
        p=Plot('mmass_{0}_puonfidsig_{1}{2}'.format(mass,icat,pfix),com='13 TeV')
        #p.cmsLabel='#bf{CMS-Totem} #it{Preliminary}'
        p.cmsLabel='#bf{CMS-TOTEM}'
        p.saveLog=True
        p.xtit='Missing mass [GeV]'
        p.ytit='Events'
        p.add(hnopu[icat],      title='No pileup',   color=ROOT.kGray, isData=False, spImpose=False, isSyst=False)
        p.add(hrpinhpur[icat],  title='Pileup+cuts', color=1,          isData=False, spImpose=True,  isSyst=False)
        p.show(outDir='./', lumi=LUMI,extraText=finalExtraText)
    

if __name__ == "__main__":
    sys.exit(main())

