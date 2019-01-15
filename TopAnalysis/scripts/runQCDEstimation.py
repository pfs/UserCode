import os
import sys
import optparse
import ROOT
import pickle
import json
from TopLJets2015.TopAnalysis.storeTools import *

"""
Get data and sum of MC from file
"""
def getTemplates(fIn,dist,tag,refName='Multijets'):
    data,refH,sumMC=None,None,None
    for key in fIn.Get(dist).GetListOfKeys():

        keyName=key.GetName()
        if 'Graph' in keyName : continue

        h=fIn.Get('%s/%s'%(dist,keyName))

        #reference
        if refName in keyName : 
            if refH:
                refH.Add(h)
            else:
                refH=h.Clone('%s_%s'%(refName,tag))
                refH.SetDirectory(0)
        #data
        elif keyName==dist:
            data=h.Clone('data_'+tag)
            data.SetDirectory(0)

        #other processes
        else:
            if sumMC:
                sumMC.Add(h)
            else:
                sumMC=h.Clone('mcsum_'+tag)
                sumMC.SetDirectory(0)

    return data,sumMC,refH


"""
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('--iso',          dest='iso',          help='plotter file with the iso selection',        default=None,       type='string')
    parser.add_option('--noniso',       dest='noniso',       help='plotter file with the non iso selection',    default=None,       type='string')
    parser.add_option('--out',          dest='outdir',       help='output directory',                           default='./',       type='string')
    parser.add_option('--sels',         dest='sels',         help='selections',                                 default='0j,1j,2j,3j,4j',       type='string')
    (opt, args) = parser.parse_args()
    #prepare output
    os.system('mkdir -p %s'%opt.outdir)

    #prepare fitter
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    #ROOT.gROOT.SetBatch(True)
    ROOT.FWLiteEnabler.enable()
    #ROOT.AutoLibraryLoader.enable()
    #ROOT.gSystem.Load('libTopLJets2015TopAnalysis.so')
    #ROOT.gROOT.LoadMacro('src/TemplatedFitTools.cc+')
    #from ROOT import TemplatedFitTools
    #fracFitter=ROOT.TemplatedFitTools()
        
    #open inputs
    fNonIso=ROOT.TFile.Open(opt.noniso)
    fIso=ROOT.TFile.Open(opt.iso)

    #perform a fit to a variable of interest
    nonIsoTemplateSF={}
    qcdNormUnc={}
    sels=opt.sels.split(',')
    for sel in sels: 

        #data in the sideband
        extraSel=''
        if sel in ['1j', '2j'] : extraSel='1t'
        if sel in ['3j','4j']  :
            extraSel='0t'
            if 'enoniso' in opt.noniso: extraSel='1t'

        postfix='_%s%s'%(sel,extraSel)
        if sel=='' : postfix=''

        dataNonIso,      sumMCNonIso, _  = getTemplates(fIn=fNonIso, dist='metpt%s'%postfix, tag='noniso')
        dataNonIsoAlt, sumMCNonIsoAlt, _ = getTemplates(fIn=fNonIso, dist='mt%s'%postfix,    tag='nonisoalt')

        #data in the signal region
        dataIso,    sumMCIso, _     = getTemplates(fIn=fIso,    dist='metpt%s'%postfix,  tag='iso')
        dataIsoAlt, sumMCIsoAlt, _  = getTemplates(fIn=fIso,    dist='mt%s'%postfix,     tag='isoalt')

        #normalized QCD template below the MT cut
        bin0           = 0
        binN           = dataNonIso.GetXaxis().FindBin(20.)
        niso           = dataIso.Integral(bin0,binN)
        nmciso         = sumMCIso.Integral(bin0,binN)
        nnoniso        = dataNonIso.Integral(bin0,binN)
        nmcnoniso      = sumMCNonIso.Integral(bin0,binN)

        #normalized QCD template above the Alt cut
        bin0         = 0
        binN         = dataNonIsoAlt.GetXaxis().FindBin(20.0)
        nisoAlt      = dataIsoAlt.Integral(bin0,binN)
        nmcisoAlt    = sumMCIsoAlt.Integral(bin0,binN)
        nnonisoAlt   = dataNonIsoAlt.Integral(bin0,binN)
        nmcnonisoAlt = sumMCNonIsoAlt.Integral(bin0,binN)

        #scale factors to apply and relative uncertainty (maximised)
        nonIsoSF=ROOT.TMath.Max(niso-nmciso,0.)/(nnoniso-nmcnoniso)
        nonIsoSFAlt=ROOT.TMath.Max(nisoAlt-nmcisoAlt,0.)/(nnonisoAlt-nmcnonisoAlt)
        unc= ROOT.TMath.Abs(nonIsoSFAlt/nonIsoSF-1) if nonIsoSF>0 else 1.0
        print sel,niso,nmciso,nnoniso,nmcnoniso

        nonIsoTemplateSF[sel]=(nonIsoSF,unc)

        #free mem
        dataNonIso.Delete()
        sumMCNonIso.Delete()
        dataIso.Delete()
        sumMCIso.Delete()

    fIso.Close()

    #produce the QCD templates
    fOut=ROOT.TFile.Open('%s/Data_QCDMultijets.root'%(opt.outdir),'RECREATE')
    for distKey in fNonIso.GetListOfKeys():
        dist=distKey.GetName()

        if 'iso' in dist or 'ratevsrun' in dist: continue
        category=dist.split('_')[-1][:2]    
        if dist.find('_')<0 : category=''
        dataNonIso, sumMCNonIso,_ = getTemplates(fIn=fNonIso, dist=dist, tag=dist)

        #do not subtract anything in the CR
        dataNonIsoUp=dataNonIso.Clone('%s_QCD%sUp'%(dist,category))
        dataNonIsoUp.SetTitle('QCD multijets - up')

        #subtract 2xMC in the CR
        dataNonIsoDown=dataNonIso.Clone('%s_QCD%sDown'%(dist,category))
        dataNonIsoDown.Add(sumMCNonIso,-2)       
        dataNonIsoDown.SetTitle('QCD multijets - down')

        #subtract 1xMC in the CR
        dataNonIso.Add(sumMCNonIso,-1)  
        dataNonIso.SetName(dist)
        dataNonIso.SetTitle('QCD multijets')

        try:
            dataNonIso.Scale( nonIsoTemplateSF[category][0] )
            totalQCD=dataNonIso.Integral()
            if totalQCD>0:
                for xbin in xrange(0,dataNonIsoUp.GetNbinsX()):
                    val=dataNonIsoUp.GetBinContent(xbin+1)
                    if val<0 : dataNonIsoUp.SetBinContent(xbin+1,1e-5)
                    val=dataNonIsoDown.GetBinContent(xbin+1)
                    if val<0 : dataNonIsoDown.SetBinContent(xbin+1,1e-5)
                dataNonIsoUp.Scale(totalQCD/dataNonIsoUp.Integral())
                dataNonIsoDown.Scale(totalQCD/dataNonIsoDown.Integral())
        except:            
            print 'unable to normalize ',dist,' with category=',category
        sumMCNonIso.Delete()

        #dump to file
        fOut.cd()
        for h in [dataNonIso,dataNonIsoUp,dataNonIsoDown]:
            h.SetDirectory(fOut)
            h.Write()

        
    #dump to file    
    print 'QCD scale factors CR->SR'
    print nonIsoTemplateSF
    cachefile=open('%s/.qcdscalefactors.pck'%opt.outdir,'w')
    pickle.dump(nonIsoTemplateSF, cachefile, pickle.HIGHEST_PROTOCOL)
    cachefile.close()

    #all done
    print 'QCD templates and systematic variations stored in %s'%fOut.GetName()
    fOut.Close()
    fNonIso.Close()

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
