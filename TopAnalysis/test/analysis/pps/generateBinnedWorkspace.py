import ROOT
import os
import sys
import argparse
import pickle
from TopLJets2015.TopAnalysis.roofitTools import showFitResult,shushRooFit

#sigma=1pb distributed accross crossing angles 
#NB this does not sum to 1 as we don't use all crossing angles in the analysis
SIGNALXSECS={120:0.269,130:0.273,140:0.143,150:0.293}
VALIDLHCXANGLES=[120,130,140,150]

def defineCsiAcceptanceAndBinning(url,mass):

    """reads the pickle file with the linear parameterization of the csi acceptances from simulation
    computes the acceptance at a given mass for each crossing angle"""

    if not url : return {}
    with open(url,'r') as f:
        csiaccParam=pickle.load(f)
        resolParam=pickle.load(f)

    csiacc={}
    for xangle in VALIDLHCXANGLES:
        csiacc[xangle]=[]
        for rp in [23,123]:
            a,b=csiaccParam[(xangle,rp)]
            csiacc[xangle].append( (a+b*mass/1000.,0.18) )
    return csiacc,resolParam


def defineProcessTemplates(histos):

    """defines the nominal template and the variations and checks fo 0's in the histograms"""

    templates=[]

    #nominal
    templates.append( histos[0] )
    nomStats=templates[-1].Integral()

    #systematic variations
    #if Up/Down already in the name store directly updating the name
    #if not, mirror the variation given 
    for i in xrange(1,len(histos)):        
        templates.append( histos[i] )
        key=templates[-1].GetName()
        if not 'Up' in key and not 'Down' in key :
            templates[-1].SetName(key+'Up')
            templates.append( histos[i].Clone(key+'Down') )
            for xbin in range(templates[0].GetNbinsX()):
                templates[-1].SetBinContent(xbin+1,2*templates[0].GetBinContent(xbin+1)-templates[-2].GetBinContent(xbin+1))
    
    #don't leave bins with 0's
    for h in templates:
        h.SetDirectory(0)
        iStats=h.Integral()
        if iStats>0: h.Scale(nomStats/iStats)
        for xbin in range(h.GetNbinsX()):
            if h.GetBinContent(xbin+1)>0: continue
            h.SetBinContent(xbin+1,1e-6)
            
    return templates

        
def fillBackgroundTemplates(opt):

    """fills the background and observed data histograms"""

    totalBkg={}
    templates=[]

    #import signal events
    data=ROOT.TChain('data')
    for f in [os.path.join(opt.input,x) for x in os.listdir(opt.input) if 'Data13TeV' in x]:
        if 'MuonEG' in f : continue
        data.AddFile(f)

    #define final preselection cuts
    cuts='xangle==%d'%opt.xangle
    if len(opt.presel) : cuts += ' && ' + opt.presel    
    if opt.csiacc:
        csiCuts ='csi1>%f && csi1<%f && '%opt.csiacc[opt.xangle][0]
        csiCuts+='csi2>%f && csi2<%f'%opt.csiacc[opt.xangle][1]
        cuts=csiCuts if len(cuts)==0 else '{0} && {1}'.format(cuts,csiCuts)

    #loop over categories build templates
    for icat in range(len(opt.categs)):

        #apply category cuts
        categCut=opt.categs[icat]
        categCut=cuts if len(categCut)==0 else '%s && %s'%(categCut,cuts)
        catName='%s_a%d_%d'%(opt.chTag,opt.xangle,icat)

        print '\t',catName,categCut

        #background modelling histos
        histos=[]
        data_obs=None
        for name,pfix in [('bkg_'+catName,'mix'),('bkg_%s_bkgShape'%catName,'mixem')]:

            templCuts=categCut.replace('csi1',pfix+'csi1')
            templCuts=templCuts.replace('csi2',pfix+'csi2')
            data.Draw('{0}mmiss >> h({1},{2},{3})'.format(pfix,opt.nbins,opt.mMin,opt.mMax),templCuts,'goff')
            h=data.GetHistogram()
            histos.append(h.Clone(name))
            histos[-1].SetDirectory(0)

            if len(histos)==1:
                totalBkg[icat]=h.Integral()
                if not opt.unblind :
                    data_obs=h.Clone('data_obs_'+catName)
                    data_obs.SetDirectory(0)

            h.Reset('ICE')
        templates += defineProcessTemplates(histos)

        #observed data in this category if unblinding
        if opt.unblind:
            data.Draw('mmiss >> h({1},{2},{3})'.format(opt.nbins,opt.mMin,opt.mMax),categCut,'goff')
            data_obs=data.GetHistogram().Clone('data_obs_'+catName)
            data_obs.SetDirectory(0)

        templates.append(data_obs)

    print '\t total background:',totalBkg
    return totalBkg,templates


def fillSignalTemplates(opt):

    """fills the signal histograms"""

    totalSig={}
    templates=[]

    #import signal events
    data=ROOT.TChain('data')
    data.AddFile(os.path.join(opt.input,opt.sig))

    #define final preselection cuts
    cuts='xangle==%d'%opt.xangle
    if len(opt.presel) : cuts += ' && ' + opt.presel
    if opt.csiacc:
        csiCuts ='csi1>%f && csi1<%f && '%opt.csiacc[opt.xangle][0]
        csiCuts+='csi2>%f && csi2<%f'%opt.csiacc[opt.xangle][1]
        cuts=csiCuts if len(cuts)==0 else '{0} && {1}'.format(cuts,csiCuts)

    #loop over categories build templates
    for icat in range(len(opt.categs)):

        #apply category cuts
        categCut=opt.categs[icat]
        categCut=cuts if len(categCut)==0 else '%s && %s'%(categCut,cuts)

        catName='%s_a%d_%d'%(opt.chTag,opt.xangle,icat)
        print '\t',catName,categCut

        #signal modelling histograms
        histos=[]
        for name,pfix in [('sig_'+catName,''),('sig_%s_sigShape'%catName,'mix')]:

            templCuts=categCut.replace('csi1',pfix+'csi1')
            templCuts=templCuts.replace('csi2',pfix+'csi2')
            wgtExpr='wgt*%f'%(SIGNALXSECS[opt.xangle]*opt.lumi)
            data.Draw('{0}mmiss >> h({1},{2},{3})'.format(pfix,opt.nbins,opt.mMin,opt.mMax),
                      '{0}*({1})'.format(wgtExpr,templCuts),
                      'goff')
            h=data.GetHistogram()
            histos.append( h.Clone(name) )         
            histos[-1].SetDirectory(0)

            if len(histos)==1:
                totalSig[icat]=h.Integral()

            h.Reset('ICE')
        templates += defineProcessTemplates(histos)
    
    print '\t total signal:',totalSig
    return totalSig,templates


def writeDataCards(opt,sigExp,bkgExp,shapesURL):

    """writes the datacard and the workspace"""

    #create a card per category
    dcList=[]
    for icat in range(len(opt.categs)):
        cat='%s_a%d_%d'%(opt.chTag,opt.xangle,icat)
        dcTxt='%s/shapes-parametric.datacard_%s.dat'%(opt.output,cat)
        dcList.append(dcTxt)
        with open(dcTxt,'w') as dc:
            dc.write('#\n')
            dc.write('# datacard was automatically generated with generateWorkspace.py\n')
            dc.write('# the options passed are printed below\n')
            dc.write('# %s\n'%opt)
            dc.write('#\n')
            dc.write('imax *\n')
            dc.write('jmax *\n')
            dc.write('kmax *\n')
            dc.write('-'*50+'\n')
            dc.write('shapes *        * {0} $PROCESS_{1} $PROCESS_$SYSTEMATIC\n'.format(shapesURL,cat))
            dc.write('shapes data_obs * {0} $PROCESS_{1}\n'.format(shapesURL,cat))
            dc.write('-'*50+'\n')
            dc.write('bin %s\n'%cat)
            dc.write('observation -1\n')
            dc.write('-'*50+'\n')
            dc.write('%15s %15s %15s\n'%('bin',cat,cat))
            dc.write('%15s %15s %15s\n'%('process','sig','bkg'))
            dc.write('%15s %15s %15s\n'%('process','0', '1'))
            dc.write('%15s %15s %15s\n'%('rate','%3.2f'%sigExp[icat], '%3.2f'%bkgExp[icat]))
            dc.write('-'*50+'\n')
            
            #float the background normalization as well as the signal
            dc.write('mu_bkg{0} rateParam {0} bkg 1\n'.format(cat))

            #uncertainties
            dc.write('lumi %8s %15s %15s\n'%('lnN','1.027','-'))
            dc.write('%s_sigShape %8s %15s %15s\n'%(cat,'shape','1','-'))
            dc.write('%s_bkgShape %8s %15s %15s\n'%(cat,'shape','-','1'))
            dc.write('{0} autoMCStats 0.0 1\n'.format(cat))
        
    print '\tshapes available @',shapesURL
    print '\tgenerated the following datacards',dcList


def main(args):

    parser = argparse.ArgumentParser(description='usage: %prog [options]')
    parser.add_argument('-i', '--input',
                        dest='input',   
                        default='/eos/cms/store/cmst3/user/psilva/ExclusiveAna/ab05162/analysis/',
                        help='input directory with the files [default: %default]')
    parser.add_argument('--xangle',
                        dest='xangle',
                        default=150,
                        type=int,
                        help='signal point [%default]')
    parser.add_argument('--sig',
                        dest='sig',
                        default='MC13TeV_ppxz_m800_x150.root',
                        help='signal point [%default]')
    parser.add_argument('--presel',
                        dest='presel', 
                        default='cat==169 && l1pt>30 && l2pt>20 && bosonpt>50',
                        help='preselection [default: %default]')
    parser.add_argument('--csiacc',
                        dest='csiacc', 
                        default='test/analysis/pps/signal_resol_acc.pck',
                        help='parametrization of the csi acceptance cuts and mass resolution [default: %default]')
    parser.add_argument('--categs',
                        dest='categs',
                        default='nvtx<20,nvtx>=20',
                        help='Sub-categories [default: %default]')
    parser.add_argument('--mix',
                        dest='mix',
                        default='mix',
                        help='Mixing values to use [default: %default]')
    parser.add_argument('--lumi',
                        dest='lumi',
                        default=37500.,
                        type=float,
                        help='integrated luminosity [default: %default]')
    parser.add_argument('--mMin',
                        dest='mMin',
                        default=0,
                        type=float,
                        help='minimum missing mass [default: %default]')
    parser.add_argument('--mMax',
                        dest='mMax',
                        default=2500,
                        type=float,
                        help='maximum missing mass [default: %default]')
    parser.add_argument('-o', '--output',
                        dest='output', 
                        default='analysis/stat',
                        help='Output directory [default: %default]')
    parser.add_argument('--unblind',
                        dest='unblind', 
                        default=False,
                        action='store_true',
                        help='Use non-mixed data in the final fit [default: %default]')
    opt=parser.parse_args(args)

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    shushRooFit()

    #parse channels from pre-selection string
    import re
    ch_dict={'169':'zmm','121':'zee','22':'g'}
    chtit_dict={'169':'Z#rightarrow#mu#mu','121':'Z#rightarrowee','22':'#gamma'}
    regex=re.compile('cat==(\d+)') 
    allCh=regex.findall(opt.presel)
    setattr(opt,'chTag','_'.join([ch_dict[x] for x in allCh]))
    setattr(opt,'chTitle',','.join([chtit_dict[x] for x in allCh]))

    #configuration
    if not opt.xangle in VALIDLHCXANGLES : 
        print 'Crossing angle',opt.xangle,'is not valid a valid one'
        return
    opt.categs=opt.categs.split(',')
    if len(opt.categs)==0 : opt.categs=[]
    midx=3 if 'gamma' in opt.sig else 2
    opt.mass=int(opt.sig.split('_')[midx].replace('m','')) #hardcoded...

    if opt.csiacc: 
        opt.csiacc,resolParam=defineCsiAcceptanceAndBinning(opt.csiacc,opt.mass)
        setattr(opt,'massResol',resolParam[opt.xangle].Eval(800)) #opt.mass))
        setattr(opt,'nbins', ROOT.TMath.FloorNint( (opt.mMax-opt.mMin)/(opt.massResol)) )

    print '[generateWorkspace]'
    print '\t signal from',opt.sig,'mass=',opt.mass
    print '\t will apply the following preselection:',opt.presel
    print '\t channel tag=',opt.chTag
    if len(opt.categs) : 
        print '\t will categorize in:',opt.categs
    if opt.csiacc:
        print '\t will apply csi acceptance cuts',opt.csiacc
        print '\t mass resolution=',opt.massResol,' GeV =>',opt.nbins,' bins for templates',
        print 'in the range (',opt.mMin,',',opt.mMax,')'

    #start the output
    os.system('mkdir -p %s'%opt.output)
    shapesURL=os.path.join(opt.output,'shapes_%d.root'%opt.xangle)
    fOut=ROOT.TFile.Open(shapesURL,'RECREATE')

    #define background templates
    print '\t filling background templates and observed data'
    bkgExp,bkgTemplates=fillBackgroundTemplates(opt)    
    for h in bkgTemplates:        
        h.SetDirectory(fOut)
        h.Write()

    #parametrize the signal
    print '\t filling signal templates'
    sigExp,sigTemplates=fillSignalTemplates(opt)
    for h in sigTemplates:
        h.SetDirectory(fOut)
        h.Write()

    #all done
    fOut.Close()


    #write summary in datacards
    print '\t writing datacard'
    writeDataCards(opt,sigExp,bkgExp,os.path.basename(shapesURL))

    print '\t all done, output can be found in',opt.output

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))



