import ROOT
import os
import sys
import argparse
import itertools
from generateBinnedWorkspace import defineProcessTemplates,smoothWithLowess
from TopLJets2015.TopAnalysis.Plot import *
import numpy as np
import pickle

def parseTitleFromCut(cut):
    title=[]
    for icut in cut.split('&&'):
        for var,token in [ ('p_{T}(V)',      'bosonpt'),
                           ('p_{T}(1)',      'l1pt'),
                           ('p_{T}(2)',      'l2pt'),
                           ('multi-multi',   'protonCat==1'),
                           ('multi-single',  'protonCat==2'),
                           ('single-multi',  'protonCat==3'),
                           ('single-single', 'protonCat==4'),
                           ('N(vtx)',        'nvtx'),
                           ('nch',           'N_{ch}(PV)'),
                           ('=',             '=='),
                           ('#geq',          '>='),
                           ('#leq',          '<='),
                           ('',              ' ')]:
            icut=icut.replace(token,var)
        if 'era' in icut:
            icut='2017%s'%str(unichr(int(icut.split('=')[1])))
        title.append(icut)
    return title

def rebinUnequalBinSize(h,newBins):
    nBin = len(newBins)
    newname = h.GetName()+'_new'
    hnew = h.Rebin(nBin-1, newname, np.array(newBins))
    return hnew


def fillShapes(inputDir,selCuts,proc='MuonEG'):

    """fills the shapes from a set of files"""

    #import signal events
    data=ROOT.TChain('data')
    for f in [os.path.join(inputDir,x) for x in os.listdir(inputDir) if proc in x]:
        data.AddFile(f)
        
    histos={}
    sf=None
    for dist,hist,title in [('mmiss',                    '(50,0,2000)',             'Missing mass [GeV]'),
                            #('(bosony>=0?mmiss:-mmiss)', '(100,-2000,2000)',        'Missing mass x sgn[y(e#mu)] [GeV]'),
                            ('csi1',                     '(25,0.025,0.18)',         '#xi(+)'),
                            ('csi2',                     '(25,0.025,0.18)',         '#xi(-)'),
                            ('mpp',                      '(25,500,2500)',           'Di-proton mass [GeV]'),
                            #('ypp',                      '(50,-1,1)',               'Di-proton rapidity'),
                            #('mmiss:mpp',                '(25,500,2500,30,0,2500)', 'Missing mass [GeV];Di-proton mass [GeV]'),
                            #('ypp:mpp',                  '(25,500,2500,50,-1,1)',   'Di-proton rapidity;Di-proton mass [GeV]'),                        
                            ]:

        histos[dist]={}

        bkgHistos=[]
        for tag,pf,mixType in [ ('data',               '',        0),   #observed data
                                ('bkg',                '',        1),   #mixed data
                                ('bkgshape',           'syst',    1),   #mixed data from e-mu
                                ('bkgsinglediffUp',    '',        2),   #mix only pos. arm
                                ('bkgsinglediffDown',  'syst',    2),   #mix only neg. arm
                                ]:
            finaldist=pf+dist
                
            #if dist=='mmiss' or dist=='(bosony>=0?mmiss:-mmiss)':
            #    from array import array
            #    newBins=range(0,2000,40)
            #    #newBins=[0,100]+range(200,1000,40)+[1000,1100,1200,1500,2000]
            #    if dist!='mmiss':
            #        newBins=range(-2000,2000,40)
            #        #nbins=len(newBins)
            #        #newBins=[-newBins[ix] for ix in range(nbins-1,0,-1)] + newBins
            #    h=ROOT.TH1F( 'hmiss',';Missing mass [GeV];Events',len(newBins)-1,array('d',newBins) )
            #    data.Draw('{0} >> hmiss'.format(finaldist),
            #              'wgt*({0} && {1}mmiss>0 && mixType=={2})'.format(selCuts,pf,mixType),
            #              'goff')
            #else:

            finalSelCuts=selCuts
            if pf=='syst' : 
                finalSelCuts=finalSelCuts.replace('protonCat','systprotonCat')

            print '{0} >> h{1}'.format(finaldist,hist)
            print 'wgt*({0} && {1}mmiss>0 && mixType=={2})'.format(finalSelCuts,pf,mixType)
            print '-'*50
            data.Draw('{0} >> h{1}'.format(finaldist,hist),
                      'wgt*({0} && {1}mmiss>0 && mixType=={2})'.format(finalSelCuts,pf,mixType),
                      'goff')
                
            h=data.GetHistogram()
            histos[dist][tag]=h.Clone('{0}_{1}_{2}_obs'.format(dist,tag,proc))
            if ';' in title:
                xtit,ytit=title.split(';')
                histos[dist][tag].GetZaxis().SetTitle('Events')
                histos[dist][tag].GetYaxis().SetTitle(ytit)
                histos[dist][tag].GetXaxis().SetTitle(xtit)
            else:
                histos[dist][tag].GetYaxis().SetTitle('Events')
                histos[dist][tag].GetXaxis().SetTitle(title)
            histos[dist][tag].SetDirectory(0)
            h.Reset('ICE')

            if tag=='data' : 
                continue

            #apply a smoothing procedure for 1D mmiss
            #if 'mmiss' in dist and not h.InheritsFrom('TH2'):
            #    rawH=histos[dist][tag]
            #    histos[dist][tag]=smoothMissingMass(rawH)
            #    histos[dist][tag].SetName(rawH.GetName())
            #    histos[dist][tag].SetDirectory(0)
            #    rawH.Delete()
            
            bkgHistos.append(histos[dist][tag])

        h.Delete()

        #scale background estimates
        for tag in histos[dist]:
            if not 'bkg' in tag: continue
            sf=histos[dist]['data'].Integral()/histos[dist][tag].Integral()            
            histos[dist][tag].Scale(sf)

        #mirror shapes
        finalTemplates=defineProcessTemplates(bkgHistos)
        for h in finalTemplates:
            histos[dist][h.GetName()]=h

    return histos



def main(args):

    parser = argparse.ArgumentParser(description='usage: %prog [options]')
    parser.add_argument('-i', '--input',
                        dest='input',   
                        default='/eos/cms/store/cmst3/user/psilva/ExclusiveAna/final/ab05162/analysis_0p04/',
                        help='input directory with the files [default: %default]')
    parser.add_argument('--selCuts',
                        dest='selCuts', 
                        default='bosonpt>40 && l1pt>30 && l2pt>20',
                        help='preselection for Z categories [default: %default]')
    parser.add_argument('--doPerNch',
                        dest='doPerNch', 
                        default=False,
                        action='store_true',
                        help='break-down per Nch [default: %default]')
    parser.add_argument('--doPerEra',
                        dest='doPerEra', 
                        default=False,
                        action='store_true',
                        help='break-down per era [default: %default]')
    parser.add_argument('--doPerPU',
                        dest='doPerPU', 
                        default=False,
                        action='store_true',
                        help='break-down per pileup category [default: %default]')
    parser.add_argument('--doPerAngle',
                        dest='doPerAngle', 
                        default=False,
                        action='store_true',
                        help='break-down per angle [default: %default]')
    parser.add_argument('--lumi',
                        dest='lumi',
                        default=37500.,
                        type=float,
                        help='integrated luminosity [default: %default]')
    parser.add_argument('-o', '--output',
                        dest='output', 
                        default='analysis/bkg',
                        help='Output directory [default: %default]')
    opt=parser.parse_args(args)

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    
    os.system('mkdir -p %s'%opt.output)

    #do all possible combinations of these categories
    baseCats = ['protonCat==%d'%(i+1) for i in range(4)]
    subCats  = ['']
    if opt.doPerEra:   subCats += ['era==%d'%int(ord(x)) for x in 'BCDEF']
    if opt.doPerPU:    subCats += ['nvtx<20','nvtx>=20']
    if opt.doPerAngle: subCats += ['xangle==%d'%i for i in [120,130,140,150]]
    if opt.doPerNch:   subCats += ['nch<15','nch>=15']
    catList = list(itertools.product(baseCats,subCats))
    with open('{}/categ_defs.pck'.format(opt.output),'wb') as cache:
        pickle.dump(catList,cache,pickle.HIGHEST_PROTOCOL)
    
    outF='%s/plotter_embkg.root'%opt.output
    os.system('rm %s'%outF)

    print '[doBackgroundValidation] with %d variations to test'%len(catList)
    print '\t output will be available in',outF
    for i,cat in enumerate(catList):
        
        selCut=''
        for c in cat:
            if len(c)==0 : continue
            selCut += '%s&&'%c
        selCut += opt.selCuts
        titleCut=parseTitleFromCut(selCut)
        print '\t',i,selCut

        try:
            data=fillShapes(inputDir=opt.input,selCuts=selCut,proc='MuonEG')
        except Exception as e:
            print '-'*50
            print 'Unable to fill shapes for',i,cat
            print e
            print '-'*50
            continue

        for dist in data:

            pname='%s_%d'%(dist,i)
            for c in [':',',','>','=','(',')','-','<','?']: pname=pname.replace(c,'')
            p=Plot(pname)
            p.doChi2=False #True
            p.nominalDistForSystsName='background'
            
            #if dist=='mmiss':
            #    newBins=range(0,1000,40)+[1000,1100,1200,1500,2000]
            #    print 'Rebinning',dist,'to',newBins
            #    for k in data[dist]:
            #        data[dist][k]=rebinUnequalBinSize(data[dist][k],newBins)

            p.add(data[dist]['data'], title='data',       color=ROOT.kBlack,  isData=True,  spImpose=False, isSyst=False)
            p.add(data[dist]['bkg'],  title='background', color=ROOT.kCyan-6, isData=False, spImpose=False, isSyst=False)

            #background systematics
            ci=1
            for syst in ['{0}_bkgshape_MuonEG_obsUp',        
                         '{0}_bkgshape_MuonEG_obsDown',
                         '{0}_bkgsinglediffUp_MuonEG_obs', 
                         '{0}_bkgsinglediffDown_MuonEG_obs',
            ]:
                ci=ci+1
                p.add(data[dist][syst.format(dist)],  
                      title=syst, 
                      color=ROOT.kCyan-6, #ci
                      isData=False, 
                      spImpose=False,  #True,
                      isSyst=True)     #False)

            #p.ratiorange=[0.78,1.22]
            p.ratiorange=[0.58,1.43]
            p.show(opt.output,opt.lumi,extraText='\\'.join(titleCut))
            p.appendTo(outF)
            p.reset()


if __name__ == "__main__":

    sys.exit(main(sys.argv[1:]))



