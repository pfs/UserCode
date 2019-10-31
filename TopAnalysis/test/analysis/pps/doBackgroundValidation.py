import ROOT
import os
import sys
import argparse
import pickle
import copy
from generateBinnedWorkspace import defineProcessTemplates,VALIDLHCXANGLES
from TopLJets2015.TopAnalysis.Plot import *
        
def fillShapes(inputDir,selCuts,proc='MuonEG'):

    """fills the shapes from a set of files"""

    #import signal events
    data=ROOT.TChain('data')
    for f in [os.path.join(inputDir,x) for x in os.listdir(inputDir) if proc in x]:
        data.AddFile(f)

    histos={}
    sf=None
    for dist,hist,title in [('mmiss',   '(30,0,2500)',    'Missing mass [GeV]'),
                            ('csi1',    '(25,0.04,0.18)', '#xi'),
                            ('csi2',    '(25,0.04,0.18)', '#xi'),                            
                            ('mpp',     '(25,500,2500)',  'Di-proton mass [GeV]'),
                            ('ypp',     '(50,-1,1)',      'Di-proton rapidity')
                            ]:

        print dist
        histos[dist]={}

        bkgHistos=[]
        for tag,pf,mixType in [ ('data',          '',    0),
                                ('bkg',           '',    1),
                                ('bkgsinglediff', 'syst',1),
                                ('bkgshape',      '',    2), 
                                ]:
            data.Draw('{0}{1} >> h{2}'.format(pf,dist,hist),
                      'wgt*({0} && {1}mmiss>0 && mixType=={2})'.format(selCuts,pf,mixType),
                      'goff')
            h=data.GetHistogram()
            histos[dist][tag]=h.Clone('{0}_{1}_{2}_obs'.format(dist,tag,proc))
            histos[dist][tag].GetYaxis().SetTitle('Events')
            histos[dist][tag].GetXaxis().SetTitle(title)
            histos[dist][tag].SetDirectory(0)
            h.Reset('ICE')
            if tag=='data' : continue
            bkgHistos.append(histos[dist][tag])

        h.Delete()

        #use low missing mass to fix possible normalization differences
        if sf is None:
            upBin=histos[dist]['data'].GetXaxis().FindBin(500.)
            sf=histos[dist]['data'].Integral(0,upBin)/histos[dist]['bkg'].Integral(0,upBin)
            print 'Scale factor:',sf

        #scale background estimates
        for tag in histos[dist]:
            if not 'bkg' in tag: continue
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
    parser.add_argument('--sig',
                        dest='sig',
                        default='{boson}_m_X_{mass}_xangle_{xangle}_2017_preTS2_opt_v1_simu_reco.root',
                        help='signal point [%default]')
    parser.add_argument('--massList',
                        dest='massList',
                        default='800,1000,1200',
                        help='signal mass list (CSV) [%default]')
    parser.add_argument('--selCuts',
                        dest='selCuts', 
                        default='bosonpt>40 && l1pt>30 && l2pt>20',
                        help='preselection for Z categories [default: %default]')
    parser.add_argument('--doPerAngle',
                        dest='doPerAngle', 
                        default=False,
                        help='do per crossing angle [default: %default]',
                        action='store_true')
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

    selCuts=[ (opt.selCuts,'','e#mu') ]
    if opt.doPerAngle:
        for xangle in VALIDLHCXANGLES:
            selCuts.append( (opt.selCuts + ' && xangle==%d'%xangle,'_%d'%xangle,'e#mu, %d #murad'%xangle) )
    for cuts,pfix,catTitle in selCuts:

        data=fillShapes(inputDir=opt.input,selCuts=cuts,proc='MuonEG')
        #sigs={}
        #for m in opt.massList.split():
        #    if not m in sigs: sigs[m]={}
        #    for xangle in VALIDLHCXANGLES:
        #        newSigs=fillShapes(inputDir=opt.input,selCuts=opt.selCuts,tag=opt.sig.format(mass=m,xangle=xangle))
        #        for dist in newSigs:
        #            if not dist in sigs[m]:
        #                sigs[m]=newSigs[dist]['data'].Clone('{0}_{1}'.format(dist,m))
        #            #FIXME scale me according to xsec
        #FIXME plot me

        for dist in data:
            p=Plot(dist+pfix)
            p.doChi2=True
            p.nominalDistForSystsName='background'

            #main distributions
            p.add(data[dist]['data'], title='data',       color=ROOT.kBlack, isData=True, spImpose=False, isSyst=False)
            p.add(data[dist]['bkg'],  title='background', color=ROOT.TColor.GetColor('#1f78b4'), isData=False, spImpose=False, isSyst=False)

            #background systematics
            for syst in ['{0}_bkgshape_MuonEG_obsUp',        
                         '{0}_bkgshape_MuonEG_obsDown',
                         '{0}_bkgsinglediff_MuonEG_obsDown', 
                         '{0}_bkgsinglediff_MuonEG_obsUp']:
                p.add(data[dist][syst.format(dist)],  
                      title=syst, 
                      color=ROOT.TColor.GetColor('#1f78b4'), 
                      isData=False, 
                      spImpose=False, 
                      isSyst=True)       

            p.ratiorange=[0.88,1.12]
            p.show(opt.output,opt.lumi,extraText=catTitle)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))



