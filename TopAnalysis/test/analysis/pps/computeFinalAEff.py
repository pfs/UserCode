import ROOT
import os
from runExclusiveAnalysis import getTracksPerRomanPot
from generateBinnedWorkspace import SIGNALXSECS,PHOTONSIGNALXSECS
from prepareOptimScanCards import OPTIMLIST
import numpy as np
import re
import pickle
import sys

import pandas as pd

_eraFracs= {'preTS2':14586.4464/41529.3,'postTS2':1.-14586.4464/41529.3}
_rgx=re.compile('.*/(.*)_m_X_(.*)_xangle_(.*)_2017_(.*).root')

def getFiducialAcceptance(args):

    """counts the number of events which are RP-fiducial"""

    idx,f=args

    #parse filename
    tags=_rgx.match(f).groups()[0:]
    boson=tags[0]
    gen_mX=float(tags[1])
    xangle=int(tags[2])
    era=tags[3]

    #open file and get tree
    fIn=ROOT.TFile.Open(f)
    tree=fIn.Get('analysis/data')
    nEntries=tree.GetEntriesFast()

    #pz weighting
    nSignalWgtSum,nFidSignalWgtSum=0.,0.
    pzwid=0.391*gen_mX+624
    for i in xrange(0,nEntries):
        tree.GetEntry(i)

        #require in-fiducial
        tracks=getTracksPerRomanPot(tree,era=era,run=-1,xangle=xangle,mcTruth=True,applyPxFid=False)
        if len(tracks[0][0])==0 or len(tracks[1][0])==0 : continue
        gencsi1=tracks[0][0][0]
        gencsi2=tracks[1][0][0]
        if gencsi1<0.02 or gencsi1>0.16 or gencsi2<0.03 or gencsi2>0.18: continue

        wpz=ROOT.TMath.Gaus(tree.gen_pzpp,0,pzwid)
        nSignalWgtSum += wpz

        if tree.evcat==22:
            if tree.bosonpt<95 : continue
        elif tree.evcat==121 or tree.cat==169:
            if tree.l1pt<30: continue
            if tree.l2pt<20: continue
            if tree.bosonpt<40 : continue # this should be the only cut removing some events
        else:
            continue

        nFidSignalWgtSum += wpz

    #all done with this file
    fIn.Close()

    #finalize computation
    fid_acc=nFidSignalWgtSum/nSignalWgtSum if nSignalWgtSum>0 else 0.
    xsec=(SIGNALXSECS[xangle] if boson=='Z' else PHOTONSIGNALXSECS[xangle])*_eraFracs[era]
    fid_xsec=xsec*fid_acc

    print('@ {} xsec={:3.3f}pb acc={:3.3f}'.format(idx,fid_xsec, fid_acc))
    with open('acc_job_%d.dat'%idx,'w') as fOut:
        fOut.write('%s %d %d %s %f %f %f'%(boson,gen_mX,xangle,era,xsec,fid_acc,fid_xsec))

def buildAcceptanceSummary():

    """build the acceptance summary looping over all the signals"""

    sig_dir='/eos/cms/store/cmst3/group/top/RunIIReReco/2017/vxsimulations_7Sep2020'
    task_list=[(i,os.path.join(sig_dir,f)) for i,f in enumerate(os.listdir(sig_dir))]
    import multiprocessing as MP
    pool = MP.Pool(20)
    pool.map( getFiducialAcceptance, task_list )

    summary='test/analysis/pps/acc_summary.dat'
    os.system("awk '{print}' acc_job_*.dat  > %s"%summary)
    os.system('rm acc_job_*.dat')

    return summary



if len(sys.argv)<2: summary=buildAcceptanceSummary()
else:               summary=sys.argv[1]

if len(sys.argv)<3:
    print 'python computeFinalAEff.py acc_summary.dat stat_analysis_info.dat'
    sys.exit(-1)


def getSelectedXsec(shapes):

    """computes the final analysis acceptance"""

    lumi={'zee':37193.,'zmm':37193.,'g':2288,}

    optim_rgx=re.compile('.*/optim_(.*)')

    summary=[]
    for s in shapes:    

        optim_idx=int(optim_rgx.match(s).group(1))
        optim_cut=OPTIMLIST[optim_idx]

        protonCat=0
        xangle=0
        nvtx=0
        if 'protonCat' in optim_cut:
            protonCat=int(re.match('.*protonCat==([0-9]*).*',optim_cut).group(1))
        if 'xangle' in optim_cut:
            xangle=int(re.match('.*xangle==([0-9]*).*',optim_cut).group(1))
        if 'nvtx<' in optim_cut:
            nvtx=-int(re.match('.*nvtx<([0-9]*).*',optim_cut).group(1))
        if 'nvtx>=' in optim_cut:
            nvtx=int(re.match('.*nvtx>=([0-9]*).*',optim_cut).group(1))
            
        for ch in [22,121,169]:
            for m in [600,660,720,780,800,840,900,960,1000,1020,1080,1140,1200,1260,1320,1380,1400,1440,1500,1560,1600,1620]:
                url=os.path.join(s,'shapes_{}.root'.format(ch))
                boson='g'
                if ch==121 : boson='zee'
                if ch==169 : boson='zmm'
                fIn=ROOT.TFile.Open(url)
                try:
                    nevts=fIn.Get('fidsig_{}_m{}'.format(boson,m)).Integral()
                except:
                    nevts=-1
                fIn.Close()

                summary.append( [boson,m,xangle,nvtx,protonCat,nevts,nevts/lumi[boson]] )
        
    return pd.DataFrame( np.array(summary), columns=['boson','mX','xangle','nvtx','protonCat','n_sel','xsec_sel'])
                

info=sys.argv[2]

print('Reading acceptance summary from',summary)
summary=pd.read_csv(summary,sep=' ',header=None, names=['boson','mX','xangle','era','xsec','fid_acc','fid_xsec'])
print(summary.head())

print('Stat analysis shapes based on',info)
shapeDirs=[]
with open(info,'r') as fIn:
    for line in fIn:
        if not '.dat' in line: continue
        shapeDirs=line.split()
        shapeDirs=[os.path.dirname(info)+'/'+os.path.dirname(s) for s in shapeDirs]
sel_summary=getSelectedXsec(shapeDirs)
print(sel_summary.head())

outf='acceptance_summary.h5'
summary.to_hdf(outf,key='acceptance',mode='w')
sel_summary.to_hdf(outf,key='selection',mode='a')
print('See summary in',outf)


