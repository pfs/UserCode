#!/usr/bin/env python

import ROOT
import copy
import optparse
import json
import sys
import os
import numpy
from array import array
import json
import random
import pickle
import re
from collections import OrderedDict,defaultdict
from TopLJets2015.TopAnalysis.HistoTool import *
from EventMixingTool import *
from EventSummary import EventSummary
from MixedEventSummary import MixedEventSummary
from PPSEfficiencyReader import PPSEfficiencyReader,isPixelFiducial,doFinalCheck2017
from TopLJets2015.TopAnalysis.myProgressBar import *

VALIDLHCXANGLES=[120,130,140,150]
DIMUONS=13*13
EMU=11*13
DIELECTRONS=11*11
SINGLEPHOTON=22
MINCSI=0.035
MIXEDRPSIG=None # this is only for a test
ALLOWPIXMULT=[1,2]
USESINGLERP=False

def isValidRunLumi(run,lumi,runLumiList):

    """checks if run is available and lumi section was certified"""

    #no run/lumi to select, all is good by default
    if not runLumiList:
        return True

    #run is available
    if not run in runLumiList: 
        return True

    for lran in runLumiList[run]:
        if lumi>=lran[0] and lumi<=lran[1]:
            return False

    #reached this far, nothing found in list
    return True

def computeCosThetaStar(lm,lp):
    dil=lm+lp
    costhetaCS = 1./dil.M()
    if dil.Pz()<0 : costhetaCS *= -1
    costhetaCS *= (lm.E() + lm.Pz()) * (lp.E() - lp.Pz()) - (lm.E() - lm.Pz()) * (lp.E() + lp.Pz());
    costhetaCS /= ROOT.TMath.Sqrt( (dil.M())**2 + (dil.Pt())**2)
    return costhetaCS


def getTracksPerRomanPot(tree,era,run,xangle,mcTruth=False,minCsi=0,orderByDecreasingCsi=True,useXY=False,applyPxFid=True):

    """
    loops over the availabe tracks in the event and groups them by roman pot id
    returns two list of protons [multi,Far(pixel), Near(strips)] for the positive and negativesides
    """

    tkPos=[[],[],[]]
    tkNeg=[[],[],[]]
    for itk in xrange(0,tree.nProtons):

        #read the reconstructed csi (or XY)
        try:

            if mcTruth:
                csi=tree.genProtonCsi[itk]
            else:
                csi=tree.protonCsi[itk]

            if useXY:
                csi=(tree.protonX[itk],tree.protonY[itk])

        except:
            csi=-99
            
            if useXY:
                csi=(-99,-99)

        if not useXY and csi<minCsi : 
            continue

        #check how the track was reconstructed
        isFar   = tree.isFarRPProton[itk]
        isMulti = tree.isMultiRPProton[itk]
        isPosRP = tree.isPosRPProton[itk]

        #fiducial cut
        if applyPxFid and (isMulti or isFar):

            sector=45 if isPosRP else 56
            passPxFid=isPixelFiducial(era,sector,
                                      tree.protonX[itk],tree.protonTX[itk],tree.protonY[itk],tree.protonTY[itk],
                                      csi if not useXY else None,xangle)
            if not passPxFid : continue

        idx = (0 if isMulti else (1 if isFar else 2))
        if mcTruth: idx=0

        if isPosRP : tkPos[idx].append(csi)
        else :  tkNeg[idx].append(csi)
    

    if orderByDecreasingCsi:   

        if useXY:
            for idx in range(3):
                tkPos[idx].sort(key=lambda x: x[0], reverse=True)
                tkNeg[idx].sort(key=lambda x: x[0], reverse=True)
        else:
            for idx in range(3):
                tkPos[idx].sort(reverse=True)
                tkNeg[idx].sort(reverse=True)


        
    #reset multi and strips
    if USESINGLERP and not mcTruth:
        tkPos[0]=[]
        tkPos[2]=[]
        tkNeg[0]=[]
        tkNeg[2]=[]

    #return all tracks
    return tkPos,tkNeg


def buildDiProton(csi_pos,csi_neg,sqrts=13000.):

    """build a diproton system from to tracks in the arms of the roman pots"""

    beamP=0.5*sqrts
    diproton=ROOT.TLorentzVector(0.,0.,beamP*(csi_pos-csi_neg),beamP*(csi_pos+csi_neg))

    return diproton


def buildMissingMassSystem(pp,boson,vis=None,sqrts=13000.):

    """builds the system corresponding to the missing mass"""

    if pp is None or boson is None:
        return None

    inPP        = ROOT.TLorentzVector(0,0,0,sqrts)
    mmassSystem = (pp-boson)
    if vis:
        mmassSystem=mmassSystem-vis

    return mmassSystem

def getDiProtonCategory(pos_protons,neg_protons,boson,allowPixMult):
            
    """combines the protons reconstructed in a diproton system and returns the event category"""

    proton_cat      = -1
    csi_pos,csi_neg = None,None
    pp              = None
    mmass           = None

    #check which csi's have been reconstructed
    nmulti_pos, nmulti_neg = len(pos_protons[0]),len(neg_protons[0])
    npix_pos,   npix_neg   = len(pos_protons[1]),len(neg_protons[1])
    if nmulti_pos==1 and nmulti_neg==1:
        proton_cat=1
        csi_pos,csi_neg=pos_protons[0][0],neg_protons[0][0]
    elif nmulti_pos==1 and (nmulti_neg==0 and npix_neg in allowPixMult):    
        proton_cat=2
        csi_pos,csi_neg=pos_protons[0][0],neg_protons[1][0]
    elif (nmulti_pos==0 and npix_pos in allowPixMult) and nmulti_neg==1:
        proton_cat=3
        csi_pos,csi_neg=pos_protons[1][0],neg_protons[0][0]
    elif (nmulti_pos==0 and npix_pos in allowPixMult) and (nmulti_neg==0 and npix_neg in allowPixMult):
        proton_cat=4
        csi_pos,csi_neg=pos_protons[1][0],neg_protons[1][0]

    #build diproton and missing mass
    if proton_cat>0:
        pp=buildDiProton(csi_pos,csi_neg)
        mmass=buildMissingMassSystem(pp,boson)

    return proton_cat,csi_pos,csi_neg,pp,mmass


def getRandomEra(isSignal,isPreTS2Signal):

    """generates a random era according to the integrated luminosity in each one"""

    r=random.random()

    cum_fracs=[('2017B',0.115),('2017C',0.348),('2017D',0.451),('2017E',0.671),('2017F',1)]
    if isSignal:
        if isPreTS2Signal:
            cum_fracs=[('2017B',0.327),('2017C',0.991),('2017D',1)]
        else:
            cum_fracs=[('2017D',0.154),('2017E',0.493),('2017F',1)]

    for era,cf in cum_fracs:
        if r>cf : continue
        return era

    return era

def isSignalFile(inFile):

    isSignal=False
    validTags=['gamma_m_X_','Z_m_X_',
               'Zee_m_X_','Zmm_m_X_',
               'SDZee','SDZmm']
    for t in validTags:
        if t in inFile:
            isSignal=True
            break
    isPreTS2Signal=True if 'preTS2' in inFile else False
    return isSignal,isPreTS2Signal

def isDYFile(inFile):
    isDY=True if 'DY50toInf' in inFile else False
    return isDY

def isPhotonSignalFile(inFile):
    isSignal=True if 'gamma_m_X_' in inFile else False
    return isSignal

def signalMassPoint(inFile):
    if 'm_X_' in inFile : 
        return float(re.search('m_X_(\d+)', inFile).group(1))
    return 0.

def isSignalFiducial(csiPos,csiNeg,gen_pzpp):

    if csiNeg<0.03 or csiPos<0.02 : return False
    if csiNeg>0.16 or csiPos>0.13 : return False
    return True

    #if gen_pzpp<-300: return False
    #if gen_pzpp>500 : return False
    #return True




def runExclusiveAnalysis(inFile,outFileName,runLumiList,effDir,ppsEffFile,maxEvents=-1,sighyp=0,mixDir=None):
    
    """event loop"""

    global MIXEDRPSIG
    global ALLOWPIXMULT

    isData=True if 'Data' in inFile else False
    era=os.path.basename(inFile).split('_')[1] if isData else None
    isDY=isDYFile(inFile)
    isSignal,isPreTS2Signal=isSignalFile(inFile)
    isFullSimSignal = False
    isSDZSim = True if 'SDZ' in inFile else False
    if isSignal:
        if 'fullsim' in inFile or isSDZSim:
            isFullSimSignal=True
    isPhotonSignal=isPhotonSignalFile(inFile)
    gen_mX=signalMassPoint(inFile) if isSignal else 0.


    #open this just once as it may be quite heavy in case it's not data or signal
    MIXEDRP=None
    if mixDir:
 
        print 'Collecting events from the mixing bank'
        mixFiles=[f for f in os.listdir(mixDir) if '.pck' in f]
        
        #open just the necessary for signal and data
        if isSignal or isData:
            allowedEras=['2017%s'%x for x in 'BCDEF']
            if isSignal: 
                try:
                    allowedAngle=re.search('xangle_(\d+)',inFile).group(1)
                except:
                    allowedAngle=re.search('_(\d+)_',inFile).group(1)
                mixFiles=[f for f in mixFiles if allowedAngle in f]
                if isPreTS2Signal : 
                    allowedEras=['2017%s'%x for x in 'BCD']
                else : 
                    allowedEras=['2017%s'%x for x in 'DEF']
            if isData:
                allowedEras=[era]
            mixFiles=[f for f in mixFiles if f.split('_')[1] in allowedEras]

        MIXEDRP=defaultdict(list)
        for f in mixFiles:
            print '\t',f
            with open(os.path.join(mixDir,f),'r') as cachefile:
                rpData=pickle.load(cachefile)
                for key in rpData:
                    MIXEDRP[key] += rpData[key]
        print '\t size of mixing bank is',sys.getsizeof(MIXEDRP),'byte'




    #bind main tree with pileup discrimination tree, if failed return
    tree=ROOT.TChain('analysis/data' if isSignal and not isFullSimSignal else 'tree')
    tree.AddFile(inFile)
    #try:
    #    pudiscr_tree=ROOT.TChain('pudiscr')
    #    baseName=os.path.basename(inFile)
    #    baseDir=os.path.dirname(inFile)
    #    pudiscr_file=os.path.join(baseDir,'pudiscr',baseName)
    #    if not os.path.isfile(pudiscr_file):
    #        raise ValueError(pudiscr_file+' does not exist')
    #    pudiscr_tree.AddFile(pudiscr_file)
    #    tree.AddFriend(pudiscr_tree)
    #    print 'Added pu tree for',inFile
    #except:
    #    #print 'Failed to add pu discrimination tree as friend for',inFile
    #    return

    
    #check if it is signal and load     
    signalPt=[]
    ppsEffReader=None
    mcEff={}
    if isSignal: 
        ppsEffReader=PPSEfficiencyReader(ppsEffFile)
        signalPt=[float(x) for x in re.findall(r'\d+', os.path.basename(inFile) )[2:]]
        for ch in ['eez','mmz','a']:
            effIn=ROOT.TFile.Open('%s/effsummary_%s_ptboson.root'%(effDir,ch))
            pname='gen%srec_ptboson_ZH#rightarrowllbb_eff'%ch
            if ch=='a': pname='genarec_ptboson_#gamma+jets_eff'
            mcEff[ch]=effIn.Get(pname)
        effIn.Close()        

    #start event mixing tool
    if MIXEDRP : print MIXEDRP.keys()
    evMixTool=EventMixingTool(mixedRP=MIXEDRP,validAngles=VALIDLHCXANGLES,usePixelOnly=USESINGLERP)
    print 'Allowed pixel multiplicity is',ALLOWPIXMULT 

    #start histograms
    ht=HistoTool()   

    if isSignal:
        ht.add(ROOT.TH2F('sighyp', ';Initial category; Final category;Events',16,0,16,16,0,16))
        for i in range(16):
            lab="|{0:04b}>".format(i)
            ht.histos['sighyp']['inc'].GetXaxis().SetBinLabel(i+1,lab)
            ht.histos['sighyp']['inc'].GetYaxis().SetBinLabel(i+1,lab)
    ht.add(ROOT.TH1F('catcount',';Proton selection category;Events',6,0,6))
    for i,c in enumerate(['inc','=2s','mm','ms','sm','ss']):
        ht.histos['catcount']['inc'].GetXaxis().SetBinLabel(i+1,c)

    #main analysis histograms
    ht.add(ROOT.TH1F('nvtx',';Vertex multiplicity;Events',50,0,100))
    ht.add(ROOT.TH1F('rho',';Fastjet #rho;Events',50,0,50))
    ht.add(ROOT.TH1F('xangle',';LHC crossing angle [#murad];Events',4,120,160))
    ht.add(ROOT.TH1F('mll',';Invariant mass [GeV];Events',50,76,106))
    ht.add(ROOT.TH1F('mll_full',';Invariant mass [GeV];Events',50,0,250))
    ht.add(ROOT.TH1F('yll',';Rapidity;Events',50,-3,3))
    ht.add(ROOT.TH1F('etall',';Pseudo-rapidity;Events',50,-6,6))
    ht.add(ROOT.TH1F('ptll',';Transverse momentum [GeV];Events',50,0,250))
    ht.add(ROOT.TH1F('ptll_high',';Transverse momentum [GeV];Events',50,50,500))
    ht.add(ROOT.TH1F('l1eta',';Pseudo-rapidiy;Events',50,0,2.5))
    ht.add(ROOT.TH1F('l1pt',';Transverse momentum [GeV];Events',50,0,250))
    ht.add(ROOT.TH1F('l2eta',';Pseudo-rapidiy;Events',50,0,2.5))
    ht.add(ROOT.TH1F('l2pt',';Transverse momentum [GeV];Events',50,0,250))
    ht.add(ROOT.TH1F('acopl',';A=1-|#Delta#phi|/#pi;Events',50,0,1))
    ht.add(ROOT.TH1F('costhetacs',';cos#theta_{CS};Events',50,-1,1))

    #pileup control
    #ht.add(ROOT.TH1F('rfc',';Random forest classifier probability;Events',50,0,1))
    for d in ['HF','HE','EE','EB']:
        ht.add(ROOT.TH1F('PFMult'+d,';PF multiplicity (%s);Events'%d,50,0,1000))
        ht.add(ROOT.TH1F('PFHt'+d,';PF HT (%s) [GeV];Events'%d,50,0,1000))
        ht.add(ROOT.TH1F('PFPz'+d,';PF P_{z} (%s) [TeV];Events'%d,50,0,40))
    ht.add(ROOT.TH1F('met',';Missing transverse energy [GeV];Events',50,0,200))
    ht.add(ROOT.TH1F('mpf',';MPF;Events',50,-5,5))
    ht.add(ROOT.TH1F('metbits',';MET filters;Events',124,0,124))
    ht.add(ROOT.TH1F('njets',';Jet multiplicity;Events',5,0,5))
    ht.add(ROOT.TH1F('zjb',';Z-jet balance [GeV];Events',50,-150,150))
    ht.add(ROOT.TH1F('zj2b',';Z-2 jets balance [GeV];Events',50,-150,150))
    ht.add(ROOT.TH1F('nch', ';Charged particle multiplicity;Events',50,0,50))
    ht.add(ROOT.TH1F('nextramu',';Additional muons ;Events',10,0,10))
    ht.add(ROOT.TH1F('extramupt',';Additional muon p_{T} [GeV] ;Events',10,0,50))
    ht.add(ROOT.TH1F('extramueta',';Additional muon pseudo-rapidty ;Events',10,0,2.5))

    #RP control
    ht.add(ROOT.TH1F('mpp',';Di-proton invariant mass [GeV];Events',50,0,3000))
    ht.add(ROOT.TH1F('pzpp',';Di-proton p_{z} [GeV];Events',50,-750,750))
    ht.add(ROOT.TH1F('ypp',';Di-proton rapidity;Events',50,-2.5,2.5))
    ht.add(ROOT.TH1F('mmass_full',';Missing mass [GeV];Events',50,-1000,3000))
    ht.add(ROOT.TH1F('mmass',';Missing mass [GeV];Events',50,0,3000))
    ht.add(ROOT.TH1F('ntk',';Track multiplicity;Events',5,0,5))
    ht.add(ROOT.TH1F('ppcount',';pp candidates;Events',3,0,3))
    ht.add(ROOT.TH1F('csi',';#xi;Events',50,0,0.3))
        
    nEntries=tree.GetEntries()              
    print '....analysing',nEntries,'in',inFile,', with output @',outFileName
    if maxEvents>0:
        nEntries=min(maxEvents,nEntries)
        print '      will process',nEntries,'events'

    #compute number of events weighted by target pz spectrum
    nSignalWgtSum=0.
    if isSignal and not isSDZSim:
        print 'Checking how many events are in the fiducial RP area...'
        for i in xrange(0,nEntries):
            tree.GetEntry(i)
            nSignalWgtSum += ROOT.TMath.Gaus(tree.gen_pzpp,0,0.391*gen_mX+624)
        print '...signal weight sum set to',nSignalWgtSum,' from ',nEntries,'raw events'
    else:
        nSignalWgtSum=1

    #start output and tree
    fOut=ROOT.TFile.Open(outFileName,'RECREATE')
    evSummary=EventSummary()
    tOut=ROOT.TTree('data','data')
    evSummary.attachToTree(tOut)

    #summary events for the mixing
    rpData={}

    #loop over events
    nfail=[0,0,0]
    for i in xrange(0,nEntries):

        tree.GetEntry(i)
        evRun=tree.run if isData else -1

        if i%500==0 : 
            drawProgressBar(float(i)/float(nEntries))

        if isSignal:
            if isPhotonSignal     and tree.evcat!=SINGLEPHOTON : continue
            if not isPhotonSignal and tree.evcat==SINGLEPHOTON : continue
        
        isOffZ=False
        if tree.mboson>101:
            if tree.evcat==DIELECTRONS or tree.evcat==DIMUONS:
                isOffZ=True

        #base event selection
        if tree.evcat==DIELECTRONS and tree.isZ: 
            evcat='ee'
        elif tree.evcat==EMU       and not tree.isSS: 
            evcat='em'
        elif tree.evcat==DIMUONS   and tree.isZ: 
            evcat='mm'
        elif isOffZ:
            evcat='offz'
        elif tree.evcat==SINGLEPHOTON and (isPhotonSignal or tree.hasATrigger) :         
            evcat="a" 
        elif tree.evcat==0 and tree.hasZBTrigger : 
            evcat=='zbias'
        else : 
            nfail[0]+=1
            continue

        #assign data-taking era and crossing angle
        evEra = era
        beamXangle = tree.beamXangle        
        if not isData:
            evEra=getRandomEra(isSignal,isPreTS2Signal)            
            if not isSignal: 
                xbin=evMixTool.getRandomLHCCrossingAngle(evEra=evEra,
                                                         evCat=SINGLEPHOTON if tree.isA else DIMUONS)
                beamXangle=VALIDLHCXANGLES[xbin]

        #check if RP is in (MC assume true by default)
        isRPIn=False if isData else True
        if isData and beamXangle in VALIDLHCXANGLES and isValidRunLumi(evRun,tree.lumi,runLumiList):
            isRPIn=True
        if not isRPIn : nfail[1]+=1

        #lepton kinematics
        l1p4=ROOT.TLorentzVector(0,0,0,0)
        l2p4=ROOT.TLorentzVector(0,0,0,0)   
        costhetacs=0
        acopl=0
        if tree.evcat!=SINGLEPHOTON and tree.evcat!=0:
            acopl=1.0-abs(ROOT.TVector2.Phi_mpi_pi(tree.l1phi-tree.l2phi))/ROOT.TMath.Pi()
            l1p4.SetPtEtaPhiM(tree.l1pt,tree.l1eta,tree.l1phi,tree.ml1)
            l2p4.SetPtEtaPhiM(tree.l2pt,tree.l2eta,tree.l2phi,tree.ml2)
            costhetacs=computeCosThetaStar(l1p4,l2p4)
            
            #force ordering by pT (before they were ordered by charge to compute costhetacs)
            if l1p4.Pt()<l2p4.Pt(): 
                l1p4,l2p4=l2p4,l1p4

        #boson kinematics
        boson=ROOT.TLorentzVector(0,0,0,0)
        boson.SetPtEtaPhiM(tree.bosonpt,tree.bosoneta,tree.bosonphi,tree.mboson)
        isZ=tree.isZ
        isA=tree.isA
       
        #for the signal-electron hypothesis this cut needs to be applied
        hasEEEBTransition=False
        if isZ:
            if abs(l1p4.Eta())>1.4442 and abs(l1p4.Eta())<1.5660 : hasEEEBTransition=True
            if abs(l2p4.Eta())>1.4442 and abs(l2p4.Eta())<1.5660 : hasEEEBTransition=True

        #PU-related variables
        #for signal most of these will be overriden by mixing       
        n_extra_mu,nvtx,nch,rho,met,njets,PFMultSumHF,PFHtSumHF,PFPzSumHF,rfc=0,0,0,0,0,0,0,0,0,0
        mpf,zjb,zj2b=0,0,0
        extra_muons=[]
        if isFullSimSignal or not isSignal:
            nvtx=tree.nvtx
            nch=tree.nchPV
            rho=tree.rho
            met=tree.met_pt
            metphi=tree.met_phi
            njets=tree.nj 
            PFMultSumHF=tree.PFMultSumHF
            PFHtSumHF=tree.PFHtSumHF
            PFPzSumHF=tree.PFPzSumHF
            #rfc=getattr(tree,'rfc_%d'%beamXangle)
            for im in range(tree.nrawmu):
                mup4=ROOT.TLorentzVector(0,0,0,0)
                mup4.SetPtEtaPhiM(tree.rawmu_pt[im],tree.rawmu_eta[im]/10.,tree.rawmu_phi[im]/10.,0.105)
                if mup4.DeltaR(l1p4)<0.05 : continue
                if mup4.DeltaR(l2p4)<0.05 : continue
                extra_muons.append( ROOT.TLorentzVector(mup4) )
            n_extra_mu=len(extra_muons)
            
            metp4=ROOT.TLorentzVector(0,0,0,0)
            metp4.SetPtEtaPhiM(met,0,metphi,0)
            mpf=1.+(metp4.Px()*boson.Px()+metp4.Py()*boson.Py())/(boson.Pt()**2+1.0e-6)            
            if njets>0: 
                zjb=tree.j1pt-boson.Pt()
                if njets>1:
                    j1p4=ROOT.TLorentzVector(0,0,0,0)
                    j1p4.SetPtEtaPhiM(tree.j1pt,tree.j1eta,tree.j1phi,tree.j1m)
                    j2p4=ROOT.TLorentzVector(0,0,0,0)
                    j2p4.SetPtEtaPhiM(tree.j2pt,tree.j2eta,tree.j2phi,tree.j2m)
                    zj2b=(j1p4+j2p4).Pt()-boson.Pt()

        #proton tracks (standard and mixed)
        ev_pos_protons,ev_neg_protons = [[],[],[]],[[],[],[]]
        ev_pos_protons_xy,ev_neg_protons_xy = [[],[],[]],[[],[],[]]

        ppsPosEff,ppsPosEffUnc=1.0,0.0
        ppsNegEff,ppsNegEffUnc=1.0,0.0
        if isSignal or (isData and isRPIn):     
            ev_pos_protons,ev_neg_protons  = getTracksPerRomanPot(tree,evEra,evRun if isData else -1,beamXangle,minCsi=MINCSI)              
            orig_ev_pos_protons = copy.deepcopy(ev_pos_protons)
            orig_ev_neg_protons = copy.deepcopy(ev_neg_protons)
            ev_pos_protons_xy,ev_neg_protons_xy  = getTracksPerRomanPot(tree,evEra,evRun if isData else -1,beamXangle,minCsi=MINCSI,useXY=True)  

        #if data and there is nothing to mix store the main characteristics of the event and continue
        if evMixTool.isIdle():
            if isData and isRPIn:
                if (isZ and tree.evcat==DIMUONS and boson.Pt()<10) or evcat=='em':
                    rpDataKey=(evEra,beamXangle,int(tree.evcat))
                    if not rpDataKey in rpData: rpData[rpDataKey]=[]
                    _,ev_pos_protons,ev_neg_protons=doFinalCheck2017(ev_pos_protons,ev_neg_protons,evRun,evEra)
                    rpData[rpDataKey].append( MixedEventSummary(puDiscr=[len(extra_muons),nvtx,rho,PFMultSumHF,PFHtSumHF,PFPzSumHF,rfc],
                                                                pos_protons=ev_pos_protons,
                                                                neg_protons=ev_neg_protons) )
            continue

        #event mixing
        mixed_pos_protons,mixed_neg_protons,mixed_pudiscr=evMixTool.getNew(evEra=evEra,
                                                                           beamXangle=beamXangle,
                                                                           isData=isData,
                                                                           validAngles=VALIDLHCXANGLES,
                                                                           mixEvCategs=[DIMUONS,EMU])

        #prepare efficiencies per arm
        ppsEff,ppsEffUnc=1.0,0.0
        if isSignal:
            rawSigHyp=0
            if len(ev_neg_protons[1])>0: rawSigHyp += 1
            if len(ev_neg_protons[0])>0: rawSigHyp += 2
            if len(ev_pos_protons[1])>0: rawSigHyp += 4
            if len(ev_pos_protons[0])>0: rawSigHyp += 8

            debugSig=False #True if 'B' in evEra else False

            if debugSig:
                print '\n'
                print '{} {:b} -> {:b}'.format(evEra,rawSigHyp,sighyp)
                print '\t(+)',ev_pos_protons,'\nt\t(-)',ev_neg_protons

            #multi-RP
            ppsMultiPosEff,ppsMultiPosEffUnc=0.0,0.0
            ppsMultiNegEff,ppsMultiNegEffUnc=0.0,0.0
            if len(ev_pos_protons[0])>0 and len(ev_pos_protons[2])>0:
                x,y=ev_pos_protons_xy[0][0] if len(ev_pos_protons_xy[0])>0 else [-99,-99]               
                ppsMultiPosEff,ppsMultiPosEffUnc=ppsEffReader.getPPSEfficiency(evEra,beamXangle,
                                                                               ev_pos_protons[2][0],
                                                                               x,
                                                                               y,
                                                                               rp=3,
                                                                               isMulti=True)
                #ppsMultiPosEff,ppsMultiPosEffUnc=1.,0.
            if len(ev_neg_protons[0])>0 and len(ev_neg_protons[2])>0:
                x,y=ev_neg_protons_xy[0][0] if len(ev_neg_protons_xy[0])>0 else [-99,-99]                
                ppsMultiNegEff,ppsMultiNegEffUnc=ppsEffReader.getPPSEfficiency(evEra,beamXangle,
                                                                               ev_neg_protons[2][0],
                                                                               x,
                                                                               y,
                                                                               rp=103,
                                                                               isMulti=True)
                #ppsMultiNegEff,ppsMultiNegEffUnc=1.,0.

            #pixels
            ppsPixelPosEff,ppsPixelPosEffUnc=0.0,0.0
            ppsPixelNegEff,ppsPixelNegEffUnc=0.0,0.0
            if len(ev_pos_protons[1])>0:
                x,y=ev_pos_protons_xy[1][0] if len(ev_pos_protons_xy[1])>0 else [-99,-99]                
                ppsPixelPosEff,ppsPixelPosEffUnc=ppsEffReader.getPPSEfficiency(evEra,beamXangle,
                                                                               ev_pos_protons[1][0],
                                                                               x,
                                                                               y,
                                                                               rp=3,
                                                                               isMulti=False)
                #ppsPixelPosEff,ppsPixelPosEffUnc=1.0,0.0

            if len(ev_neg_protons[1])>0:
                x,y=ev_neg_protons_xy[1][0] if len(ev_neg_protons_xy[1])>0 else [-99,-99]                
                ppsPixelNegEff,ppsPixelNegEffUnc=ppsEffReader.getPPSEfficiency(evEra,beamXangle,
                                                                               ev_neg_protons[1][0],
                                                                               x,
                                                                               y,
                                                                               rp=103,
                                                                               isMulti=False)
                #ppsPixelNegEff,ppsPixelNegEffUnc=1.0,0.0
                
            if debugSig: 
                print '\tEfficiencies'
                print '\t\t[+multi]',ppsMultiPosEff,'+/-',ppsMultiPosEffUnc
                print '\t\t[-multi]',ppsMultiNegEff,'+/-',ppsMultiNegEffUnc
                print '\t\t[+px]',ppsPixelPosEff,'+/-',ppsPixelPosEffUnc
                print '\t\t[-px]',ppsPixelNegEff,'+/-',ppsPixelNegEffUnc

            #assign the final list of reconstructed protons depending on how the sighyp is requested
            ev_pos_protons,ev_neg_protons,ppsEff,ppsEffUnc = ppsEffReader.getProjectedFinalState( ev_pos_protons, ppsMultiPosEff, ppsMultiPosEffUnc, ppsPixelPosEff, ppsPixelPosEffUnc,
                                                                                                  ev_neg_protons, ppsMultiNegEff, ppsMultiNegEffUnc, ppsPixelNegEff, ppsPixelNegEffUnc,
                                                                                                  sighyp,evRun,evEra)
            
            if debugSig:
                print '\t[+]',ev_pos_protons
                print '\t[-]',ev_neg_protons
                print '\tEfficiency:',ppsEff,'+/-',ppsEffUnc
                print '='*100
            #mixed_pos_protons={DIMUONS:ev_pos_protons,EMU:ev_pos_protons}
            #mixed_neg_protons={DIMUONS:ev_neg_protons,EMU:ev_neg_protons}
            mixed_pos_protons, mixed_neg_protons = evMixTool.mergeWithMixedEvent(ev_pos_protons, 
                                                                                 mixed_pos_protons,
                                                                                 ev_neg_protons,
                                                                                 mixed_neg_protons)


            orig_mixed_pos_protons, orig_mixed_neg_protons = evMixTool.mergeWithMixedEvent(orig_ev_pos_protons, 
                                                                                           mixed_pos_protons,
                                                                                           orig_ev_neg_protons,
                                                                                           mixed_neg_protons)


            #control before and after projection
            ht.fill((rawSigHyp,sighyp,1.0),    'sighyp',  ['raw'])
            ht.fill((rawSigHyp,sighyp,ppsEff), 'sighyp',  ['wgt'])

            n_extra_mu,nvtx,rho,PFMultSumHF,PFHtSumHF,PFPzSumHF,rfc = mixed_pudiscr[DIMUONS]

        #kinematics using RP tracks
        pos_protons = ev_pos_protons if isData else mixed_pos_protons[DIMUONS]
        neg_protons = ev_neg_protons if isData else mixed_neg_protons[DIMUONS]        
        _,pos_protons,neg_protons=doFinalCheck2017(pos_protons,neg_protons,evRun,evEra)
        proton_cat,csi_pos,csi_neg,ppSystem,mmassSystem = getDiProtonCategory(pos_protons,neg_protons,boson,ALLOWPIXMULT)
     
        #compare categorization with fully exclusive selection of pixels
        ht.fill((0,ppsEff),'catcount',['inc'])
        if len(pos_protons[1]) in ALLOWPIXMULT and len(neg_protons[1]) in ALLOWPIXMULT : 
            ht.fill((1,ppsEff),'catcount',['inc'])
        ht.fill((proton_cat+1,ppsEff),'catcount',['inc'])
        if isSignal:
            ht.fill((0,1.),'catcount',['single'])
            if len(orig_mixed_pos_protons[DIMUONS][1]) in ALLOWPIXMULT and len(orig_mixed_neg_protons[DIMUONS][1]) in ALLOWPIXMULT:
                ht.fill((1,1.),'catcount',['single'])


        #event categories
        cats=[]            
        cats.append(evcat)
        if isRPIn:
            cats += [evcat+'rpin']
            if proton_cat>0:
                cats += [evcat+'rpinhpur']

        #fill control plots (for signal correct wgt by ee efficiency curve and duplicate for mm channel)
        wgt        = tree.evwgt
        finalPlots = [[wgt,cats]]        
        gen_pzpp   = 0
        gen_pzwgt  = [1.,1.,1.]
        gen_csiPos = 0.
        gen_csiNeg = 0.
        if isSignal:
            true_pos_protons,true_neg_protons = getTracksPerRomanPot(tree,evEra,evRun if isData else -1,beamXangle,True)
            if len(true_pos_protons[0])>0 : gen_csiPos=true_pos_protons[0][0]
            if len(true_neg_protons[0])>0 : gen_csiNeg=true_neg_protons[0][0]
            
            gen_pzpp     = tree.gen_pzpp
            if not isSDZSim:
                pzwid        = 0.391*gen_mX+624
                gen_pzwgt[0] = ROOT.TMath.Gaus(gen_pzpp,0,pzwid)
                gen_pzwgt[1] = ROOT.TMath.Gaus(gen_pzpp,0,pzwid*1.1)/gen_pzwgt[0]
                gen_pzwgt[2] = ROOT.TMath.Gaus(gen_pzpp,0,pzwid*0.9)/gen_pzwgt[0]

            #use the sum of pz weighted events as normalization factor
            if not isFullSimSignal and not isSDZSim:
                if isZ:
                    finalPlots=[ [wgt*ppsEff*gen_pzwgt[0]*mcEff['eez'].Eval(boson.Pt())/nSignalWgtSum , cats],
                                 [wgt*ppsEff*gen_pzwgt[0]*mcEff['mmz'].Eval(boson.Pt())/nSignalWgtSum, [c.replace(evcat,'mm') for c in cats if c[0:2]=='ee']] ]

                    #reject Z->ee if one electron in the transition
                    if hasEEEBTransition:
                        finalPlots[0][0]=0.

                elif isPhotonSignal:
                    finalPlots=[ [wgt*ppsEff*gen_pzwgt[0]*mcEff['a'].Eval(boson.Pt())/nSignalWgtSum, cats] ]
            else:
                finalPlots=[ [wgt*ppsEff*gen_pzwgt[0]/nSignalWgtSum, cats] ]

        for pwgt,pcats in finalPlots:   

            #fill plots only with fiducial signal contribution
            if isSignal and not isSignalFiducial(gen_csiPos,gen_csiNeg,tree.gen_pzpp): continue

            #boson kinematics
            ht.fill((l1p4.Pt(),pwgt),             'l1pt',         pcats)
            ht.fill((l2p4.Pt(),pwgt),             'l2pt',         pcats)
            ht.fill((abs(l1p4.Eta()),pwgt),       'l1eta',        pcats)
            ht.fill((abs(l2p4.Eta()),pwgt),       'l2eta',        pcats)
            ht.fill((acopl,pwgt),                 'acopl',        pcats)
            ht.fill((boson.M(),pwgt),             'mll',          pcats)
            ht.fill((boson.M(),pwgt),             'mll_full',     pcats)
            ht.fill((boson.Rapidity(),pwgt),      'yll',          pcats)
            ht.fill((boson.Eta(),pwgt),           'etall',        pcats)
            ht.fill((boson.Pt(),pwgt),            'ptll',         pcats)
            ht.fill((boson.Pt(),pwgt),            'ptll_high',    pcats)
            ht.fill((costhetacs,pwgt),            'costhetacs',   pcats)
            
            #pileup related
            ht.fill((beamXangle,pwgt),            'xangle', pcats)
            ht.fill((nvtx,pwgt),                  'nvtx',   pcats)
            ht.fill((rho,pwgt),                   'rho',    pcats)
            ht.fill((met,pwgt),                   'met',    pcats)
            ht.fill((mpf,pwgt),                   'mpf',    pcats)
            ht.fill((njets,pwgt),                 'njets',  pcats)
            if njets>0: ht.fill((zjb,pwgt),       'zjb',    pcats)
            if njets>1: ht.fill((zj2b,pwgt),      'zj2b',   pcats)
            ht.fill((nch,pwgt),                   'nch',    pcats) 
            #ht.fill((getattr(tree,'rfc_%d'%beamXangle),pwgt), 'rfc',         pcats)
            ht.fill((PFMultSumHF,pwgt),     'PFMultHF',    pcats)
            ht.fill((PFHtSumHF,pwgt),       'PFHtHF',      pcats)
            ht.fill((PFPzSumHF/1.e3,pwgt),  'PFPZHF',      pcats)
            ht.fill((n_extra_mu,pwgt), 'nextramu', pcats)
            if isFullSimSignal or not isSignal:
                ht.fill((tree.metfilters,pwgt), 'metbits', pcats)
                for sd in ['HE','EE','EB']:
                    ht.fill((getattr(tree,'PFMultSum'+sd),pwgt),    'PFMult'+sd, pcats)
                    ht.fill((getattr(tree,'PFHtSum'+sd),pwgt),      'PFHt'+sd,   pcats)
                    ht.fill((getattr(tree,'PFPzSum'+sd)/1.e3,pwgt), 'PFPZ'+sd,   pcats)
                for mp4 in extra_muons:
                    ht.fill((mp4.Pt(),pwgt), 'extramupt', pcats)
                    ht.fill((abs(mp4.Eta()),pwgt), 'extramueta', pcats)

            #proton counting and kinematics
            for ip in range(3):
                for irp,rpside in [(0,'%dpos'%ip),(1,'%dneg'%ip)]:
                    csiColl=pos_protons[ip] if irp==0 else neg_protons[ip]
                    ht.fill((len(csiColl),pwgt), 'ntk', pcats,rpside)
                    for csi in csiColl:
                        ht.fill((csi,pwgt), 'csi', pcats,rpside)                        

            #diproton kinematics
            if proton_cat<0:
                ht.fill((0,pwgt), 'ppcount', pcats)
            else:
                ht.fill((1,pwgt),                   'ppcount', pcats)
                ht.fill((ppSystem.M(),pwgt),        'mpp',     pcats)
                ht.fill((ppSystem.Pz(),pwgt),       'pzpp',    pcats)
                ht.fill((ppSystem.Rapidity(),pwgt), 'ypp',     pcats)                    
                mmass=mmassSystem.M()
                ht.fill((mmass,pwgt), 'mmass_full', pcats)
                ht.fill((mmass,pwgt), 'mmass_full', pcats, '%d'%proton_cat)
                if mmass>0:
                    ht.fill((mmass,pwgt), 'mmass',  pcats)
                    ht.fill((mmass,pwgt), 'mmass',  pcats, '%d'%proton_cat)

            #signal characteristics in the absense of pileup and other effects
            if isSignal:
                nopu_proton_cat,nopu_csi_pos,nopu_csi_neg,nopu_ppSystem,nopu_mmassSystem = getDiProtonCategory(ev_pos_protons,ev_neg_protons,boson,ALLOWPIXMULT)
                if nopu_proton_cat>0:
                    nopu_mmass = nopu_mmassSystem.M()
                    ht.fill((nopu_ppSystem.M(),pwgt),  'mpp',         pcats, 'nopu')
                    ht.fill((nopu_mmass,pwgt),         'mmass_full',  pcats, 'nopu')
                    ht.fill((nopu_mmass,pwgt),         'mmass_full',  pcats, '%dnopu'%nopu_proton_cat)
                    if nopu_mmass>0:
                        ht.fill((nopu_mmass,pwgt), 'mmass',  pcats, 'nopu')
                        ht.fill((nopu_mmass,pwgt), 'mmass',  pcats, '%dnopu'%nopu_proton_cat)

        if not isData and not isSignal and not isDY : continue

        #save the event summary for the statistical analysis        
        nMixTries=100 if isData else 1
        for itry in range(2*nMixTries+1):

            itry_wgt=wgt
            
            #nominal 
            if itry==0:
                mixType           = 0 
                i_pos_protons     = copy.deepcopy(pos_protons)
                i_neg_protons     = copy.deepcopy(neg_protons)

                #shift csi by 1%
                i_pos_protons_syst=[]
                i_neg_protons_syst=[]
                for ialgo in range(3):
                    i_pos_protons_syst.append( [1.01*x for x in pos_protons[ialgo]] )
                    i_neg_protons_syst.append( [1.01*x for x in neg_protons[ialgo]] )

            else:
                #get a new event to mix
                i_mixed_pos_protons, i_mixed_neg_protons, i_mixed_pudiscr = evMixTool.getNew(evEra=evEra,
                                                                                             beamXangle=beamXangle,
                                                                                             isData=isData,
                                                                                             validAngles=VALIDLHCXANGLES,
                                                                                             mixEvCategs=[DIMUONS,EMU])

                #FIXME this is broken in this new version
                #if MIXEDRPSIG:
                #    sigCsi=random.choice( MIXEDRPSIG[beamXangle] )
                #    for mixEvCat in mixed_far_rptks:
                #        tksPos=mixed_far_rptks[mixEvCat][0]+[sigCsi[0]]
                #        shuffle(tksPos)
                #        tksNeg=mixed_far_rptks[mixEvCat][1]+[sigCsi[1]]
                #        shuffle(tksNeg)
                #        mixed_far_rptks[mixEvCat]=(tksPos,tksNeg)

                #merge signal protons with pileup protons for first attempt
                if isSignal and itry==1:
                    i_mixed_pos_protons, i_mixed_neg_protons = evMixTool.mergeWithMixedEvent(ev_pos_protons, 
                                                                                             i_mixed_pos_protons,
                                                                                             ev_neg_protons,
                                                                                             i_mixed_neg_protons)
                    if not isFullSimSignal:
                        n_extra_mu,nvtx,rho,PFMultSumHF,PFHtSumHF,PFPzSumHF,rfc = i_mixed_pudiscr[DIMUONS]

                itry_wgt = wgt/float(nMixTries)

                if itry<=nMixTries or isSignal:
                    mixType           = 1
                    if isSignal: mixType=itry
                    i_pos_protons      = i_mixed_pos_protons[DIMUONS]
                    i_neg_protons      = i_mixed_neg_protons[DIMUONS]
                    i_pos_protons_syst = i_mixed_pos_protons[EMU]
                    i_neg_protons_syst = i_mixed_neg_protons[EMU]
                else:
                    mixType            = 2
                    i_pos_protons      = i_mixed_pos_protons[DIMUONS]
                    i_neg_protons      = copy.deepcopy(neg_protons)
                    i_pos_protons_syst = copy.deepcopy(pos_protons)
                    i_neg_protons_syst = i_mixed_neg_protons[DIMUONS]                   

            _,i_pos_protons,i_neg_protons=doFinalCheck2017(i_pos_protons,i_neg_protons,evRun,evEra)
            i_proton_cat,     i_csi_pos,     i_csi_neg,     i_ppSystem,     i_mmassSystem      = getDiProtonCategory(i_pos_protons,     i_neg_protons,      boson,ALLOWPIXMULT)
            _,i_pos_protons_syst,i_neg_protons_syst=doFinalCheck2017(i_pos_protons_syst,i_neg_protons_syst,evRun,evEra)
            i_proton_cat_syst,i_csi_pos_syst,i_csi_neg_syst,i_ppSystem_syst,i_mmassSystem_syst = getDiProtonCategory(i_pos_protons_syst,i_neg_protons_syst, boson,ALLOWPIXMULT)
            
            #if itry>nMixTries:
            #    print itry,mixType
            #    print '\t',i_pos_protons,i_pos_protons_syst
            #    print '\t--->',i_proton_cat,     i_csi_pos,     i_csi_neg
            #    print '\t',i_neg_protons,i_neg_protons_syst
            #    print '\t--->',i_proton_cat_syst,i_csi_pos_syst,i_csi_neg_syst

            #passAtLeastOneSelection=(i_proton_cat>0 or i_proton_cat_syst>0)
            passAtLeastOneSelection=(isSignal and mixType==1) or (i_proton_cat>0 or i_proton_cat_syst>0)

            #start event summary
            evSummary.reset()
            evSummary.sighyp[0]=int(sighyp)
            if isData:
                evSummary.run[0]=int(evRun)
                evSummary.event[0]=long(tree.event)
                evSummary.lumi[0]=int(tree.lumi)

            evSummary.era[0]=int(ord(evEra[-1]))            
            evSummary.cat[0]=int(tree.evcat)
            evSummary.isOffZ[0]=int(isOffZ)
            evSummary.wgt[0]=itry_wgt
            evSummary.xangle[0]=int(beamXangle)
            evSummary.l1pt[0]=l1p4.Pt()
            evSummary.l1eta[0]=l1p4.Eta()
            evSummary.l2pt[0]=l2p4.Pt()
            evSummary.l2eta[0]=l2p4.Eta()
            evSummary.bosonm[0]=boson.M()
            evSummary.bosonpt[0]=boson.Pt()
            evSummary.bosoneta[0]=boson.Eta()
            evSummary.bosony[0]=boson.Rapidity()
            evSummary.acopl[0]=acopl
            evSummary.costhetacs[0]=costhetacs            
            evSummary.njets[0]=int(njets)
            evSummary.mpf[0]=mpf
            evSummary.zjb[0]=zjb
            evSummary.zj2b[0]=zj2b
            evSummary.nch[0]=int(nch)
            evSummary.nvtx[0]=int(nvtx)
            evSummary.rho[0]=rho
            evSummary.PFMultSumHF[0]=int(PFMultSumHF)
            evSummary.PFHtSumHF[0]=PFHtSumHF
            evSummary.PFPzSumHF[0]=PFPzSumHF
            evSummary.rfc[0]=rfc
            evSummary.gen_pzpp[0]=gen_pzpp
            evSummary.gen_pzwgtUp[0]=gen_pzwgt[1]
            evSummary.gen_pzwgtDown[0]=gen_pzwgt[2]
            evSummary.gencsi1[0]=gen_csiPos
            evSummary.gencsi2[0]=gen_csiNeg

            #vary boson energy scale
            if i_ppSystem:
                boson_up=boson*1.03
                evSummary.mmissvup[0]= buildMissingMassSystem(i_ppSystem,boson_up).M()
                boson_dn=boson*0.97
                evSummary.mmissvdn[0]= buildMissingMassSystem(i_ppSystem,boson_dn).M() 

            evSummary.mixType[0]=mixType
            evSummary.protonCat[0]=i_proton_cat
            if i_proton_cat>0:
                evSummary.csi1   [0]= i_csi_pos
                evSummary.csi2   [0]= i_csi_neg
                evSummary.mpp    [0]= i_ppSystem.M()
                evSummary.ypp    [0]= i_ppSystem.Rapidity()
                evSummary.pzpp   [0]= i_ppSystem.Pz()
                evSummary.mmiss  [0]= i_mmassSystem.M()
                evSummary.ymmiss [0]= i_mmassSystem.Rapidity()
                evSummary.ppsEff[0]=ppsEff
                evSummary.ppsEffUnc[0]=ppsEffUnc                

            evSummary.systprotonCat[0]=i_proton_cat_syst
            if i_proton_cat_syst>0:
                evSummary.systcsi1   [0]= i_csi_pos_syst
                evSummary.systcsi2   [0]= i_csi_neg_syst
                evSummary.systmpp    [0]= i_ppSystem_syst.M()
                evSummary.systypp    [0]= i_ppSystem_syst.Rapidity()
                evSummary.systpzpp   [0]= i_ppSystem_syst.Pz()
                evSummary.systmmiss  [0]= i_mmassSystem_syst.M()
                evSummary.systymmiss [0]= i_mmassSystem_syst.Rapidity()
                evSummary.systppsEff[0]=ppsEff
                evSummary.systppsEffUnc[0]=ppsEffUnc                

            #if no selection passes the cuts ignore its summary
            if not passAtLeastOneSelection: continue

            #for signal update the event weight for ee/mm/photon hypothesis
            if isData or isDY:
                tOut.Fill()

            elif isFullSimSignal:
                origWgt          = evSummary.wgt[0]
                evSummary.wgt[0] = origWgt*gen_pzwgt[0]/nSignalWgtSum
                tOut.Fill()

            elif isSignal:

                origWgt=evSummary.wgt[0]

                if isZ:
                    #add a copy for ee
                    if not hasEEEBTransition:                        
                        evSummary.cat[0]=DIELECTRONS
                        evSummary.wgt[0]=origWgt*gen_pzwgt[0]*mcEff['eez'].Eval(boson.Pt())/nSignalWgtSum
                        tOut.Fill()
                
                    #add a copy for mm
                    evSummary.cat[0]=DIMUONS
                    evSummary.wgt[0]=origWgt*gen_pzwgt[0]*mcEff['mmz'].Eval(boson.Pt())/nSignalWgtSum
                    tOut.Fill()
                
                if isA:
                    #add a copy for the photon
                    evSummary.cat[0]=SINGLEPHOTON
                    evSummary.wgt[0]=origWgt*gen_pzwgt[0]*mcEff['a'].Eval(boson.Pt())/nSignalWgtSum
                    tOut.Fill()

    #dump events for the mixing
    nSelRPData=sum([len(rpData[x]) for x in rpData])
    if nSelRPData>0:
        rpDataOut=outFileName.replace('.root','.pck')                
        print 'Saving',nSelRPData,'events for mixing in',rpDataOut
        with open(rpDataOut,'w') as cachefile:
            pickle.dump(rpData,cachefile, pickle.HIGHEST_PROTOCOL)        

    #if there was no mixing don't do anything else
    if not MIXEDRP: return 

    #save results
    fOut.cd()
    tOut.Write()
    ht.writeToFile(fOut)
    fOut.Close()


def runExclusiveAnalysisPacked(args):

    """wrapper for parallel execution"""

    try:
        runExclusiveAnalysis(*args)
    except Exception as e:
        print 50*'<'
        print "  Problem with", args[1], "continuing without"
        print e
        print 50*'<'
        return False
    
def runAnalysisTasks(opt):

    """create analysis tasks"""

    #read samples to process
    task_dict={}    
    for jsonF in opt.json.split(','):
        with open(jsonF,'r') as cachefile:
            samples=json.load(cachefile,  encoding='utf-8', object_pairs_hook=OrderedDict).items()
            for x in samples:
                task_dict[str(x[0])]=[]

    #group chunks matching the same name   
    for file_path in os.listdir(opt.input):
        
        if opt.only:
            if not file_path in opt.only: 
                continue

        file_name,ext=os.path.splitext(file_path)
        if ext != '.root' : continue

        #check if file tag is already in the list of samples to process
        isSignal,_=isSignalFile(file_name)
        lastTkn=file_name.rfind('_')
        tag=file_name[0:lastTkn]
        if isSignal:
            if not 'fullsim' in file_name and not 'SDZ' in file_name:
                tag=file_name
            print 'This file will be treated as signal with tag=',tag

        if not tag in task_dict: continue
        task_dict[tag].append( os.path.join(opt.input,file_path) )

    #parse json file with list of run/lumi sections
    with open(opt.RPout,'r') as cachefile:        
        runLumi=json.load(cachefile,  encoding='utf-8', object_pairs_hook=OrderedDict).items()
        runLumi={int(x[0]):x[1] for x in runLumi}

    global ALLOWPIXMULT
    global USESINGLERP
    if opt.allowPix:
        print 'Updating the number of allowed pixels'
        ALLOWPIXMULT=[int(x) for x in opt.allowPix.split(',')]
        if sum([x<0 for x in ALLOWPIXMULT]):
            ALLOWPIXMULT=[abs(x) for x in ALLOWPIXMULT]
            USESINGLERP=True
        print 'Use singleRP=',USESINGLERP, 'with pix multiplicity in=',ALLOWPIXMULT

    #used to test inclusion of potential signal (test case from conveners request)
    global MIXEDRPSIG
    if opt.mixSignal:

        print 'Going to add mixed signal here from %s'%opt.mixSignal
        print 'Use this at your own risk...'

        MIXEDRPSIG={}
        for a in VALIDLHCXANGLES:
            MIXEDRPSIG[a]=[]
            fIn=ROOT.TFile.Open(opt.mixSignal.format(a))
            data=fIn.Get('data')
            for i in range(data.GetEntriesFast()):
                data.GetEntry(i)
                if data.mixType!=0 : continue
                MIXEDRPSIG[a].append( (data.csi1,data.csi2) )
            fIn.Close()
            print '\t',a,'murad has',len(MIXEDRPSIG[a]),'events'
    
    #create the tasks and submit them
    import multiprocessing as MP
    pool = MP.Pool(opt.jobs)
    task_list=[]
    mergeList={}
    for x in task_dict.keys():
        
        isData     = True if 'Data' in x else False
        isSignal,_ = isSignalFile(x)        
        if opt.step==0 and not isData : continue
        runLumiList=runLumi if isData else None
        for f in task_dict[x]:
            if isSignal:
                sigOut='%s/Chunks/%s'%(opt.output,os.path.basename(f))
                mergeList[sigOut]=[]
                for sighyp in range(16):
                    fOut=sigOut.replace('.root','_%d.root'%sighyp)
                    mergeList[sigOut].append(fOut)
                    task_list.append( (f,fOut,runLumiList,opt.effDir,opt.ppsEffFile,opt.maxEvents,sighyp,opt.mix) )
            else:
                fOut='%s/Chunks/%s'%(opt.output,os.path.basename(f))
                task_list.append( (f,fOut,runLumiList,opt.effDir,opt.ppsEffFile,opt.maxEvents,0,opt.mix) )

    pool.map(runExclusiveAnalysisPacked,task_list)

    #signals need to be merged
    for fOut,fList in mergeList.items():
        fListStr=' '.join(fList)
        os.system('hadd -f -k %s %s'%(fOut,fListStr))
        os.system('rm %s'%fListStr)



def main():
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--input',
                      dest='input',   
                      default='/eos/cms/store/cmst3/user/psilva/ExclusiveAna/Chunks',
                      help='input directory with the files [default: %default]')
    parser.add_option('--jobs',
                      dest='jobs', 
                      default=2,
                      type=int,
                      help='# of jobs to process in parallel the trees [default: %default]')
    parser.add_option('--maxEvents',
                      dest='maxEvents', 
                      default=-1,
                      type=int,
                      help='# of events to process [default: %default]')
    parser.add_option('--allowPix',
                      dest='allowPix', 
                      default='1,2',
                      help='number of pixel protons allowed [default: %default]')
    parser.add_option('--json',
                      dest='json', 
                      default='pps_samples.json',
                      type='string',
                      help='json with the files to process')
    parser.add_option('--only',
                      dest='only', 
                      default=None,
                      type='string',
                      help='only this file (or CSV list of files)')
    parser.add_option('--RPout',
                      dest='RPout', 
                      default='golden_noRP.json',
                      type='string',
                      help='json with the runs/lumi sections in which RP are out')
    parser.add_option('--effDir',
                      dest='effDir', 
                      default='/eos/cms//store/cmst3/user/psilva/ExclusiveAna/final/ab05162/plots/',
                      type='string',
                      help='directory with efficiency files for signal weighting')
    parser.add_option('--ppsEffFile',
                      dest='ppsEffFile', 
                      default='${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/pps/PreliminaryEfficiencies_July132020_1D2DMultiTrack.root,${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/pps/pixelEfficiencies_multiRP.root,${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/pps/pixelEfficiencies_radiation.root',
                      type='string',
                      help='file with PPS reconstructed efficiency')
    parser.add_option('--step',
                      dest='step', 
                      default=0,
                      type=int,
                      help='analysis step: 0 - prepare event mixing bank; 1 - analysis')
    parser.add_option('--mix',
                      dest='mix',
                      default=None,
                      type='string',
                      help='bank of events to use for the mixing')
    parser.add_option('--mixSignal',
                      dest='mixSignal',
                      default=None,
                      type='string',
                      help='this is just for a test : what if signal protons are added on top of the pileup protons?')
    parser.add_option('-o', '--output',
                      dest='output', 
                      default='analysis',
                      help='Output directory [default: %default]')
    (opt, args) = parser.parse_args()
    
    #transform CSV to list
    if opt.only: opt.only=[os.path.basename(x) for x in opt.only.split(',')]

    os.system('mkdir -p %s/Chunks' % opt.output)
    runAnalysisTasks(opt)
        

if __name__ == "__main__":
    sys.exit(main())
