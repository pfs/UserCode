#! /usr/bin/env python
import ROOT
import copy
import optparse
import json
import sys
import os
import numpy
import array
import json
import random
import pickle
import array
import re
from random import shuffle
from collections import OrderedDict
from TopLJets2015.TopAnalysis.HistoTool import *
from mixedEvent import *

VALIDLHCXANGLES=[120,130,140,150]
DIMUONS=13*13
EMU=11*13

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

def getTracksPerRomanPot(tree,mcTruth=False,usefarRP=True):

    """loops over the availabe tracks in the event and groups them by roman pot id"""

    tkPos=[]
    tkNeg=[]
    for itk in xrange(0,tree.nRPtk):
        rpid=tree.RPid[itk]
        try:
            csi=tree.RPfarcsi[itk] if usefarRP else tree.RPnearcsi[itk]
        except:            
            if notMCTruth:
                csi=tree.RPfarx[itk] if useFarRP else tree.RPnearx[itk]
            else:
                csi=tree.RPtruecsi[itk] 
        if rpid==23  and csi>0: tkPos.append(csi)
        if rpid==123 and csi>0: tkNeg.append(csi)
        
    return (tkPos,tkNeg)

def buildDiproton(rptks,sqrts=13000.):

    """build a diproton system from to tracks in the arms of the roman pots"""

    if not rptks : return None
    if len(rptks[0])==0 or len(rptks[1])==0 : return None
    beamP=0.5*sqrts
    csiL,csiR=rptks[0][0],rptks[1][0]
    diproton=ROOT.TLorentzVector(0.,0.,beamP*(csiL-csiR),beamP*(csiL+csiR))
    return diproton

def getRandomEra():

    """generates a random era according to the integrated luminosity in each one"""

    r=random.random()
    if r<0.115   : return '2017B'
    elif r<0.348 : return '2017C'
    elif r<0.451 : return '2017D'
    elif r<0.671 : return '2017E'
    return '2017F'

def runExclusiveAnalysis(inFile,outFileName,runLumiList,mixFile):
    
    """event loop"""

    isData=True if 'Data' in inFile else False
    isSignal=True if 'MC13TeV_ppxz_' in inFile else False

    #bind main tree with pileup discrimination tree, if failed return
    tree=ROOT.TChain('analysis/data' if isSignal else 'tree')
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


    #identify data-taking era
    era=None
    if isData:
        era=os.path.basename(inFile).split('_')[1]
    
    #check if it is signal and load     
    signalPt=[]
    mcEff={}
    if isSignal: 
        signalPt=[float(x) for x in re.findall(r'\d+', os.path.basename(inFile) )[2:]]
        for ch in ['eez','mmz','a']:
            effIn=ROOT.TFile.Open('$CMSSW_BASE/src/TopLJets2015/TopAnalysis/plots/effsummary_%s_ptboson.root'%ch)
            pname='gen%srec_ptboson_ZH#rightarrowllbb_eff'%ch
            if ch=='a': pname='genarec_ptboson_EWK #gammajj_eff'
            mcEff[ch]=effIn.Get(pname)
            effIn.Close()        

    #filter events to mix according to tag if needed
    mixedRP=None
    xangleRelFracs={}
    try:
        print 'Analysing mixfile',mixFile
        with open(mixFile,'r') as cachefile:
            mixedRP=pickle.load(cachefile)        

        #build the list of probabilities for the crossing angles in each era
        for key in mixedRP:
            mix_era,mix_xangle,mix_evcat=key
            if not mix_xangle in VALIDLHCXANGLES: continue
            n=len(mixedRP[key])
            xangleKey=(mix_era,mix_evcat)
            if not xangleKey in xangleRelFracs:
                xangleRelFracs[xangleKey]=ROOT.TH1F('xanglefrac_%s_%d'%xangleKey,'',len(VALIDLHCXANGLES),0,len(VALIDLHCXANGLES))
            xbin=(mix_xangle-120)/10+1
            xangleRelFracs[xangleKey].SetBinContent(xbin,n)
        for xangleKey in xangleRelFracs:
            xangleRelFracs[xangleKey].Scale(1./xangleRelFracs[xangleKey].Integral())
    except Exception as e:
        if mixFile : print e
        pass
    
    #start histograms
    ht=HistoTool()   

    #main analysis histograms
    ht.add(ROOT.TH1F('mll',';Dilepton invariant mass [GeV];Events',50,20,250))
    ht.add(ROOT.TH1F('yll',';Dilepton rapidity;Events',50,0,5))
    ht.add(ROOT.TH1F('ptll',';Dilepton transverse momentum [GeV];Events',50,0,250))
    ht.add(ROOT.TH1F('l1eta',';Lepton pseudo-rapidiy;Events',50,0,2.5))
    ht.add(ROOT.TH1F('l1pt',';Lepton transverse momentum [GeV];Events',50,0,250))
    ht.add(ROOT.TH1F('l2eta',';Lepton pseudo-rapidiy;Events',50,0,2.5))
    ht.add(ROOT.TH1F('l2pt',';Lepton transverse momentum [GeV];Events',50,0,250))
    ht.add(ROOT.TH1F('acopl',';A=1-|#Delta#phi|/#pi;Events',50,0,1))
    ht.add(ROOT.TH1F('xangle',';LHC crossing angle [#murad];Events',4,120,160))
    ht.add(ROOT.TH1F('mpp',';Di-proton invariant mass [GeV];Events',50,0,3000))
    ht.add(ROOT.TH1F('ypp',';Di-proton rapidity;Events',50,0,2))
    ht.add(ROOT.TH2F('mpp2d',';Far di-proton invariant mass [GeV];Near di-proton invariant mass [GeV];Events',50,0,3000,50,0,3000))
    ht.add(ROOT.TH2F('ypp2d',';Far di-proton rapidity;Near di-proton rapidity;Events',50,0,2,50,0,2))
    ht.add(ROOT.TH1F('mmass',';Missing mass [GeV];Events',50,0,3000))
    ht.add(ROOT.TH1F('ntk',';Track multiplicity;Events',5,0,5))
    ht.add(ROOT.TH1F('ppcount',';pp candidates;Events',3,0,3))
    ht.add(ROOT.TH1F('csi',';#xi;Events',50,0,0.3))
    ht.add(ROOT.TH2F('csi2d',';#xi(far);#xi(near);Events',50,0,0.3,50,0,0.3))

    #pileup control
    ht.add(ROOT.TH1F('nvtx',';Vertex multiplicity;Events',50,0,100))
    ht.add(ROOT.TH1F('rho',';Fastjet #rho;Events',50,0,30))
    #ht.add(ROOT.TH1F('rfc',';Random forest classifier probability;Events',50,0,1))
    for d in ['HF','HE','EE','EB']:
        ht.add(ROOT.TH1F('PFMult'+d,';PF multiplicity (%s);Events'%d,50,0,1000))
        ht.add(ROOT.TH1F('PFHt'+d,';PF HT (%s) [GeV];Events'%d,50,0,1000))
        ht.add(ROOT.TH1F('PFPz'+d,';PF P_{z} (%s) [TeV];Events'%d,50,0,40))
    ht.add(ROOT.TH1F('met',';Missing transverse energy [GeV];Events',50,0,200))
    ht.add(ROOT.TH1F('metbits',';MET filters;Events',124,0,124))
    ht.add(ROOT.TH1F('njets',';Jet multiplicity;Events',5,0,5))
    ht.add(ROOT.TH1F('nch', ';Charged particle multiplicity;Events',50,0,50))
    ht.add(ROOT.TH1F('nextramu',';Additional muons ;Events',10,0,10))
    ht.add(ROOT.TH1F('extramupt',';Additional muon p_{T} [GeV] ;Events',10,0,50))
    ht.add(ROOT.TH1F('extramueta',';Additional muon pseudo-rapidty ;Events',10,0,2.5))
        
    nEntries=tree.GetEntries()          
    print '....analysing',nEntries,'in',inFile,', with output @',outFileName
    if mixedRP : print '    events mixed with',mixFile

    #loop over events
    rpData={}
    selEvents=[]
    summaryVars='cat:wgt:xangle:'
    summaryVars+='l1pt:l1eta:l2pt:l2eta:acopl:bosonpt:'
    summaryVars+='nch:nvtx:rho:PFMultSumHF:PFHtSumHF:PFPzSumHF:rfc:'
    summaryVars+='csi1:csi2:mpp:mmiss:'
    summaryVars+='mixcsi1:mixcsi2:mixmpp:mixmmiss:'
    summaryVars+='mixemcsi1:mixemcsi2:mixemmpp:mixemmmiss'
    summaryVars=summaryVars.split(':')
    for i in xrange(0,nEntries):

        tree.GetEntry(i)

        #init "golden selected events for final statistcal analysis"
        goldenSel=None

        if i%1000==0 : sys.stdout.write('\r [ %d/100 ] done' %(int(float(100.*i)/float(nEntries))))
    
        #base event selection
        if tree.evcat==11*11 and not tree.isSS : evcat='ee'
        elif tree.evcat==EMU and not tree.isSS : evcat='em'
        elif tree.evcat==DIMUONS and not tree.isSS : evcat='mm'
        elif tree.evcat==22 :
            if isSignal: 
                evcat="a" 
            else:
                if tree.hasATrigger : evcat="a"
        elif tree.evcat==0  and tree.hasZBTrigger : evcat=='zbias'
        else : continue
        
        #assign data-taking era and crossing angle
        evEra = era
        beamXangle = tree.beamXangle
        if not isData:
            evEra=getRandomEra()
            if not isSignal: 
                xangleKey=(evEra,DIMUONS)
                xbin=ROOT.TMath.FloorNint(xangleRelFracs[xangleKey].GetRandom())
                beamXangle=VALIDLHCXANGLES[xbin]

        #data specific filters
        if isData:

            #reject invalid beam crossing angles
            if not beamXangle in VALIDLHCXANGLES : continue

            #check RPs are in
            if not isValidRunLumi(tree.run,tree.lumi,runLumiList): continue




        #lepton kinematics
        l1p4=ROOT.TLorentzVector(0,0,0,0)
        l2p4=ROOT.TLorentzVector(0,0,0,0)        
        acopl=0
        if tree.evcat!=22 and tree.evcat!=0:
            acopl=1.0-abs(ROOT.TVector2.Phi_mpi_pi(tree.l1phi-tree.l2phi))/ROOT.TMath.Pi()
            l1p4.SetPtEtaPhiM(tree.l1pt,tree.l1eta,tree.l1phi,tree.ml1)
            l2p4.SetPtEtaPhiM(tree.l2pt,tree.l2eta,tree.l2phi,tree.ml2)
            if l1p4.Pt()<l2p4.Pt(): l1p4,l2p4=l2p4,l1p4
            #if abs(l1p4.Eta())>2.1 : continue

        #boson kinematics
        boson=ROOT.TLorentzVector(0,0,0,0)
        boson.SetPtEtaPhiM(tree.bosonpt,tree.bosoneta,tree.bosonphi,tree.mboson)
        isZ=tree.isZ
        isA=tree.isA
        isHighPtZ=(boson.Pt()>50)
        isLowPtZ=(boson.Pt()<10)

        #PU-related variables
        #for signal most of these will be overriden by mixing       
        n_extra_mu,nvtx,nch,rho,met,njets,PFMultSumHF,PFHtSumHF,PFPzSumHF,rfc=0,0,0,0,0,0,0,0,0,0
        extra_muons=[]
        if not isSignal:
            nvtx=tree.nvtx
            nch=tree.nchPV
            rho=tree.rho
            met=tree.met_pt
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
            
        #proton tracks (standard and mixed)
        far_rptks, near_rptks = None, None
        if isSignal or isData: 
            far_rptks  = getTracksPerRomanPot(tree)
            near_rptks = getTracksPerRomanPot(tree,False,False)            


        #if data and there is nothing to mix store the main characteristics of the event and continue
        if isData and not mixedRP:

            #for Z->mm use only 10% otherwise we have way too many events to do this efficiently
            if (isZ and tree.evcat==DIMUONS and isLowPtZ and random.random()<0.1) or evcat=='em':

                rpDataKey=(era,beamXangle,int(tree.evcat))
                if not rpDataKey in rpData: rpData[rpDataKey]=[]
                rpData[rpDataKey].append( MixedEvent(beamXangle,
                                                     [len(extra_muons),nvtx,rho,PFMultSumHF,PFHtSumHF,PFPzSumHF,rfc],
                                                     far_rptks,
                                                     near_rptks
                                                     ) )
            continue

        #do the mixing
        mixed_far_rptks,mixed_far_1rptk = {}, {}
        try:
            for mixEvCat in [DIMUONS,EMU]:
                mixedEvKey=(evEra,beamXangle,mixEvCat)
                mixedEv=random.choice( mixedRP[mixedEvKey] )
                mixed_far_rptks[mixEvCat]=mixedEv.far_rptks
                if far_rptks and mixed_far_rptks[mixEvCat]:
                    if random.random()<0.5:
                        mixed_far_1rptk[mixEvCat] = (mixed_far_rptks[mixEvCat][0],far_rptks[1])
                    else:
                        mixed_far_1rptk[mixEvCat] = (far_rptks[0],mixed_far_rptks[mixEvCat][1])

                #for signal add the tracks to the simulated ones
                #assign pileup characteristics from the Z->mm low pT events
                if isSignal:
                    if mixEvCat==DIMUONS:
                        n_extra_mu,nvtx,rho,PFMultSumHF,PFHtSumHF,PFPzSumHF,rfc=mixedEv.puDiscr
                    tksPos=mixed_far_rptks[mixEvCat][0]+far_rptks[0]
                    shuffle(tksPos)
                    tksNeg=mixed_far_rptks[mixEvCat][1]+far_rptks[1]
                    shuffle(tksNeg)
                    mixed_far_rptks[mixEvCat]=(tksPos,tksNeg)

        except Exception as e:
            print e  
            print evEra,beamXangle,mixEvCat,'->',mixedEvKey
            pass

        #prepare the combinations of protons with central detector to fill histograms for
        mon_tasks=[]
        if isData:
            mon_tasks.append( (far_rptks,near_rptks,'') )
            if DIMUONS in mixed_far_1rptk:
                mon_tasks.append( (mixed_far_1rptk[DIMUONS],None,'_mix1') )
                mon_tasks.append( (mixed_far_rptks[DIMUONS],None,'_mix2') )
            if EMU in mixed_far_1rptk:
                mon_tasks.append( (mixed_far_1rptk[EMU],None,'_mixem1') )
                mon_tasks.append( (mixed_far_rptks[EMU],None,'_mixem2') )
        else:
            if DIMUONS in mixed_far_rptks:
                mon_tasks.append( (mixed_far_rptks[DIMUONS],None,'') )            
            if EMU in mixed_far_rptks:
                mon_tasks.append( (mixed_far_rptks[EMU],None,'_mixem2') )            
            if isSignal:
                mon_tasks.append( (far_rptks,near_rptks,'nopu') )
            
        #fill the histograms
        wgt=tree.evwgt
        hasAHighPurSelection=False
        for protons,near_protons, pfix in mon_tasks:

            #no calibration, not worth it...
            if not beamXangle in VALIDLHCXANGLES: continue
    
            #high purity selection for proton tracks
            highPur = True if protons and len(protons[0])==1 and len(protons[1])==1 else False  
            if highPur:
                if protons[0][0]<0.02 or protons[0][0]>0.18 : highPur=False
                if protons[1][0]<0.02 or protons[1][0]>0.18 : highPur=False

            #check if Z and di-proton combination is consistent with elastic scattering
            pp=buildDiproton(protons)
            near_pp=buildDiproton(near_protons) if near_protons else None
            isElasticLike=False
            mmass=0
            if pp:
                isElasticLike=(13000.-boson.E()-pp.E()>0)
                inPP=ROOT.TLorentzVector(0,0,0,13000.)
                if isElasticLike: 
                    mmass=(pp-boson).M()

            #categories to fill
            cats=[]            
            cats.append(evcat)
            if isZ : 
                cats.append(evcat+'Z')
                if isHighPtZ : cats.append(evcat+'hptZ')
                if isLowPtZ  : cats.append(evcat+'lptZ')
            if isElasticLike and highPur :
                ppCats=[c+'hpur' for c in cats]
                cats += ppCats
            beamAngleCats=[c+'%d'%beamXangle for c in cats]
            cats += beamAngleCats
                

            #save he basic info on golden events
            saveGoldenSel=False
            if isSignal and pfix in ['','_mixem2']:
                saveGoldenSel=True
            if isData and (isZ or isA):
                if pfix in ['','_mix2','_mixem2']:
                    saveGoldenSel=True                

            if saveGoldenSel:                    
                if not goldenSel:
                    goldenSel=[tree.evcat,wgt,beamXangle,
                               l1p4.Pt(),l1p4.Eta(),l2p4.Pt(),l2p4.Eta(),acopl,boson.Pt(),
                               nch,nvtx,rho,PFMultSumHF,PFHtSumHF,PFPzSumHF,rfc]
                task_protonInfo=[0,0,0,0]
                if isElasticLike and highPur:
                    hasAHighPurSelection=True
                    task_protonInfo=[protons[0][0],protons[1][0],pp.M(),mmass]
                        
                goldenSel += task_protonInfo


            #final plots (for signal correct wgt by efficiency curve and duplicate for mm channel)
            finalPlots=[[wgt,cats]]
            if isSignal:
                #signal has been pre-selected in the fiducial phase space so nEntries is 
                #the effective number of events (or sum of weights) generated
                if isZ:
                    finalPlots=[ [wgt*mcEff['eez'].Eval(boson.Pt())/nEntries, cats],
                                 [wgt*mcEff['mmz'].Eval(boson.Pt())/nEntries, [c.replace(evcat,'mm') for c in cats if c[0:2]=='ee']] ]
                else:
                    finalPlots=[ [wgt*mcEff['a'].Eval(boson.Pt())/nEntries, cats] ]

            for pwgt,pcats in finalPlots:   

                #boson kinematics
                ht.fill((l1p4.Pt(),pwgt),             'l1pt',   pcats,pfix)
                ht.fill((l2p4.Pt(),pwgt),             'l2pt',   pcats,pfix)
                ht.fill((abs(l1p4.Eta()),pwgt),       'l1eta',  pcats,pfix)
                ht.fill((abs(l2p4.Eta()),pwgt),       'l2eta',  pcats,pfix)
                ht.fill((acopl,pwgt),                 'acopl',  pcats,pfix) 
                ht.fill((boson.M(),pwgt),             'mll',    pcats,pfix)
                ht.fill((abs(boson.Rapidity()),pwgt), 'yll',    pcats,pfix)
                ht.fill((boson.Pt(),pwgt),            'ptll',   pcats,pfix)
                
                #pileup related
                ht.fill((beamXangle,pwgt),            'xangle', pcats,pfix)
                ht.fill((nvtx,pwgt),                  'nvtx',   pcats,pfix)
                ht.fill((rho,pwgt),                   'rho',    pcats,pfix)
                ht.fill((met,pwgt),                   'met',    pcats,pfix)
                ht.fill((njets,pwgt),                 'njets',  pcats,pfix)
                ht.fill((nch,pwgt),                   'nch',    pcats,pfix) 
                #ht.fill((getattr(tree,'rfc_%d'%beamXangle),pwgt), 'rfc',         pcats,pfix)
                ht.fill((PFMultSumHF,pwgt),     'PFMultHF',    pcats,pfix)
                ht.fill((PFHtSumHF,pwgt),       'PFHtHF',      pcats,pfix)
                ht.fill((PFPzSumHF/1.e3,pwgt),  'PFPZHF',      pcats,pfix)
                ht.fill((n_extra_mu,pwgt), 'nextramu', pcats, pfix)
                if not isSignal:
                    ht.fill((tree.metfilters,pwgt), 'metbits', pcats,pfix)
                    for sd in ['HE','EE','EB']:
                        ht.fill((getattr(tree,'PFMultSum'+sd),pwgt),    'PFMult'+sd, pcats, pfix)
                        ht.fill((getattr(tree,'PFHtSum'+sd),pwgt),      'PFHt'+sd,   pcats, pfix)
                        ht.fill((getattr(tree,'PFPzSum'+sd)/1.e3,pwgt), 'PFPZ'+sd,   pcats, pfix)
                    for mp4 in extra_muons:
                        ht.fill((mp4.Pt(),pwgt), 'extramupt', pcats,pfix)
                        ht.fill((abs(mp4.Eta()),pwgt), 'extramueta', pcats,pfix)

                #proton counting and kinematics
                for irp,rpside in [(0,'pos'),(1,'neg')]:
                    ht.fill((len(protons[irp]),pwgt), 'ntk',    pcats,rpside+pfix)
                    for csi in protons[irp]:
                        ht.fill((csi,pwgt), 'csi', pcats,rpside+pfix)                        
                        if not near_protons : continue
                        if len(near_protons[irp])==0 : continue
                        csi_near=near_protons[irp][0]
                        ht.fill((csi,csi_near,pwgt), 'csi2d', pcats,rpside+pfix)

                #diproton kinematics
                if not pp: 
                    ht.fill((0,pwgt), 'ppcount', pcats,pfix)
                    continue        
                ht.fill((1,pwgt),                  'ppcount',   pcats, pfix)
                ht.fill((pp.M(),pwgt),             'mpp',       pcats, pfix)
                ht.fill((abs(pp.Rapidity()),pwgt), 'ypp',       pcats, pfix)                
                if near_pp:
                    ht.fill((2,pwgt),                                          'ppcount', pcats, pfix)
                    ht.fill((pp.M(),near_pp.M(),pwgt),                         'mpp2d',   pcats, pfix)
                    ht.fill((abs(pp.Rapidity()),abs(near_pp.Rapidity()),pwgt), 'ypp2d',   pcats, pfix)

                #the final variable
                ht.fill((mmass,pwgt),              'mmass',     pcats,pfix)                    

        #select events
        if not hasAHighPurSelection: continue
        if not goldenSel: continue

        #fill missing variables with 0.'s (signal)
        nVarsMissed=len(summaryVars)-len(goldenSel)
        if nVarsMissed>0: goldenSel += [0.]*nVarsMissed

        if isSignal:
                
            if isZ:
                #add a copy for ee
                eeGoldenSel=copy.copy(goldenSel)
                eeGoldenSel[0]=11*11
                eeGoldenSel[1]=goldenSel[1]*mcEff['eez'].Eval(boson.Pt())/nEntries
                selEvents.append(eeGoldenSel)
                
                #add a copy for mm
                mmGoldenSel=copy.copy(goldenSel)
                mmGoldenSel[0]=DIMUONS
                mmGoldenSel[1]=goldenSel[1]*mcEff['mmz'].Eval(boson.Pt())/nEntries
                selEvents.append(mmGoldenSel)
                
            if isA:
                
                #add a copy for the photon
                aGoldenSel=copy.copy(goldenSel)
                aGoldenSel[0]=22
                aGoldenSel[1]=goldenSel[1]*mcEff['a'].Eval(boson.Pt())/nEntries
                selEvents.append(aGoldenSel)
                    
        else:
            
            selEvents.append(goldenSel)
            

    #dump events for the mixing
    nSelRPData=sum([len(rpData[x]) for x in rpData])
    if nSelRPData>0:
        rpDataOut=outFileName.replace('.root','.pck')        
        print 'Saving',nSelRPData,'events for mixing in',rpDataOut
        with open(rpDataOut,'w') as cachefile:
            pickle.dump(rpData,cachefile, pickle.HIGHEST_PROTOCOL)        

    #if there was no mixing don't do anything else
    if not mixFile: return 

    #save results
    ht.writeToFile(outFileName)

    #dump events for fitting
    nSelEvents=len(selEvents)
    if nSelEvents>0:
        print 'Adding',nSelEvents,'selected events to',outFileName
        fOut=ROOT.TFile.Open(outFileName,'UPDATE')
        fOut.cd()
        t=ROOT.TNtuple('data','data',':'.join(summaryVars))
        for v in selEvents : t.Fill(array.array("f",v))
        t.Write()
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
    with open(opt.json,'r') as cachefile:
        samples=json.load(cachefile,  encoding='utf-8', object_pairs_hook=OrderedDict).items()
        for x in samples:
            task_dict[x[0]]=[]

    #group chunks matching the same name
    for file_path in os.listdir(opt.input):
        
        if opt.only and file_path!=os.path.basename(opt.only): continue

        file_name,ext=os.path.splitext(file_path)
        if ext != '.root' : continue
        
        #check if file tag is already in the list of samples to process
        lastTkn=file_name.rfind('_')
        tag=file_name[0:lastTkn]

        if not tag in task_dict: continue
        task_dict[tag].append( os.path.join(opt.input,file_path) )

    #parse json file with list of run/lumi sections
    with open(opt.RPout,'r') as cachefile:        
        runLumi=json.load(cachefile,  encoding='utf-8', object_pairs_hook=OrderedDict).items()
        runLumi={int(x[0]):x[1] for x in runLumi}

    #create the tasks and submit them
    import multiprocessing as MP
    pool = MP.Pool(opt.jobs)
    task_list=[]
    for x in task_dict.keys():
        
        isData = True if 'Data' in x else False
        if opt.step==0 and not isData : continue
        runLumiList=runLumi if isData else None
        for f in task_dict[x]:
            fOut='%s/Chunks/%s'%(opt.output,os.path.basename(f))
            task_list.append( (f,fOut,runLumiList,opt.mix) )

    pool.map(runExclusiveAnalysisPacked,task_list)


def main():
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--input',
                      dest='input',   
                      default='/eos/cms/store/cmst3/user/psilva/ExclusiveAna/Chunks',
                      help='input directory with the files [default: %default]')
    parser.add_option('--jobs',
                      dest='jobs', 
                      default=8,
                      type=int,
                      help='# of jobs to process in parallel the trees [default: %default]')
    parser.add_option('--json',
                      dest='json', 
                      default='pps_samples.json',
                      type='string',
                      help='json with the files to process')
    parser.add_option('--only',
                      dest='only', 
                      default=None,
                      type='string',
                      help='only this file')
    parser.add_option('--RPout',
                      dest='RPout', 
                      default='golden_noRP.json',
                      type='string',
                      help='json with the runs/lumi sections in which RP are out')
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
    parser.add_option('-o', '--output',
                      dest='output', 
                      default='analysis',
                      help='Output directory [default: %default]')
    (opt, args) = parser.parse_args()
    
    os.system('mkdir -p %s/Chunks' % opt.output)
    runAnalysisTasks(opt)
        

if __name__ == "__main__":
    sys.exit(main())
