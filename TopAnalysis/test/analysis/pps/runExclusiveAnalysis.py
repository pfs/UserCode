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
        if rpid==23  : tkPos.append(csi)
        if rpid==123 : tkNeg.append(csi)
        
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

    #identify data-taking era
    era=None
    isData=True if 'Data' in inFile else False
    if isData:
        era=os.path.basename(inFile).split('_')[1]
    
    #check if it is signal and load 
    isSignal=True if 'MC13TeV_2017_PPZX_' in inFile else False
    signalPt=[]
    mcEff={}
    if isSignal: 
        signalPt=[float(x) for x in re.findall(r'\d+', os.path.basename(inFile) )[2:]]
        for ch in ['ee','mm']:
            effIn=ROOT.TFile.Open('$CMSSW_BASE/src/TopLJets2015/TopAnalysis/plots/effsummary_%sz_ptboson.root'%ch)
            mcEff[ch]=effIn.Get('gen%sz2trec_ptboson_ZH#rightarrowllbb_eff'%ch)
            effIn.Close()

    #filter events to mix according to tag if needed
    mixedRP=None
    try:
        with open(mixFile,'r') as cachefile:
            mixedRP=pickle.load(cachefile)
    except:
        pass

    
    #start histograms
    ht=HistoTool()
    ht.add(ROOT.TH1F('nvtx',';Vertex multiplicity;Events',50,0,100))
    ht.add(ROOT.TH1F('rho',';Fastjet #rho;Events',50,0,30))
    ht.add(ROOT.TH1F('met',';Missing transverse energy [GeV];Events',50,0,200))
    ht.add(ROOT.TH1F('metbits',';MET filters;Events',124,0,124))
    ht.add(ROOT.TH1F('njets',';Jet multiplicity;Events',5,0,5))
    ht.add(ROOT.TH1F('nch', ';Charged particle multiplicity;Events',50,0,50))
    ht.add(ROOT.TH1F('acopl',';A=1-|#Delta#phi|/#pi;Events',50,0,1))
    ht.add(ROOT.TH1F('l1eta',';Lepton pseudo-rapidiy;Events',50,0,2.5))
    ht.add(ROOT.TH1F('l1pt',';Lepton transverse momentum [GeV];Events',50,0,250))
    ht.add(ROOT.TH1F('l2eta',';Lepton pseudo-rapidiy;Events',50,0,2.5))
    ht.add(ROOT.TH1F('l2pt',';Lepton transverse momentum [GeV];Events',50,0,250))
    ht.add(ROOT.TH1F('mll',';Dilepton invariant mass [GeV];Events',50,20,250))
    ht.add(ROOT.TH1F('yll',';Dilepton rapidity;Events',50,0,5))
    ht.add(ROOT.TH1F('ptll',';Dilepton transverse momentum [GeV];Events',50,0,250))
    ht.add(ROOT.TH1F('minenfwd',';min(E_{+},E_{-}) [GeV];Events',20,0,300))
    ht.add(ROOT.TH1F('maxenfwd',';max(E_{+},E_{-}) [GeV];Events',20,0,300))
    ht.add(ROOT.TH1F('deltaenfwd',';|E_{+}-E_{-}| [GeV];Events',20,0,300))
    ht.add(ROOT.TH1F('sgny',';y x sgn(LRG) ;Events',20,-3,3))
    ht.add(ROOT.TH1F('nextramu',';Additional muons ;Events',10,0,10))
    ht.add(ROOT.TH1F('extramupt',';Additional muon p_{T} [GeV] ;Events',10,0,50))
    ht.add(ROOT.TH1F('extramueta',';Additional muon pseudo-rapidty ;Events',10,0,2.5))
    ht.add(ROOT.TH1F('xangle',';LHC crossing angle [#murad];Events',4,120,160))
    ht.add(ROOT.TH1F('mpp',';Di-proton invariant mass [GeV];Events',50,0,3000))
    ht.add(ROOT.TH1F('ypp',';Di-proton rapidity;Events',50,0,2))
    ht.add(ROOT.TH1F('mmass',';Missing mass [GeV];Events',50,0,3000))
    ht.add(ROOT.TH1F('ntk',';Track multiplicity;Events',5,0,5))
    ht.add(ROOT.TH1F('csi',';#xi;Events',50,0,0.3))

    #start analysis
    tree=ROOT.TChain('analysis/data' if isSignal else 'tree')
    tree.AddFile(inFile)
    nEntries=tree.GetEntries()          
    print '....analysing',nEntries,'in',inFile,', with output @',outFileName
    if mixedRP : print '    events mixed with',mixFile

    #loop over events
    rpData={era:[]}
    selEvents=[]
    summaryVars='cat:wgt:nvtx:nch:xangle:l1pt:l1eta:l2pt:l2eta:acopl:bosonpt:mpp:mmiss:mpp2:mmiss2'.split(':')
    for i in xrange(0,nEntries):

        tree.GetEntry(i)
        
        if i%1000==0 : sys.stdout.write('\r [ %d/100 ] done' %(int(float(100.*i)/float(nEntries))))
    
        #base event selection
        if tree.evcat==11*11   : evcat='ee'
        elif tree.evcat==11*13 : evcat='em'
        elif tree.evcat==13*13 : evcat='mm'
        else : continue
        if tree.isSS : continue        
        if isData and not isValidRunLumi(tree.run,tree.lumi,runLumiList): continue

        wgt=tree.evwgt
        nvtx=tree.nvtx
        nch=tree.nch
        rho=tree.rho
        met=tree.met_pt
        njets=0 if isSignal else tree.nj 

        #acoplanarity
        acopl=1.0-abs(ROOT.TVector2.Phi_mpi_pi(tree.l1phi-tree.l2phi))/ROOT.TMath.Pi()

        l1p4=ROOT.TLorentzVector(0,0,0,0)
        l1p4.SetPtEtaPhiM(tree.l1pt,tree.l1eta,tree.l1phi,tree.ml1)
        l2p4=ROOT.TLorentzVector(0,0,0,0)
        l2p4.SetPtEtaPhiM(tree.l2pt,tree.l2eta,tree.l2phi,tree.ml2)
        if l1p4.Pt()<l2p4.Pt(): l1p4,l2p4=l2p4,l1p4
        if abs(l1p4.Eta())>2.1 : continue

        #boson kinematics
        boson=ROOT.TLorentzVector(0,0,0,0)
        boson.SetPtEtaPhiM(tree.bosonpt,tree.bosoneta,tree.bosonphi,tree.mboson)
        isZ=tree.isZ
        isHighPt=(boson.Pt()>50)

        #possible diffractive-sensitive variabes
        en_posRG=tree.jsumposhfen
        en_negRG=tree.jsumneghfen
        extra_muons=[]
        for im in range(tree.nrawmu):
            mup4=ROOT.TLorentzVector(0,0,0,0)
            mup4.SetPtEtaPhiM(tree.rawmu_pt[im],tree.rawmu_eta[im]/10.,tree.rawmu_phi[im]/10.,0.105)
            if mup4.DeltaR(l1p4)<0.05 : continue
            if mup4.DeltaR(l2p4)<0.05 : continue
            extra_muons.append( ROOT.TLorentzVector(mup4) )
            
        #proton tracks (standard and mixed)
        rptks, near_rptks = None, None
        beamXangle       = tree.beamXangle
        if isSignal: 
            beamXangle=signalPt[0]
            rptks = getTracksPerRomanPot(tree)
            far_rptks = getTracksPerRomanPot(tree,False,False)
            rptks = ( [x/0.0964 for x in rptks[0]], [x/0.06159 for x in rptks[1]] )
            far_rptks = ( [x/0.0964 for x in far_rptks[0]], [x/0.06159 for x in far_rptks[1]] )
        if isData:
            rptks = getTracksPerRomanPot(tree) 
            far_rptks = getTracksPerRomanPot(tree,False,False) 

        mixed_beamXangle = None
        mixed_rptks      = None
        mixed_1rptk      = None
        if not isData:
            evEra=getRandomEra()
        else:
            evEra=era
            if not mixedRP:
                if evcat=='em' and tree.bosonpt>50:
                    rpData[era].append( MixedEvent(beamXangle,
                                                   nvtx,
                                                   rho,
                                                   rptks,
                                                   far_rptks,
                                                   extra_muons,
                                                   en_posRG,
                                                   en_negRG
                                                   ) )
                continue
        try:
            mixedEv=random.choice( mixedRP[evEra] )
            mixed_beamXangle=mixedEv.beamXangle
            mixed_rptks=mixedEv.far_rptks
            if rptks and mixed_rptks:
                if random.random()<0.5:
                    mixed_1rptk = (mixed_rptks[0],rptks[1])
                else:
                    mixed_1rptk = (rptks[0],mixed_rptks[1])
        except :
            pass

        #fill histograms
        goldenSel=None
        mon_tasks=[]
        if isData:
            mon_tasks.append( (rptks,beamXangle,'') )
            mon_tasks.append( (mixed_1rptk,mixed_beamXangle,'_mix1') )
            mon_tasks.append( (mixed_rptks,mixed_beamXangle,'_mix2') )
        else:
            if isSignal:
                #embed pileup to signal
                tksPos=mixed_rptks[0]+rptks[0]
                shuffle(tksPos)
                tksNeg=mixed_rptks[1]+rptks[1]
                shuffle(tksNeg)
                mixed_rptks=(tksPos,tksNeg)
                mixed_beamXangle=beamXangle
                mon_tasks.append( (rptks,beamXangle,'nopu') )
            mon_tasks.append( (mixed_rptks,mixed_beamXangle,'') )

        for protons, xangle, pfix in mon_tasks:

            #no calibration, not worth it...
            if not xangle in VALIDLHCXANGLES: continue
    
            #high purity selection for proton tracks
            highPur = True if protons and len(protons[0])==1 and len(protons[1])==1 else False  
            if highPur:
                if protons[0][0]<0.05 or protons[0][0]>0.20 : highPur=False
                if protons[1][0]<0.05 or protons[1][0]>0.20 : highPur=False

            #the famous cut
            noExtraMu = True
            if len(extra_muons)>1: noExtraMu=False

            #check if Z and di-proton combination is consistent with elastic scattering
            pp=buildDiproton(protons)
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
            if isZ              : cats.append(evcat+'Z')
            if isHighPt         : cats.append(evcat+'hpt')
            if isZ and isHighPt : cats.append(evcat+'hptZ')
            if isElasticLike    :
                ppCats=[c+'elpp' for c in cats]
                if highPur: 
                    ppCats+=[c+'elpphighPur' for c in cats]
                if isZ:
                    cats += ppCats + [c+'%d'%xangle for c in ppCats]
            if noExtraMu:
                extrmucats=[c+'noextramu' for c in cats]
                cats += extrmucats
            if nvtx<2 and len(protons[0])+len(protons[1])==1:
                if (en_negRG==0 and en_posRG>0) or (en_negRG>0 and en_posRG==0):
                    diffCats=[c+'diff' for c in cats]
                    cats += diffCats

            if len(pfix)!=0:
                cats=[c for c in cats if 'elpp' in c]
                
            if (isData and 'mix' in pfix) or (isSignal and pfix==''):
                if isZ and isElasticLike and highPur:
                    if goldenSel:
                        goldenSel += [pp.M(),mmass]
                    else:
                        goldenSel=[tree.evcat,wgt,nvtx,nch,xangle,l1p4.Pt(),l1p4.Eta(),l2p4.Pt(),l2p4.Eta(),acopl,boson.Pt(),pp.M(),mmass]

            #final plots (for signal correct wgt by efficiency curve and duplicate for mm channel)
            finalPlots=[[wgt,cats]]
            if isSignal:
                finalPlots=[ [wgt*mcEff['ee'].Eval(boson.Pt())/nEntries, cats],
                             [wgt*mcEff['mm'].Eval(boson.Pt())/nEntries, [c.replace(evcat,'mm') for c in cats if c[0:2]=='ee']] ]

            for pwgt,pcats in finalPlots:
                
                ht.fill((nvtx,pwgt),                  'nvtx',   pcats,pfix)
                ht.fill((rho,pwgt),                   'rho',   pcats,pfix)
                ht.fill((met,pwgt),                   'met',   pcats,pfix)
                ht.fill((tree.metfilters,pwgt),       'metbits',   pcats,pfix)
                ht.fill((njets,pwgt),                 'njets',  pcats,pfix)
                ht.fill((nch,pwgt),                   'nch',    pcats,pfix)                                
                ht.fill((l1p4.Pt(),pwgt),             'l1pt',   pcats,pfix)
                ht.fill((l2p4.Pt(),pwgt),             'l2pt',   pcats,pfix)
                ht.fill((abs(l1p4.Eta()),pwgt),       'l1eta',  pcats,pfix)
                ht.fill((abs(l2p4.Eta()),pwgt),       'l2eta',  pcats,pfix)
                ht.fill((acopl,pwgt),                 'acopl',  pcats,pfix) 
                ht.fill((boson.M(),pwgt),             'mll',    pcats,pfix)
                ht.fill((abs(boson.Rapidity()),pwgt), 'yll',    pcats,pfix)
                ht.fill((boson.Pt(),pwgt),            'ptll',   pcats,pfix)
                ht.fill((xangle,pwgt),                'xangle',    pcats,pfix)
                ht.fill((min(en_posRG,en_negRG),pwgt),'minenfwd',    pcats,pfix)
                ht.fill((max(en_posRG,en_negRG),pwgt),'maxenfwd',    pcats,pfix)                
                ht.fill((abs(en_posRG-en_negRG),pwgt),'deltaenfwd',    pcats,pfix)                
                sgnY=boson.Rapidity() if en_posRG<en_negRG else -boson.Rapidity()
                ht.fill((sgnY,pwgt),   'sgny', pcats,pfix)
                ht.fill((len(extra_muons),pwgt), 'nextramu', pcats, pfix)
                for mp4 in extra_muons:
                    ht.fill((mp4.Pt(),pwgt), 'extramupt', pcats,pfix)
                    ht.fill((abs(mp4.Eta()),pwgt), 'extramueta', pcats,pfix)

                for irp,rpside in [(0,'pos'),(1,'neg')]:
                    ht.fill((len(protons[irp]),pwgt), 'ntk',    pcats,rpside+pfix)
                if not pp: continue        
                ht.fill((pp.M(),pwgt),             'mpp',       pcats,pfix)
                ht.fill((abs(pp.Rapidity()),pwgt), 'ypp',       pcats,pfix)
                ht.fill((mmass,pwgt),              'mmass',     pcats,pfix)
                for irp,rpside in [(0,'pos'),(1,'neg')]:
                    ht.fill((len(protons[irp]),pwgt), 'ntk',    pcats,rpside+pfix)
                    for csi in protons[irp]:
                        ht.fill((csi,pwgt), 'csi',             pcats,rpside+pfix)

        #select events
        if goldenSel:
            nVarsMissed=len(summaryVars)-len(goldenSel)
            if nVarsMissed>0: goldenSel += [0.]*nVarsMissed

            if isSignal:
                                
                #add a copy for ee
                eeGoldenSel=copy.copy(goldenSel)
                eeGoldenSel[0]=11*11
                eeGoldenSel[1]=goldenSel[1]*mcEff['ee'].Eval(boson.Pt())/nEntries
                selEvents.append(eeGoldenSel)
                
                #add a copy for mm
                mmGoldenSel=copy.copy(goldenSel)
                mmGoldenSel[0]=13*13
                mmGoldenSel[1]=goldenSel[1]*mcEff['mm'].Eval(boson.Pt())/nEntries
                selEvents.append(mmGoldenSel)

            else:

                selEvents.append(goldenSel)
            

    #dump events for the mixing
    nSelRPData=len(rpData[era])
    if nSelRPData:
        rpDataOut=outFileName.replace('.root','.pck')
        print 'Saving',nSelRPData,'events for mixing in',rpDataOut
        with open(rpDataOut,'w') as cachefile:
            pickle.dump(rpData,cachefile, pickle.HIGHEST_PROTOCOL)        

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
            if opt.only and not opt.only in x[0]: continue
            task_dict[x[0]]=[]

    #group chunks matching the same name
    for file_path in os.listdir(opt.input):
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
                      help='only this process')
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
