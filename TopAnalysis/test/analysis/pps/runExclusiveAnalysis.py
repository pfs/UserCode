
#! /usr/bin/env python
import ROOT
import optparse
import json
import sys
import os
import numpy
import array
import json
import random
import pickle
from collections import OrderedDict

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

def getTracksPerRomanPot(tree):

    """loops over the availabe tracks in the event and groups them by roman pot id"""

    tkPos=[]
    tkNeg=[]
    try:
        for itk in xrange(0,tree.nRPtk):
            rpid=tree.RPid[itk]
            csi=tree.RPfarcsi[itk]
            if rpid==23  : tkPos.append(csi)
            if rpid==123 : tkNeg.append(csi)
    except:
        pass

    return (tkPos,tkNeg)

def buildDiproton(rptks,sqrts=13000.):

    """build a diproton system from to tracks in the arms of the roman pots"""

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


def createEventMixBank(fileList,outdir,runLumiList):
    
    """build a list of events for the mixing"""

    print '... @ createEventMixBank'

    rpData={}
    for tag in fileList:        
        rpData[tag]=[]
        tree=ROOT.TChain('tree')
        for f in fileList[tag]: tree.AddFile(f)        
        nEntries=tree.GetEntries()

        print 'Starting',tag,'with',nEntries
        for i in xrange(0,nEntries):
            if i%1000==0 : sys.stdout.write('\r [ %d/100 ] done' %(int(float(100.*i)/float(nEntries))))
            tree.GetEntry(i)
            if not isValidRunLumi(tree.run,tree.lumi,runLumiList): continue
            tkPos,tkNeg=getTracksPerRomanPot(tree)
            rpData[tag].append( (tree.beamXangle,tkPos,tkNeg) )

        random.shuffle(rpData[tag])

    with open(os.path.join(outdir,'evmix.pck'),'w') as cachefile:
        pickle.dump(rpData,cachefile, pickle.HIGHEST_PROTOCOL)


def runExclusiveAnalysis(inFile,outFileName,runLumiList,mixFile,mixSel):
    
    """event loop"""

    tree=ROOT.TChain('tree')
    tree.AddFile(inFile)
    nEntries=tree.GetEntries()  

    with open(mixFile,'r') as cachefile:
        mixedRP=pickle.load(cachefile)

    #start histograms
    histos={}
    histos['nvtx']        = {'inc':ROOT.TH1F('nvtx',';Vertex multiplicity;Events',50,0,100)}
    histos['mll']         = {'inc':ROOT.TH1F('mll',';Dilepton invariant mass [GeV];Events',50,20,250)}
    histos['ptll']        = {'inc':ROOT.TH1F('ptll',';Dilepton transverse momentum [GeV];Events',50,0,250)}
    histos['xangle']      = {'inc':ROOT.TH1F('xangle',';LHC crossing angle [#murad];Events',4,120,160)}
    histos['xanglevsnvtx']= {'inc':ROOT.TH2F('xanglevsnvtx',';LHC crossing angle [#murad];Vertex multiplicity;Events',4,120,160,5,0,100)}
    histos['mpp']         = {'inc':ROOT.TH1F('mpp',';Di-proton invariant mass [GeV];Events',100,0,2000)}
    histos['mmass']       = {'inc':ROOT.TH1F('mmass',';Missing mass [GeV];Events',100,0,2000)}
    histos['mmassvsnvtx'] = {'inc':ROOT.TH2F('mmassvsnvtx',';Missing mass [GeV];Vertex multiplicity;Events',100,0,2000,5,0,100)}
    histos['ntk']         = {'inc':ROOT.TH1F('ntk',';Track multiplicity;Events',5,0,5)}
    histos['csi']         = {'inc':ROOT.TH1F('csi',';#xi;Events',50,0,0.3)}
    for name in histos:
        for cat in histos[name]:
            histos[name][cat].SetDirectory(0)
            histos[name][cat].Sumw2()
            
    def fillHisto(val,key):
        name,cat=key
        if not cat in histos[name]:
            histos[name][cat]=histos[name]['inc'].Clone('%s_%s'%key)
            histos[name][cat].SetDirectory(0)
            histos[name][cat].Reset('ICE')
        histos[name][cat].Fill(*val)
        

    print '....analysing',nEntries,'in',inFile,', with output @',outFileName
    print '    events mixed with',mixSel if mixSel else 'lumi-weighted data','from',mixFile
    isData=True if 'Data' in inFile else False

    #loop over events in the tree and fill histos
    irand=0
    for i in xrange(0,nEntries):

        tree.GetEntry(i)

        if not isValidRunLumi(tree.run,tree.lumi,runLumiList): continue
        
        if i%1000==0 : sys.stdout.write('\r [ %d/100 ] done' %(int(float(100.*i)/float(nEntries))))
        
        if tree.isSS : continue
        dilId=abs(int(tree.l1id*tree.l2id))
        if dilId==11*11 : dilcat='ee'
        if dilId==11*13 : dilcat='em'
        if dilId==13*13 : dilcat='mm'

        wgt=tree.evwgt
        nvtx=tree.nvtx

        beamXangle=tree.beamXangle
        rptks   = getTracksPerRomanPot(tree)
        eraTag  = mixSel if mixSel else getRandomEra()
        mixedEv = random.choice( mixedRP[eraTag] )
        mixed_beamXangle=mixedEv[0]
        mixed_rptks=(mixedEv[1],mixedEv[2])

        #filter based on the beam crossing angle and high purity of the event
        passBeamXangle = True if beamXangle in [120,130,140] else False
        highPur=True if len(rptks[0])==1 and len(rptks[1])==1 else False        
        mix_passBeamXangle = True if mixed_beamXangle in [120,130,140] else False
        mix_highPur=True if len(mixed_rptks[0])==1 and len(mixed_rptks[1])==1 else False
        
        ll=ROOT.TLorentzVector(0,0,0,0)
        ll.SetPtEtaPhiM(tree.llpt,tree.lleta,tree.llphi,tree.mll)

        pp=buildDiproton(rptks)
        mixed_pp=buildDiproton(mixed_rptks)

        mmass,mixed_mmass=0.,0.
        if pp:        mmass=(ll-pp).M()
        if mixed_pp : mixed_mmass=(ll-mixed_pp).M()

        cats=[dilcat]
        if tree.isZ                 : cats.append(dilcat+'Z')
        if ll.Pt()>50               : cats.append(dilcat+'hpt')
        if tree.isZ and ll.Pt()>50  : cats.append(dilcat+'hptZ')
        #if highPur                  : cats += [x+'highpur' for x in cats]

        blind=False
        if tree.isZ and ll.Pt()>50:
            if pp and (mmass>1000 and mmass<2000):
                blind=True

        #fill histograms
        for cat in cats:
            
            fillHisto(val=(nvtx,wgt),             key=('nvtx',cat))
            fillHisto(val=(ll.M(),wgt),           key=('mll',cat))
            fillHisto(val=(ll.Pt(),wgt),          key=('ptll',cat)) 
       
            fillHisto(val=(beamXangle,wgt),            key=('xangle',cat))
            fillHisto(val=(beamXangle,nvtx,wgt),       key=('xanglevsnvtx',cat))
            fillHisto(val=(mixed_beamXangle,wgt),      key=('xangle',cat+'_mix'))
            fillHisto(val=(mixed_beamXangle,nvtx,wgt), key=('xanglevsnvtx',cat+'_mix'))

            if pp and passBeamXangle:
                if not (blind and 'hpt' in cat):
                    fillHisto(val=(pp.M(),wgt),     key=('mpp',cat))
                    fillHisto(val=(mmass,wgt),      key=('mmass',cat))
                    fillHisto(val=(mmass,nvtx,wgt), key=('mmassvsnvtx',cat))
                    
                    for irp,rpside in [(0,'pos'),(1,'neg')]:
                        fillHisto(val=(len(rptks[irp]),wgt), key=('ntk',cat+'_'+rpside))
                        for csi in rptks[irp]:
                            fillHisto(val=(csi,wgt), key=('csi',cat+'_'+rpside))

            if mixed_pp and mix_passBeamXangle:
                fillHisto(val=(mixed_pp.M(),wgt),     key=('mpp',cat+'_mix'))
                fillHisto(val=(mixed_mmass,wgt),      key=('mmass',cat+'_mix'))
                fillHisto(val=(mixed_mmass,nvtx,wgt), key=('mmassvsnvtx',cat+'_mix'))

                for irp,rpside in [(0,'pos'),(1,'neg')]:
                    fillHisto(val=(len(rptks[irp]),wgt), key=('ntk',cat+'_'+rpside))
                    for csi in rptks[irp]:
                        fillHisto(val=(csi,wgt), key=('csi',cat+'_'+rpside))

                    fillHisto(val=(len(mixed_rptks[irp]),wgt), key=('ntk',cat+'_mix'+rpside))
                    for csi in mixed_rptks[irp]:
                        fillHisto(val=(csi,wgt), key=('csi',cat+'_mix'+rpside))

    #save results
    fOut=ROOT.TFile.Open(outFileName,'RECREATE')
    fOut.cd()
    for name in histos:
        for cat in histos[name]:
            if histos[name][cat].GetEntries()==0 : 
                histos[name][cat].Delete()
                continue
            histos[name][cat].SetDirectory(fOut)
            histos[name][cat].Write()
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

    #prepare event mixing 
    if opt.step==0:
        dataPerEra={}
        for x in task_dict.keys():
            if 'MC' in x : continue
            if not 'MuonEG' in x : continue
            era=x.split('_')[1]
            if not era in dataPerEra: dataPerEra[era]=[]
            dataPerEra[era]+=task_dict[x]
        createEventMixBank(dataPerEra,opt.output,runLumi)

    #run the analysis
    elif opt.step==1:

        import multiprocessing as MP
        pool = MP.Pool(opt.jobs)
        task_list=[]
        for x in task_dict.keys():
            runLumiList=None
            mixSel=None
            if 'Data' in x: 
                runLumiList=runLumi        
                mixSel=x.split('_')[1]
            for f in task_dict[x]:
                fOut='%s/Chunks/%s'%(opt.output,os.path.basename(f))
                task_list.append( (f,fOut,runLumiList,opt.mix,mixSel) )

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
                      default='analysis/evmix.pck',
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
