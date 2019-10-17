import ROOT
import os
import generateBinnedWorkspace
import numpy as np
import sys
import pickle
       
MASSPOINTS=[780,800,840,900,960,1000,1020,1080,1140,1200,1260,1320,1380,1400,1440,1500,1560,1600,1620]

MASSPOINTS=[800,900,1000,1080,1200,1320,1400,1500,1600]

def generateDataCards(baseDir,l1pt,l2pt,bosonpt,categs,runGen):

    """generates datacards and returns the collection of directories with datacards to combine and fit"""

    dataCards={}

    os.system('mkdir -p %s'%baseDir)
    errorList=[]
    for m in MASSPOINTS:
        dataCards[m]={}
        try:
            #ee/mumu cases
            dataCards[m]['Z']={'out':'%s/m%d/Z'%(baseDir,m),'datacards':[]}
            for cat in [121,169]:
                jobDir='%s/m%d/cat%d'%(baseDir,m,cat)
                if runGen:
                    for xangle in generateBinnedWorkspace.VALIDLHCXANGLES:
                        generateBinnedWorkspace.main(['--xangle', '%d'%xangle,
                                                      '--sig',    'Z_m_X_{0}_xangle_{1}_2017_preTS2_opt_v1_simu_reco.root'.format(m,xangle),
                                                      '--presel', 'cat=={0} && l1pt>{1} && l2pt>{2} && bosonpt>{3}'.format(cat,l1pt,l2pt,bosonpt),
                                                      '--categs', categs,
                                                      '--csiacc', '{cmssw}/src/TopLJets2015/TopAnalysis/test/analysis/pps/csiaccparam_Z.pck'.format(cmssw=os.environ['CMSSW_BASE']),
                                                      '-o',       jobDir
                                                  ])
                dataCards[m]['Z']['datacards'].append( jobDir )
                
            #photon case
            continue
            dataCards[m]['gamma']={'out':'%s/m%d/gamma'%(baseDir,m),'datacards':[]}
            jobDir='%s/m%d/cat22'%(baseDir,m)
            if runGen:
                for xangle in generateBinnedWorkspace.VALIDLHCXANGLES:
                    generateBinnedWorkspace.main(['--xangle', '%d'%xangle,
                                                  '--lumi',   '2288',
                                                  '--sig',    'gamma_m_X_{0}_xangle_{1}_2017_preTS2_opt_v1_simu_reco.root'.format(m,xangle),
                                                  '--presel', 'cat==22 && bosonpt>{0} && bosonpt<200'.format(bosonpt),
                                                  '--csiacc', '{cmssw}/src/TopLJets2015/TopAnalysis/test/analysis/pps/csiaccparam_gamma.pck'.format(cmssw=os.environ['CMSSW_BASE']),
                                                  '--categs', categs,
                                                  '-o',       jobDir
                                              ])
                dataCards[m]['gamma'].append( jobDir )
        except Exception as e:
            errorList.append( (m,e) )

    if len(errorList)>0:
        for m,e in errorList:
            print m
            print e
        raw_input()
        

    return dataCards

def runFits(fit_desc):

    results=[]

    workDir=fit_desc['out']
    inputDirs=fit_desc['datacards']

    #dump a script to steer combine
    os.system('mkdir -p %s'%workDir)
    steerScript=os.path.join(workDir,'runstat.sh')
    with open(steerScript,'w') as sh: 

        sh.write('#!/bin/bash\n')

        #combine individual cards and create workspace
        dc=[]
        for d in inputDirs:
            dc+=[os.path.join(os.path.abspath(d),f) for f in os.listdir(d) if '.dat' in f and f!='combined_card.dat']
        combStr=' '.join( ['cat%d=%s'%(i,dc[i]) for i in range(len(dc))] )   
        sh.write('combineCards.py {0} > combined_card.dat\n'.format(combStr))
        sh.write('text2workspace.py combined_card.dat -o workspace.root\n')    
    
        #run expected limits
        sh.write('combine workspace.root -M AsymptoticLimits -t -1\n')
        sh.write('combine workspace.root -M Significance -t -1 --expectSignal=1\n')
    
    os.system('cd {0} && sh runstat.sh && cd -'.format(workDir))
    
    #parse results
    fIn=ROOT.TFile.Open(os.path.join(workDir,'higgsCombineTest.AsymptoticLimits.mH120.root'))
    limit=fIn.Get('limit')
    for ientry in range(5):
        limit.GetEntry(ientry)
        results.append(limit.limit)
    fIn=ROOT.TFile.Open(os.path.join(workDir,'higgsCombineTest.Significance.mH120.root'))
    limit=fIn.Get('limit')
    limit.GetEntry(0)
    results.append(limit.limit)

    return results
    

runGen=False #True

#baseDir='analysis/nicolafit'
#l1pt=30
#l2pt=20
#bosonpt=50

baseDir='analysis/finalfit'
l1pt=30
l2pt=20
bosonpt=50
categs='nvtx<20,nvtx>=20'

dataCards=generateDataCards(baseDir,l1pt,l2pt,bosonpt,categs,runGen)

results=[]
for m in dataCards:
    for proc in dataCards[m]:
        iresults=[m,22 if proc=='gamma' else 23]
        if len(dataCards[m][proc]['datacards'])==0 : continue
        iresults += runFits(dataCards[m][proc])
        results.append(iresults)

with open('{0}/results.pck'.format(baseDir),'w') as cache:
    pickle.dump(results,cache,pickle.HIGHEST_PROTOCOL) 
        
