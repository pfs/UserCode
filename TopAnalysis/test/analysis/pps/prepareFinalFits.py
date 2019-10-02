import ROOT
import os
import generateBinnedWorkspace
import numpy as np
import sys
import pickle
       
MASSPOINTS=[780,800,840,900,960,1000,1020,1080,1140,1200,1260,1320,1380,1400,1440,1500,1560,1600,1620]

def generateDataCards(baseDir,l1pt,l2pt,bosonpt,categs,runGen):

    """generates datacards and returns the collection of directories with datacards to combine and fit"""

    dataCards={}

    os.system('mkdir -p %s'%baseDir)
    errorList=[]
    for m in MASSPOINTS:
        try:
            dataCards[m]={'Z':[],'gamma':[]}

            #ee/mumu cases
            for cat in [121,169]:
                jobDir='%s/m%d/cat%d'%(baseDir,m,cat)
                if runGen:
                    for xangle in generateBinnedWorkspace.VALIDLHCXANGLES:
                        generateBinnedWorkspace.main(['--xangle', '%d'%xangle,
                                                      '--sig',    'Z_m_X_{0}_xangle_{1}_2017_preTS2_opt_v1_simu_reco.root'.format(m,xangle),
                                                      '--presel', 'cat=={0} && l1pt>{1} && l2pt>{2} && bosonpt>{3}'.format(cat,l1pt,l2pt,bosonpt),
                                                      '--categs', categs,
                                                      '-o',       jobDir
                                                  ])
                dataCards[m]['Z'].append( jobDir )
                
            #photon case
            continue
            jobDir='%s/m%d/cat22'%(baseDir,m)
            if runGen:
                for xangle in generateBinnedWorkspace.VALIDLHCXANGLES:
                    generateBinnedWorkspace.main(['--xangle', '%d'%xangle,
                                                  '--lumi',   '2288',
                                                  '--sig',    'gamma_m_X_{0}_xangle_{1}_2017_preTS2_opt_v1_simu_reco.root'.format(m,xangle),
                                                  '--presel', 'cat==22 && bosonpt>{0} && bosonpt<200'.format(bosonpt),
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

def runFits(inputDirs):


    results=[]

    #combine all available datacards
    dc=[]
    for d in inputDirs:
        dc+=[os.path.join(d,f) for f in os.listdir(d) if '.dat' in f and not '_1.dat' in f]
    combStr=' '.join( ['cat%d=%s'%(i,dc[i]) for i in range(len(dc))] )   
    os.system('combineCards.py {0} > combined_card.dat\n'.format(combStr))
    os.system('text2workspace.py combined_card.dat -o workspace.root\n')    
    
    #run limits and read quantiles
    os.system('combine workspace.root -M AsymptoticLimits -t -1\n')
    fIn=ROOT.TFile.Open('higgsCombineTest.AsymptoticLimits.mH120.root')
    limit=fIn.Get('limit')
    for ientry in range(5):
        limit.GetEntry(ientry)
        results.append(limit.limit)

    #run significance
    os.system('combine workspace.root -M Significance -t -1 --expectSignal=1\n')
    fIn=ROOT.TFile.Open('higgsCombineTest.Significance.mH120.root')
    limit=fIn.Get('limit')
    limit.GetEntry(0)
    results.append(limit.limit)

    return results
    

runGen=True

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
        iresults += runFits(dataCards[m][proc])
        results.append(iresults)

print results

with open('{0}/results.pck'.format(baseDir),'w') as cache:
    pickle.dump(results,cache,pickle.HIGHEST_PROTOCOL) 
        
