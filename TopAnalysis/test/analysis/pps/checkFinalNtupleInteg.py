#!/usr/bin/env python

import sys
import os
import ROOT

def getEntries(url,tname):
    if not os.path.isfile(url): 
        raise NameError('No file for %s'%tname)

    f=ROOT.TFile.Open(url)
    if f.IsZombie():
        raise RuntimeError('Zombie file for %s'%tname)
    n=0
    #n=f.Get(tname).GetEntriesFast()
    f.Close()

    return n

def runAnaPacked(args):
    url,tag,step,stepDir,mix_file=args

    os.system('sh test/analysis/pps/wrapAnalysis.sh {step} {url}/{stepDir} {url}/Chunks {tag} {mix_file}'.format(step=step,
                                                                                                                 url=url,
                                                                                                                 tag=tag,
                                                                                                                 stepDir=stepDir,
                                                                                                                 mix_file=mix_file))


url=sys.argv[1]
runMode=int(sys.argv[2])
step=int(sys.argv[3])
stepDir='analysis' if step==1 else 'mixing'
mix_file=url+'/mixing/mixbank.pck' if step==1 else ''

toCheck=[]
nTot=0
for f in os.listdir(url+'/Chunks'):

    #if not 'Data' in f : continue
    if step==0 :
        if not ('SingleMuon' in f or 'MuonEG' in f or 'DoubleMuon' in f) :
            continue
            
    nTot+=1

    fullf=os.path.join(url+'/%s/Chunks'%stepDir,f)

    if step==0:
        if not os.path.isfile(fullf.replace('.root','.pck')): 
            print 'Missing in action',fullf
            toCheck.append(f)
    else:
        if not os.path.isfile(fullf) : 
            print 'Missing in action',fullf
            toCheck.append(f)
        else:
            try:
                nentries=getEntries(fullf,'data')                    
            except Exception as e:
                print 'Corrupted or empty file for',f,e
                toCheck.append(f)

if runMode==0:
    for x in toCheck:
        print x
elif runMode==1:
    print 'Running locally',len(toCheck),'/',nTot,'jobs'
    import multiprocessing as MP
    pool = MP.Pool(8)
    pool.map( runAnaPacked,[ (url,x,step,stepDir,mix_file) for x in toCheck ])
elif runMode==2:
    print 'Submitting',len(toCheck),'/',nTot,'jobs'
    cmssw_base=os.environ['CMSSW_BASE']
    with open('zxana_recover.sub','w') as cache:
        cache.write("executable  = %s/src/TopLJets2015/TopAnalysis/test/analysis/pps/wrapAnalysis.sh\n"%cmssw_base)
        cache.write("output      = zxana_recover.out\n")
        cache.write("error       = zxana_recover.err\n")
        cache.write("log         = zxana_recover.log\n")
        for x in toCheck:
            cache.write("arguments   = {step} {predout} {predin} {filein} {mix_file}\n".format(step=step,
                                                                                               predout=url+'/'+stepDir,
                                                                                               predin=url+'/Chunks',
                                                                                               filein=x,
                                                                                               mix_file=mix_file))
            cache.write("queue 1\n")
        os.system('condor_submit zxana_recover.sub')
