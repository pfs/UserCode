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


inurl=sys.argv[1]
outurl=sys.argv[2]
runMode=int(sys.argv[3])
step=int(sys.argv[4])
mix_file=''
if len(sys.argv)>5:
    mix_file=sys.argv[5]
stepDir='analysis' if step==1 else 'mixing'

toCheck=[]
nTot=0
for f in os.listdir(inurl):

    #if not 'Data' in f : continue
    if step==0 :
        if not ('SingleMuon' in f or 'MuonEG' in f or 'DoubleMuon' in f) :
            continue
            
    nTot+=1

    fullf=os.path.join(outurl,f)

    if step==0:
        fullf=fullf.replace('.root','.pck')
        if not os.path.isfile(fullf): 
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


print '-'*50
print 'These outputs are missing'
print [x.replace('.root','') for x in toCheck]
print '-'*50


outurl=outurl.replace('/Chunks','')
if runMode==1:
    print 'Running locally',len(toCheck),'/',nTot,'jobs'
    os.system('sh test/analysis/pps/wrapAnalysis.sh {cmssw} {step} {outurl} {inurl} {tag} {mix_file}'.format(cmssw=os.environ['CMSSW_BASE'],
                                                                                                             step=step,
                                                                                                             outurl=outurl,
                                                                                                             inurl=inurl,
                                                                                                             tag=','.join(toCheck),
                                                                                                             stepDir=stepDir,
                                                                                                             mix_file=mix_file))

elif runMode==2:
    print 'Submitting',len(toCheck),'/',nTot,'jobs'
    cmssw_base=os.environ['CMSSW_BASE']
    with open('zxana_recover.sub','w') as cache:
        cache.write("executable  = %s/src/TopLJets2015/TopAnalysis/test/analysis/pps/wrapAnalysis.sh\n"%cmssw_base)
        cache.write("output      = zxana_recover.out\n")
        cache.write("error       = zxana_recover.err\n")
        cache.write("log         = zxana_recover.log\n")
        cache.write("+JobFlavour = \"tomorrow\"\n")
        cache.write("request_cpus = 4\n")
        for x in toCheck:
            cache.write("arguments   = {cmssw} {step} {predout} {predin} {filein} {mix_file}\n".format(cmssw=os.environ['CMSSW_BASE'],
                                                                                                       step=step,
                                                                                                       predout=outurl,
                                                                                                       predin=inurl,
                                                                                                       filein=x,
                                                                                                       mix_file=mix_file))
            cache.write("queue 1\n")
        os.system('condor_submit zxana_recover.sub')
