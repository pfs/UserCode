import ROOT
import sys
import os
import pickle
import re
from prepareOptimScanCards import OPTIMLIST

def readLimitsFrom(url,getObs=False):

    """parses the r95 limits from the tree"""

    fIn=ROOT.TFile.Open(url)
    try:
        t=fIn.Get('limit')
        vals=[]
        for i in range(t.GetEntriesFast()):
            t.GetEntry(i)
            if getObs and t.quantileExpected<0:
                vals=[t.limit]*5
            elif not getObs and t.quantileExpected>0:
                vals.append(t.limit)
        fIn.Close()

        if len(vals)<5:
            raise Exception('Could not recover 5 limit quantiles')

    except Exception as e:
        print(e,'@',url)
        vals=[999.]*5

    return [vals[2],vals[3]-vals[2],vals[1]-vals[2],vals[4]-vals[2],vals[0]-vals[2]]


def readSignificanceFrom(url):

    """parses the significance value from combine tree and converts it to a p-value"""

    fIn=ROOT.TFile.Open(url)
    try:
        t=fIn.Get('limit')
        t.GetEntry(0)
        vals=[t.limit,ROOT.RooStats.SignificanceToPValue(t.limit)]
        fIn.Close()
    except:
        vals=[0,0.5]

    return vals


def main():

    baseDir=sys.argv[1]

    #scan all results in directory
    results=[]
    toCheck=[]
    for ana in os.listdir(baseDir):
        
        #if 'optim_' in ana: continue
        anaDir=os.path.join(baseDir,ana)
        if not os.path.isdir(anaDir) : continue
        print anaDir
        limitFiles=[os.path.join(anaDir,x) for x in os.listdir(anaDir) if 'X.obs.AsymptoticLimits' in x]
        for f in limitFiles:
            ch        = re.search('PP([egmz]+)X', f).group(1)
            mass      = int(re.search('.mH(\d+)', f).group(1))
            iresults  = [ana,ch,mass] 
            if ana=='inc_xangle_nvtx' and ch=='z': print '\t',mass,
            iresults += readLimitsFrom(f)
            iresults += readLimitsFrom(f,True)[0:1]
            sigFile='X.Significance' if not '/obs/' in anaDir else 'X.obs.Significance'
            iresults += readSignificanceFrom(f.replace('X.obs.AsymptoticLimits',sigFile))
            if iresults[3]<900:
                results.append( iresults )
            else:
                toCheck.append( (ch,mass,anaDir) )

    print ''
    #save summary in a pandas dataformat
    import pandas as pd
    columns=['ana','channel','mass','r95','drup68','drdn68','drup95','drdn95','r95obs','sig','pval']
    df=pd.DataFrame(data=results, columns=columns)
    df.to_hdf('%s/summary.h5'%baseDir,key='scan')

    print df.head()
    print 'All results are available in %s/summary.h5'%baseDir

    if len(toCheck)>0:
        print 'Recovering %d missing jobs'%len(toCheck)
        #run with condor
        with open('%s/zxstatana_resub_run.sub'%baseDir,'w') as condor:
            condor.write("executable  = $(odir)/statAnaJob.sh\n")
            condor.write("arguments   = $(mass) $(channel)\n")
            condor.write("output       = zxstatana_run.out\n")
            condor.write("error        = zxstatana_run.err\n")
            condor.write("log          = zxstatana_run.log\n")
            condor.write("+JobFlavour = \"longlunch\"\n")
            condor.write("request_cpus = 4\n")        
            for ch,m,odir in toCheck:
                condor.write("odir =%s\n"%os.path.abspath(odir))
                condor.write("mass = %s\n"%m)
                condor.write("channel = %s\n"%ch)
                condor.write("queue 1\n")
        #os.system('condor_submit %s/zxstatana_resub_run.sub'%baseDir)


if __name__ == "__main__":
    main()
