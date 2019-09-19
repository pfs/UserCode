import ROOT
import sys
import os
from subprocess import check_output

def modification_date(f):
    t = os.path.getmtime(f)
    return datetime.datetime.fromtimestamp(t)

def checkIntegrity(f,tname):

    errorCode=None
    nentries=None

    #try to open the file and do a basic check
    try:
        inF=ROOT.TFile.Open(f)
        keys=inF.GetListOfKeys()
        if inF.IsZombie():
            errorCode=1
        if keys.GetSize()==0:
            errorCode=2
    except:
        errorCode=0
        return (errorCode,nentries)
        
    #read tree if available
    if tname:
        try:
            t=inF.Get(tname)
            nentries=t.GetEntriesFast()
        except: 
            errorCode=3

    #close file
    inF.Close()
    
    tstamp=None
    try:
        tsamp=modification_date(f.replace('root://eoscms//','/eos/cms'))
    except:
        pass

    return (errorCode,nentries,tstamp)




#read files
FARMDIR=sys.argv[1]
f=open(os.path.join(FARMDIR,'checkIntegList.dat'))
files = [line.rstrip('\n').split() for line in f]
f.close()
print len(files),'files to check'

itname=sys.argv[2] if len(sys.argv)>2 else None
otname=sys.argv[3] if len(sys.argv)>3 else None

runLocal=False
toRun=[]

with open(os.path.join(FARMDIR,'localanalysis_integ_report.dat'),'w') as r:
    for i,o in files:

        if i.find('root://eoscms//')>=0:
            i=i.replace('/eos/cms/','')

        ires=checkIntegrity(i,itname)
        ores=checkIntegrity('root://eoscms//'+o.replace('/eos/cms/',''),otname)

        msg='IS OK' if ires[0] is None and ores[0] is None else 'NOTOK'
        ijob='/'.join(i.split('/')[-2:])
        r.write('[{0:6s}] job with input {1} input={2} output={3}\n'.format(msg,ijob,ires,ores))
        
        if msg=='NOTOK':            
            shScript=check_output(["grep -ir {0} {1}/*.sh".format(ijob,FARMDIR)],shell=True).split(':')[0]
            if runLocal:
                print 'Running locally the job with input',ijob                
                os.system('sh {0}'.format(shScript))
            else:
                base=os.path.basename(shScript)
                shScript=os.path.splitext(base)[0]
                toRun += [shScript]

print 'Report in localanalysis_integ_report.dat'
if not runLocal and len(toRun)>0:
    
    print 'Resubmitting',len(toRun),'jobs for',FARMDIR
    OpSysAndVer = str(os.system('cat /etc/redhat-release'))
    if 'SLC' in OpSysAndVer:
        OpSysAndVer = "SLCern6"
    else:
        OpSysAndVer = "CentOS7"

    with open('%s/condor.sub'%FARMDIR,'r') as f:
        condorLines = f.readlines()

    with open ('%s/condor_resub.sub'%FARMDIR,'w') as condor:
        for l in condorLines:
            if 'cfgFile=' in l:
                cfg=l.split('=')[1].rstrip('\n')
                if not cfg in toRun : 
                    continue
                condor.write(l)
                condor.write('queue 1\n')
            elif 'queue' in l : 
                continue
            else:
                condor.write(l)

    os.system('condor_submit %s/condor_resub.sub'%FARMDIR)


    


