import ROOT
import sys
import os
from subprocess import check_output
from checkLocalAnalysisInteg import checkIntegrity

def main():

    #read files
    FARMDIR=sys.argv[1]
    f=open(os.path.join(FARMDIR,'condor.sub'))
    tasks = [line.rstrip('\n').split()[2:] for line in f if 'arguments' in line]
    f.close()
    print len(tasks),'tasks to check'

    for t in tasks:
        outf=t[0]
        ores=checkIntegrity(outf,'analysis/data')
        if ores[0] is None: continue
        print 'Will redo',outf,'out of ',len(t)-1,'inputs'
        args=' '.join(t)
        os.system('sh %s/custom_merge.sh %s'%(FARMDIR,args))

if __name__ == "__main__":
    sys.exit(main())

