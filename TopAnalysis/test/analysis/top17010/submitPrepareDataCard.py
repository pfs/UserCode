import os
import sys
import optparse
import ROOT

def getScanAnchors(opt):

    """ checks the directories available for likelihood scan anchors """

    scanAnchors=[]

    #list all alternatives (convention is they are found in the scenario* sub-directories)
    scenarioDir=os.path.dirname(opt.templ)
    for sd in os.listdir(scenarioDir):
        if opt.nom and sd==opt.nom:
            scanAnchors.append(('nom',os.path.join(scenarioDir,opt.nom)))
        elif 'scenario' in sd: 
            scenarioURL=os.path.join(scenarioDir,sd)
            root=[f for f in os.listdir(scenarioURL) if  '.root' in f]
            if len(root)==0 : continue        
            scanAnchors.append((sd,os.path.join(scenarioURL,root[0])))

    return scanAnchors

def getSignals(opt):

    """opens the nominal and syst plotter and finds all t#bar{t} plots available"""

    signals=[]
    plotDir=os.path.join(os.path.dirname(opt.templ),'plots')
    plotters=['plotter.root'] #,'syst_plotter.root']
    for i in range(len(plotters)): 
        f=plotters[i]
        f=os.path.join(plotDir,f)
        if not os.path.isfile(f): continue

        if i==0:
            signals.append('data,%s,$(dist)/$(dist)'%f)

        fIn=ROOT.TFile.Open(f)
        for key in fIn.Get('em_mlb').GetListOfKeys():
            keyname=key.GetName()
            pos=keyname.find('t#bar{t}')
            if pos<0: continue
            keyname=keyname[pos:]
            signals.append( 'sig,%s,$(dist)/$(dist)_%s'%(f,keyname) )
    return signals
    
def generateJobs(scanAnchors,signals,opt):

    """ loops over distributions to generate the jobs """

    cmssw_base=os.environ['CMSSW_BASE']

    condor=open('datacard_condor.sub','w')
    condor.write('executable  = %s/src/TopLJets2015/TopAnalysis/scripts/wrapLocalAnalysisRun.sh\n'%cmssw_base)
    condor.write('output      = datacard_condor.out\n')
    condor.write('error       = datacard_condor.err\n')
    condor.write('log         = datacard_condor.log\n')
    condor.write('+JobFlavour = "workday"\n')

    dataStr=''
    for s in signals:
        dataStr='dataDef=%s %s'%(s,dataStr)
    
    script='{0}/src/TopLJets2015/TopAnalysis/test/analysis/top17010/prepareDataCard.py'.format(cmssw_base)
    condor.write('arguments = %s %s -d $(dist) -t %s %s -s $(anchor) -o %s --systs %s\n'%(cmssw_base,script,opt.templ,dataStr,opt.outDir,opt.systs))
    condor.write('queue dist, anchor from (\n')
    ijobs=0
    for d in opt.dists.split(','):
        for a in scanAnchors:
            condor.write('%s %s\n'%(d,','.join(a)))
            ijobs+=1     
    condor.write(')\n')

    condor.close()

    return ijobs


def main():

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('--systs',          
                      dest='systs',       
                      help='description of the systematics [%default]',
                      default='test/analysis/top17010/systs_dict.json',
                      type='string')
    parser.add_option('--nom',
                      dest='nom',
                      help='nominal simlation file [%default]',  
                      default=None,
                      type='string')
    parser.add_option('-t', '--templ',          
                      dest='templ',
                      help='template directory [%default]',  
                      default='/eos/cms/store/cmst3/group/top/TOP17010/0c522df/templates',
                      type='string')
    parser.add_option('-o', '--outdir',          
                      dest='outDir',
                      help='output directory [%default]',  
                      default='datacard',
                      type='string')
    parser.add_option('--dists',          
                      dest='dists',       
                      help='CSV list of distributions [%default]',
                      default='em_mlb',
                      type='string')
    (opt, args) = parser.parse_args()

    #build the list of scan anchors
    scanAnchors=getScanAnchors(opt)
    print '%d scan likelihood scan anchors have been found'%len(scanAnchors)

    #build the list of data/signals
    signals=getSignals(opt)
    print '%d potential signals have been found'%len(signals)

    #generate the jobs
    njobs=generateJobs(scanAnchors,signals,opt)
    print '%d jobs will be submitted to create the likelihood scan datacards - see datacard_condor.sub'%njobs

    os.system('condor_submit datacard_condor.sub')

if __name__ == "__main__":
    sys.exit(main())
