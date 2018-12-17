#!/usr/bin/env python

import os,sys
import json
import commands
import optparse
from TopLJets2015.TopAnalysis.miniAnalyzer_cfi import ANALYSISTRIGGERLISTS
import ROOT

"""
runs lumicalc
"""
def getLumiTable(args):

    trig,opt=args      
    if len(trig)==0: trig='inc'

    #run brilcalc
    extraOpt='--hltpath %s*'%trig if trig!='inc' else ''    
    rawLumiTable=commands.getstatusoutput('brilcalc lumi -b "STABLE BEAMS" -u /pb -i %s --normtag %s %s'%(opt.lumiMask,opt.normtag,extraOpt))[1].split('\n')

    #parse lumi output file
    valList=[]
    for line in rawLumiTable:
        tokens=line.split()
        try:
            run,lumisec=tokens[1].split(':')
            lumi=tokens[-2]
            valList.append( (int(run),float(lumi)) )
        except:
            pass

    #fill histogram
    nRuns=len(valList)
    try:
        f=ROOT.TFile.Open('%s/lumisec_inc.root'%opt.out)
        h=f.Get('lumisec_inc').Clone('lumisec_%s'%trig)
        h.SetDirectory(0)
        f.Close()
        h.Reset('ICE')
        for i in xrange(0,nRuns):
            run,lumi=valList[i]
            binLabel='%d'%run
            for xbin in xrange(1,h.GetNbinsX()+1):
                ilabel=h.GetXaxis().GetBinLabel(xbin)
                if ilabel!=binLabel : continue
                h.SetBinContent(xbin,lumi)
                break
    except:
        h=ROOT.TH1F('lumisec_%s'%trig,'lumisec;Run number;Lumi recorded (pb);',nRuns,0,nRuns)
        for i in xrange(0,nRuns):
            run,lumi=valList[i]
            h.GetXaxis().SetBinLabel(i+1,'%d'%run)
            h.SetBinContent(i+1,lumi)

    #save to file
    print trig,' lumi=',h.Integral(),'/pb'
    f=ROOT.TFile.Open('%s/lumisec_%s.root'%(opt.out,trig),'RECREATE')
    h.SetDirectory(f)
    h.Write()
    f.Close()



"""
steer the script
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-o', '--out',      dest='out'   ,      help='output dir',         default='data/era2017/',    type='string')
    parser.add_option(      '--normtag',  dest='normtag'   ,  help='norm tag',           default='/cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json',    type='string')
    parser.add_option('-y', '--year',     dest='year'   ,     help='year [%default]',    default=2017, type=int)
    parser.add_option('-l', '--lumi',     dest='lumiMask',    help='json with list of good lumis', default='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt')
    (opt, args) = parser.parse_args()
    

    #inclusive (run separetely to create template histogram)
    getLumiTable( ('',opt) )

    #per trigger path
    task_list=[]
    triggerPaths=ANALYSISTRIGGERLISTS[opt.year]
    for trig in triggerPaths:
        task_list.append( (trig,opt) )
    from multiprocessing import Pool
    pool=Pool(8)
    pool.map(getLumiTable, task_list)


    os.system('hadd -f %s/lumisec.root %s/lumisec_*.root'%(opt.out,opt.out))
    os.system('rm %s/lumisec_*.root'%opt.out)


"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
