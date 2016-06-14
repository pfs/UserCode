import optparse
import os,sys
import json
import commands
import ROOT
from SimGeneral.MixingModule.mix_2015_25ns_FallMC_matchData_PoissonOOTPU_cfi import *

"""
steer the script
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('--json',      dest='inJson'  ,      help='json file with processed runs',      default=None,    type='string')
    parser.add_option('--mbXsec',    dest='mbXsec'  ,      help='minimum bias cross section to use',  default=69000,   type=float)
    parser.add_option('--puJson',    dest='puJson'  ,      help='pileup json file',      
                      default='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/PileUp/pileup_latest.txt',    type='string')
    (opt, args) = parser.parse_args()
    
    #simulated pileup
    NPUBINS=len(mix.input.nbPileupEvents.probValue)
    simPuH=ROOT.TH1F('simPuH','',NPUBINS,float(0),float(NPUBINS))
    for xbin in xrange(0,NPUBINS):
        probVal=mix.input.nbPileupEvents.probValue[xbin]
        simPuH.SetBinContent(xbin,probVal)
    simPuH.Scale(1./simPuH.Integral())

    #compute pileup in data assuming different xsec
    puDist=[]
    puWgts=[]
    MINBIASXSEC={'nom':opt.mbXsec,'up':opt.mbXsec*1.05,'down':opt.mbXsec*0.95}
    for scenario in MINBIASXSEC:
        print scenario, 'xsec=',MINBIASXSEC[scenario]
        cmd='pileupCalc.py -i %s --inputLumiJSON %s --calcMode true --minBiasXsec %f --maxPileupBin %d --numPileupBins %s Pileup.root'%(opt.inJson,opt.puJson,MINBIASXSEC[scenario],NPUBINS,NPUBINS)
        commands.getstatusoutput(cmd)

        fIn=ROOT.TFile.Open('Pileup.root')
        pileupH=fIn.Get('pileup')
        pileupH.Scale(1./pileupH.Integral())
        puDist.append( ROOT.TGraph(pileupH) )
        puDist[-1].SetName('pu_'+scenario)

        pileupH.Divide(simPuH)
        puWgts.append( ROOT.TGraph(pileupH) )
        puWgts[-1].SetName('puwgts_'+scenario)
        fIn.Close()
    commands.getstatusoutput('rm Pileup.root')

    #save pileup weights to file
    fOut=ROOT.TFile.Open('$CMSSW_BASE/src/TopLJets2015/TopAnalysis/data/pileupWgts.root','RECREATE')
    for gr in puWgts: gr.Write()
    for gr in puDist: gr.Write()
    simPuH.Write()
    fOut.Close()

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
