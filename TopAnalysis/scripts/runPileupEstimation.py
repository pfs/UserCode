import optparse
import os,sys
import json
import commands
import ROOT
from SimGeneral.MixingModule.mix_2016_25ns_Moriond17MC_PoissonOOTPU_cfi import *

"""
steer the script
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('--json',      dest='inJson'  ,      help='json file with processed runs [%default]',      
                      default='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt',    
                      type='string')
    parser.add_option('--out',       dest='output'  ,      help='output file with weights',           default=None,    type='string')
    parser.add_option('--mbXsec',    dest='mbXsec'  ,      help='minimum bias cross section to use [%default]',  
                      default=69200,   type=float)
    parser.add_option('--puJson',    dest='puJson'  ,      help='pileup json file [%default]',      
                      default='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt',
                      type='string')
    (opt, args) = parser.parse_args()
    
    #simulated pileup
    NPUBINS=len(mix.input.nbPileupEvents.probValue)
    MAXPU=NPUBINS    
    simPuH=ROOT.TH1F('simPuH','',NPUBINS,float(0),MAXPU)
    for xbin in xrange(0,len(mix.input.nbPileupEvents.probValue)):
        probVal=mix.input.nbPileupEvents.probValue[xbin]
        simPuH.SetBinContent(xbin,probVal)
    simPuH.Scale(1./simPuH.Integral())

    #compute pileup in data assuming different xsec
    puDist=[]
    puWgts=[]
    MINBIASXSEC={'nom':opt.mbXsec,'up':opt.mbXsec*1.05,'down':opt.mbXsec*0.95}
    for scenario in MINBIASXSEC:
        print scenario, 'xsec=',MINBIASXSEC[scenario]
        cmd='pileupCalc.py -i %s --inputLumiJSON %s --calcMode true --minBiasXsec %f --maxPileupBin %d --numPileupBins %s Pileup.root'%(opt.inJson,opt.puJson,MINBIASXSEC[scenario],MAXPU,NPUBINS)
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
    fOut=ROOT.TFile.Open(opt.output,'RECREATE')
    for gr in puWgts: gr.Write()
    for gr in puDist: gr.Write()
    simPuH.Write()
    fOut.Close()

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
