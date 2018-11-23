import os
import sys
import optparse
import ROOT
import pickle
import json
from TopLJets2015.TopAnalysis.storeTools import *
from runQCDEstimation import getTemplates
from TopLJets2015.TopAnalysis.rounding import *

"""
implementation of RinRout a la TOP-16-005
"""
def doRinRout(nin,nout,ndyin,ndyout,ch1,ch2):
    
    rinll    = nin[ch1][0]/nin[ch2][0]
    rinllUnc = rinll*ROOT.TMath.Sqrt(1./nin[ch1][0]+1./nin[ch2][0])

    kll    = ROOT.TMath.Sqrt(rinll)
    kllUnc = 0.5*kll*rinllUnc/rinll
    
    routin    = ndyout[ch1][0]/ndyin[ch1][0]
    routinUnc = routin*ROOT.TMath.Sqrt( (ndyout[ch1][1]/ndyout[ch1][0])**2 
                                        + (ndyin[ch1][1]/ndyin[ch1][0])**2 )

    nllout    = routin*(nin[ch1][0]-0.5*nin['EM'][0]*kll)
    nlloutUnc = ROOT.TMath.Sqrt( (routinUnc*(nin[ch1][0]-0.5*nin['EM'][0]*kll))**2
                                 + (routin**2)*nin[ch1][0]
                                 + nin['EM'][0]*(0.5*routin*kll)**2
                                 + (0.5*routin*nin['EM'][0]*kllUnc)**2 )
    
    sf    = nllout/ndyout[ch1][0]
    sfUnc = sf*ROOT.TMath.Sqrt( (nlloutUnc/nllout)**2+(ndyout[ch1][1]/ndyout[ch1][0])**2 )

    #print report
    print '-'*20,ch1,'-'*20
    print 'R_{out/nin}(MC)=',toLatexRounded(routin,routinUnc)
    print 'k_{\\ell,\\ell}=',toLatexRounded(kll,kllUnc)
    print 'N_{\\rm in}(data)=',nin[ch1][0]
    print 'N_{\\rm out}=',toLatexRounded(nllout,nlloutUnc)
    print 'SF_{\\rm DY}=',toLatexRounded(sf,sfUnc)
    print '-'*50

    return (sf,sfUnc)
    

"""
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('--in',           dest='input',   help='plotter file with the dilepton selection', default=None,       type='string')
    parser.add_option('--categs',       dest='categs',  help='categories to apply Rin/Rout',             default='1b,2b',    type='string')
    parser.add_option('--out',          dest='outdir',  help='output directory',                         default='./',       type='string')
    (opt, args) = parser.parse_args()

    #prepare output
    os.system('mkdir -p %s'%opt.outdir)

    #prepare fitter
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(True)
    ROOT.FWLiteEnabler.enable()
            
    #open inputs
    fIn=ROOT.TFile.Open(opt.input)
    dySF={}
    for sel in opt.categs.split(','):
        #compute integrals
        nin,nout,ndyin,ndyout={},{},{},{}
        intErr=ROOT.Double(0)
        for ch in ['EE','MM','EM']:
            try:
                data, sumMC, refH  = getTemplates(fIn=fIn, dist='%s%s_mll'%(ch,sel), tag=ch+sel,refName='DY')
            except:
                data, sumMC, refH  = getTemplates(fIn=fIn, dist='mll_%s%s'%(ch,sel), tag=ch+sel,refName='DY')
            binMin,binMax=data.GetXaxis().FindBin(20),data.GetNbinsX()+1
            bin1,bin2=data.GetXaxis().FindBin(76),data.GetXaxis().FindBin(106)

            nin[ch]=(data.Integral(bin1,bin2),0.)
            nout[ch]=(data.Integral(binMin,binMax)-nin[ch][0],0.)

            dyCts=refH.IntegralAndError(bin1,bin2,intErr)
            ndyin[ch]=(dyCts,intErr)

            dyCts=refH.IntegralAndError(binMin,binMax,intErr)-dyCts
            ndyout[ch]=(dyCts,ROOT.TMath.Sqrt(intErr**2-ndyin[ch][1]**2))

        #table a la TOP-16-015 
        print '*'*20,sel,'*'*20
        dySF['EE'+sel]=doRinRout(nin,nout,ndyin,ndyout,'EE','MM')
        dySF['MM'+sel]=doRinRout(nin,nout,ndyin,ndyout,'MM','EE')

        sfProd=dySF['EE'+sel][0]*dySF['MM'+sel][0]
        sfProdUnc=ROOT.TMath.Sqrt((dySF['EE'+sel][0]*dySF['MM'+sel][1])**2+(dySF['MM'+sel][0]*dySF['EE'+sel][1])**2)        
        dySF['EM'+sel]=(ROOT.TMath.Sqrt(sfProd),0.5*sfProdUnc/ROOT.TMath.Sqrt(sfProd))
        print 'SF_{e\\mu}=',toLatexRounded(dySF['EM'+sel][0],dySF['EM'+sel][1])
        print '*'*50

    fIn.Close()

    #dump to file    
    cachefile=open('%s/.dyscalefactors.pck'%opt.outdir,'w')
    print 'DY scale factors saved in ',cachefile.name
    pickle.dump(dySF, cachefile, pickle.HIGHEST_PROTOCOL)
    cachefile.close()



"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
