#!/usr/bin/env python

import ROOT
import sys
import os
import pickle


NUISGROUPS={'trigsel'    : ['eetrig','emtrig','mmtrig','esel','msel','l1prefire'],
            'lepen'      : ["messtat","meszpt","mesewk","mesdm","eesstat","eesgain","eessyst","eessigma","eessphi","eessrho","eesscalet"],
            'btag'       : ["beffhf","befflf"],
            'jer'        : ['JER',"JERstat","JERJEC","JERPU","JERPLI","JERptCut","JERtrunc","JERpTdep","JERSTmFE"],
            'bfrag'      : ["bfrag","slepbr"],
            'toppt'      : ['toppt'],
            'jes'        : ['pileup',"AbsoluteStatJEC","AbsoluteScaleJEC","AbsoluteMPFBiasJEC","FragmentationJEC","SinglePionECALJEC","SinglePionHCALJEC","FlavorPureGluonJEC","FlavorPureQuarkJEC","FlavorPureCharmJEC","FlavorPureBottomJEC","TimePtEtaJEC","RelativeJEREC1JEC","RelativeJEREC2JEC","RelativeJERHFJEC","RelativePtBBJEC","RelativePtEC1JEC","RelativePtEC2JEC","RelativePtHFJEC","RelativeBalJEC","RelativeFSRJEC","RelativeStatFSRJEC","RelativeStatECJEC","RelativeStatHFJEC","PileUpDataMCJEC","PileUpPtRefJEC","PileUpPtBBJEC","PileUpPtEC1JEC","PileUpPtEC2JEC","PileUpPtHFJEC"],
            'qcdscale'   : ["muR","muF","combMuRmuF"],
            'pdf'        : ["PDFenv","PDFaS"],
            'fsr'        : ['fsrg2ggmuR',"fsrg2qqmuR","fsrq2qgmuR","fsrx2xgmuR",'fsrg2ggcNS',"fsrg2qqcNS","fsrq2qgcNS","fsrx2xgcNS"],
            'isr'        : ['isrg2ggmuR',"isrg2qqmuR","isrq2qgmuR","isrx2xgmuR",'isrg2ggcNS',"isrg2qqcNS","isrq2qgcNS","isrx2xgcNS"],
            'hdamp'      : ['hdamp'],
            'cr'         : ['UE','CRerd','CRqcd','CRgmove'],
            'tw'         : ['mtoptw','drdstw'],
            'allsoftqcd' : ['CRerd','CRqcd','CRgmove','FSR'],
            }

def main():

    """loops over the nuisance groups found in a datacard and fixes all nuisances in the group repeating the fit"""

    url=sys.argv[1]

    #read datacard
    with open(url,'r') as f:
        dc=[l.split() for l in f]

    #create the datacards and run 
    nllVals=[]
    for group in NUISGROUPS:
        nuisList=NUISGROUPS[group]

        #skip nuisances to be fixed
        with open('datacard_fixed.dat','w') as f:
            for l_tkns in dc:
                nTkns=len(l_tkns)
                toWrite=True 
                if nTkns>0:

                    for nuisTag in nuisList:
                        if nuisTag in l_tkns[0]: 
                            toWrite=False
                            break

                    if nTkns>1 and l_tkns[1]=='group' : 
                        toWrite=False

                    if l_tkns[0]=='kmax'  : 
                        l_tkns[1]='*'

                if not toWrite : 
                    continue

                f.write(' '.join(l_tkns)+'\n')

        #run the fit and save the likelihood
        try:
            os.system('rm fitresults_fixed.root')
            os.system('text2hdf5.py datacard_fixed.dat')
            os.system('combinetf.py datacard_fixed.dat.hdf5 -o fitresults_fixed.root')
            inF=ROOT.TFile.Open('fitresults_fixed.root')
            tree=inF.Get('fitresults')
            tree.GetEntry(0)
            nll=tree.nllvalfull
            inF.Close()
        except:
            nll=None
        nllVals.append( (group,nll) )

    #remove temporary files
    os.system('rm datacard_fixed.dat')
    os.system('rm fitresults_fixed.rot')

    #dump results to pickle file
    with open('fitresults_fixedgroups.pck','w') as cache:
        pickle.dump(nllVals,cache,pickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":
    sys.exit(main())
