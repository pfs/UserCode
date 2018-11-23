#!/usr/bin/env python

import sys
import json
import ROOT
import re
from array import array

def getRange(tkn):
    return [float(x) for x in re.findall(r"[-+]?\d*\.\d+|\d+",tkn)]

url=sys.argv[1]
fOut=ROOT.TFile.Open(url.replace('.json','.root'),'RECREATE')

with open(url,'r') as f:
    results = json.load(f)

    for i in ['Loose','Medium','Tight']:
        key='NUM_%sID_DEN_genTracks'%i
        if not key in results:
            if i=='Loose' : key='NUM_LooseRelIso_DEN_LooseID'
            if i=='Medium': key='NUM_TightRelIso_DEN_MediumID'
            if i=='Tight':  key='NUM_TightRelIso_DEN_TightIDandIPCut'

        effVals=[]
        for etaTkn in results[key]['abseta_pt']:

            ieffVals=[]
            for ptTkn in results[key]['abseta_pt'][etaTkn]:
                ieffVals.append( getRange(ptTkn)+[results[key]['abseta_pt'][etaTkn][ptTkn]["value"],results[key]['abseta_pt'][etaTkn][ptTkn]["error"]] )
            ieffVals.sort(key=lambda x : x[0])
            effVals.append( (getRange(etaTkn),ieffVals) )

        effVals.sort(key=lambda x : x[0][0])

        etaVals=[]
        for i in xrange(0,len(effVals)):
            etaVals.append( effVals[i][0][0] )
        etaVals.append(effVals[-1][0][1])

        ptVals=[]
        for i in xrange(0,len(effVals[0][1])):
            ptVals.append( effVals[0][1][i][0] )
        ptVals.append( effVals[0][1][-1][1] )

        fOut.cd()
        sfH=ROOT.TH2F(key,key,len(etaVals)-1,array('d',etaVals),len(ptVals)-1,array('d',ptVals))
        for i in xrange(0,len(effVals)):
            for j in xrange(0,len(effVals[i][1])):
                sf,sfUnc=effVals[i][1][j][2:4]
                sfH.SetBinContent(i+1,j+1,sf)
                sfH.SetBinError(i+1,j+1,sfUnc)
        sfH.Write()

fOut.Close()
        


