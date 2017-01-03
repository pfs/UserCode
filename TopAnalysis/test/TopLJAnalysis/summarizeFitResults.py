#!/usr/bin/env python

import ROOT
import sys
import os

def parseFitResults(inputDir,toParse):
    results={}

    #open the MultiDimFit results, load fit snapshot and access central values and uncertainties 
    fitResults = [f for f in os.listdir(inputDir) if 'datacard' in f]
    for res in fitResults:
        _,var,bin = res.split('_')
        if not var in results: results[var]={}
        if not bin in results[var]: results[var][bin]={}
        for dataType in toParse:
            results[var][bin][dataType]=[0,0,0,0,0]
            for uncType in ['','_stat']:
                fIn=ROOT.TFile.Open('%s/%s/mlfit_%s%s.root'%(inputDir,res,dataType,uncType))
                fit_s=fIn.Get('fit_s')
                r=fit_s.floatParsFinal().find("r")
                if len(uncType)==0 :
                    results[var][bin][dataType][0]=r.getVal()
                    results[var][bin][dataType][3]=r.getErrorHi()
                    results[var][bin][dataType][4]=r.getErrorLo()
                else:
                    results[var][bin][dataType][1]=r.getErrorHi()
                    results[var][bin][dataType][2]=r.getErrorLo()
    return results


"""
"""
def main():
    inputDir=sys.argv[1]
    toParse=['exp','obs']
    if len(sys.argv)>2:
        toParse=sys.argv[2].split(',')
    results=parseFitResults(inputDir,toParse)
    for var in results:   
        for bin in results[var]:
            for dataType in results[var][bin]:
                val,statHi,statLo,totalHi,totalLo=results[var][bin][dataType]
                systHi=ROOT.TMath.Sqrt(totalHi**2-(statHi**2))
                systLo=ROOT.TMath.Sqrt(totalLo**2-(statLo**2))
                perc=50*(abs(totalHi)+abs(totalLo))                
                resultStr='$%3.2f^{+%3.2f}_{%3.2f}(stat)^{+%3.2f}_{%3.2f}(syst) (%3.0f%%)$'%(val,statHi,statLo,systHi,systLo,perc)
                print '%12s & %12s & %12s & %s'%(var,bin,dataType,resultStr)

if __name__ == "__main__":
    main()
