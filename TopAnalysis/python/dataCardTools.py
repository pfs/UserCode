#!/usr/bin/env python

import optparse
import os,sys
import commands
from array import array
import ROOT

"""
Prepares the fitting script
"""
def prepareFitScript(datacard,POIs,unblind=False):
    baseDir=os.path.dirname(datacard)
    fitScriptUrl='%s/runCombine.sh'%baseDir
    fitScript=open(fitScriptUrl,'w')

    fitScript.write('cd %s\n'%baseDir)
    
    fitScript.write('\n# convert datacard to workspace\n')
    fitScript.write('echo \"Converting datacard to workspace\"\n')
    fitScript.write('text2workspace.py %s -m 0 -o workspace.root\n' % os.path.basename(datacard))

    fitScript.write('\n# likelihood scans\n')
    for parameter in POIs:

        minParVal=0.8 if parameter=='r' else -2.0
        maxParVal=1.2 if parameter=='r' else +2.0
        rangeOpt='--setPhysicsModelParameterRanges %s=%f,%f'%(parameter,minParVal,maxParVal)

        poiOpt='' if parameter=='r' else '--redefineSignalPOIs %s'%parameter

        if parameter=='r':
            fitScript.write('\n## max likelihood fit\n')
            fitScript.write('echo \"Running MaxLikelihoodFit for r\"\n')
            fitScript.write('combine workspace.root -M MaxLikelihoodFit -t -1 --expectSignal=1 -m 0\n')
            fitScript.write('mv mlfit.root mlfit_exp.root\n')
            
            fitScript.write('\n## impacts\n')
            fitScript.write('echo \"To compute expected impacts re-run runCombine.sh uncommenting the appropriate lines\n')
            fitScript.write('#combineTool.py -M Impacts -d workspace.root -m 0  -t -1 --expectSignal=1 --doInitialFit\n')
            fitScript.write('#combineTool.py -M Impacts -d workspace.root -m 0  -t -1 --expectSignal=1 --doFits\n')
            fitScript.write('#combineTool.py -M Impacts -d workspace.root -m 0 -o impacts.json\n')
            
            fitScript.write('\n## toys\n')
            fitScript.write('echo \n"To run toys  re-run runCombine.sh uncommenting the appropriate lines\n')
            fitScript.write('#rscan=(0.9 1.0 1.1)\n')
            fitScript.write('#for r in ${rscan[@]}; do\n')
            fitScript.write('#\t combine workspace.root -M MaxLikelihoodFit -t 100 --expectSignal=${r} -m ${r} --toysFrequentist --noErrors --minos none;\n')
            fitScript.write('#done\n\n')

            if unblind:
                fitScript.write('combine workspace.root -M MaxLikelihoodFit -m 0\n')
                fitScript.write('mv mlfit.root mlfit_obs.root\n')
                            
        fitScript.write('\n## function of %s\n'%parameter)
        fitScript.write('echo \"Running likelihood scan for %s\"\n'%parameter)
        fitScript.write('combine workspace.root -M MultiDimFit -P %s -t -1 --expectSignal=1 --algo=grid --points=50 %s %s -m 0\n'%(parameter,rangeOpt,poiOpt))
        fitScript.write('mv higgsCombineTest.MultiDimFit.mH0.root exp_plr_scan_%s.root\n'%parameter)
        fitScript.write('combine workspace.root -M MultiDimFit -P %s -t -1 --expectSignal=1 --algo=grid --points=50 %s %s -m 0 -S 0\n'%(parameter,rangeOpt,poiOpt))
        fitScript.write('mv higgsCombineTest.MultiDimFit.mH0.root exp_plr_scan_stat_%s.root\n'%parameter)
        if unblind:
            fitScript.write('combine workspace.root -M MultiDimFit -P %s --algo=grid --points=50 %s %s -m 0\n'%(parameter,rangeOpt,poiOpt))
            fitScript.write('mv higgsCombineTest.MultiDimFit.mH0.root obs_plr_scan_%s.root\n'%parameter)
            fitScript.write('combine workspace.root -M MultiDimFit -P %s --algo=grid --points=50 %s %s -m 0 -S 0\n'%(parameter,rangeOpt,poiOpt))
            fitScript.write('mv higgsCombineTest.MultiDimFit.mH0.root obs_plr_scan_stat_%s.root\n'%parameter)

    fitScript.write('\n# 2D likelihood scans\n')
    for i in xrange(0,len(POIs)):
        for j in xrange(i+1,len(POIs)):
            minPiVal=0.8 if POIs[i]=='r' else -2.0
            maxPiVal=1.2 if POIs[i]=='r' else +2.0
            minPjVal=0.8 if POIs[j]=='r' else -2.0
            maxPjVal=1.2 if POIs[j]=='r' else +2.0
            rangeOpt='--setPhysicsModelParameterRanges %s=%f,%f:%s=%f,%f'%(POIs[i],minPiVal,maxPiVal,POIs[j],minPjVal,maxPjVal)
            
            fitScript.write('## function of %s,%s\n'%(POIs[i],POIs[j]))
            fitScript.write('echo \"Running 2D likelihood scan for %s vs %s\"\n'%(POIs[i],POIs[j]))
            fitScript.write('combine workspace.root -M MultiDimFit --redefineSignalPOIs %s,%s -P %s -P %s  -t -1 --expectSignal=1 --algo=grid --points=1000 %s -m 0\n'%(POIs[i],POIs[j],POIs[i],POIs[j],rangeOpt)) 
            fitScript.write('mv higgsCombineTest.MultiDimFit.mH0.root  exp_plr_scan_%svs%s.root\n'%(POIs[i],POIs[j]))
            if unblind:
                fitScript.write('combine workspace.root -M MultiDimFit --redefineSignalPOIs %s,%s -P %s -P %s --algo=grid --points=1000  %s -m 0\n'%(POIs[i],POIs[j],POIs[i],POIs[j],rangeOpt)) 
                fitScript.write('mv higgsCombineTest.MultiDimFit.mH0.root obs_plr_scan_%svs%s.root\n'%(POIs[i],POIs[j]))

    fitScript.write('\ncd -\n')
    
    fitScript.close()
    return fitScriptUrl

"""
test if variation is significant enough i.e. if sum_{bins} |var-nom| > tolerance
"""
def acceptVariationForDataCard(nomH,upH,downH,tol=1e-2):
    diffUp,diffDown=0,0
    for xbin in xrange(1,nomH.GetNbinsX()):
        val,valUp,valDown=nomH.GetBinContent(xbin),upH.GetBinContent(xbin),downH.GetBinContent(xbin)
        diffUp+=ROOT.TMath.Abs(valUp-val)
        diffDown+=ROOT.TMath.Abs(valDown-val)
    accept = True if (diffUp>tol or diffDown>tol) else False
    return accept


"""
in case of multiple signals, remove others
"""
def filterShapeList(exp,signalList,rawSignalList,signalSub=None):
    newExp={}
    for key in exp:

        matchFound=False
        for rawKey in rawSignalList:
            if rawKey in key:
                matchFound=True
                try:
                    exp[key].Add(exp[signalSub],-1)
                    print 'Subtracted',signalSub,'from',key
                except:
                    pass
                        
        if matchFound and not key in signalList : continue

        newExp[key]=exp[key]
    return newExp


"""
get distributions from file
"""
def getDistsFrom(directory,keyFilter=''):
    obs=None
    exp={}
    dirName=directory.GetName()
    for key in directory.GetListOfKeys():
        if len(keyFilter)>0 and key.GetName()!='%s_%s'%(dirName,keyFilter) : continue
        obj=directory.Get(key.GetName())        
        if not obj.InheritsFrom('TH1') : continue
        if obj.GetName()==dirName : 
            obs=obj.Clone('data_obs')
            obs.SetDirectory(0)
        else : 
            newName=obj.GetName().split(dirName+'_')[-1]
            for token in ['+','-','*',' ','#','{','(',')','}','@']:
                newName=newName.replace(token,'')
            newName=newName.replace('=','eq')
            newName=newName.replace('.','p')
            exp[newName]=obj.Clone(newName+'_'+dirName)
            exp[newName].SetDirectory(0)
            if exp[newName].InheritsFrom('TH2'):
                for xbin in xrange(1,exp[newName].GetNbinsX()+1):
                    for ybin in xrange(1,exp[newName].GetNbinsY()+1):
                        binContent=exp[newName].GetBinContent(xbin,ybin)
                        if binContent>0: continue
                        newBinContent=ROOT.TMath.Max(ROOT.TMath.Abs(binContent),1e-3)
                        exp[newName].SetBinContent(xbin,ybin,newBinContent)
                        exp[newName].SetBinError(xbin,ybin,newBinContent)
            else:
                for xbin in xrange(1,exp[newName].GetNbinsX()+1):
                    binContent=exp[newName].GetBinContent(xbin)
                    if binContent>0: continue
                    newBinContent=ROOT.TMath.Max(ROOT.TMath.Abs(binContent),1e-3)
                    exp[newName].SetBinContent(xbin,newBinContent)
                    exp[newName].SetBinError(xbin,newBinContent)

    return obs,exp

"""
save distributions to file
"""
def saveToShapesFile(outFile,shapeColl,directory='',rebin=0):
    fOut=ROOT.TFile.Open(outFile,'UPDATE')
    if len(directory)==0:
        fOut.cd()     
    else:
        outDir=fOut.mkdir(directory)
        outDir.cd()
    for key in shapeColl:
        #remove bin labels  
        shapeColl[key].GetXaxis().Clear()

        #convert to TH1D (projections are TH1D)
        #if not shapeColl[key].InheritsFrom('TH1D') :            
        #    h=ROOT.TH1D()
        #    shapeColl[key].Copy(h)
        #    shapeColl[key]=h
        
        #reinforce that all histos have the same type
        #projections from 2D histos will be TH1D will nominal may be TH1F...
            
        h2write=ROOT.TH1F(key,
                          '%s;%s;%s'%(shapeColl[key].GetTitle(),shapeColl[key].GetXaxis().GetTitle(),shapeColl[key].GetYaxis().GetTitle()),
                          shapeColl[key].GetNbinsX(),
                          shapeColl[key].GetXaxis().GetXmin(),
                          shapeColl[key].GetXaxis().GetXmax())
        for xbin in xrange(1,shapeColl[key].GetNbinsX()+1):
            h2write.SetBinContent(xbin,shapeColl[key].GetBinContent(xbin))

        #rebin final shape, if required
        #h2write=shapeColl[key].Clone()
        if rebin!=0: h2write.Rebin(rebin)

        #save to file
        h2write.SetDirectory(outDir)        
        h2write.Write(key,ROOT.TObject.kOverwrite)

    #all done here
    fOut.Close()

"""
"""
def replaceBadCharacters(inputStr):
    newStr = inputStr
    for token in ['+','-','*',' ','#','{','(',')','}','@']:
        newStr=newStr.replace(token,'')
    return newStr
