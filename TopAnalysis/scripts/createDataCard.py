import os
import sys
import optparse
import ROOT
import commands
import getpass
import pickle

<<<<<<< HEAD
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
def filterShapeList(exp,signalList,rawSignalList):
    newExp={}
    for key in exp:

        matchFound=False
        for rawKey in rawSignalList:
            if rawKey in key:
                matchFound=True
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
            exp[newName]=obj.Clone(newName)
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
        if not shapeColl[key].InheritsFrom('TH1D') :
            h=ROOT.TH1D()
            shapeColl[key].Copy(h)
            shapeColl[key]=h

        #rebin final shape, if required
        if rebin!=0: shapeColl[key].Rebin(rebin)

        #write
        shapeColl[key].Write(key,ROOT.TObject.kOverwrite)
    fOut.Close()

=======
from TopLJets2015.TopAnalysis.dataCardTools import *
from TopLJets2015.TopAnalysis.xsecSystSpecs import *
>>>>>>> 9058e3898a0b728430a2f63a9c409bf903ecdb62

"""
steer the script
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--input',          dest='input',       help='input plotter',                                      default=None,          type='string')
    parser.add_option(      '--systInput',      dest='systInput',   help='input plotter for systs from alternative samples',   default=None,          type='string')
    parser.add_option('-d', '--dist',           dest='dist',        help='distribution',                                       default='nbtags',      type='string')
    parser.add_option('-s', '--signal',         dest='signal',      help='signal (csv)',                                       default='tbart',       type='string')
    parser.add_option('--signalSub',            dest='signalSub',   help='subtract this component out of signal',              default=None,          type='string')
    parser.add_option(      '--specs',          dest='specs',       help='specifications [%default]',                          default='TOP-16-006',       type='string')
    parser.add_option('-m', '--mass',           dest='mass',        help='signal mass',                                        default=0,             type=float)
    parser.add_option('-c', '--cat',            dest='cat',         help='categories (csv)',                                   default='1j,2j,3j,4j', type='string')
    parser.add_option('-q', '--qcd',            dest='qcd',         help='qcd normalization file',                             default=None,          type='string')
    parser.add_option('-w', '--wjets',          dest='wjets',       help='wjets normalization file',                           default=None,          type='string')
    parser.add_option('--addBinByBin',          dest='addBinByBin', help='add bin-by-bin stat uncertainties @ threshold',      default=-1,            type=float)
    parser.add_option(      '--rebin',          dest='rebin',       help='histogram rebin factor',                             default=0,             type=int)
    parser.add_option('-o', '--output',         dest='output',      help='output directory',                                   default='datacards',   type='string')
    parser.add_option(      '--rebin',          dest='rebin',       help='histogram rebin factor',                             default=0,             type=int)
    (opt, args) = parser.parse_args()

    rawSignalList=opt.signal.split(',')
    signalList=rawSignalList[:]
    if opt.mass!=0:
        for i in xrange(0,len(rawSignalList)):
            signalList[i]='%sm%3.0f'%(rawSignalList[i],10*opt.mass)

    catList=opt.cat.split(',')

    #read qcd normalization
    qcdNorm=None
    if opt.qcd:
        cache=open(opt.qcd,'r')
        qcdNorm=pickle.load(cache)


    #read wjets normalization
    wjetsNorm=None
    if opt.wjets:
        cache=open(opt.wjets,'r')
        wjetsNorm=pickle.load(cache)

    #prepare output directory 
    os.system('mkdir -p %s'%opt.output)

    anCat=''
    for subDir in opt.input.split('/'):
        if 'analysis_' not in subDir: continue
        anCat=subDir.replace('analysis_','')

    #get data and nominal expectations
    fIn=ROOT.TFile.Open(opt.input)
    systfIn=None
    if opt.systInput:
        systfIn=ROOT.TFile.Open(opt.systInput)

    #loop over categories
    for cat in catList:

        print 'Initiating %s datacard for %s'%(opt.dist,cat)

<<<<<<< HEAD
        #nominal expectations
=======
        #get syst specifications
        rateSysts,sampleSysts=xsecSystSpecs(opt.specs)

        #nomimal expectations
>>>>>>> 9058e3898a0b728430a2f63a9c409bf903ecdb62
        obs,exp=getDistsFrom(directory=fIn.Get('%s_%s'%(opt.dist,cat)))
        exp=filterShapeList(exp,signalList,rawSignalList,opt.signalSub)

        #prepare output ROOT file
        outFile='%s/shapes_%s.root'%(opt.output,cat)
        fOut=ROOT.TFile.Open(outFile,'RECREATE')
        fOut.Close()
        
        #start the datacard
        datacard=open('%s/datacard_%s.dat'%(opt.output,cat),'w')
        datacard.write('#\n')
        datacard.write('# Generated by %s with git hash %s for analysis category %s\n' % (getpass.getuser(),
                                                                                          commands.getstatusoutput('git log --pretty=format:\'%h\' -n 1')[1],
                                                                                          cat) )

        datacard.write('#\n')
        datacard.write('imax *\n')
        datacard.write('jmax *\n')
        datacard.write('kmax *\n')
        datacard.write('-'*50+'\n')
        datacard.write('shapes *        * shapes_%s.root nom/$PROCESS $SYSTEMATIC/$PROCESS\n' % cat)
        datacard.write('-'*50+'\n')
        datacard.write('bin 1\n')
        datacard.write('observation %3.1f\n' % obs.Integral())
        datacard.write('-'*50+'\n')

        #expectations
        datacard.write('\t\t\t %15s'%'bin')
        for i in xrange(0,len(exp)): datacard.write('%15s'%'1')
        datacard.write('\n')
        datacard.write('\t\t\t %15s '%'process')
        for sig in signalList: datacard.write('%15s'%sig)
        for proc in exp: 
            if proc in signalList: continue
            datacard.write('%15s'%proc)
        datacard.write('\n')
        datacard.write('\t\t\t %15s'%'process')
        for i in xrange(0,len(signalList)) : datacard.write('%15s'%str(i+1-len(signalList)))
        i=0
        for proc in exp: 
            if proc in signalList: continue
            i=i+1
            datacard.write('%15s'%str(i))
        datacard.write('\n')
        datacard.write('\t\t\t %15s'%'rate')
        for sig in signalList: datacard.write('%15s'%('%3.2f'%exp[sig].Integral()))
        for proc in exp: 
            if proc in signalList: continue
            datacard.write('%15s'%('%3.2f'%exp[proc].Integral()))
        datacard.write('\n')
        datacard.write('-'*50+'\n')

        nomShapes=exp.copy()
        nomShapes['data_obs']=obs
<<<<<<< HEAD
        saveToShapesFile(outFile,nomShapes,'nom',rebin=opt.rebin)

        #experimental systematics
        _,expVarShapes=getDistsFrom(directory=fIn.Get('%sshapes_%s_exp'%(opt.dist,cat)))
        expVarShapes=filterShapeList(expVarShapes,signalList,rawSignalList)
        nExpSysts=expVarShapes[signalList[0]].GetNbinsY()/2
        shapeOnlyExpSysts=['Pileup']
        for isyst in xrange(1,nExpSysts+1):
            
            #test which variations are significant
            bwList={}
            ybin=2*(isyst-1)+1
            systVar=''
            upShapes,downShapes={},{}
=======
        saveToShapesFile(outFile,nomShapes,'nom',opt.rebin)

        #MC stats systematics for bins with large stat uncertainty
        if opt.addBinByBin>0:
>>>>>>> 9058e3898a0b728430a2f63a9c409bf903ecdb62
            for proc in exp:

                #qcd is handled separately
                if proc=='Multijetsdata' : continue

                finalNomShape=exp[proc].Clone('tmp')
                if opt.rebin>0 : finalNomShape.Rebin(opt.rebin)
                for xbin in xrange(1,finalNomShape.GetXaxis().GetNbins()+1):
                    val,unc=finalNomShape.GetBinContent(xbin),finalNomShape.GetBinError(xbin)
                    if val==0 : continue
                    if ROOT.TMath.Abs(unc/val)<opt.addBinByBin: continue

                    binShapes={}
                    systVar='%sbin%d%s'%(proc,xbin,cat)

<<<<<<< HEAD
                    #set as shape only if in the list
                    if systVar in shapeOnlyExpSysts : 
                        downShapeH.Scale(exp[proc].Integral()/downShapeH.Integral())
                        upShapeH.Scale(exp[proc].Integral()/upShapeH.Integral())

                    bwList[proc]     = acceptVariationForDataCard(nomH=exp[proc], upH=upShapeH, downH=downShapeH)
                    if not bwList[proc]: continue
=======
                    binShapes[proc]=finalNomShape.Clone('%sUp'%systVar)
                    binShapes[proc].SetBinContent(xbin,val+unc)
                    saveToShapesFile(outFile,binShapes,binShapes[proc].GetName(),False)
>>>>>>> 9058e3898a0b728430a2f63a9c409bf903ecdb62

                    binShapes[proc]=finalNomShape.Clone('%sDown'%systVar)
                    binShapes[proc].SetBinContent(xbin,ROOT.TMath.Max(val-unc,1e-3))
                    saveToShapesFile(outFile,binShapes,binShapes[proc].GetName(),False)
                    
                    #write to datacard
                    datacard.write('%32s shape'%systVar)        
                    for sig in signalList: 
                        if proc==sig:
                            datacard.write('%15s'%'1') 
                        else:
                            datacard.write('%15s'%'-')
                    for iproc in exp: 
                        if iproc in signalList: continue
                        if iproc==proc:
                            datacard.write('%15s'%'1')
                        else:
                            datacard.write('%15s'%'-')
                    datacard.write('\n')

                finalNomShape.Delete()
                    

<<<<<<< HEAD
            #export to shapes file
            saveToShapesFile(outFile,downShapes,systVar+'Down',rebin=opt.rebin)
            saveToShapesFile(outFile,upShapes,systVar+'Up',rebin=opt.rebin)
=======
        #experimental systematics
        try:
            _,expVarShapes=getDistsFrom(directory=fIn.Get('%sshapes_%s_exp'%(opt.dist,cat)))
            expVarShapes=filterShapeList(expVarShapes,signalList,rawSignalList)
            nExpSysts=expVarShapes[signalList[0]].GetNbinsY()/2
>>>>>>> 9058e3898a0b728430a2f63a9c409bf903ecdb62

            for isyst in xrange(1,nExpSysts+1):

                #test which variations are significant
                bwList={}
                ybin=2*(isyst-1)+1
                systVar=''
                upShapes,downShapes={},{}
                iRateVars={}
                for proc in exp:

                    nbins=exp[proc].GetNbinsX()
                    n=exp[proc].Integral(0,nbins+1)
                    if n==0 : continue

                    if proc in expVarShapes:
                    
                        systVarDown = expVarShapes[proc].GetYaxis().GetBinLabel(ybin)
                        systVarUp   = expVarShapes[proc].GetYaxis().GetBinLabel(ybin+1)
                        systVar     = systVarUp[:-2]
                    
                        downShapeH  = expVarShapes[proc].ProjectionX('%s%dDown'%(proc,isyst), ybin,   ybin)
                        upShapeH    = expVarShapes[proc].ProjectionX('%s%dUp'%(proc,isyst),   ybin+1, ybin+1)

                        #normalize varied shapes to nominal expectations
                        nUp=upShapeH.Integral(0,nbins+1)
                        nDn=downShapeH.Integral(0,nbins+1)
                        if nUp>0: upShapeH.Scale(n/nUp)
                        if nDn>0: downShapeH.Scale(n/nDn)

                        #save a rate systematic from the variation of the yields
                        iRateVars[proc]=(nUp/n,nDn/n)
                        if iRateVars[proc][0]<1.001 and iRateVars[proc][1]<1.001 : del iRateVars[proc]

                        bwList[proc]     = acceptVariationForDataCard(nomH=exp[proc], upH=upShapeH, downH=downShapeH)
                        if not bwList[proc]: continue
                        
                        downShapes[proc] = downShapeH
                        upShapes[proc]   = upShapeH
                    
                #test if at least one process has been white listed
                if len(upShapes)+len(downShapes)==0 and len(iRateVars)==0:
                    print '\t skipping',systVar,'for %s (no significant shapes of rate variations found)'%(anCat+cat)
                    continue
                
                if 'eeff' in systVar or 'meff' in systVar:
                    print 'will decouple %s for category %s%s'%(systVar,cat,anCat)
                    systVar += cat+anCat
            
                #export to shapes file     
                if len(upShapes)+len(downShapes)>0:
                    
                    saveToShapesFile(outFile,downShapes,systVar+'Down',opt.rebin)
                    saveToShapesFile(outFile,upShapes,systVar+'Up',opt.rebin)
                    
                    #shape variations
                    datacard.write('%32s shape'%systVar)        
                    for sig in signalList: 
                        if sig in bwList and bwList[sig]:
                            datacard.write('%15s'%'1') 
                        else:
                            datacard.write('%15s'%'-')
                    for proc in exp: 
                        if proc in signalList: continue
                        if proc in bwList and bwList[proc] :
                            datacard.write('%15s'%'1')
                        else:
                            datacard.write('%15s'%'-')
                    datacard.write('\n')
                    

                #rate variations written separately
                if len(iRateVars)>0:
                    
                    datacard.write('%32s %8s'%(systVar+'Rate','lnN'))
                    for sig in signalList:
                        if sig in iRateVars :                            
                            datacard.write('%15s'%'%3.3f/%3.3f'%(iRateVars[sig][0],iRateVars[sig][1]))
                        else:
                            datacard.write('%15s'%'-')
                    for proc in exp:
                        if proc in  signalList: continue
                        if proc in iRateVars :                            
                            datacard.write('%15s'%'%3.3f/%3.3f'%(iRateVars[proc][0],iRateVars[proc][1]))
                        else:
                            datacard.write('%15s'%'-')
                    datacard.write('\n')
                    
        except:
            pass

        #rate systematics for QCD
        try:
            jetCat=cat[:-2] if cat.endswith('t') else cat
            rateSysts.append( ('MultiJetsNorm%s%s'%(cat,anCat),  ROOT.TMath.Min(1+2*ROOT.TMath.Abs(qcdNorm[jetCat][1]),1.5),                       'lnN',    ['Multijetsdata']    ,[]) )
            #rateSysts.append( ('Wnorm%s'%jetCat,          1+ROOT.TMath.Abs(1-wjetsNorm[jetCat][0]), 'lnU',    ['Wl','Wc','Wb']     ,[]) )
        except:
            pass

        for syst,val,pdf,whiteList,blackList in rateSysts:

            datacard.write('%32s %8s'%(syst,pdf))
            entryTxt=''
            try:
                entryTxt='%15s'%('%3.3f/%3.3f'%(ROOT.TMath.Max(val[0],0.01),val[1]))
            except:
                entryTxt='%15s'%('%3.3f'%val)

            for sig in signalList: 
                if (len(whiteList)==0 and not sig in blackList) or sig in whiteList:
                    datacard.write(entryTxt)
                else:
                    datacard.write('%15s'%'-')
            for proc in exp: 
                if proc in signalList: continue
                if (len(whiteList)==0 and not proc in blackList) or proc in whiteList:
                    datacard.write(entryTxt)
                else:
                    datacard.write('%15s'%'-')
            datacard.write('\n')

        #generator level systematics 
        if systfIn is None: continue
<<<<<<< HEAD
        sampleSysts=[

            #ttbar modelling
            ('Mtop',            {'tbart'         : ['tbartm=169.5','tbartm=175.5'],  'tW':['tWm=169.5','tWm=175.5'] },                True ,  True, False),
            ('ttPartonShower',  {'tbart'         : ['tbartscaledown','tbartscaleup']},                                                False , True, False),            
            ('NLOgenerator',    {'tbart'         : ['tbartaMCNLO']},                                                                  False,  True, False),
            ('Hadronizer',      {'tbart'         : ['tbartHerwig']},                                                                  False , True, True),

            #tWinterference
            ('tWttinterf',       {'tW'            : ['tWDS']},                                                                        False , True, True),            

            #QCD SCALES
            ('tWscale',         {'tW'            : ['tWscaledown','tWscaleup']},                                                      False , True, False),            

            #Madgraph W+jets
            #('wFactScale',           { 'Wl': ['id3mur1muf0.5','id2mur1muf2'], 
            #                           'Wc': ['id3mur1muf0.5','id2mur1muf2'], 
            #                           'Wb': ['id3mur1muf0.5','id2mur1muf2'] },  False, False),
            #('wRenScale',            { 'Wl': ['id7mur0.5muf1','id4mur2muf1'],  
            #                           'Wc': ['id7mur0.5muf1','id4mur2muf1'],  
            #                           'Wb': ['id7mur0.5muf1','id4mur2muf1'] },  False, False),
            #('wCombScale',           { 'Wl': ['id9mur0.5muf0.5','id5mur2muf2'], 
            #                           'Wc': ['id9mur0.5muf0.5','id5mur2muf2'],
            #                           'Wb': ['id9mur0.5muf0.5','id5mur2muf2'] },  False, False),

            #amc@NLO W+jets
            ('wFactScale',           { 'W': ['id1003muR0.10000E+01muF0.50000E+00','id1002muR0.10000E+01muF0.20000E+01'] },    False, False, False),
            ('wRenScale',            { 'W': ['id1007muR0.50000E+00muF0.10000E+01','id1004muR0.20000E+01muF0.10000E+01'] },    False, False, False),
            ('wCombScale',           { 'W': ['id1009muR0.50000E+00muF0.50000E+00','id1005muR0.20000E+01muF0.20000E+01'] },  False, False, False),
            
            #ttbar Powheg
            ('ttFactScale',          { 'tbart': ['muR1muF0.5hdampmt172.5',  'muR1muF2hdampmt172.5'] },     True , False, False),
            ('ttRenScale',           { 'tbart': ['muR0.5muF1hdampmt172.5',  'muR2muF1hdampmt172.5'] },     True , False, False),
            ('ttCombScale',          { 'tbart': ['muR0.5muF0.5hdampmt172.5','muR2muF2hdampmt172.5'] },     True , False, False),
            ]
=======
>>>>>>> 9058e3898a0b728430a2f63a9c409bf903ecdb62

        _,genVarShapes = getDistsFrom(directory=fIn.Get('%sshapes_%s_gen'%(opt.dist,cat)))
        genVarShapes=filterShapeList(genVarShapes,signalList,rawSignalList)
        _,altExp       = getDistsFrom(directory=systfIn.Get('%s_%s'%(opt.dist,cat)))        
        if signalList[0]!=rawSignalList[0]:
            altExp=filterShapeList(altExp,signalList,rawSignalList)
        for systVar, procsToApply, normalize, useAltShape, projectRelToNom in sampleSysts:

            #prepare shapes and check if variation is significant
            downShapes, upShapes = {}, {}
            
            for iproc in procsToApply:

                nomH=exp[iproc]

                #check which shape to use
                if useAltShape:

                    #get directly from another file
                    downH  = altExp[ procsToApply[iproc][0] ]
                    if len( procsToApply[iproc] ) > 1 :
                        upH    = altExp[ procsToApply[iproc][1] ]
                    else:
                        #if only one variation is available, mirror it
                        upH = downH.Clone( '%s%sUp'%(iproc,systVar) )
                        for xbin in xrange(1,upH.GetNbinsX()+1):
                            diff=upH.GetBinContent(xbin)-nomH.GetBinContent(xbin)
                            upH.SetBinContent(xbin,nomH.GetBinContent(xbin)-diff)
                else:

                    #project from 2D histo (re-weighted from nominal sample)
                    ybinUp, ybinDown = -1, -1
                    for ybin in xrange(1,genVarShapes[ iproc ].GetNbinsY()+1):
                        label = genVarShapes[ iproc ].GetYaxis().GetBinLabel(ybin)
                        if procsToApply[iproc][0] in label and ybinDown<0 : ybinDown=ybin
                        if procsToApply[iproc][1] in label and ybinUp<0   : ybinUp=ybin

                    downH = genVarShapes[ iproc ].ProjectionX('%s%sDown'%(iproc,systVar), ybinDown, ybinDown)
                    upH   = genVarShapes[ iproc ].ProjectionX('%s%sUp'%(iproc,systVar),   ybinUp,   ybinUp)

                # use do down/up x nom to generate the variation, then mirror it
                if projectRelToNom:
                    ratioH=downH.Clone()
                    ratioH.Divide(upH)
                    for xbin in xrange(1,nomH.GetNbinsX()+1):
                        nomVal=nomH.GetBinContent(xbin)
                        varVal = ratioH.GetBinContent(xbin) * nomVal
                        upH.SetBinContent(xbin, varVal)
                        varVal = varVal- nomVal
                        downH.SetBinContent(xbin, nomVal-varVal)

                #normalize (shape only variation is considered)
                if normalize : downH.Scale( nomH.Integral()/downH.Integral() ) 
                if normalize : upH.Scale( nomH.Integral()/upH.Integral() )

                #check if variation is meaningful
                accept = acceptVariationForDataCard(nomH=nomH, upH=upH, downH=downH)
                if not accept : continue
                
                #save
                downShapes[iproc]=downH
                upShapes[iproc]=upH

            #check if something has been accepted
            if len(upShapes)==0 : continue

            #export to shapes file
<<<<<<< HEAD
            saveToShapesFile(outFile,downShapes,systVar+'Down',rebin=opt.rebin)
            saveToShapesFile(outFile,upShapes,systVar+'Up',rebin=opt.rebin)
=======
            saveToShapesFile(outFile,downShapes,systVar+'Down',opt.rebin)
            saveToShapesFile(outFile,upShapes,systVar+'Up',opt.rebin)
>>>>>>> 9058e3898a0b728430a2f63a9c409bf903ecdb62

            #write to datacard
            datacard.write('%32s shape'%systVar)
            for sig in signalList: 
                if sig in procsToApply and sig in upShapes:
                    datacard.write('%15s'%'1')
                else:
                    datacard.write('%15s'%'-')
            for proc in exp: 
                if proc in signalList: continue
                if proc in procsToApply and proc in upShapes:
                    datacard.write('%15s'%'1')
                else:
                    datacard.write('%15s'%'-')
            datacard.write('\n')

        #
        # QCD shapes
        # 
        systName='MultiJetsShape%s%s'%(cat,anCat)
        qcdExp=exp['Multijetsdata'].Integral()
        if qcdExp>0 :
            datacard.write('%32s shape'%systName)
            for sig in signalList: 
                if (len(whiteList)==0 and not sig in blackList) or sig in whiteList:
                    datacard.write(entryTxt)
                else:
                    datacard.write('%15s'%'-')
            for proc in exp: 
                if proc in signalList: continue
                if proc=='Multijetsdata':
                    datacard.write('%15s'%'1')
                else:
                    datacard.write('%15s'%'-')
            datacard.write('\n')

            _,qcdShapesUp = getDistsFrom(directory=fIn.Get('%s_%s_QCD%sUp'%(opt.dist,cat,jetCat)))
            qcdShapesUp['Multijetsdata'].Scale( qcdExp/qcdShapesUp['Multijetsdata'].Integral() ) 
<<<<<<< HEAD
            saveToShapesFile(outFile,qcdShapesUp,systName+'Up',rebin=opt.rebin)

            _,qcdShapesDown = getDistsFrom(directory=fIn.Get('%s_%s_QCD%sDown'%(opt.dist,cat,jetCat)))
            qcdShapesDown['Multijetsdata'].Scale( qcdExp/qcdShapesDown['Multijetsdata'].Integral() ) 
            saveToShapesFile(outFile,qcdShapesDown,systName+'Down',rebin=opt.rebin)
=======
            saveToShapesFile(outFile,qcdShapesUp,systName+'Up',opt.rebin)

            _,qcdShapesDown = getDistsFrom(directory=fIn.Get('%s_%s_QCD%sDown'%(opt.dist,cat,jetCat)))
            qcdShapesDown['Multijetsdata'].Scale( qcdExp/qcdShapesDown['Multijetsdata'].Integral() ) 
            saveToShapesFile(outFile,qcdShapesDown,systName+'Down',opt.rebin)
>>>>>>> 9058e3898a0b728430a2f63a9c409bf903ecdb62

        #all done
        datacard.close()


"""                                                                                                                                                                                                               
for execution from another script                                                                                                                                                                           
"""
if __name__ == "__main__":
    sys.exit(main())
