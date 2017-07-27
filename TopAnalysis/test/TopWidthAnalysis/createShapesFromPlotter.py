import os
import sys
import optparse
import ROOT
import commands
import getpass
import numpy


ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
#ROOT.gROOT.SetBatch()

def replaceBadCharacters(inputStr):
    newStr = inputStr
    for token in ['+','-','*',' ','#','{','(',')','}','@']:
        newStr=newStr.replace(token,'')
    return newStr

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
def getDistsFrom(directory,addList=None,addTitle='',keyFilter=''):
    obs=None
    exp={}
    dirName=directory.GetName()
    addedHistProduced=False
    for key in directory.GetListOfKeys():
        if len(keyFilter)>0 and key.GetName()!='%s_%s'%(dirName,keyFilter) : continue
        obj=directory.Get(key.GetName())
        if not obj.InheritsFrom('TH1') : continue
        inAddList = False;
        if obj.GetName()==dirName :
            obs=obj.Clone('data_obs')
            obs.SetDirectory(0)
        else :
            newName=obj.GetName().split(dirName+'_')[-1]
            for token in ['+','-','*',' ','#','{','(',')','}','@']:
                newName=newName.replace(token,'')
            if addList and newName in addList :
                newName='added' if addTitle == '' else addTitle
                if addedHistProduced :
                    exp[newName].Add(obj)
                    addedHistProduced = True;
                    inAddList=True
                else :
                    exp[newName] = obj.Clone(newName)
            else :
                exp[newName]=obj.Clone(newName)
                exp[newName].SetDirectory(0)

            #newName = ('added' if addTitle == '' else addTitle) if inAddList and addedHistProduced else newName
            if obj.InheritsFrom('TH2'):
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
merge many directories with the same directory base-name
  '.' --> 'p'
"""
def getMergedDists(fIn,basedir='',addList=None,addTitle='',
                   widList=None,nomWid='',sigList=None,keyFilter=''):
    obs=None
    exp={}
    fIn.Get._creates=True
    for wid in widList :
        dirName='%s%s'%(basedir,wid)
        #print dirName
        directory = fIn.Get(dirName)
        if wid == nomWid :
            obs,mergeExp=getDistsFrom(directory=directory,addList=addList,addTitle=addTitle,keyFilter=keyFilter)
            for key in mergeExp :
                newName=key
                if key in sigList :
                    newName=('%s%s'%(key,wid.replace('.','p')))
                exp[newName]=mergeExp[key].Clone(newName)
        else :
            _,mergeExp=getDistsFrom(directory=directory,addList=addList,addTitle=addTitle,keyFilter=keyFilter)
            for key in mergeExp :
                if key in sigList :
                    newName=('%s%s'%(key,wid.replace('.','p')))
                    exp[newName]=mergeExp[key].Clone(newName)
                else : continue
    return obs,exp


"""
save distributions to file
"""
def saveToShapesFile(outFile,shapeColl,directory=''):
    fOut=ROOT.TFile.Open(outFile,'UPDATE')
    if len(directory)==0:
        fOut.cd()
    else:
        if not fOut.Get(directory):
            fOut.mkdir(directory)
        outDir=fOut.Get(directory)
        outDir.cd()
    for key in shapeColl:
        #remove bin labels
        shapeColl[key].GetXaxis().Clear()

        #convert to TH1D (projections are TH1D)
        if not shapeColl[key].InheritsFrom('TH1D') :
            h=ROOT.TH1D()
            shapeColl[key].Copy(h)
            shapeColl[key]=h

        shapeColl[key].Write(key,ROOT.TObject.kOverwrite)
    fOut.Close()

"""
make an MC truth dataset (sigs+bkgs)
"""
def makeMCTruthHist(hypothesis,sigList,dists,infile=None,dirhandle="",filtername="",exthypo=""):
    outputHist=None
    firstLoop=True
    print "Producing MC Truth data for hypothesis %s"%hypothesis

    tdists=dists.copy()

    #replace with whatever we find in the provided dir
    if dirhandle != "" and (infile is not None) and filtername != "" and exthypo!="":
        _,pdists=getMergedDists(infile,dirhandle,None,"",[exthypo.replace('p','.')],"1.0w",sigList)
        for key in pdists :
            if filtername not in key : continue
            tdists[key.replace(filtername,'')+hypothesis] = pdists[key].Clone()
            tdists[key.replace(filtername,'')+hypothesis].SetDirectory(0)

    for sig in sigList :
        sigHist=tdists["%s%s"%(sig,hypothesis)].Clone()
        if firstLoop :
            outputHist=sigHist.Clone()
            firstLoop=False
        else :
            outputHist.Add(sigHist)
    for dist in tdists :
        isSig=False
        for sig in sigList :
            if sig in dist and 'V' not in dist: isSig=True
        if isSig : continue
        bkgHist=tdists[dist].Clone()
        outputHist.Add(bkgHist)
        outputHist.SetName("data_obs")

    outputHist.SetDirectory(0)

    return outputHist

"""
steer the script
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--input',     dest='input',           help='input plotter',                            default=None,               type='string')
    parser.add_option(      '--systInput', dest='systInput',       help='input plotter for systs from alt samples', default=None,               type='string')
    parser.add_option('-s', '--signal',    dest='signal',          help='signal (csv)',                             default='tbart,Singletop',  type='string')
    parser.add_option('-c', '--cat',       dest='cat',             help='categories (csv)',                         default='1b,2b',            type='string')
    parser.add_option('-o', '--output',    dest='output',          help='output directory',                         default='datacards',        type='string')
    parser.add_option('-n', '--outname',   dest='outname',         help='output file name',                         default='shapes',           type='string')
    parser.add_option(      '--addSigs',   dest='addSigs',         help='signal processes to add',                  default=False,        action='store_true')
    parser.add_option(      '--lfs',       dest='lfsInput',        help='lepton final states to consider',          default='EE,EM,MM',         type='string')
    parser.add_option(      '--lbCat',     dest='lbCat',           help='pt categories to consider',                default='highpt,lowpt',     type='string')
    parser.add_option(      '--truth',     dest='truthDataset',    help='make data out of MC truth',                default='',                 type='string')
    parser.add_option(      '--trx4',      dest='truthExtDataset', help='make data out of 4xMC truth',              default='',                 type='string')
    parser.add_option(      '--trx4Width', dest='truthExtWidth',   help='4xMC truth has this width',                default='',                 type='string')
    parser.add_option(      '--trx4Input', dest='truthExtfIn'  ,   help='4xMC truth has this infile',               default='',                 type='string')
    parser.add_option(      '--min',       dest='minIndex',        help='start at this card',                       default=-5,                 type=int)
    parser.add_option(      '--max',       dest='maxIndex',        help='end at this card',                         default=10000,              type=int)
    parser.add_option('-d', '--dists',     dest='distList',        help='distribution',                             default='incmlb',           type='string')
    parser.add_option('-w', '--wids',      dest='widList',         help='signal widths',
            default='0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.5,4.0',
            type='string')
    (opt, args) = parser.parse_args()

    # parse the dists to consider
    distList = opt.distList.split(',')

    # parse the channels, pt(l,b), b categories to consider
    lfsList=opt.lfsInput.split(',')
    lbCatList=opt.lbCat.split(',')
    catList=opt.cat.split(',')

    # what are our signal processes?
    rawSignalList=opt.signal.split(',')
    signalList = { 'top' } if opt.addSigs else rawSignalList[:]

    # what are our widths?
    rawWidList = opt.widList.split(',')

    # widths with '...w'
    widList = rawWidList[:]
    for i in xrange(0,len(rawWidList)) :
        widList[i] = ('%sw'%(rawWidList[i]))

    # widths with '.p.w'
    modWidList=widList[:]
    for i in xrange(0,len(widList)) :
        modWidList[i]=widList[i].replace('.','p')

    # helper nominal widths
    nomWid="1.0w"
    modNomWid="1p0w"

    # helper functions
    def modWidIn(tproc) :
        for tmwid in modWidList :
            if tmwid in proc: return True
        return False

    # HARDCODE systematics
    tWRateSystList=['Singletop']
    ttRateSystList=['tbart']
    tWIsSig=('Singletop'    in rawSignalList)
    ttIsSig=('tbart' in rawSignalList)
    for wid in modWidList :
        if tWIsSig : tWRateSystList += ["Singletop%s"%wid]
        if ttIsSig : ttRateSystList += ["tbart%s"%wid]

    rateSysts=[
          ('lumi_13TeV',       1.062,    'lnN',    []                   ,['Multijetsdata']),
          ('DYnorm_th',        1.30,     'lnN',    ['DY']  ,[]),
          ('Wnorm_th',         1.30,     'lnN',    ['W']   ,[]),
          ('tWnorm_th',        1.054,    'lnN',    tWRateSystList,[]),
          ('tnorm_th',         1.044,    'lnN',    ['tch']              ,[]),
          ('VVnorm_th',        1.20,     'lnN',    ['Multiboson']       ,[]),
          ('tbartVnorm_th',    1.30,     'lnN',    ['tbartV']           ,[])
          #('sel',              1.02,     'lnN',    tWRateSystList+ttRateSystList+['DY','W','tch','Multiboson','tbartV'], []),
    ]

    Mtop       =['tbartm=171.5','tbartm=173.5']#,'Singletopm=171.5','Singletopm=173.5']
    #ttParton_tt=['tbartscaledown','tbartscaleup']
    #ttParton_tW=['tWscaledown','tWscaleup']
    #ME={    "muF": ['%sgen%i'%(sig,ind) for sig in ['tbart'] for ind in [3,4]],
    #        "muR": ['%sgen%i'%(sig,ind) for sig in ['tbart'] for ind in [5,8]],
    #        "tot": ['%sgen%i'%(sig,ind) for sig in ['tbart'] for ind in [6,10]]
    #        }
    #tWinterf=['tWDS']
    #amcNLO  =['tbartaMCNLO']
    #Herwig  =['tbartHerwig']
    ISR     =['tbartisrdn','tbartisrup']#,'Singletopisrdn','Singletopisrup']
    FSR     =['tbartfsrdn','tbartfsrup']#,'Singletopfsrdn','Singletopfsrup']
    hdamp   =['tbarthdampdn','tbarthdampup']#,'Singletophdampdn','Singletophdampup']
    undEve  =['tbartUEdn','tbartUEup']#,'SingletopUEdn','SingletopUEup']
    #cflip   =['tbartcflip']
    #noSC    =['tbartnoSC']

    systSignalList=Mtop+ISR+FSR+hdamp+undEve

    MtopMap={}
    #ttPartonMap={}
    #twPartonMap={}
    #tWinterfMap={}
    #amcNLOMap={}
    #HerwigMap={}
    ISRMap={}
    FSRMap={}
    hdampMap={}
    UEMap={}
    #cflMap={}
    #noSCMap={}

    for wid in modWidList :
        MtopMap[    'tbart%s'%wid]=['tbartm=171.5%s'    %(wid),'tbartm=173.5%s'%(wid)]
        #MtopMap['Singletop%s'%wid]=['Singletopm=171.5%s'%(wid),'Singletopm=173.5%s'%(wid)]

        #ttPartonMap['tbart%s'%wid]=['tbartscaledown%s'%(wid),'tbartscaleup%s'%(wid)]
        #twPartonMap['tW%s'   %wid]=['tWscaledown%s'   %(wid),'tWscaleup%s'   %(wid)]

        #amcNLOMap['tbart%s'%wid]=['tbartaMCNLO%s'%(wid)]
        #HerwigMap['tbart%s'%wid]=['tbartHerwig%s'%(wid)]
        ISRMap['tbart%s'%wid]=['tbartisrdn%s'%(wid),'tbartisrup%s'%(wid)]
        #ISRMap['Singletop%s'%wid]=['Singletopisrdn%s'%(wid),'Singletopisrup%s'%(wid)]

        FSRMap['tbart%s'%wid]=['tbartfsrdn%s'%(wid),'tbartfsrup%s'%(wid)]
        #FSRMap['Singletop%s'%wid]=['Singletopfsrdn%s'%(wid),'Singletopfsrup%s'%(wid)]

        hdampMap['tbart%s'%wid]=['tbarthdampdn%s'%(wid),'tbarthdampup%s'%(wid)]
        #hdampMap['Singletop%s'%wid]=['Singletophdampdn%s'%(wid),'Singletophdampup%s'%(wid)]

        UEMap['tbart%s'%wid]=['tbartUEdn%s'%(wid),'tbartUEup%s'%(wid)]
        #UEMap['Singletop%s'%wid]=['SingletopUEdn%s'%(wid),'SingletopUEup%s'%(wid)]

        #tWinterfMap['tW%s'   %wid]=['tWDS%s'   %(wid)]

    # systVar, procsToApply, normalize, useAltShape, projectRelToNom
    sampleSysts=[
          #ttbar modelling
          ('Mtop'          , MtopMap    , True , True, False),
          #('ttPartonShower', ttPartonMap, False, True, False),
          #('tWQCDScale'    , twPartonMap, False, True, False),
          #('Herwig'     , HerwigMap , False, True, True),
          #('amcnlo'     , amcNLOMap , False, True, True),
          ('ISR'        , ISRMap    , False, True, False),
          ('FSR'        , FSRMap    , False, True, False),
          ('hdamp'      , hdampMap  , False, True, False),
          ('UE'         , UEMap     , False, True, False),
          #('cflip'      , cflMap    , False, True, True),
          #('noSC'       , noSCMap   , False, True, True)
          #tWinterference
          #('tWttinterf'    , tWinterfMap, False, True, True),
    ]

    # systVar, normalize, useAltShape, projectRelToNom
    genSysts=[
        ('jes',         True,   False,  False),
        ('les',         False,  False,  False),
        ('ltag',        False,  False,  False),
        ('trig',        False,  False,  False),
        ('sel',         False,  False,  False),
        ('toppt',       True,   False,  False),
        ('jer',         False,  False,  False),
        ('btag',        False,  False,  False),
        ('pu',          True,   False,  False)
        #('MEqcdscale',  True,   False,  False),
        #('PDF',         False,  False,  False)
        #('nloproddec',  True,   False,  True)
    ]

    # prepare output directory
    os.system('mkdir -p %s'%opt.output)

    # keep track of analysis name
    anCat=''
    for subDir in opt.input.split('/'):
        if 'analysis_' not in subDir: continue
        anCat=subDir.replace('analysis_','')

    # get data and nominal expectations
    fIn=ROOT.TFile.Open(opt.input)
    fIn.Get._creates=True
    systfIn=None
    if opt.systInput:
        systfIn=ROOT.TFile.Open(opt.systInput)
        systfIn.Get._creates=True

    # prepare output ROOT file
    outFile='%s/%s.root'%(opt.output,opt.outname)
    fOut=ROOT.TFile.Open(outFile,'UPDATE')
    fOut.Close()

    # keep track of progress
    cardIndex=-1
    numCards =(len(rawWidList))*len(lfsList)*len(lbCatList)*len(catList)*len(distList)

    #loop over categories
    for (lbCat,lfs,cat,dist) in [(a,b,c,d) for a in lbCatList
            for b in lfsList
            for c in catList
            for d in distList] :

        # keep track of indices
        if cardIndex >= opt.maxIndex : break
        if cardIndex+len(rawWidList) < opt.minIndex:
            cardIndex+=len(rawWidList)
            continue

        obs=ROOT.TH1F('','',100,0,100)
        exp={}
        #nominal expectations
        if opt.addSigs :
            obs,exp=getMergedDists(fIn,('%s%s%s_%s_'%(lbCat,lfs,cat,dist)),rawSignalList,'top',widList,nomWid,signalList)
        else :
            obs,exp=getMergedDists(fIn,('%s%s%s_%s_'%(lbCat,lfs,cat,dist)),None,'',widList,nomWid,signalList)

        # inject MC truth if so desired
        if not opt.truthDataset=="" :
            obs=makeMCTruthHist(opt.truthDataset,signalList,exp,
                    (ROOT.TFile(opt.truthExtfIn) if opt.truthExtfIn != "" else None),
                    ('%s%s%s_%s_'%(lbCat,lfs,cat,dist)),
                    opt.truthExtDataset,opt.truthExtWidth)
        #exp=filterShapeList(exp,signalList,rawSignalList)

        nomShapes=exp.copy()
        nomShapes['data_obs']=obs

        for tSig in signalList :
            nomShapes['%sNOM'%(tSig)] = exp['%s%s'%(tSig,modNomWid)].Clone()

        saveToShapesFile(outFile,nomShapes,('%s%s%s_%s'%(lbCat,lfs,cat,dist)))

        #loop over categories, widths
        for wid in modWidList:

            # skip cards if not within index
            cardIndex+=1
            if cardIndex < opt.minIndex :
                print 'Skipping %s datacards for %s%s%s_%s \t\t [%i/%i]'%(dist,lbCat,lfs,cat,wid,cardIndex,numCards)
                continue
            if cardIndex > opt.maxIndex :
                break
            print 'Initiating %s datacard for %s%s%s_%s \t\t [%i/%i]'%(dist,lbCat,lfs,cat,wid,cardIndex,numCards)


            #start the datacard
            datacard=open('%s/datacard__%s_%s%s%s_%s.dat'%(opt.output,wid,lbCat,lfs,cat,dist),'w')
            datacard.write('#\n')
            datacard.write('# Generated by %s with git hash %s for analysis category %s%s%s_%s\n' % (getpass.getuser(),
                commands.getstatusoutput('git log --pretty=format:\'%h\' -n 1')[1],
                lbCat, lfs, cat, wid) )

            datacard.write('#\n')
            datacard.write('imax *\n')
            datacard.write('jmax *\n')
            datacard.write('kmax *\n')
            datacard.write('-'*50+'\n')
            datacard.write('shapes *        * shapes.root %s%s%s_%s/$PROCESS %s%s%s_%s_$SYSTEMATIC/$PROCESS\n'%(lbCat,lfs,cat,dist,lbCat,lfs,cat,dist))
            datacard.write('-'*50+'\n')
            datacard.write('bin 1\n')
            datacard.write('observation %3.1f\n' % obs.Integral())
            datacard.write('-'*50+'\n')

            #expectations
            datacard.write('\t\t\t %32s'%'bin')
            for i in xrange(0,len(exp)+(2 if opt.addSigs else 0)-(len(modWidList)-1)*len(signalList)+len(signalList)):
                datacard.write('%15s'%'1')
            datacard.write('\n')
            datacard.write('\t\t\t %32s'%'process')
            for sig in signalList: datacard.write('%15s'%('%s%s'%(sig,modNomWid)))
            if modNomWid == wid :
                for sig in signalList: datacard.write('%15s'%('%sNOM'%(sig)))
            else :
                for sig in signalList: datacard.write('%15s'%('%s%s'%(sig,wid)))
            for proc in exp:
                isSig=modWidIn(proc)
                if isSig : continue
                datacard.write('%15s'%proc)
            datacard.write('\n')
            datacard.write('\t\t\t %32s'%'process')
            for i in xrange(0,2*len(signalList)) : datacard.write('%15s'%str(i+1-2*len(signalList)))
            i=0
            for proc in exp:
                isSig=modWidIn(proc)
                if isSig : continue
                i=i+1
                datacard.write('%15s'%str(i))
            datacard.write('\n')
            datacard.write('\t\t\t %32s'%'rate')
            for sig in signalList :
                datacard.write('%15s'%('%3.2f'%(exp['%s%s'%(sig,modNomWid)].Integral())))
            for sig in signalList :
                datacard.write('%15s'%('%3.2f'%(exp['%s%s'%(sig,wid)].Integral())))
            for proc in exp:
                isSig=modWidIn(proc)
                if isSig : continue
                datacard.write('%15s'%('%3.2f'%exp[proc].Integral()))
            datacard.write('\n')
            datacard.write('-'*50+'\n')


            #rate systematics
            try:
                jetCat=cat[:-2] if cat.endswith('t') else cat
                rateSysts.append( ('MultiJetsNorm%s%s'%(jetCat,anCat),
                                    1+qcdNorm[jetCat][1],
                                    'lnN',
                                    ['Multijetsdata'],
                                    []) )
            except:
                pass

            for syst,val,pdf,whiteList,blackList in rateSysts:
                if syst in ['sel','DYnorm_th'] : syst="%s%s"%(syst,lfs)


                datacard.write('%32s %8s'%(syst,pdf))
                entryTxt=''
                try:
                    entryTxt='%15s'%('%3.3f/%3.3f'%(ROOT.TMath.Max(val[0],0.01),val[1]))
                except:
                    entryTxt='%15s'%('%3.3f'%val)

                for sig in signalList:
                    newSig='%s%s'%(sig,modNomWid)
                    if (len(whiteList)==0 and not newSig in blackList) or newSig in whiteList:
                        datacard.write(entryTxt)
                    else:
                        datacard.write('%15s'%'-')
                for sig in signalList:
                    newSig='%s%s'%(sig,wid)
                    if (len(whiteList)==0 and not newSig in blackList) or newSig in whiteList:
                        datacard.write(entryTxt)
                    else:
                        datacard.write('%15s'%'-')
                for proc in exp:
                    isSig=modWidIn(proc)
                    if isSig : continue
                    if (len(whiteList)==0 and not proc in blackList) or proc in whiteList:
                        datacard.write(entryTxt)
                    else:
                        datacard.write('%15s'%'-')
                datacard.write('\n')

            #generator level systematics
            if systfIn is None :
                datacard.close()
                continue


            # sample systematics
            _,genVarShapes = getMergedDists(fIn,
                    ('%s%s%s_%s_'%(lbCat,lfs,cat,dist)),
                    (rawSignalList if opt.addSigs else None),
                    ('top' if opt.addSigs else ''),
                    widList,nomWid,systSignalList)
            _,altExp       = getMergedDists(systfIn,
                    ('%s%s%s_%s_'%(lbCat,lfs,cat,dist)),
                    (rawSignalList if opt.addSigs else None),
                    ('top' if opt.addSigs else ''),
                    widList,nomWid,systSignalList)

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
                                if nomH.GetBinContent(xbin)-diff >= 0 :
                                    upH.SetBinContent(xbin,nomH.GetBinContent(xbin)-diff)
                                else :
                                    upH.SetBinContent(xbin,upH.GetBinContent(xbin-1))
                    else:

                        #project from 2D histo (re-weighted from nominal sample)
                        ybinUp, ybinDown = -1, -1
                        for ybin in xrange(1,genVarShapes[ iproc ].GetNbinsY()+1):
                            label = genVarShapes[ iproc ].GetYaxis().GetBinLabel(ybin)
                            if procsToApply[iproc][0] in label : ybinDown=ybin
                            if procsToApply[iproc][1] in label : ybinUp=ybin

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
                            if upH.GetBinContent(xbin) <= 0 :
                                upH.SetBinContent(xbin,1e-010)
                            varVal = varVal- nomVal
                            downH.SetBinContent(xbin, nomVal-varVal)
                            if downH.GetBinContent(xbin) <= 0 :
                                downH.SetBinContent(xbin,1e-010)

                    #normalize (shape only variation is considered)
                    if normalize : downH.Scale( nomH.Integral()/downH.Integral() )
                    if normalize : upH.Scale( nomH.Integral()/upH.Integral() )

                    #check if variation is meaningful
                    accept = acceptVariationForDataCard(nomH=nomH, upH=upH, downH=downH)
                    if not accept : continue

                    #save
                    downShapes[iproc]=downH
                    upShapes[iproc]=upH

                    if modNomWid in iproc:
                        for tSig in signalList :
                            downShapes['%sNOM'%tSig]=downH
                            upShapes['%sNOM'%tSig]=upH

                #check if something has been accepted
                if len(upShapes)==0 : continue

                #export to shapes file
                saveToShapesFile(outFile,downShapes,lbCat+lfs+cat+"_"+dist+"_"+systVar+'Down')
                saveToShapesFile(outFile,upShapes,lbCat+lfs+cat+"_"+dist+"_"+systVar+'Up')

                #write to datacard
                datacard.write('%26s shape'%systVar)
                for sig in signalList:
                    if ("%s%s"%(sig,wid)) in procsToApply and ("%s%s"%(sig,wid)) in upShapes:
                        if "Mtop" in systVar:
                            datacard.write('%15s'%'0.5')
                        else:
                            datacard.write('%15s'%'1')
                    else:
                        datacard.write('%15s'%'-')
                for sig in signalList:
                    if ("%s%s"%(sig,wid)) in procsToApply and ("%s%s"%(sig,wid)) in upShapes:
                        if "Mtop" in systVar:
                            datacard.write('%15s'%'0.5')
                        else:
                            datacard.write('%15s'%'1')
                    else:
                        datacard.write('%15s'%'-')
                for proc in exp:
                    isSig=modWidIn(proc)
                    if isSig : continue
                    if proc in procsToApply and proc in upShapes:
                        if "Mtop" in systVar:
                            datacard.write('%15s'%'0.5')
                        else:
                            datacard.write('%15s'%'1')
                    else:
                        datacard.write('%15s'%'-')
                datacard.write('\n')


            #
            for systVar, normalize, useAltShape, projectRelToNom in genSysts:

                genVarShapesUp,genVarShapesDn=None,None

                if 'PDF'==systVar:

                    #fill in 2D histos y-axis is the replica number
                    genVarShapesDn2D={}
                    nPdfReplicas=100
                    startIdx=11
                    for igen in xrange(startIdx,nPdfReplicas+startIdx+1):
                        dirForSyst=('gen%d%s%s%s_%s_'%(igen,lbCat,lfs,cat,dist))
                        _,iGenVarShapesDn = getMergedDists((systfIn if useAltShape else fIn),
                                                           dirForSyst,
                                                           (rawSignalList if opt.addSigs else None),
                                                           ('top' if opt.addSigs else ''),
                                                           widList,nomWid,signalList)
                        for proc in iGenVarShapesDn:
                            nbins=iGenVarShapesDn[proc].GetXaxis().GetNbins()
                            if not proc in genVarShapesDn2D:
                                name=iGenVarShapesDn[proc].GetName()+'2D'
                                title=iGenVarShapesDn[proc].GetTitle()
                                xmin=iGenVarShapesDn[proc].GetXaxis().GetXmin()
                                xmax=iGenVarShapesDn[proc].GetXaxis().GetXmax()
                                genVarShapesDn2D[proc]=ROOT.TH2D(name,title,nbins,xmin,xmax,nPdfReplicas,0,nPdfReplicas)
                                genVarShapesDn2D[proc].SetDirectory(0)
                            for xbin in xrange(1,nbins+1):
                                genVarShapesDn2D[proc].SetBinContent(xbin,igen-startIdx+1,iGenVarShapesDn[proc].GetBinContent(xbin))

                    #fill the down shapes with the bin contents being the mean-RMS of the replicas
                    genVarShapesDn={}
                    for proc in genVarShapesDn2D:
                        #down shape
                        name=genVarShapesDn2D[proc].GetName()[:-2]+"PDF"
                        genVarShapesDn[proc]=genVarShapesDn2D[proc].ProjectionX(name,1,1)
                        genVarShapesDn[proc].SetDirectory(0)
                        genVarShapesDn[proc].Reset('ICE')

                        #project to compute the mean and RMS
                        for xbin in xrange(1,genVarShapesDn2D[proc].GetXaxis().GetNbins()+1):
                            tmp=genVarShapesDn2D[proc].ProjectionY("tmp",xbin,xbin)

                            tmean=numpy.zeros(tmp.GetXaxis().GetNbins())
                            for txbin in xrange(1,tmp.GetXaxis().GetNbins()+1) :
                                tmean[txbin-1]=tmp.GetBinContent(txbin)

                            mean=numpy.mean(tmean)
                            rms=numpy.std(tmean)

                            tmp.Delete()
                            genVarShapesDn[proc].SetBinContent(xbin,mean-rms)

                        #all done, can remove this histo from memory
                        genVarShapesDn2D[proc].Delete()

                else:

                    genVarShapesUp,genVarShapesDn={},{}
                    #construct an envelope from the 6 QCD scale variations, then mirror it
                    if 'MEqcdscale' in systVar:
                        for igen in [3,5,6,4,8,10]:
                            dirForSyst=('gen%d%s%s%s_%s_'%(igen,lbCat,lfs,cat,dist))
                            _,igenVarShapes = getMergedDists((systfIn if useAltShape else fIn),
                                                             dirForSyst,
                                                             (rawSignalList if opt.addSigs else None),
                                                             ('top' if opt.addSigs else ''),
                                                             widList,nomWid,signalList)
                            for proc in igenVarShapes:
                                if proc not in genVarShapesDn:
                                    name=igenVarShapes[proc].GetName()+'envdn'
                                    genVarShapesDn[proc]=exp[proc].Clone(name)
                                    genVarShapesDn[proc].SetDirectory(0)
                                    name=igenVarShapes[proc].GetName()+'envup'
                                    genVarShapesUp[proc]=exp[proc].Clone(name)
                                    genVarShapesUp[proc].SetDirectory(0)
                                for xbin in xrange(1,igenVarShapes[proc].GetNbinsX()+1):
                                    yvar = igenVarShapes[proc].GetBinContent(xbin)
                                    ydn  = genVarShapesDn[proc].GetBinContent(xbin)
                                    yup  = genVarShapesUp[proc].GetBinContent(xbin)
                                    if yvar>yup:
                                        genVarShapesUp[proc].SetBinContent(xbin,yvar)
                                    elif yvar<ydn:
                                        genVarShapesDn[proc].SetBinContent(xbin,yvar)

                    #proceed as usual
                    else:
                        dirForSyst=('%s%s%s%s%s_%s_'%(systVar,'up' if not projectRelToNom else '', lbCat,lfs,cat,dist))
                        _,genVarShapesUp = getMergedDists((systfIn if useAltShape else fIn),
                                                        dirForSyst,
                                                        (rawSignalList if opt.addSigs else None),
                                                        ('top' if opt.addSigs else ''),
                                                        widList,nomWid,signalList)

                        dirForSyst=('%s%s%s%s%s_%s_'%(systVar,'dn' if not projectRelToNom else '', lbCat,lfs,cat,dist))
                        _,genVarShapesDn = getMergedDists((systfIn if useAltShape else fIn),
                                                         dirForSyst,
                                                        (rawSignalList if opt.addSigs else None),
                                                        ('top' if opt.addSigs else ''),
                                                        widList,nomWid,signalList)

                #prepare shapes and check if variation is significant
                downShapes, upShapes = {}, {}

                for iproc in exp:
                    nomH=exp[iproc]

                    if not iproc in genVarShapesDn :
                        continue
                    #get directly from another file
                    downH  = genVarShapesDn[iproc]
                    if genVarShapesUp is not None and iproc in genVarShapesUp :
                        upH = genVarShapesUp[iproc]
                    else:
                        #if only one variation is available, mirror it
                        upH = downH.Clone( '%s%sUp'%(iproc,systVar) )
                        for xbin in xrange(1,upH.GetNbinsX()+1):
                            diff=upH.GetBinContent(xbin)-nomH.GetBinContent(xbin)
                            upH.SetBinContent(xbin,nomH.GetBinContent(xbin)-diff)

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

                    if modNomWid in iproc:
                        for tSig in signalList :
                            downShapes['%sNOM'%tSig]=downH
                            upShapes['%sNOM'%tSig]=upH

                #check if something has been accepted
                if len(upShapes)==0 : continue

                #export to shapes file
                if systVar in ["sel","les","trig"] :
                    systVar=systVar+lfs
                extraRateSyst={}
                if systVar in ['ISR','FSR','hdamp','UE','btag','btag%s'%lfs,'ltag','les','les%s'%lfs,'jer','pu','sel','sel%s'%lfs,'toppt','jes','trig%s'%lfs] :
                    for proc in downShapes:
                        nomproc = proc if "NOM" not in proc else proc.replace('NOM','1p0w')
                        nevtsnm=exp[nomproc].Integral()

                        nevtsup,nevtsdn=upShapes[proc].Integral(),downShapes[proc].Integral()
                        if nevtsdn==0 : continue

                        # compute rate
                        rateVarInProc=1.0+0.5*ROOT.TMath.Abs(1.0-nevtsup/nevtsdn)
                        extraRateSyst[proc]=rateVarInProc

                        # normalize variations to nominal integral
                        upShapes[  proc].Scale(nevtsnm/nevtsup)
                        downShapes[proc].Scale(nevtsnm/nevtsdn)

                saveToShapesFile(outFile,downShapes,lbCat+lfs+cat+"_"+dist+"_"+systVar+'Down')
                saveToShapesFile(outFile,upShapes,lbCat+lfs+cat+"_"+dist+"_"+systVar+'Up')

                #write to datacard
                datacard.write('%26s shape'%systVar)
                for sig in signalList:
                    if ("%s%s"%(sig,modNomWid)) in exp and ("%s%s"%(sig,modNomWid)) in upShapes:
                        datacard.write('%15s'%'1')
                    else:
                        datacard.write('%15s'%'-')
                for sig in signalList:
                    if ("%s%s"%(sig,wid)) in exp and ("%s%s"%(sig,wid)) in upShapes:
                        datacard.write('%15s'%'1')
                    else:
                        datacard.write('%15s'%'-')
                for proc in exp:
                    isSig=modWidIn(proc)
                    if isSig : continue
                    if proc in upShapes:
                        datacard.write('%15s'%'1')
                    else:
                        datacard.write('%15s'%'-')
                datacard.write('\n')

                #add extra rate systs
                if len(extraRateSyst)>0:
                    #if "trigMM" in systVar or "selMM" in systVar : continue
                    datacard.write('%26s lnN'%(systVar+'rate'))
                    for sig in signalList:
                        if ("%s%s"%(sig,modNomWid)) in exp and ("%s%s"%(sig,modNomWid)) in extraRateSyst:
                            datacard.write('%15s'%('%3.3f'%extraRateSyst[sig+modNomWid]))
                        else:
                            datacard.write('%15s'%'-')
                    for sig in signalList:
                        if ("%s%s"%(sig,wid)) in exp and ("%s%s"%(sig,wid)) in extraRateSyst:
                            datacard.write('%15s'%('%3.3f'%extraRateSyst[sig+wid]))
                        else:
                            datacard.write('%15s'%'-')
                    for proc in exp:
                        isSig=modWidIn(proc)
                        if isSig : continue
                        if proc in extraRateSyst:
                            datacard.write('%15s'%('%3.3f'%extraRateSyst[proc]))
                        else:
                            datacard.write('%15s'%'-')
                    datacard.write('\n')


            #all done
            datacard.close()


"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
