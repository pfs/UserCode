import pprint
import os
import sys
import optparse
import ROOT
import commands
import getpass
import pickle
import numpy
import gc


ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
ROOT.gROOT.SetBatch()

def replaceBadCharacters(inputStr):
    newStr = inputStr
    for token in ['+','-','*',' ','#','{','(',')','}','@']:
        newStr=newStr.replace(token,'')
    return newStr

"""
get distributions from file
"""
def getDistsFrom(directory,keyFilter=[]) :
    exp={}
    dirName=directory.GetName()
    for key in directory.GetListOfKeys():
        obj=directory.Get(key.GetName())
        if not obj.InheritsFrom('TH1') : continue
        newName=replaceBadCharacters(obj.GetName().split(dirName+'_')[-1])
        for flt in keyFilter :
            newName=newName.replace(flt,'')

        exp[newName]=obj.Clone(newName)
        #exp[newName].SetDirectory(0)

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
    return exp

"""
"""
def doLinearRegression(widValArray,binArray,nomVal) :
    x = numpy.vstack([numpy.array([wid - nomVal for wid in widValArray]),
            numpy.ones(len(widValArray))]).T
    y = numpy.array(binArray)
    m,b = numpy.linalg.lstsq(x,y)[0]
    return m,b

"""
"""
def doTF2Fit(histo, minX, nomX, maxX, minY, nomY, maxY) :
    formula = "[0]*(x-%3.3f)+[1]*(y-%3.3f)+[2]*(x-%3.3f)*(y-%3.3f)+[3]"
    formula = formula%(nomX,nomY,nomX,nomY)
    fit=ROOT.TF2("",formula,minX,maxX,minY,maxY)
    histo.Fit(fit)
    return [fit.GetParameter(0),fit.GetParameter(1),fit.GetParameter(2),fit.GetParameter(3)]


"""
steer the script
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--input',     dest='input',     help='input shapes file',               default=None,            type='string')
    parser.add_option('-d', '--dists',     dest='distList',  help='distribution',                    default='incmlb',        type='string')
    parser.add_option('-s', '--signal',    dest='signal',    help='signal (csv)',                    default='tbart,tW',      type='string')
    parser.add_option('-c', '--cat',       dest='cat',       help='categories (csv)',                default='1b,2b',         type='string')
    parser.add_option('-o', '--output',    dest='output',    help='output directory',                default='datacards',     type='string')
    parser.add_option('-n', '--outname',   dest='outname',   help='output file name',                default='shapes',        type='string')
    parser.add_option(      '--lfs',       dest='lfsInput',  help='lepton final states to consider', default='EE,EM,MM',      type='string')
    parser.add_option(      '--lbCat',     dest='lbCat',     help='pt categories to consider',       default='highpt,lowpt',  type='string')
    parser.add_option(      '--nomW',      dest='nomWid',    help='nominal width [GeV]',             default='1.324',         type='string')
    parser.add_option(      '--nomM',      dest='nomMas',    help='nominal mass  [GeV]',             default='172.5',         type='string')
    parser.add_option(      '--widUnits',  dest='widUnits',  help='use if the input list has units', default=False,     action='store_true')
    parser.add_option('-m', '--makeSplit', dest='makeSplit', help='make split shapes file?',         default=False,     action='store_true')
    parser.add_option(      '--nocards',   dest='nocards',   help='do not produce datacards',        default=False,     action='store_true')
    parser.add_option('-k', '--bkgs',      dest='bkgList',   help='backgrounds (csv)',               default='W,Multiboson,tbartV,DY,tch',
            type='string')
    parser.add_option(      '--mass',      dest='mass',      help='masses to consider',              default='169.5,172.5,175.5',
            type='string')
    parser.add_option('-w', '--wids',      dest='widList',   help='signal widths',
            default='0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.5,4.0',
            type='string')
    (opt, args) = parser.parse_args()


    # parse the channels, pt(l,b), b categories to consider
    sigList  =opt.signal.split(',')
    bkgList  =opt.bkgList.split(',')
    distList =opt.distList.split(',')
    lfsList  =opt.lfsInput.split(',')
    lbCatList=opt.lbCat.split(',')
    catList  =opt.cat.split(',')
    massList =opt.mass.split(',')

    # what are our signal processes?
    rawSignalList=opt.signal.split(',')

    # what are our widths?
    widList = opt.widList.split(',')
    for i in xrange(0,len(widList)) :
        widList[i] = ('%sw'%(widList[i]))
    modNomWid="1p0w"

    modWidList=widList[:]
    for i in xrange(0,len(widList)) :
        modWidList[i]=widList[i].replace('.','p')

    # helpers for coefficient computation
    widValArray=sorted([float(wid.replace('p','.').replace('w','')) for wid in widList])
    masValArray=sorted([float(mas) for mas in massList])
    nomWid=float(opt.nomWid)
    nomMas=float(opt.nomMas)

    # HARDCODE systematics
    tWRateSystList=['tW']
    ttRateSystList=['tbart']
    tWIsSig=('tW'    in rawSignalList)
    ttIsSig=('tbart' in rawSignalList)
    for wid in modWidList :
        if tWIsSig : tWRateSystList += ["tW%s"%wid]
        if ttIsSig : ttRateSystList += ["tbart%s"%wid]

    rateSysts=[
          ('lumi_13TeV',       1.062,    'lnN',    []                   ,['Multijetsdata']),
          ('DYnorm_th',        1.30,     'lnN',    ['DY']  ,[]),
          ('Wnorm_th',         1.30,     'lnN',    ['W']   ,[]),
          ('tWnorm_th',        1.054,    'lnN',    tWRateSystList,[]),
          ('tnorm_th',         1.044,    'lnN',    ['tch']              ,[]),
          ('VVnorm_th',        1.20,     'lnN',    ['Multiboson']       ,[]),
          ('tbartVnorm_th',    1.30,     'lnN',    ['tbartV']           ,[]),
          ('sel',              1.02,     'lnN',    tWRateSystList+ttRateSystList+['DY','W','tch','Multiboson','tbartV'], []),
    ]

    ttParton_tt=['tbartscaledown','tbartscaleup']
    ttParton_tW=['tWscaledown','tWscaleup']
    tWinterf=['tWDS']
    amcNLO  =['tbartamcnloFxFx']
    Herwig  =['tbartHerwig']

    systSignalList=ttParton_tt+amcNLO+Herwig+tWinterf+ttParton_tW

    MtopMap=['m=175.5','m=169.5']
    ttPartonMap={}
    twPartonMap={}
    tWinterfMap={}
    amcNLOMap={}
    HerwigMap={}

    for wid in modWidList :
        ttPartonMap['tbart%s'%wid]=['tbartscaledown%s'%(wid),'tbartscaleup%s'%(wid)]
        twPartonMap['tW%s'   %wid]=['tWscaledown%s'   %(wid),'tWscaleup%s'   %(wid)]
        amcNLOMap['tbart%s'%wid]=['tbartamcnloFxFx%s'%(wid)]
        HerwigMap['tbart%s'%wid]=['tbartHerwig%s'%(wid)]
        tWinterfMap['tW%s'   %wid]=['tWDS%s'   %(wid)]

    sampleSysts=[
          #ttbar modelling
          ('ttPartonShower', ttPartonMap, False, True, False),
          ('tWQCDScale'    , twPartonMap, False, True, False),
          ('Herwig'        , HerwigMap  , False, True, True),
          ('amcnloFxFx'    , amcNLOMap  , False, True, True),
          #tWinterference
          ('tWttinterf'    , tWinterfMap, False, True, True),
    ]

    genSysts=[
        ('jes',   True, False,False,False),
        ('les',   False,False,False,False),
        ('ltag',  False,False,False,False),
        ('trig',  False,False,False,True),
        #('sel',   False,False,False,True),
        ('toppt', True, False,False,False),
        ('jer',   False,False,False,False),
        ('btag',  False,False,False,False),
        ('pu',    True, False,False,False),
        ('MEqcdscale',  True,False,False,False),
        ('PDF',   False,False,False,False)
    ]

    # prepare output directory
    os.system('mkdir -p %s'%opt.output)

    # get data and nominal expectations
    fIn=ROOT.TFile.Open(opt.input,"READ")

    # keep track of progress
    cardIndex=0
    numCards =len(lfsList)*len(lbCatList)*len(catList)*len(distList)*30

    # prepare array for picklefile
    binCoeffs={}

    #loop over lepton final states
    for (lbCat,lfs,cat,dist) in [(a,b,c,d) for a in lbCatList
            for b in lfsList
            for c in catList
            for d in distList] :
        if opt.nocards : break;
        category="%s%s%s"%(lbCat,lfs,cat)
        exp=getDistsFrom(fIn.Get("%s_%s"%(category,dist)))
        obs=exp["data_obs"].Clone()
        for bin in range(1,31) :
            #loop over categories, widths
            cardIndex+=1
            print 'Initiating %s datacard for %s \t\t [%i/%i]'%(dist,category,cardIndex,numCards)

            #get distributions for this category
            mUp,mDn = None,None
            if bin==1 :
                mUp=getDistsFrom(fIn.Get("%s_%s_MtopUp"%(category,dist)),keyFilter=MtopMap)
                mDn=getDistsFrom(fIn.Get("%s_%s_MtopDown"%(category,dist)),keyFilter=MtopMap)

            #prep coefficient storage
            if len(binCoeffs)<1:
                for ibin in range(1,obs.GetNbinsX()+1):
                    binCoeffs[ibin]={}

            #compute coefficients
            for sig,ibin in [(a,b)
                    for a in sigList
                    for b in range(1,obs.GetNbinsX()+1)] :
                if bin!= 1 : break
                histo=ROOT.TH2D("","",
                        len(widValArray),widValArray[0],widValArray[-1],
                        len(masValArray),masValArray[0],masValArray[-1])
                for wid in widValArray :
                    histo.Fill(wid,masValArray[0] ,mDn[sig+("%2.1fw"%(wid)).replace('.','p')].GetBinContent(ibin))
                    histo.Fill(wid,nomMas         ,exp[sig+("%2.1fw"%(wid)).replace('.','p')].GetBinContent(ibin))
                    histo.Fill(wid,masValArray[-1],mUp[sig+("%2.1fw"%(wid)).replace('.','p')].GetBinContent(ibin))

                binCoeffs[ibin][category+"_"+sig]=doTF2Fit(histo,
                        widValArray[0],nomWid,widValArray[-1],
                        masValArray[0],nomMas,masValArray[-1])

                #binArray=[]
                #mUpArray=[]
                #mDnArray=[]
                #mBnArray=[]

                ## collect bin information
                #for wid in widValArray :
                #    binArray += [exp[sig+("%2.1fw"%(wid)).replace('.','p')].GetBinContent(ibin)]
                #    mUpArray += [mUp[sig+("%2.1fw"%(wid)).replace('.','p')].GetBinContent(ibin)]
                #    mDnArray += [mDn[sig+("%2.1fw"%(wid)).replace('.','p')].GetBinContent(ibin)]

                #for mas in masValArray :
                #    mBnArray += [exp[sig+modNomWid].GetBinContent(ibin)]

                # get linear terms for mass, width
                #a,d1=doLinearRegression(widValArray,binArray,nomWid)
                #b,d2=doLinearRegression(masValArray,mBnArray,nomMas)

                # get corner term (d/dM)(d/dW)
                #c,d3= doLinearRegression(widValArray,mUpArray,nomWid)
                #c2,d4=doLinearRegression(widValArray,mDnArray,nomWid)
                #c-=c2
                #c/=(masValArray[-1]-masValArray[0])

                #d=d1
                #binCoeffs[ibin][category+"_"+sig]=[a,b,c,d]

            #start the datacard
            datacard=open('%s/datacard__%s_bin%i_%s_%s.dat'%(opt.output,modNomWid,bin,category,dist),'w')
            datacard.write('#\n')
            datacard.write('# Generated by %s for analysis category %s_%s\n' % (getpass.getuser(),
                category, wid) )

            datacard.write('#\n')
            datacard.write('imax *\n')
            datacard.write('jmax *\n')
            datacard.write('kmax *\n')
            datacard.write('-'*50+'\n')
            datacard.write('shapes *        * shapes.root $CHANNEL/$PROCESS $CHANNEL_$SYSTEMATIC/$PROCESS\n')
            datacard.write('-'*50+'\n')
            datacard.write('bin bin%i_%s_%s\n'%(bin,category,dist))
            datacard.write('observation %3.1f\n' % obs.GetBinContent(bin))
            datacard.write('-'*50+'\n')

            if "data_obs" in exp : del exp["data_obs"]

            #expectations
            datacard.write('\t\t\t %32s'%'bin')
            for i in xrange(0,len(sigList)+len(bkgList)):
                datacard.write('%30s'%('bin%i_%s_%s'%(bin,category,dist)))
            datacard.write('\n')
            datacard.write('\t\t\t %32s'%'process')
            for sig in sigList: datacard.write('%15s'%('%s%s'%(sig,modNomWid)))
            for proc in exp:
                isSig=False
                for sig in modWidList :
                    if sig in proc and 'V' not in proc: isSig=True
                if isSig : continue
                datacard.write('%15s'%proc)
            datacard.write('\n')
            datacard.write('\t\t\t %32s'%'process')
            for i in xrange(0,len(sigList)) : datacard.write('%15s'%str(i+1-len(sigList)))
            i=0
            for proc in exp:
                isSig=False
                for sig in modWidList :
                    if sig in proc and 'V' not in proc: isSig=True
                if isSig : continue
                i=i+1
                datacard.write('%15s'%str(i))
            datacard.write('\n')
            datacard.write('\t\t\t %32s'%'rate')
            for sig in sigList :
                datacard.write('%15s'%('%3.2f'%(exp['%s%s'%(sig,modNomWid)].GetBinContent(bin))))
            for proc in exp:
                isSig=False
                for sig in modWidList :
                    if sig in proc and 'V' not in proc: isSig=True
                if isSig : continue
                datacard.write('%15s'%('%3.2f'%exp[proc].GetBinContent(bin)))
            datacard.write('\n')
            datacard.write('-'*50+'\n')


            for syst,val,pdf,whiteList,blackList in rateSysts:
                if syst in ['sel','DYnorm_th'] : syst="%s%s"%(syst,lfs)

                datacard.write('%32s %8s'%(syst,pdf))
                entryTxt=''
                try:
                    entryTxt='%15s'%('%3.3f/%3.3f'%(ROOT.TMath.Max(val[0],0.01),val[1]))
                except:
                    entryTxt='%15s'%('%3.3f'%val)

                for sig in sigList :
                    newSig='%s%s'%(sig,modNomWid)
                    if (len(whiteList)==0 and not newSig in blackList) or newSig in whiteList:
                        datacard.write(entryTxt)
                    else:
                        datacard.write('%15s'%'-')
                for proc in exp:
                    isSig=False
                    for sig in modWidList :
                        if sig in proc and 'V' not in proc: isSig=True
                    if isSig : continue
                    if (len(whiteList)==0 and not proc in blackList) or proc in whiteList:
                        datacard.write(entryTxt)
                    else:
                        datacard.write('%15s'%'-')
                datacard.write('\n')

            # sample systematics
            for systVar, procsToApply, normalize, useAltShape, projectRelToNom in sampleSysts:
                upShapes=getDistsFrom(fIn.Get("%s_%s_%sUp"%(category,dist,systVar)))

                #write to datacard
                datacard.write('%26s shape'%systVar)
                for sig in sigList :
                    isInShapes=False
                    for upS in upShapes :
                        isInShapes=isInShapes or (sig in upS)
                    if isInShapes:
                        datacard.write('%15s'%'1')
                    else:
                        datacard.write('%15s'%'-')
                for proc in exp:
                    isSig=False
                    for sig in modWidList :
                        if sig in proc and 'V' not in proc: isSig=True
                    if isSig : continue
                    isInShapes=False
                    for upS in upShapes :
                        sigIn=False
                        for sig in sigList :
                            sigIn=sigIn or (sig in upS)
                        if sigIn : continue
                        isInShapes=isInShapes or (proc in upS)
                    if isInShapes:
                        datacard.write('%15s'%'1')
                    else:
                        datacard.write('%15s'%'-')
                datacard.write('\n')
                del upShapes

            # sample systematics
            for systVar, normalize, useAltShape, projectRelToNom, chanSep in genSysts:
                chan = (lfs if chanSep else '')
                upShapes=getDistsFrom(fIn.Get("%s_%s_%s%sUp"%(category,dist,systVar,chan)))

                #write to datacard
                datacard.write('%26s shape'%(systVar+chan))
                for sig in sigList :
                    if ("%s%s"%(sig,modNomWid)) in exp and ("%s%s"%(sig,modNomWid)) in upShapes:
                        datacard.write('%15s'%'1')
                    else:
                        datacard.write('%15s'%'-')
                for proc in exp:
                    isSig=False
                    for sig in modWidList :
                        if sig in proc and 'V' not in proc: isSig=True
                    if isSig : continue
                    if proc in upShapes:
                        datacard.write('%15s'%'1')
                    else:
                        datacard.write('%15s'%'-')
                datacard.write('\n')
                del upShapes


            #all done
            datacard.close()
            del datacard
            del mDn
            del mUp
            gc.collect()
        del exp
        del obs
        gc.collect()

    # save to picklefile
    if not opt.nocards:
        pp=pprint.PrettyPrinter(indent=2)
        pp.pprint(binCoeffs)
        with open('%s/coeff_cachefile.pck'%opt.output,'w') as cachefile:
            pickle.dump(binCoeffs, cachefile, pickle.HIGHEST_PROTOCOL)

    if not opt.makeSplit : return
    fOut=ROOT.TFile("%s/shapes_split.root"%opt.output,"RECREATE")

    for key in fIn.GetListOfKeys() :
        cdir=fIn.Get(key.GetName())
        dirName=cdir.GetName();
        formatName="bin%i_"+dirName

        print formatName

        fOut.cd()
        for ibin in range(1,31) :
            fOut.mkdir(formatName%ibin)
        for histKey in cdir.GetListOfKeys() :
            hist=cdir.Get(histKey.GetName());
            #hist.SetDirectory(0)
            for ibin in range(1,31) :
                fOut.cd(formatName%ibin)
                ohist=ROOT.TH1F(histKey.GetName(),"",1,
                        hist.GetBinLowEdge(ibin),
                        hist.GetBinLowEdge(ibin)+hist.GetBinWidth(ibin))
                ohist.Fill(hist.GetBinCenter(ibin),hist.GetBinContent(ibin))
                ohist.Write()
                del ohist
            del hist
            #gc.collect()

    fOut.Close();


"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
