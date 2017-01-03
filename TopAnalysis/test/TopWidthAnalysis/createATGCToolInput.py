import pprint
import os
import sys
import optparse
import ROOT
import commands
import getpass
import gc

from math import sqrt, pow

ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
ROOT.gROOT.SetBatch()

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
    parser.add_option('-n', '--chNameAdd', dest='chNAdd',    help='name to append to channels',      default='',              type='string')
    parser.add_option(      '--lfs',       dest='lfsInput',  help='lepton final states to consider', default='EE,EM,MM',      type='string')
    parser.add_option(      '--lbCat',     dest='lbCat',     help='pt categories to consider',       default='highpt,lowpt',  type='string')
    parser.add_option(      '--injWid',    dest='injWid',    help='injected width [GeV]',            default='',              type='string')
    parser.add_option(      '--injMas',    dest='injMas',    help='injected mass  [GeV]',            default='',              type='string')
    parser.add_option(      '--nomW',      dest='nomWid',    help='nominal width [GeV]',             default='1.324',         type='string')
    parser.add_option(      '--nomM',      dest='nomMas',    help='nominal mass  [GeV]',             default='172.5',         type='string')
    parser.add_option(      '--widUnits',  dest='widUnits',  help='use if the input list has units', default=True,     action='store_true')
    parser.add_option(      '--masUnits',  dest='masUnits',  help='use if the input list has units', default=True,      action='store_true')
    parser.add_option('-m', '--makeSplit', dest='makeSplit', help='make split shapes file?',         default=False,     action='store_true')
    parser.add_option(      '--nocards',   dest='nocards',   help='do not produce datacards',        default=False,     action='store_true')
    parser.add_option('-k', '--bkgs',      dest='bkgList',   help='backgrounds (csv)',               default='W,Multiboson,tbartV,DY,tch',
            type='string')
    parser.add_option(      '--mass',      dest='mass',      help='masses to consider',              default='169.5,172.5,175.5',
            type='string')
    parser.add_option('-w', '--wids',      dest='widList',   help='signal widths',
            default='0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0',#,3.5,4.0',
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

    rateSysts=[
          ('lumi_13TeV',       1.062,    'lnN',    ['tbart','tW','DY','W','tch','Multiboson','tbartV'], False),
          ('DYnorm_th',        1.30,     'lnN',    ['DY']               ,True),
          ('Wnorm_th',         1.30,     'lnN',    ['W']                ,False),
          ('tWnorm_th',        1.054,    'lnN',    ['tW']               ,False),
          ('tnorm_th',         1.044,    'lnN',    ['tch']              ,False),
          ('VVnorm_th',        1.20,     'lnN',    ['Multiboson']       ,False),
          ('tbartVnorm_th',    1.30,     'lnN',    ['tbartV']           ,False),
          ('sel',              1.02,     'lnN',    ['tbart','tW','DY','W','tch','Multiboson','tbartV'], True),
    ]

    # prepare output directory
    os.system('mkdir -p %s'%opt.output)
    chNAdd=opt.chNAdd

    # get data and nominal expectations
    fIn=ROOT.TFile.Open(opt.input,"READ")

    # keep track of progress
    cardIndex=0
    numCards =len(lfsList)*len(lbCatList)*len(catList)*len(distList)*30

    def replaceAllSystNames(tkey) :
        key=tkey
        for genName,_,_,_,_ in genSysts :
            key.replace("_"+genName+"Up","")
            key.replace("_"+genName+"Down","")
        for sampleName,_,_,_,_ in sampleSysts :
            key.replace("_"+sampleName+"Up","")
            key.replace("_"+sampleName+"Down","")
        return key

    def isSig(proc) :
        for sig,wid in [(a,b) for a in sigList for b in modWidList] :
            if sig==proc or (sig in proc and proc !='tbartV'): return True
        return False

    def makeMCTemplate(channel,injWid,injMas) :
        modInjWid=(injWid+"w").replace('.','p')
        tChannel=channel
        hist=None

        if   injMas=="175.5" : tChannel+="_MtopUp"
        elif injMas=="169.5" : tChannel+="_MtopDown"

        cdir=fIn.Get(tChannel)

        for histKey in cdir.GetListOfKeys() :
            if isSig(histKey.GetName()) and modInjWid not in histKey.GetName() :
                continue
            thist=cdir.Get(histKey.GetName()).Clone()
            if hist is None : hist=thist.Clone("data_obs")
            else : hist.Add(thist)

        if hist is None : print "ERROR: data_obs is None for channel,",tChannel
        return hist

    #############################
    # Get input settings card   #
    #############################
    datacard=open('%s/config_all%s.dat'%(opt.output,chNAdd),'w')
    datacard.write('[Global]\n')
    datacard.write('model=par1par2_TH2\n')
    datacard.write('par1Name = mt\n')
    datacard.write('par1Low = %6.6f\n' %(min(masValArray)))
    datacard.write('par1High = %6.6f\n'%(max(masValArray)))
    datacard.write('par2Name = gam\n')
    datacard.write('par2Low = %6.6f\n' %(min(widValArray)))
    datacard.write('par2High = %6.6f\n'%(max(widValArray)))

    channelList= [a+b+c+"_"+d for a in lbCatList
            for b in lfsList
            for c in catList
            for d in distList]

    #get the log-normal systs and add them in
    nRateSysts=len(rateSysts)
    # TODO: separate lnN uncertainties by ch
    #nRateSysts+=len([x for _,_,_,_,x in rateSysts if x])*(len(lfsList)-1)

    iRateSyst=1
    datacard.write('NlnN=%i\n'%(nRateSysts))
    for syst,val,pdf,whitelist,sepByCh in rateSysts :
        if opt.nocards : break;
        sigWhiteList=[x for x in whitelist if isSig(x)]
        bkgWhiteList=[x for x in whitelist if not isSig(x)]

        sigValStr=''
        sigNameStr=''
        for cat in channelList :
            if len(sigWhiteList) == 0 : break
            if sigValStr != '' : sigValStr+=','
            if sigNameStr != '' : sigNameStr+=','
            sigNameStr+=(cat+chNAdd+"_signal")
            newRate=0

            cDir=fIn.Get(cat)
            totalNorm=0
            for sig in sigList :
                integral=cDir.Get(sig+"1p0w").Integral()
                totalNorm+=integral
                if sig in sigWhiteList :
                    newRate+=pow((val-1)*integral,2)

            newRate=sqrt(newRate)
            newRate/=totalNorm
            newRate+=1

            sigValStr+='%1.2f'%newRate

        datacard.write('lnN%i_name=%s\n'%(iRateSyst,syst))
        datacard.write('lnN%i_value=%s'%(iRateSyst,
            ','.join([("%1.2f"%val)]*len(channelList)*len(bkgWhiteList))))
        if len(bkgWhiteList) > 0 and len(sigWhiteList) > 0 : datacard.write(',')
        datacard.write(sigValStr)

        datacard.write('\nlnN%i_for='%iRateSyst)
        datacard.write(','.join([x+chNAdd+"_"+y for x in channelList for y in bkgWhiteList]))
        if len(bkgWhiteList) > 0 and len(sigWhiteList) > 0 : datacard.write(',')
        datacard.write(sigNameStr)
        datacard.write('\n')

        iRateSyst+=1

    datacard.write('\n')

    #loop over categories
    for (lbCat,lfs,cat,dist) in [(a,b,c,d) for a in lbCatList
            for b in lfsList
            for c in catList
            for d in distList] :
        if opt.nocards : break;
        category="%s%s%s"%(lbCat,lfs,cat)
        mainDirName="%s_%s"%(category,dist)

        fOut = ROOT.TFile("%s/%s%s.root"%(opt.output,mainDirName,chNAdd),"RECREATE")

        datacard.write('\n[%s%s]\n'%(mainDirName,chNAdd))

        # Get SM/backgrounds
        # (add signals)
        bkgMap={}
        for histKey in fIn.Get(mainDirName).GetListOfKeys() :
            cName=histKey.GetName()

            # special handling for data
            if cName=="data_obs" and opt.injWid!="" and opt.injMas!="":
                bkgMap[cName] = {}
                bkgMap[cName][cName] = makeMCTemplate(mainDirName,opt.injWid,opt.injMas)
                continue

            # special handling for signals
            tisSig = isSig(cName)
            if tisSig and '1p0w' not in cName: continue
            if tisSig and "diboson" not in bkgMap : bkgMap["diboson"]={}


            # FIXME: Is it possible to get rid of this annoying hardcoding?
            # add signals or just save bkg histogram
            tcName=cName
            if tisSig and "diboson" not in bkgMap["diboson"] :
                bkgMap["diboson"]["diboson"] = fIn.Get(mainDirName+"/"+cName).Clone("diboson")
                tcName=cName.replace('1p0w','')
            elif tisSig :
                bkgMap["diboson"]["diboson"].Add(fIn.Get(mainDirName+"/"+cName).Clone())
                tcName=cName.replace('1p0w','')

            if tcName not in bkgMap : bkgMap[tcName]={}
            bkgMap[tcName][tcName] = fIn.Get(mainDirName+"/"+cName).Clone()

        # Get SM/bkg systematics
        for tDirKey in fIn.GetListOfKeys() :
            if "Mtop" in tDirKey.GetName() : continue
            if "Down" in tDirKey.GetName() : continue
            if category not in tDirKey.GetName() : continue
            if tDirKey.GetName() == mainDirName  : continue
            systName=tDirKey.GetName().replace(category+"_","")
            for histKey in fIn.Get(tDirKey.GetName()).GetListOfKeys() :
                cName=histKey.GetName()
                tisSig=isSig(cName)
                if tisSig and '1p0w' not in cName: continue
                elif tisSig : cName=cName.replace('1p0w','')

                upName=cName+"_"+systName
                dnName=upName.replace("Up","Down")

                if cName not in bkgMap : bkgMap[cName]={}

                bkgMap[cName][upName]=fIn.Get(tDirKey.GetName()+"/"+histKey.GetName()).Clone(upName)
                bkgMap[cName][dnName]=fIn.Get(tDirKey.GetName().replace("Up","Down")+"/"+histKey.GetName()).Clone(dnName)

        pp= pprint.PrettyPrinter(indent=2)
        #pp.pprint(bkgMap)

        # Fix signal uncertainties (add signals in quadrature)
        for sig,tsig in [(a,b) for a in sigList for b in sigList] :
            if sig not in bkgMap : continue
            if tsig not in bkgMap : continue
            if sig==tsig : continue


            #print sig," ",tsig
            # loop over signals
            for histName in bkgMap[sig] :
                if histName == sig : continue
                if 'Down' in histName or 'Up' not in histName : continue
                #print ' -- ',histName

                histNameDn=histName.replace('Up','Down')
                ctHistNameUp=histName.replace(sig,tsig)
                ctSigsNameUp=histName.replace(sig,"diboson")
                ctHistNameDn=ctHistNameUp.replace('Up','Down')
                ctSigsNameDn=ctSigsNameUp.replace('Up','Down')

                # make the corresponding histogram in the combined signal set
                if histName in bkgMap["diboson"] : continue
                bkgMap["diboson"][ctSigsNameUp]=bkgMap[sig][histName  ].Clone(ctSigsNameUp)
                bkgMap["diboson"][ctSigsNameDn]=bkgMap[sig][histNameDn].Clone(ctSigsNameDn)
                bkgMap["diboson"][ctSigsNameUp].Add(bkgMap[sig][sig],-1)
                bkgMap["diboson"][ctSigsNameDn].Add(bkgMap[sig][sig],-1)

                # get the upward and downward variations - add appropriately
                # (downward becomes deeper, upward becomes steeper)
                for thistName in bkgMap[tsig] :
                    if thistName != ctHistNameUp : continue
                    for ibin in range(1,31) :
                        ntBinContent=bkgMap[tsig][tsig].GetBinContent(ibin)
                        tBinContentUp=bkgMap[tsig][ctHistNameUp].GetBinContent(ibin)
                        tBinContentDn=bkgMap[tsig][ctHistNameDn].GetBinContent(ibin)

                        tBinContentUp-=ntBinContent
                        tBinContentDn-=ntBinContent

                        cBinContentUp=bkgMap["diboson"][ctSigsNameUp].GetBinContent(ibin)
                        cBinContentDn=bkgMap["diboson"][ctSigsNameDn].GetBinContent(ibin)

                        if cBinContentUp == cBinContentDn : cBinContentDn *= -1
                        if tBinContentUp == tBinContentDn : tBinContentDn *= -1

                        # swap bin contents to match variations
                        newContentsUp=0
                        newContentsDn=0
                        upVar=max(tBinContentUp,tBinContentDn)
                        dnVar=min(tBinContentUp,tBinContentDn)

                        #print " --- ",tBinContentUp," ",tBinContentDn," "

                        if cBinContentUp >= 0 and cBinContentDn < 0:
                            newContentsUp=sqrt(pow(cBinContentUp,2)+pow(upVar,2))
                            newContentsDn=sqrt(pow(cBinContentDn,2)+pow(dnVar,2))
                        elif cBinContentUp <= 0 and cBinContentDn >= 0:
                            newContentsUp=sqrt(pow(cBinContentUp,2)+pow(dnVar,2))
                            newContentsDn=sqrt(pow(cBinContentDn,2)+pow(upVar,2))
                        else :
                            print "ERROR: Problem with shape uncertainty ",ctSigsNameUp
                            exit(0)

                        bkgMap["diboson"][ctSigsNameUp].SetBinContent(ibin,newContentsUp)
                        bkgMap["diboson"][ctSigsNameDn].SetBinContent(ibin,newContentsDn)

        # add back in the nominal information
        for hist in bkgMap["diboson"] :
            if hist=="diboson" : continue
            bkgMap["diboson"][hist].Add(bkgMap["diboson"]["diboson"])

        # Write bkgMap to file
        fOut.cd()
        for process in bkgMap :
            if isSig(process) : continue
            for hist in bkgMap[process] :
                bkgMap[process][hist].Write()

        # Write in the shape systematics for signal/bkg
        nDataInMap= 1 if "data_obs" in bkgMap else 0
        nSigInMap = len({k: v for k,v in bkgMap.iteritems() if isSig(v)})+1

        datacard.write("Nbkg=%i\n"%(len(bkgMap)-nDataInMap-nSigInMap))
        iBkg=1
        for nBkg,histList in bkgMap.iteritems() :
            if "data_obs" in nBkg : continue
            if "diboson" in nBkg : continue
            if isSig(nBkg) : continue
            datacard.write("bkg%i_name=%s\n"%(iBkg,nBkg))
            bkgStr=",".join([x.replace('Up','')
                for x in histList if len(x.split('_'))>1 and 'Down' not in x])
            if bkgStr != '' :
                datacard.write("bkg%i_shape_syst=%s\n"%(iBkg,bkgStr))
            iBkg+=1;

        #all done
        fOut.Close();
        del bkgMap
        gc.collect()

    datacard.close()

    #############################
    # Get signal model TH2Ds    #
    #############################

    if not opt.makeSplit : return
    unitWidList=[float(x)*(1 if opt.widUnits else nomWid) for x in widValArray]
    unitMasList=[float(x)*(1 if opt.masUnits else nomMas) for x in masValArray]

    for (lbCat,lfs,cat,dist) in [(a,b,c,d) for a in lbCatList
            for b in lfsList
            for c in catList
            for d in distList] :
        category="%s%s%s"%(lbCat,lfs,cat)
        mainDirName="%s_%s"%(category,dist)

        fOut = ROOT.TFile("%s/signal_proc_%s%s.root"%(opt.output,mainDirName,chNAdd),"RECREATE")

        hists={}
        # Get full set of nominal information
        for key in fIn.GetListOfKeys() :
            cdir=fIn.Get(key.GetName())
            dirName=cdir.GetName();

            # FIXME: Don't want syst uncertainties to be TH2s
            if dirName != mainDirName : continue

            formatName="bin_content_par1_par2_%i"

            fOut.cd()
            for histKey in cdir.GetListOfKeys() :
                if "NOM" in histKey.GetName() : continue
                hist=cdir.Get(histKey.GetName());
                sig=False
                cWid=False

                for tSig in sigList :
                    if tSig not in histKey.GetName() : continue
                    if 'tbartV' in histKey.GetName() : continue
                    sig=tSig
                    cWid=histKey.GetName().replace(sig,'').replace('p','.').replace('w','')
                    cWid=float(cWid)
                    break

                for ibin in range(1,31) :
                    if formatName%(ibin) not in hists :
                        hists[formatName%(ibin)]=ROOT.TH2D(formatName%(ibin),"",
                                3,  min(unitMasList)-1.5, max(unitMasList)+1.5,
                                15, min(unitWidList)-0.1, max(unitWidList)+0.1)
                    hists[formatName%(ibin)].Fill(nomMas,cWid,hist.GetBinContent(ibin))
                del hist
                gc.collect()

        # Get up and down mass variations
        for key in fIn.GetListOfKeys() :
            cdir=fIn.Get(key.GetName())
            dirName=cdir.GetName()

            cMass=False

            if category not in dirName : continue
            if "Mtop" not in dirName : continue
            if "MtopUp" in dirName : cMass = max(unitMasList)
            if "MtopDown" in dirName : cMass = min(unitMasList)
            if not opt.masUnits : cMass *= nomMas

            #print dirName
            formatName="bin_content_par1_par2_%i"

            fOut.cd()
            for histKey in cdir.GetListOfKeys() :
                if "NOM" in histKey.GetName() : continue
                hist=cdir.Get(histKey.GetName());
                sig=False
                cWid=False

                for tSig in sigList :
                    if tSig not in histKey.GetName() : continue
                    if 'tbartV' in histKey.GetName() : continue
                    sig=tSig
                    cWid=histKey.GetName().replace(sig,'').replace('p','.').replace('w','')
                    cWid=float(cWid)
                    if not opt.widUnits : cWid *= nomWid
                    break

                for ibin in range(1,31) :
                    if formatName%(ibin) not in hists or not cMass :
                        print "ERROR: Histogram %s not found and cMass is "%(formatName%(ibin)),cMass
                        exit(0)
                    #print "Filling ",formatName%ibin," at ",cMass," ",cWid," ",hist.GetBinContent(ibin)
                    hists[formatName%(ibin)].Fill(cMass,cWid,hist.GetBinContent(ibin))
                del hist
                gc.collect()

        # Normalize all templates relative to SM simulation
        for key,hist in hists.iteritems():
            nominalKey=replaceAllSystNames(key)
            nMasBin=hists[nominalKey].GetXaxis().FindBin(nomMas)
            nWidBin=hists[nominalKey].GetYaxis().FindBin(nomWid)
            nominalVal=hists[nominalKey].GetBinContent(nMasBin,nWidBin)

            hists[key].Scale(1/nominalVal)

        # Write all to outfile
        fOut.cd()
        for hist in hists:
            hists[hist].Write();
        fOut.Close();


"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
