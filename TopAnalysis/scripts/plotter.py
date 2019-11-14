import optparse
import os,sys
import json
import pickle
from collections import OrderedDict

from TopLJets2015.TopAnalysis.Plot import *


"""
steer the script
"""
	
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option(     '--mcUnc',        dest='mcUnc'  ,     help='common MC related uncertainty (e.g. lumi)',        default=0,              type=float)
    parser.add_option(     '--com',          dest='com'  ,       help='center of mass energy',                            default='13 TeV',       type='string')
    parser.add_option(     '--rawYields',    dest='rawYields',   help='do not scale by lumi, xsec, etc.',                 default=False, action='store_true')
    parser.add_option(     '--blined',       dest='blined',      help='make blind plotting',                              default=False, action='store_true')
    parser.add_option('-j', '--json',        dest='json'  ,      help='json with list of files',        default=None,              type='string')
    parser.add_option( '--systJson',         dest='systJson',    help='json with list of systematics', default=None, type='string')
    parser.add_option(      '--signalJson',  dest='signalJson',  help='signal json list',               default=None,              type='string')
    parser.add_option('-i', '--inDir',       dest='inDir' ,      help='input directory',                default=None,              type='string')
    parser.add_option('-O', '--outDir',      dest='outDir' ,     help='output directory',                default=None,              type='string')
    parser.add_option('-o', '--outName',     dest='outName' ,    help='name of the output file',        default='plotter.root',    type='string')
    parser.add_option(      '--noStack',     dest='noStack',     help='don\'t stack distributions',     default=False,             action='store_true')
    parser.add_option(      '--normToData',  dest='normToData',  help='normalize to data yields',       default=False,             action='store_true')
    parser.add_option(      '--binWid',      dest='binWid',      help='divide by bin width',            default=False,             action='store_true')
    parser.add_option(      '--saveLog',     dest='saveLog' ,    help='save log versions of the plots', default=False,             action='store_true')
    parser.add_option(      '--silent',      dest='silent' ,     help='only dump to ROOT file',         default=False,             action='store_true')
    parser.add_option(      '--ratioRange',  dest='ratioRange' , help='ratio range',                    default="0.7,1.3",         type='string')
    parser.add_option(      '--onlyData',    dest='onlyData' ,   help='only plots containing data',     default=False,             action='store_true')
    parser.add_option(      '--saveTeX',     dest='saveTeX' ,    help='save as tex file as well',       default=False,             action='store_true')
    parser.add_option(      '--rebin',       dest='rebin',       help='rebin factor',                   default=1,                 type=int)
    parser.add_option('-l', '--lumi',        dest='lumi' ,       help='lumi to print out, if == 1 draw normalized',              default=12900,              type=float)
    parser.add_option(      '--lumiSpecs',   dest='lumiSpecs',   help='lumi specifications for some channels [tag:lumi,tag2:lumi2,...]', default=None,       type=str)
    parser.add_option(      '--only',        dest='only',        help='plot only these (csv)',          default='',                type='string')
    parser.add_option(      '--strictOnly',  dest='strictOnly',  help='strict matching for only plots', default=False, action='store_true')
    parser.add_option(      '--skip',        dest='skip',        help='skip these samples (csv)',       default='MC13TeV_TTJets_cflip',                type='string')
    parser.add_option(      '--rawList',     dest='rawList',     help='don\'t scale these samples',     default='',                type='string')
    parser.add_option(      '--puNormSF',    dest='puNormSF',    help='Use this histogram to correct pu weight normalization', default=None, type='string')
    parser.add_option(      '--procSF',      dest='procSF',      help='Use this to scale a given process component e.g. "W":.wjetscalefactors.pck,"DY":dyscalefactors.pck', default=None, type='string')
    (opt, args) = parser.parse_args()

    opt.ratioRange=[float(x) for x in opt.ratioRange.split(',')]

    #read lists of samples
    samplesList=[]
    jsonList = opt.json.split(',')
    for jsonPath in jsonList:
        jsonFile = open(jsonPath,'r')
        samplesList += json.load(jsonFile, encoding='utf-8', object_pairs_hook=OrderedDict).items()
        jsonFile.close()

    #read lists of syst samples
    systSamplesList=[]
    if opt.systJson:
        systJsonList = opt.systJson.split(',')
        for jsonPath in systJsonList:
            jsonFile = open(jsonPath,'r')
            systSamplesList += json.load(jsonFile,encoding='utf-8').items()
            jsonFile.close()

    #read list of signal samples
    signalSamplesList=None
    try:
        jsonFile = open(opt.signalJson,'r')
        signalSamplesList=json.load(jsonFile, encoding='utf-8', object_pairs_hook=OrderedDict).items()
        print signalSamplesList
        jsonFile.close()
    except Exception as e:
        print e
        pass
    
    skipList=opt.skip.split(',')



    #lumi specifications per tag
    lumiSpecs={}
    if opt.lumiSpecs:
        for spec in opt.lumiSpecs.split(','):
            tag,lumi=spec.split(':')
            lumiSpecs[tag]=float(lumi)
        print lumiSpecs

    #proc SF
    procSF={}
    if opt.procSF:
        procList=opt.procSF.split(',')
        for newProc in procList:
            proc,cacheUrl=newProc.split(':')
            if os.path.isfile(cacheUrl) : 
                cache=open(cacheUrl,'r')
                procSF[proc]=pickle.load(cache)
                cache.close()
                print 'Scale factors added for',proc
            else:
                try:
                    procSF[proc]={'':(float(cacheUrl),0)}
                    print 'Scale factors added for',proc
                except:
                    pass


    onlyList=opt.only.split(',')
    rawList=opt.rawList.split(',')

    #read plots 
    plots=OrderedDict()

    report=''
    for slist,isSignal,isSyst in [ (samplesList,False,False),(signalSamplesList,True,False),(systSamplesList,False,True) ]:
        if slist is None: continue
        for tag,sample in slist:
            print "tag: %s, sample: %s" %(tag,sample)
#            if isSyst and not 't#bar{t}' in sample[3] : continue
            if tag in skipList:
                print("SKIPPED "+tag)
                continue
            skip = False
            for sl in skipList:
                if ROOT.TString(tag).Contains(sl):
                    skip = True
                    break
            if skip:
                print("SKIPPED "+tag)
                continue
            xsec=sample[0]
            isData=sample[1]
            doFlavourSplitting=sample[6]
            subProcs=[(tag,sample[3],sample[4])]
            if doFlavourSplitting:
                subProcs=[]
                for flav in [(1,sample[3]+'+l'),(4,sample[3]+'+c'),(5,sample[3]+'+b',sample[4])]:
                    subProcs.append(('%d_%s'%(flav[0],tag),flav[1],sample[4]+3*len(subProcs)))
            for sp in subProcs:
                print '%s/%s.root' % ( opt.inDir, sp[0]) 
                fIn=ROOT.TFile.Open('%s/%s.root' % ( opt.inDir, sp[0]) )
                if not fIn : continue

                #fix pileup weighting normalization
                puNormSF=1
                if opt.puNormSF and not isData:
                    puCorrH=fIn.Get(opt.puNormSF)
                    if tag not in rawList:
                        try:
                            nonWgt=puCorrH.GetBinContent(1)
                            wgt=puCorrH.GetBinContent(2)
                            if wgt>0 :
                                puNormSF=nonWgt/wgt
                                if puNormSF>1.3 or puNormSF<0.7 : 
                                    puNormSF=1
                                    report += '%s wasn\'t be scaled as too large SF was found (probably low stats)\n' % sp[0]
                                else :
                                    report += '%s was scaled by %3.3f for pileup normalization\n' % (sp[0],puNormSF)
                        except:
                            print 'Check pu weight control histo',opt.puNormSF,'for',sp[0]

                for tkey in fIn.GetListOfKeys():
                    keyIsSyst=False
                    try:
                        key=tkey.GetName()

                        #filter plots using a selection list
                        keep=False if len(onlyList)>0 else True
                        for pname in onlyList:
                            if opt.strictOnly and pname!=key: continue
                            if not opt.strictOnly and not pname in key: continue
                            keep=True
                            break
                        if not keep: continue
                        histos = []
                        obj=fIn.Get(key)
                        if obj.InheritsFrom('TH2'):
                            if key[-5:]=='_syst':
                                if sample[3]=='t#bar{t}':
                                    keyIsSyst=True
                                    key = key[:-5]
                                    for ybin in xrange(1,obj.GetNbinsY()+1):
                                        for xbin in xrange(0,obj.GetNbinsX()+2):
                                            if math.isnan(obj.GetBinContent(xbin, ybin)):
                                                obj.SetBinContent(xbin, ybin, 0)
                                                obj.SetBinError(xbin, ybin, 0)
                                        weighthist = obj.ProjectionX('_px'+str(ybin), ybin, ybin)
                                        weighthist.SetTitle(sp[1]+' weight '+str(ybin))                                  
                                        weighthist.Draw()
                                        if (weighthist.Integral() > 0): histos.append(weighthist)
                                else:
                                    continue
                            else:
                                histos.append(obj)
                                histos[-1].SetTitle(sp[1])
                        else:
                            histos.append(obj)
                            histos[-1].SetTitle(sp[1])

                        for hist in histos:

                            if "vbfmva" in hist.GetName() and isData and opt.blined:

                                tmpBin = hist.GetXaxis().FindBin(0.2)
                                for iBin in range(tmpBin,hist.GetXaxis().GetNbins()):
                                    hist.SetBinContent(iBin, 0.0000001)
                                
                            if not isData and not '(data)' in sp[1]: 

                                #check if a special scale factor needs to be applied
                                sfVal=1.0                                                 
                                for procToScale in procSF:
                                    if sp[1]==procToScale:                                        
                                        for pcat in procSF[procToScale]:                                    
                                            if pcat not in key: continue
                                            sfVal=procSF[procToScale][pcat][0]
                                            break

                                #scale by lumi
                                lumi=opt.lumi
                                for lSpec in lumiSpecs:                                    
                                    if not lSpec in key.split('_'): continue
                                    lumi=lumiSpecs[lSpec]
                                    break
                                if not opt.rawYields and not tag in rawList:
                                    hist.Scale(xsec*lumi*puNormSF*sfVal)
                                    print (xsec*lumi)
                                    print xsec
             
                            #rebin if needed
                            if opt.rebin>1:  hist.Rebin(opt.rebin)

                            #create new plot if needed
                            if not key in plots : 
                                plots[key]=Plot(key,com=opt.com)
                                plots[key].ratiorange=opt.ratioRange

                            #add process to plot
                            plots[key].add(h=hist,
                                           title=hist.GetTitle(),
                                           color=sp[2],
                                           isData=sample[1],
                                           spImpose=isSignal,
                                           isSyst=(isSyst or keyIsSyst),
                                           doDivideByBinWidth=opt.binWid)
                            
                    except Exception as e:
                        print e
                        pass

    #show plots
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gROOT.SetBatch(True)
    if (not opt.outDir): outDir = opt.inDir+'/plots/'+opt.outName.split('.root')[0]
    else:                outDir = opt.outDir
    os.system('mkdir -p %s' % outDir)
    os.system('rm %s/%s'%(outDir,opt.outName))
    for p in plots:
        plots[p].mcUnc=opt.mcUnc
        if opt.saveLog    : plots[p].savelog=True
        skipPlot=False
        if opt.onlyData and plots[p].dataH is None: skipPlot=True 
        if opt.silent : skipPlot=True
        lumi=opt.lumi
        for lSpec in lumiSpecs:
            if not lSpec in p.split('_'): continue
            lumi=lumiSpecs[lSpec]
            break

        #continue
        if opt.normToData: plots[p].normToData()
        if not skipPlot: plots[p].show(outDir=outDir,lumi=lumi,noStack=opt.noStack,saveTeX=opt.saveTeX)
        plots[p].appendTo('%s/../%s'%(outDir,opt.outName))
        plots[p].reset()

    print '-'*50
    print 'Plots and summary ROOT file can be found in %s' % outDir
    if len(report) : print report
    print '-'*50

        
"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())

