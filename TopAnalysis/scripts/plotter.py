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
    parser.add_option(     '--mcUnc',        dest='mcUnc'  ,      help='common MC related uncertainty (e.g. lumi)',        default=0,              type=float)
    parser.add_option(     '--com',          dest='com'  ,        help='center of mass energy',                            default='13 TeV',       type='string')
    parser.add_option('-j', '--json',        dest='json'  ,      help='json with list of files',        default=None,              type='string')
    parser.add_option(      '--signalJson',  dest='signalJson',  help='signal json list',               default=None,              type='string')
    parser.add_option('-i', '--inDir',       dest='inDir' ,      help='input directory',                default=None,              type='string')
    parser.add_option('-O', '--outDir',      dest='outDir' ,     help='output directory',                default=None,              type='string')
    parser.add_option('-o', '--outName',     dest='outName' ,    help='name of the output file',        default='plotter.root',    type='string')
    parser.add_option(      '--noStack',     dest='noStack',     help='don\'t stack distributions',     default=False,             action='store_true')
    parser.add_option(      '--saveLog',     dest='saveLog' ,    help='save log versions of the plots', default=False,             action='store_true')
    parser.add_option(      '--silent',      dest='silent' ,     help='only dump to ROOT file',         default=False,             action='store_true')
    parser.add_option(      '--onlyData',    dest='onlyData' ,   help='only plots containing data',     default=False,             action='store_true')
    parser.add_option(      '--saveTeX',     dest='saveTeX' ,    help='save as tex file as well',       default=False,             action='store_true')
    parser.add_option(      '--rebin',       dest='rebin',       help='rebin factor',                   default=1,                 type=int)
    parser.add_option('-l', '--lumi',        dest='lumi' ,       help='lumi to print out',              default=12900,              type=float)
    parser.add_option(      '--lumiSpecs',   dest='lumiSpecs',   help='lumi specifications for some channels [tag:lumi,tag2:lumi2,...]', default=None,       type=str)
    parser.add_option(      '--only',        dest='only',        help='plot only these (csv)',          default='',                type='string')
    parser.add_option(      '--puNormSF',    dest='puNormSF',    help='Use this histogram to correct pu weight normalization', default=None, type='string')
    parser.add_option(      '--procSF',      dest='procSF',      help='Use this to scale a given process component e.g. "W":.wjetscalefactors.pck,"DY":dyscalefactors.pck', default=None, type='string')
    (opt, args) = parser.parse_args()

    #read list of samples
    jsonFile = open(opt.json,'r')
    samplesList=json.load(jsonFile, encoding='utf-8', object_pairs_hook=OrderedDict).items()
    jsonFile.close()

    #read list of signal samples
    signalSamplesList=None
    try:
        jsonFile = open(opt.signalJson,'r')
        signalSamplesList=json.load(jsonFile, encoding='utf-8', object_pairs_hook=OrderedDict).items()
        jsonFile.close()
    except:
        pass

    #lumi specifications per tag
    lumiSpecs={}
    if opt.lumiSpecs:
        for spec in opt.lumiSpecs.split(','):
            tag,lumi=spec.split(':')
            lumiSpecs[tag]=float(lumi)

    #proc SF
    procSF={}
    if opt.procSF:
        procList=opt.procSF.split(',')
        for newProc in procList:
            proc,cacheUrl=newProc.split(':')
            if not os.path.isfile(cacheUrl) : continue
            cache=open(cacheUrl,'r')
            procSF[proc]=pickle.load(cache)
            cache.close()
            print 'Scale factors added for',proc

    onlyList=opt.only.split(',')

    #read plots 
    plots=OrderedDict()

    report=''
    for slist,isSignal in [ (samplesList,False),(signalSamplesList,True) ]:
        if slist is None: continue
        for tag,sample in slist: 
            xsec=sample[0]
            isData=sample[1]
            doFlavourSplitting=sample[6]
            subProcs=[(tag,sample[3],sample[4])]
            if doFlavourSplitting:
                subProcs=[]
                for flav in [(1,sample[3]+'+l'),(4,sample[3]+'+c'),(5,sample[3]+'+b',sample[4])]:
                    subProcs.append(('%d_%s'%(flav[0],tag),flav[1],sample[4]+3*len(subProcs)))
            for sp in subProcs:

                fIn=ROOT.TFile.Open('%s/%s.root' % ( opt.inDir, sp[0]) )
                if not fIn : continue

                #fix pileup weighting normalization
                puNormSF=1
                if opt.puNormSF and not isData:
                    puCorrH=fIn.Get(opt.puNormSF)
                    nonWgt=puCorrH.GetBinContent(1)
                    wgt=puCorrH.GetBinContent(2)
                    if wgt>0 :
                        puNormSF=nonWgt/wgt
                        if puNormSF>1.3 or puNormSF<0.7 : 
                            puNormSF=1
                            report += '%s wasn\'t be scaled as too large SF was found (probably low stats)\n' % sp[0]
                        else :
                            report += '%s was scaled by %3.3f for pileup normalization\n' % (sp[0],puNormSF)

                try:
                    for tkey in fIn.GetListOfKeys():
                        key=tkey.GetName()
                        keep=False if len(onlyList)>0 else True
                        for pname in onlyList: 
                            if pname in key: keep=True
                        if not keep: continue
                        obj=fIn.Get(key)
                        if not obj.InheritsFrom('TH1') : continue
                        if not obj.InheritsFrom('TH2') : fixExtremities(obj, False, False)

                        if not isData and not '(data)' in sp[1]: 
                            sfVal=1.0
                            for procToScale in procSF:
                                if sp[1]==procToScale:
                                    #if procToScale in sp[1]:
                                    for pcat in procSF[procToScale]:                                    
                                        if pcat not in key: continue
                                        sfVal=procSF[procToScale][pcat][0]
                                     #print '[Warning] Applying scale factor for ',sp[1],key,sfVal
                           
                            #scale by lumi
                            lumi=opt.lumi
                            for tag in lumiSpecs:
                                if not tag in key: continue
                                lumi=lumiSpecs[tag]
                                #print '[Warning] applying',lumi,'for',tag
                                break
                                        
                            obj.Scale(xsec*lumi*puNormSF*sfVal)                    
                        if opt.rebin>1:  obj.Rebin(opt.rebin)
                        if not key in plots : plots[key]=Plot(key,com=opt.com)
                        plots[key].add(h=obj,title=sp[1],color=sp[2],isData=sample[1],spImpose=isSignal)
                except:
                    pass

    #show plots
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gROOT.SetBatch(True)
    if (not opt.outDir): outDir = opt.inDir+'/plots'
    else:                outDir = opt.outDir
    os.system('mkdir -p %s' % outDir)
    os.system('rm %s/%s'%(outDir,opt.outName))
    for p in plots : 
        plots[p].mcUnc=opt.mcUnc
        if opt.saveLog    : plots[p].savelog=True
        skipPlot=False
        if opt.onlyData and plots[p].dataH is None: skipPlot=True 
        if opt.silent : skipPlot=True
        if not skipPlot : plots[p].show(outDir=outDir,lumi=opt.lumi,noStack=opt.noStack,saveTeX=opt.saveTeX)
        plots[p].appendTo('%s/%s'%(outDir,opt.outName))
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

