import optparse
import ROOT
import sys
from TopLJets2015.TopAnalysis.Plot import *

def normalize(h):
    try:
        if h.Integral()<=0 : return
        h.Scale(1./h.Integral(0,h.GetNbinsX()+1))
    except:
        h=None
        pass

def getPlotsIn(inF,dirname):
    """gets the data and mc sums in a given directory"""
    data,mc=None,None
    for key in inF.Get(dirname).GetListOfKeys():
        name=key.GetName()
        if 'Graph' in name : continue
        if name==dirname:
            data=key.ReadObj()
            data.SetDirectory(0)
        else:
            h=key.ReadObj()
            if not mc:
                mc=h.Clone(dirname+'_mc')
                mc.SetDirectory(0)
                mc.Reset('ICE')
            mc.Add(h)
    return data,mc

def computeVBFRatios(inputList,triggerBased=False):
    """opens a plotter and computes the gamma/Z ratios for data and MC
       for different distributions"""

    ratios={}

    for t,f in inputList:

        inF=ROOT.TFile.Open(f)
        for key in inF.GetListOfKeys():
            name=key.GetName()

            doRatioTo=None
            ratioTitle=''
            if not triggerBased:
                if 'A_' in name:
                    doRatioTo=('A_','MM_')
                    ratioTitle='#gamma / Z#rightarrow #mu#mu ratio'
            else:
                if 'HighPtVBFA_' in name:
                    doRatioTo=('HighPtVBFA_','HighPtOfflineVBFA_')
                    ratioTitle='High p_{T} VBF #gamma / High p_{T} #gamma'

            if not doRatioTo: continue
            try:
                dataNum,mcNum=getPlotsIn(inF,name)
                dataDen,mcDen=getPlotsIn(inF,name.replace(doRatioTo[0],doRatioTo[1]))
            except:
                print 'Failed at',name,'for',t
                continue

            if triggerBased:
                for h in [dataNum,mcNum,dataDen,mcDen] : normalize(h)

            if not name in ratios: ratios[name]=[]

            if dataNum : 
                dataNum.Divide(dataDen)
                dataNum.GetYaxis().SetTitle(ratioTitle)
                if len(t)>0 : dataNum.SetTitle('Data (%s)'%t)
                ratios[name].append(dataNum)
            if mcNum   : 
                mcNum.Divide(mcDen)
                mcNum.GetYaxis().SetTitle(ratioTitle)
                if len(t)>0 : mcNum.SetTitle('MC (%s)'%t)
                ratios[name].append(mcNum)

    inF.Close()

    return ratios

def showRatios(ratios,outUrl):
    """shows the ratio plots"""

    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gROOT.SetBatch(True)

    os.system('rm %s'%outUrl)
    outDir=os.path.dirname(outUrl)

    COLOURS=[1,'#f4a582','#bababa','#abdda4']
    for key in ratios:
        p=Plot(key+'_ratio',com='13 TeV')
        p.doPoissonErrorBars=False
        p.ratiorange=(0.68,1.34)
        p.spimposeWithErrors=True
        ic=0
        for h in ratios[key]:            
            p.add(h=h,title=h.GetTitle(),color=COLOURS[ic],isData=False,spImpose=True,isSyst=False)
            ic+=1
        p.show(outDir=outDir,lumi=41400,noStack=False,saveTeX=False,noRatio=True)
        p.appendTo(outUrl)
        p.reset()

def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i',  dest='input',  help='input plotters (CSV list) [%default]',  default=None,        type='string')
    parser.add_option('-o',  dest='output', help='output plotter [%default]', default='ratio_plotter.root',    type='string')
    parser.add_option('--titles',  dest='titles', help='titles (CSV list) [%default]',    default=None,        type='string')
    parser.add_option('-t',  dest='triggerBased', help='trigger-based ratio [%default]',  default=False,       action='store_true')
    (opt, args) = parser.parse_args()

    inputList=list(zip(opt.titles.split(','),opt.input.split(',')))
    ratios=computeVBFRatios(inputList,opt.triggerBased)
    showRatios(ratios,opt.output)

if __name__ == "__main__":
    sys.exit(main())
