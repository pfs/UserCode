from TopLJets2015.TopAnalysis.Plot import *
import sys

COLORS={
    "t#bar{t}":0,
    "Single top":"#91bfdb",
    "W":"#fee090",
    "DY":"#fc8d59",
    "Multiboson":"#d73027",
    "Data":1
}

def transformToCount(h):

    """ transform an histogram to a single bin counting experiment """

    hcount=ROOT.TH1F(h.GetName()+'_count','',1,0,1)
    hcount.Sumw2()
    err=ROOT.Double(0)    
    hcount.SetBinContent(1,h.IntegralAndError(1,h.GetNbinsX(),err))
    hcount.SetBinError(1,err)
    return hcount
  

def doPlot(plotName,chList,extraText,url,outpName,countOnly=False):

    """ do plot """

    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gROOT.SetBatch(True)

    #open the file
    inF = ROOT.TFile.Open(url)

    plotsPerProc={}
    systList={'toppt':False,
              'l1prefire':True,
              'pu':True,
              #'btagjes':True,
              'btaghf':True,
              'btaglf':True,
              #'btaghfstats1':True,
              #'btaghfstats2':True,
              #'btaglfstats1':True,
              #'btaglfstats2':True,
              #'btagcferr':True,
              'JER':True,
              'bfrag':False,
              'slbr':False,
              "AbsoluteStatJEC"     : True,
              "AbsoluteScaleJEC"    : True,
              "AbsoluteMPFBiasJEC"  : True,
              "FragmentationJEC"    : True,
              "SinglePionECALJEC"   : True,
              "SinglePionHCALJEC"   : True,
              "FlavorPureGluonJEC"  : True,
              "FlavorPureQuarkJEC"  : True,
              "FlavorPureCharmJEC"  : True,
              "FlavorPureBottomJEC" : True,
              "TimePtEtaJEC"        : True,
              "RelativeJEREC1JEC"   : True,
              "RelativeJEREC2JEC"   : True,
              "RelativeJERHFJEC"    : True,
              "RelativePtBBJEC"     : True,
              "RelativePtEC1JEC"    : True,
              "RelativePtEC2JEC"    : True,
              "RelativePtHFJEC"     : True,
              "RelativeBalJEC"      : True,
              "RelativeFSRJEC"      : True,
              "RelativeStatFSRJEC"  : True,
              "RelativeStatECJEC"   : True,
              "RelativeStatHFJEC"   : True,
              "PileUpDataMCJEC"     : True,
              "PileUpPtRefJEC"      : True,
              "PileUpPtBBJEC"       : True,
              "PileUpPtEC1JEC"      : True,
              "PileUpPtEC2JEC"      : True,
              "PileUpPtHFJEC"       : True,
              'muR':True,
              'muF':True,
              'muRmuF':True
              }
    for i in xrange(0,102): systList['PDF%03d'%(i)]=False
    altShapes={}
    scaleVars={}

    for ch in chList :
        pName='%s_%s'%(ch,plotName)        
        pdir=inF.GetDirectory(pName)
        for key in pdir.GetListOfKeys():
            keyName=key.GetName()
            if 'Graph' in keyName : continue

            h=key.ReadObj()
            print keyName,h
            if 'incmlb' in plotName:
                h.GetXaxis().SetRangeUser(20,h.GetXaxis().GetXmax())
            if 'drlb' in plotName:
                h.GetXaxis().SetRangeUser(0.4,5.0)
            if countOnly:  h=transformToCount(h)

            #do systs, if available
            if keyName==pName+'_t#bar{t}':                
                for syst in ['th','exp']:
                    systH=inF.Get('{0}_{1}/{0}_{1}_t#bar{{t}}'.format(pName,syst))
                    try:
                        for ybin in xrange(1,systH.GetNbinsY()+1):

                            #if syst is interesting project it
                            syst=systH.GetYaxis().GetBinLabel(ybin)
                            keep=None
                            for s in systList: 
                                if syst.find(s)==0: 
                                    keep=systList[s]
                                    break
                            if keep is None: continue
                            py=systH.ProjectionX('px',ybin,ybin)
                            if countOnly: py=transformToCount(py)

                            #store effect on normalization
                            totalS=py.Integral()
                            if keep: 
                                if not syst in scaleVars: scaleVars[syst]=0
                                scaleVars[syst]+=totalS

                            #store effect on shape
                            if not syst in altShapes: 
                                altShapes[syst]=py.Clone('%suptotal'%syst)
                                altShapes[syst].Reset('ICE')
                                altShapes[syst].SetDirectory(0)
                            altShapes[syst].Add(py)

                            #all done
                            py.Delete()

                    except Exception,e:
                        print e
                        pass

            title=keyName.replace(pName+'_','')
            if title==keyName : title='Data'
            if not title in plotsPerProc:
                plotsPerProc[title]=h.Clone(title)
                plotsPerProc[title].Reset('ICE')
                plotsPerProc[title].SetTitle(title)
                plotsPerProc[title].SetDirectory(0)
                plotsPerProc[title].SetFillColor(h.GetFillColor())
                plotsPerProc[title].SetLineColor(h.GetLineColor())
                plotsPerProc[title].SetMarkerColor(h.GetMarkerColor())
            plotsPerProc[title].Add(h)
 
    #finalize systematics
    sigH=plotsPerProc['t#bar{t}']
    nbins=sigH.GetNbinsX()+1
    totalExp=sigH.Integral()
    scaleUnc=[0,0]
    scaleVars['lumiup']=(1+.025)*totalExp
    scaleVars['lumidn']=(1-.025)*totalExp
    scaleVars['xsecup']=(1+.051)*totalExp
    scaleVars['xsecdn']=(1-.051)*totalExp
    scaleUnc=[0,0]
    for key in scaleVars:
        sf=scaleVars[key]/totalExp
        scaleUnc[1 if sf>1 else 0] += (1-sf)**2
    scaleUnc = [ROOT.TMath.Sqrt(x) for x in scaleUnc]
    for xbin in xrange(1,nbins+1):
        cts=sigH.GetBinContent(xbin)
        unc=sigH.GetBinError(xbin)
        totalUnc=ROOT.TMath.Sqrt(unc**2+(cts*0.5*(scaleUnc[0]+scaleUnc[1]))**2)        
        sigH.SetBinError(xbin,totalUnc)

    relShapeGr=None
    if len(altShapes):
        systShape=[[0.]*nbins,[0.]*nbins]
        for key in altShapes:
            sf=totalExp/altShapes[key].Integral()
            altShapes[key].Scale(sf)
            for xbin in xrange(1,nbins+1):
                diff=altShapes[key].GetBinContent(xbin)-sigH.GetBinContent(xbin)
                systShape[1 if diff>0 else 0][xbin-1] += diff**2

        totalMCShape=sigH.Clone('totalttshape')
        for xbin in xrange(1,nbins+1):
            cts=sigH.GetBinContent(xbin)+0.5*(ROOT.TMath.Sqrt(systShape[0][xbin-1])+ROOT.TMath.Sqrt(systShape[1][xbin-1]))        
            totalMCShape.SetBinContent(xbin,cts)
            totalMCShape.SetBinError(xbin,0.)
        totalMCShapeGr=ROOT.TGraphErrors()
        
        #normalize by integral
        totalMCShape.Scale(totalExp/totalMCShape.Integral())
        #normalize to first bin with non-zero counts
        #for xbin in xrange(1,nbins+1):
        #    nom=sigH.GetBinContent(xbin)
        #    var=totalMCShape.GetBinContent(xbin)
        #    if nom==0 or var==0 : continue
        #    totalMCShape.Scale(nom/var)
        #    break

        relShapeGr=ROOT.TGraphErrors()
        for xbin in xrange(1,nbins+1):
            xcen=sigH.GetXaxis().GetBinCenter(xbin)
            xwid=sigH.GetXaxis().GetBinWidth(xbin)
            nom=sigH.GetBinContent(xbin)
            var=totalMCShape.GetBinContent(xbin)
            if nom==0 : continue
            r=abs(1-var/nom)
            np=relShapeGr.GetN()
            relShapeGr.SetPoint(np,xcen,1)
            relShapeGr.SetPointError(np,0.5*xwid,r)
        relShapeGr.SetFillStyle(3354) #3001)
        relShapeGr.SetFillColor(ROOT.kGray+3) #ROOT.kRed)
        relShapeGr.SetLineWidth(2)
        totalMCShape.Delete()
    
    #show
    finalOutpName=outpName
    if countOnly : finalOutpName += '_count'
    plot=Plot(finalOutpName)
    plot.savelog=True
    plot.wideCanvas=False
    plot.doMCOverData = False
    plot.ratioFrameFill=3444
    plot.ratioFrameColor=1
    #plot.ratiorange=(0.655,1.385)
    doDivideByBinWidth=False
    if 'mlb' in outpName or 'ptlb' in outpName : doDivideByBinWidth=True
    if relShapeGr : plot.relShapeGr=relShapeGr
    plot.plotformats=['pdf','png']
    for key in ['Data','t#bar{t}','Single top','W','DY','Multiboson']:
        if not key in plotsPerProc : continue
        isData=True if 'Data' in plotsPerProc[key].GetTitle() else False
        color=COLORS[plotsPerProc[key].GetTitle()]
        if key=='DY': plotsPerProc[key].Scale(0.83)
        plot.add(plotsPerProc[key],
                 plotsPerProc[key].GetTitle(),
                 color,
                 isData,
                 False,
                 False,
                 doDivideByBinWidth
             )
    plot.finalize()
    plot.mcUnc=0.0

    totalMC=sigH.Clone('tmptotal')
    totalMC.Reset('ICE')
    for h in plot.mc: totalMC.Add(plot.mc[h])
    plot.normUncGr=ROOT.TGraphErrors(totalMC)
    plot.normUncGr.SetFillStyle(3354)
    plot.normUncGr.SetFillColor(ROOT.kGray+3)
    plot.normUncGr.SetMarkerStyle(1)
    plot.normUncGr.SetLineColor(1)
    plot.normUncGr.SetName("normuncgr")
    plot.normUncGr.SetTitle('Stat #oplus syst')
    totalMC.Delete()
    plot.show(outDir="plots/",lumi=35882,extraText=extraText,saveTeX=countOnly)

def main():

    os.system('mkdir -p plots')

    plots=sys.argv[1].split(',')
    chList="eehighpt1b,eehighpt2b,eelowpt1b,eelowpt2b"
    chList+=",mmhighpt1b,mmhighpt2b,mmlowpt1b,mmlowpt2b"
    chList+=",emhighpt1b,emhighpt2b,emlowpt1b,emlowpt2b"
    if len(sys.argv)>2: chList=sys.argv[2]
    
    plotter='root://eoscms//eos/cms/store/cmst3/group/top/TOP17010/final/0c522df/plots/plotter.root'
    if len(sys.argv)>3: plotter=sys.argv[3]

    countOnly=False
    if len(sys.argv)>4: countOnly=True if sys.argv[4]=="True" else False

    for poutp in plots:
        p,outp=poutp.split(':')
        print p,outp
        extraText=''
        if 'em' in chList : extraText+='e#mu, '
        if 'ee' in chList : extraText+='ee, '
        if 'mm' in chList : extraText+='#mu#mu, '
        extraText=extraText[0:-2]+'\\'
        if 'lowpt' in chList  and not 'highpt' in chList: extraText += 'p_{T}(lepton,jet)<100 GeV\\'
        if 'highpt' in chList and not 'lowpt'  in chList: extraText += 'p_{T}(lepton,jet)#geq 100 GeV\\'
        if '1b' in chList and not '2b' in chList : extraText += '=1b-tag'
        if '2b' in chList and not '1b' in chList : extraText += '#geq2 b-tags'
        if '1b' in chList and     '2b' in chList : extraText += '#geq1 b-tags'
        doPlot(p,chList.split(','),extraText,plotter,outp,countOnly=countOnly)
        

if __name__ == "__main__":
    sys.exit(main())
