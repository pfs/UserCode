#!/usr/bin/env python

import ROOT
from UEAnalysisHandler import VARS,EVAXES
from collections import defaultdict

FILLS=[1001,3004,3002,1001]
COLORS=[ROOT.kAzure+4, ROOT.kMagenta, ROOT.kGreen+3,  ROOT.kMagenta+2, ROOT.kMagenta-9,ROOT.kRed+1,ROOT.kAzure+7, ROOT.kBlue-7,ROOT.kGray,ROOT.kGray]
MARKERS=[20,22,24,27] #,23,33,20,32,24]
OBSERVABLES=['chmult','sphericity','C','D','aplanarity','chavgpt','chavgpz','chflux','chfluxz']
OBSRANGES={'sphericity':(5e-2,5),
           'aplanarity':(5e-3,30),
           'C':(5e-2,5),
           'D':(5e-3,10),
           'chmult':(1e-4,0.5),
           'chavgpt':(2e-3,2),
           'chavgpz':(2e-3,2),
           'chflux':(5e-5,5e-2),
           'chfluxz':(5e-5,5e-2)}
RATIORANGES={'sphericity':(0.8,1.27),
           'aplanarity':(0.8,1.47),
           'C':(0.8,1.27),
           'D':(0.8,1.47),
           'chmult':(0.5,2.17),
           'chavgpt':(0.6,1.67),
           'chavgpz':(0.7,1.47),
           'chflux':(0.5,1.77),
           'chfluxz':(0.5,1.77)}

SLICES=['nj','chmult',None] #'ptll','nj','chmult','ptttbar','ptll',None]

"""
"""
def shiftGraph(gr,shiftx,shifty):
    x=ROOT.Double(0)
    y=ROOT.Double(0)
    for i in xrange(0,gr.GetN()):
        gr.GetPoint(i,x,y)
        exh,exl=gr.GetErrorXhigh(i),gr.GetErrorXlow(i)
        eyh,eyl=gr.GetErrorYhigh(i),gr.GetErrorYlow(i)
        gr.SetPoint(i,x+shiftx,y*shifty)
        gr.SetPointError(i,exl,exh,eyl*shifty,eyh*shifty)
    return gr

"""
"""
class SimpleUEPlot:
    def __init__(self,h,hname,trueAxis,addXunc):

        self.h=h.Clone(hname)
        self.h.SetDirectory(0)
        self.uncH={}
        self.trueAxis=trueAxis        
        self.addXunc=addXunc
        self.xmean=0
        self.xmeanUnc=0
        self.ymin=1.0e9
        self.ymax=-1.0e9
        self.generateGraph()
        self.garbage=[self.h,self.gr]

    def clear(self):
        for o in self.garbage:
            try:
                o.Delete()
            except:
                pass

    def subtractContribution(self,bkg):
        self.h.Add(bkg,-1)
        self.generateGraph()

    def generateGraph(self):

        norm=self.h.Integral()
        if norm==0:
            self.gr=None
            return self.gr

        self.gr=ROOT.TGraphAsymmErrors()
        self.gr.SetName('%s_gr'%self.h.GetName())
        self.xmean=0
        self.xmeanUnc=0
        for xbin in xrange(1,self.trueAxis.GetNbins()+1):
            cen,hwid=self.trueAxis.GetBinCenter(xbin),self.trueAxis.GetBinWidth(xbin)*0.5
            cts=self.h.GetBinContent(xbin)
            ctsUnc=self.h.GetBinError(xbin)

            self.xmean    += cts*cen
            self.xmeanUnc += ctsUnc*cen

            #include uncertainties (add envelopes in quadrature)
            errLo,errHi=self.h.GetBinError(xbin),self.h.GetBinError(xbin)
            for syst in self.uncH:
                minLo,maxHi=9999999999,-9999999999
                for k in xrange(0,len(self.uncH[syst])):
                    dSignal=expSysts[syst][k].GetBinContent(xbin)
                    if dSignal<minLo : minLo=dSignal
                    if dSignal>maxHi : maxHi=dSignal
                if minLo<0 and maxHi>0:
                    errLo += minLo**2
                    errHi += maxHi**2
                else:
                    env=max(abs(minLo),abs(maxHi))
                    errLo += env**2
                    errHi += env**2
            errLo=ROOT.TMath.Sqrt(errLo)
            errHi=ROOT.TMath.Sqrt(errHi)

            np=self.gr.GetN()
            normCts=cts/(norm*hwid)
            self.gr.SetPoint(np, cen, normCts )
            self.gr.SetPointError(np,
                                  hwid if self.addXunc else 0,
                                  hwid if self.addXunc else 0,
                                  ROOT.TMath.Sqrt(errLo)/(norm*hwid),
                                  ROOT.TMath.Sqrt(errHi)/(norm*hwid))

            self.ymin=min(normCts,self.ymin)
            self.ymax=min(normCts,self.ymax)


        self.xmean    /= norm
        self.xmeanUnc /= norm

        return self.gr
    
"""
"""
class UESummaryPlotInfo:

    """
    """
    def __init__(self,obsAxis,sliceAxis,tag='',isReco=True):
        
        self.obsAxis=obsAxis
        self.sliceAxis=sliceAxis
        self.obs=obsAxis.GetName().split('_')[0]
        self.sliceVar=sliceAxis.GetName().split('_')[0] if sliceAxis else None
        self.tag=''

        self.data = []
        self.signal = []
        self.expSysts = []
        self.signalVars = []

    """
    """
    def addData(self,h,islice):
        self.data.append( (islice,
                           SimpleUEPlot(h,'data_%s_%s'%(str(islice),self.tag),self.obsAxis,False) ) )
                          
    """
    """
    def addSignal(self,h,islice):
        self.signal.append( (islice,
                             SimpleUEPlot(h,'signal_%s_%s'%(str(islice),self.tag),self.obsAxis,True) ) )

    """
    """
    def addExperimentalUncs(self,expSysts,islice):
        self.expSysts.append( (islice,defaultdict(list) ) )
        for x in expSysts:
            for h in expSysts[x]:
                self.expSysts[-1][1][x].append( SimpleUEPlot(h,'signalexp_%s_%s'%(h.GetName(),self.tag),self.obsAxis,True) )

    """
    """
    def addSignalVariations(self,signalVars,islice):
        self.signalVars.append( (islice,{}) )
        for compSet in signalVars:
            self.signalVars[-1][1][compSet]=defaultdict(list)
            for compVarName in signalVars[compSet]:
                for h in signalVars[compSet][compVarName]:
                    self.signalVars[-1][1][compSet][compVarName].append(  SimpleUEPlot(h,'signalth_%s_%s'%(h.GetName(),self.tag),self.obsAxis,True) )
        
    """
    subtracts the background
    """
    def subtractBackground(self,bkg,islice):
        for i,d in self.data:
            if i!=islice: continue
            d.subtractContribution(bkg)


    """
    """
    def showProfile(self,outDir):

        #start the canvas
        c=ROOT.TCanvas('c','c',500,500)
        c.SetTopMargin(0.05)
        c.SetRightMargin(0.0)
        c.SetLeftMargin(0.0)
        c.SetBottomMargin(0.0)

        c.cd()
        p1=ROOT.TPad('p1','p1',0,0,0.8,0.95)
        p1.SetTopMargin(0.01)
        p1.SetRightMargin(0.01)
        p1.SetLeftMargin(0.12)
        p1.SetBottomMargin(0.1)
        p1.Draw()

        c.cd()
        p2=ROOT.TPad('p2','p2',0.8,0,1,0.95)
        p2.SetTopMargin(0.01)
        p2.SetRightMargin(0.05)
        p2.SetLeftMargin(0.01)
        p2.SetBottomMargin(0.1)
        p2.Draw()

        #check ranges and prepare frames
        meanMin,meanMax=1.0e9,-1.0e9
        for _,p in self.data:
            meanMin=min(meanMin,p.xmean)
            meanMax=max(meanMax,p.xmean)
        frame=ROOT.TH1F('frame','frame',self.sliceAxis.GetNbins(),self.sliceAxis.GetXbins().GetArray())
        frame.GetYaxis().SetRangeUser(meanMin*0.5,meanMax*1.5)
        incframe=ROOT.TH1F('incframe','incframe',1,0.5,1.5)
        incframe.GetYaxis().SetRangeUser(meanMin*0.5,meanMax*1.5)
                
        #signal
        mcAvgProfile=ROOT.TGraphAsymmErrors()
        ci=ROOT.TColor.GetColor('#3182bd')
        mcAvgProfile.SetMarkerColor(ci)
        mcAvgProfile.SetLineColor(ci)
        mcAvgProfile.SetFillColor(ci)
        mcAvgProfile.SetFillStyle(1001)
        mcProfile=mcAvgProfile.Clone()
        
        ci=ROOT.TColor.GetColor('#9ecae1')
        mcAvgExpProfile=mcAvgProfile.Clone()
        mcAvgExpProfile.SetFillColor(ci)
        mcExpProfile=mcAvgExpProfile.Clone()

        ci=ROOT.TColor.GetColor('#deebf7')
        mcAvgTotalProfile=mcAvgProfile.Clone()
        mcAvgTotalProfile.SetFillColor(ci)
        mcTotalProfile=mcAvgTotalProfile.Clone()

        for islice,p in self.signal:
            mean,meanUnc=p.xmean,p.xmeanUnc

            #experimental uncertainties
            expUnc=meanUnc**2
            for i,d in self.expSysts:
                if i!=islice: continue   
                for x in d:
                    iexpUnc=0
                    for h in d[x]:
                        if h.gr is None : continue
                        iexpUnc=max(iexpUnc,abs(h.xmean-mean))
                    expUnc += iexpUnc**2
            expUnc=ROOT.TMath.Sqrt(expUnc)

            #total uncertainties
            totalUnc=expUnc**2
            for i,d, in self.signalVars:
                if i!=islice: continue  
                for compSet in d:
                    if compSet!=MAINMC[0] : continue
                    for x in d[compSet]:
                        ithUnc=0
                        for h in d[compSet][x]:
                            if h.gr is None : continue
                            ithUnc=max(ithUnc,abs(h.xmean-mean))
                        totalUnc += ithUnc**2
            totalUnc=ROOT.TMath.Sqrt(totalUnc)

            #add to plots
            if islice==0:
                mcAvgProfile.SetPoint(0,1,mean)
                mcAvgProfile.SetPointError(0,0.5,0.5,meanUnc,meanUnc)
                mcExpProfile.SetPoint(0,1,mean)
                mcExpProfile.SetPointError(0,0.5,0.5,expUnc,expUnc)
                mcAvgTotalProfile.SetPoint(0,1,mean)
                mcAvgTotalProfile.SetPointError(0,0.5,0.5,totalUnc,totalUnc)
            else:
                xcen=self.sliceAxis.GetBinCenter(islice)
                xwid=self.sliceAxis.GetBinWidth(islice)
                mcProfile.SetPoint(islice-1,xcen,mean)
                mcProfile.SetPointError(islice-1,0.5*xwid,0.5*xwid,meanUnc,meanUnc)
                mcExpProfile.SetPoint(islice-1,xcen,mean)
                mcExpProfile.SetPointError(islice-1,0.5*xwid,0.5*xwid,expUnc,expUnc)
                mcTotalProfile.SetPoint(islice-1,xcen,mean)
                mcTotalProfile.SetPointError(islice-1,0.5*xwid,0.5*xwid,totalUnc,totalUnc)

        #data
        dataAvgProfile=ROOT.TGraphAsymmErrors()
        dataAvgProfile.SetMarkerStyle(20)
        dataAvgProfile.SetMarkerColor(1)
        dataAvgProfile.SetLineColor(1)
        dataProfile=dataAvgProfile.Clone()
        for islice,p in self.data:
            mean,meanUnc=p.xmean,p.xmeanUnc
            if islice==0:
                dataAvgProfile.SetPoint(0,1,mean)
                dataAvgProfile.SetPointError(0,0,0,meanUnc,meanUnc)
            else:
                xcen=self.sliceAxis.GetBinCenter(islice)
                dataProfile.SetPoint(islice-1,xcen,mean)
                dataProfile.SetPointError(islice-1,0,0,meanUnc,meanUnc)


        #inclusive frame
        p2.cd()
        incframe.Draw()
        incframe.GetXaxis().SetNdivisions(0)
        incframe.GetXaxis().SetTitleSize(0.15)
        incframe.GetXaxis().SetTitle('Inclusive')
        incframe.GetXaxis().SetTitleOffset(0.25)
        mcAvgTotalProfile.Draw('2')
        mcAvgExpProfile.Draw('2')
        mcAvgProfile.Draw('2')
        dataAvgProfile.Draw('ep')

        #legend for uncertainties
        lines=[]
        for gr in [mcAvgTotalProfile,mcAvgExpProfile,mcAvgProfile]:
            lines.append(ROOT.TLine())
            dx=0.35
            if len(lines)>1: dx=0.24
            if len(lines)>2: dx=0.1
            lines[-1].SetLineWidth(5)
            lines[-1].SetLineColor(gr.GetFillColor())
            lines[-1].DrawLineNDC(0.6-dx,0.93,0.6+0.9*dx,0.93)
        tex=ROOT.TLatex()
        tex.SetTextFont(42)
        tex.SetTextSize(0.1)
        tex.SetNDC()
        tex.DrawLatex(0.1,0.95,'#it{Total Exp. Stat.}')

        p2.RedrawAxis()
                
        #differential frame
        p1.cd()
        frame.Draw()
        frame.GetYaxis().SetTitle('<%s>'%VARS[self.obs][0])
        frame.GetYaxis().SetTitleOffset(1.3)
        frame.GetYaxis().SetTitleSize(0.04)
        frame.GetYaxis().SetLabelSize(0.04)
        frame.GetXaxis().SetTitle(VARS[self.sliceVar][0])
        frame.GetXaxis().SetTitleSize(0.04)
        frame.GetXaxis().SetLabelSize(0.04)
        mcTotalProfile.Draw('2')
        mcExpProfile.Draw('2')
        mcProfile.Draw('2')
        dataProfile.Draw('ep')

        #add legend
        leg=ROOT.TLegend(0.15,0.9,0.94,0.75)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.04)
        leg.AddEntry(dataProfile,'Data','ep')
        leg.AddEntry(mcTotalProfile,'#splitline{Powheg+Pythia8}{CUETP8M2T4}','f')
        leg.Draw()
        
        #the header
        tex2=ROOT.TLatex()
        tex2.SetTextFont(42)
        tex2.SetTextSize(0.045)
        tex2.SetNDC()
        tex2.DrawLatex(0.15,0.93,'#bf{CMS} #it{preliminary}')

        p1.RedrawAxis()

        #the lumi/sqrts
        c.cd()
        tex3=ROOT.TLatex()
        tex3.SetTextFont(42)
        tex3.SetTextSize(0.04)
        tex3.SetNDC()
        tex3.DrawLatex(0.7,0.96,'#scale[0.8]{35.9 fb^{-1} (#sqrt{s}=13 TeV)}')

        # all done
        c.Modified()
        c.Update()
        outName='profile%s_%s_%s'%(self.tag,self.obs,self.sliceVar)
        for ext in ['pdf','png']:
            c.SaveAs('%s/%s.%s'%(outDir,outName,ext))

        p1.Delete()
        p2.Delete()
        c.Delete()



    """
    """
    def showMain(self,outDir):

        c=ROOT.TCanvas('c','c',500,500)
        c.SetTopMargin(0.05)
        c.SetRightMargin(0.02)
        c.SetLeftMargin(0.12)
        c.SetBottomMargin(0.1)
        c.SetLogy()

        frame=ROOT.TH1F('frame','frame',1,self.obsAxis.GetXmin(),self.obsAxis.GetXmax())
        frame.Draw()
        frame.GetYaxis().SetRangeUser(OBSRANGES[self.obs][0],OBSRANGES[self.obs][1])
        frame.GetXaxis().SetLabelSize(0.035)
        frame.GetXaxis().SetTitleSize(0.04)
        frame.GetXaxis().SetTitle(VARS[self.obs][0])
        frame.GetYaxis().SetLabelSize(0.035)
        frame.GetYaxis().SetTitleSize(0.04)
        frame.GetYaxis().SetTitle('1/N dN/d%s'%VARS[self.obs][0])
        frame.GetYaxis().SetTitleOffset(1.3)
        frame.Draw()

        leg=ROOT.TLegend(0.68,0.88,0.95,0.88-(len(self.data)+1)*0.05)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.035)

        #signal
        mcGrs=[]
        for islice,p in self.signal:
            mcGrs.append( p.gr.Clone('signal_%d'%islice) )
            if islice>0:
                shiftGraph( mcGrs[-1],shiftx=0,shifty=ROOT.TMath.Exp(-islice))
            mcGrs[-1].Draw('2')
            mcGrs[-1].SetMarkerStyle(1) 
            ci=ROOT.kAzure+7
            mcGrs[-1].SetMarkerColor(ci)
            mcGrs[-1].SetLineColor(ci)
            mcGrs[-1].SetFillColor(ci)
            mcGrs[-1].SetFillStyle(1001)
      
        #data
        dataGrs=[]
        for islice,p in self.data:
            dataGrs.append( p.gr.Clone('data_%d'%islice) )
            dataGrs[-1].SetTitle('inclusive') 
            if islice>0:
                shiftGraph( dataGrs[-1],shiftx=0,shifty=ROOT.TMath.Exp(-islice))
                dataGrs[-1].SetTitle('%d#leq#scale[0.8]{%s}<%d #scale[0.8]{(#timese^{-%d})}'%(self.sliceAxis.GetBinLowEdge(islice),
                                                                VARS[self.sliceVar][0],
                                                                self.sliceAxis.GetBinUpEdge(islice),
                                                                islice) )
            dataGrs[-1].Draw('p')
            dataGrs[-1].SetMarkerStyle(MARKERS[islice%len(MARKERS)])
            ci=COLORS[islice/len(MARKERS)]
            dataGrs[-1].SetMarkerColor(ci)
            dataGrs[-1].SetLineColor(ci)
            leg.AddEntry(dataGrs[-1],dataGrs[-1].GetTitle(),'ep')

        leg.AddEntry(mcGrs[0],'#scale[0.8]{#splitline{Powheg+Pythia8}{CUETP8M2T4}}','f')
        leg.Draw()
      

        tex=ROOT.TLatex()
        tex.SetTextFont(42)
        tex.SetTextSize(0.04)
        tex.SetNDC()
        tex.DrawLatex(0.7,0.9,'#bf{CMS} #it{preliminary}')
        tex.DrawLatex(0.72,0.96,'#scale[0.8]{35.9 fb^{-1} (#sqrt{s}=13 TeV)}')
        

        outName='main%s_%s_%s'%(self.tag,self.obs,self.sliceVar)
        for ext in ['pdf','png']:
            c.SaveAs('%s/%s.%s'%(outDir,outName,ext))

        c.Delete()




    """
    """
    def showRatio(self,outDir):
        for i in xrange(0,10):
            rpads,rframes=[],[]
            for j in xrange(1,npads):
                c.cd()
                rheight=200./cheight
                rpads.append( ROOT.TPad('p2','p2',0.0,(j-1)*rheight,1.0,j*rheight) )
                rpads[-1].Draw()
                rpads[-1].SetRightMargin(0.02)
                rpads[-1].SetLeftMargin(0.12)
                rpads[-1].SetBottomMargin(0.2 if j==1 else 0.02)
                rpads[-1].SetTopMargin(0.01)
                rpads[-1].Draw()
                rpads[-1].cd()
                rframes.append( ROOT.TH1F('rframe%d'%j,'rframe%d'%j,1,self.obsAxis.GetXmin(),self.obsAxis.GetXmax()) )
                rframes[-1].GetYaxis().SetTitleSize(0.09)
                rframes[-1].GetYaxis().SetLabelSize(0.09)
                rframes[-1].GetYaxis().SetTitle('MC/Data')
                rframes[-1].GetYaxis().SetTitleOffset(0.7)
                rframes[-1].GetXaxis().SetTitleSize(0.09 if j==1 else 0.0)
                rframes[-1].GetXaxis().SetLabelSize(0.09 if j==1 else 0.0)
                rframes[-1].GetXaxis().SetTitle(VARS[self.obs][0])
                rframes[-1].GetYaxis().SetRangeUser(RATIORANGES[self.obs][0],RATIORANGES[self.obs][1])
                rframes[-1].GetYaxis().SetNdivisions(5+100*5)
                rframes[-1].Draw()


            c.cd()
            c.Modified()
            c.Update()
            raw_input()


#
#
#
#    for islice in xrange(1,nslices+1):
#
#        idataGr=ROOT.TGraphErrors();
#        idataGr.SetMarkerStyle(20+(islice-1))
#        idataGr.SetMarkerColor(1)
#        idataGr.SetLineColor(1)
#        idataGr.SetName('data_%d'%islice)
#        title='inc'
#        if sliceAxis:
#            title='[%d,%d]'%(int(sliceAxis.GetBinLowEdge(islice)),int(sliceAxis.GetBinUpEdge(islice)))
#        idataGr.SetTitle(title)
#
#        isignalGr=ROOT.TGraphAsymmErrors();
#        isignalGr.SetFillStyle(1001)
#        ci=ROOT.kAzure+7
#        isignalGr.SetFillColor(ci)
#        isignalGr.SetMarkerColor(ci)
#        isignalGr.SetLineColor(ci)
#        isignalGr.SetMarkerStyle(1)
#
#        isignalRatioGr=isignalGr.Clone('nominalratio')
#
#
#        ratiosGr=ROOT.TMultiGraph()
#        #ratiosLeg=ROOT.TLegend(0.7,0.2,0.95,0.6)
#        ratiosLeg=ROOT.TLegend(0.1,0.85,0.98,0.9)
#        ratiosLeg.SetNColumns( len(signalVars)+1 )
#
#        #build final distribution
#        nobsBins=obsAxis.GetNbins()
#        for xbin in xrange(1,nobsBins+1):
#            cen,hwid=obsAxis.GetBinCenter(xbin),obsAxis.GetBinWidth(xbin)*0.5
#            
#            rawBin=xbin
#            if sliceAxis: rawBin+=nobsBins*(islice-1)
#            
#            #experimental systematics (add envelopes in quadrature)
#            errLo,errHi=0,0
#            for syst in expSysts:
#                minLo,maxHi=9999999999,-9999999999
#                for k in xrange(0,len(expSysts[syst])):
#                    dSignal=expSysts[syst][k].GetBinContent(rawBin)
#                    if dSignal<minLo : minLo=dSignal
#                    if dSignal>maxHi : maxHi=dSignal
#                if minLo<0 and maxHi>0:
#                    errLo += minLo**2
#                    errHi += maxHi**2
#                else:
#                    env=max(abs(minLo),abs(maxHi))
#                    errLo += env**2
#                    errHi += env**2
#            errLo=ROOT.TMath.Sqrt(errLo)
#            errHi=ROOT.TMath.Sqrt(errHi)
#
#            np=idataGr.GetN()
#
#            dataCts,dataUnc=data.GetBinContent(rawBin),data.GetBinError(rawBin)
#            idataGr.SetPoint       (np,   cen,        dataCts/(2*hwid))
#            idataGr.SetPointError  (np,   0,          dataUnc/(2*hwid))
#
#            signalCts=signal.GetBinContent(rawBin)
#            isignalGr.SetPoint     (np,   cen,        signalCts/(2*hwid))
#            isignalGr.SetPointError(np,   hwid, hwid, errLo/(2*hwid), errHi/(2*hwid))
#
#            if dataCts>0 and signalCts>0:
#                iratioVal=signalCts/dataCts
#                iratioUncLo=iratioVal*ROOT.TMath.Sqrt(pow(errLo/signalCts,2)+pow(dataUnc/dataCts,2))
#                iratioUncHi=iratioVal*ROOT.TMath.Sqrt(pow(errHi/signalCts,2)+pow(dataUnc/dataCts,2))
#                isignalRatioGr.SetPoint     (np,   cen,        iratioVal)
#                isignalRatioGr.SetPointError(np,   hwid, hwid, iratioUncLo, iratioUncHi)
#
#        ratiosGr.Add(isignalRatioGr,'2')
#        ratiosLeg.AddEntry(isignalRatioGr,'PW+PY8 CUETP8M2T4','f')
#
#        #dataAvg=averageDistribution(idataGr,obsAxis)
#        #mcAvg={'nominal':averageDistribution(isignalGr,obsAxis)}
#
#        #variations to be compared
#        ivar=0
#        for var in signalVars:
#            ivar+=1
#            isignalVarGr=ROOT.TGraphAsymmErrors();
#            isignalVarGr.SetFillStyle(1001)
#            ci=ROOT.TColor.GetColor(COLORS[ivar-1])
#            isignalVarGr.SetFillColor(ci)
#            isignalVarGr.SetMarkerColor(ci)
#            isignalVarGr.SetLineColor(ci)
#            isignalVarGr.SetMarkerStyle(MARKERS[ivar-1])
#            isignalVarGr.SetMarkerSize(0.8)
#
#            for xbin in xrange(1,nobsBins+1):
#                cen,hwid=obsAxis.GetBinLowEdge(xbin),obsAxis.GetBinWidth(xbin)*0.5
#                cen += 2*hwid*ivar/(len(signalVars)+2)
#
#                rawBin=xbin
#                if sliceAxis: rawBin+=nobsBins*(islice-1)
#
#                minR2data,maxR2data=200,-200
#                for varH in signalVars[var]:
#                    r2data=varH.GetBinContent(rawBin)
#                    minR2data=min(r2data,minR2data)
#                    maxR2data=max(r2data,maxR2data)
#                
#                y=0.5*(minR2data+maxR2data)
#                errHi=max(maxR2data-y,1e-4) #just for display purposes
#                errLo=max(y-minR2data,1e-4)
#                            
#                np=isignalVarGr.GetN()
#                isignalVarGr.SetPoint(np,cen,y)
#                isignalVarGr.SetPointError(np,0,0,errLo,errHi)
#            
#            ratiosGr.Add(isignalVarGr,'pZ')
#            ratiosLeg.AddEntry(isignalVarGr,var,'ep')
#
#
#        

#        ratiosGr.Draw('p')
#        l=ROOT.TLine(frameratio.GetXaxis().GetXmin(),1,frameratio.GetXaxis().GetXmax(),1)
#        l.SetLineColor(ROOT.kBlue)
#        l.Draw()
#        ratiosLeg.SetBorderSize(0)
#        ratiosLeg.SetFillStyle(0)
#        ratiosLeg.SetTextFont(42)
#        ratiosLeg.SetTextSize(0.06)
#        ratiosLeg.Draw()
#        
#        c.cd()
#        c.Modified()
#        c.Update()
#        outname=obs
#        if nslices>1: outname += '%s_%d'%(sliceVar,islice)
#        for ext in ['png','pdf']:
#            c.SaveAs('~/www/TopUE_ReReco2016/%s_smeared.%s'%(outname,ext))
#

