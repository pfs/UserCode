import ROOT
import math
import os,sys

from TopLJets2015.TopAnalysis.rounding import *

"""
increments the first and the last bin to show the under- and over-flows
"""
def fixExtremities(h,addOverflow=True,addUnderflow=True):
    if addUnderflow :
        fbin  = h.GetBinContent(0) + h.GetBinContent(1)
	fbine = ROOT.TMath.Sqrt(h.GetBinError(0)*h.GetBinError(0) + h.GetBinError(1)*h.GetBinError(1))
	h.SetBinContent(1,fbin)
	h.SetBinError(1,fbine)
	h.SetBinContent(0,0)
	h.SetBinError(0,0)
    if addOverflow:
        nbins = h.GetNbinsX();
	fbin  = h.GetBinContent(nbins) + h.GetBinContent(nbins+1)
	fbine = ROOT.TMath.Sqrt(h.GetBinError(nbins)*h.GetBinError(nbins)  + h.GetBinError(nbins+1)*h.GetBinError(nbins+1))
	h.SetBinContent(nbins,fbin)
	h.SetBinError(nbins,fbine)
	h.SetBinContent(nbins+1,0)
	h.SetBinError(nbins+1,0)



"""
A wrapper to store data and MC histograms for comparison
"""
class Plot(object):

    def __init__(self,name,com='13 TeV'):
        self.name = name
        self.com=com
        self.wideCanvas = True if 'ratevsrun' in self.name else False
        self.mc = {}
        self.spimpose={}
        self.dataH = None
        self.data = None
        self._garbageList = []
        self.plotformats = ['pdf','png']
        self.savelog = False
        self.doChi2 = False
        self.ratiorange = (0.76,1.24)
        self.frameMin=0.1
        self.frameMax=1.45
        self.mcUnc=0

    def add(self, h, title, color, isData,spImpose):

        if 'ratevsrun' in self.name and not isData: return

        h.SetTitle(title)
        if isData:
            try:
                self.dataH.Add(h)
            except:
                self.dataH=h
                self.dataH.SetDirectory(0)
                self.dataH.SetMarkerStyle(20)
                self.dataH.SetMarkerSize(1.4)
                self.dataH.SetMarkerColor(color)
                self.dataH.SetLineColor(ROOT.kBlack)
                self.dataH.SetLineWidth(2)
                self.dataH.SetFillColor(0)
                self.dataH.SetFillStyle(0)
                self._garbageList.append(h)
        else:
            try:
                if spImpose : self.spimpose[title].Add(h)
                else        : self.mc[title].Add(h)
            except:
                h.SetName('%s_%s' % (h.GetName(), title ) )
                h.SetDirectory(0)
                h.SetMarkerStyle(1)
                h.SetMarkerColor(color)
                if spImpose : 
                    self.spimpose[title]=h
                    h.SetFillStyle(0)
                    h.SetLineColor(color)
                    h.SetLineWidth(2)
                else : 
                    h.SetLineColor(ROOT.kBlack)
                    h.SetLineWidth(1)
                    h.SetFillColor(color)
                    h.SetFillStyle(1001)
                    self.mc[title]=h
                self._garbageList.append(h)

    def finalize(self):
        self.data = convertToPoissonErrorGr(self.dataH)

    def appendTo(self,outUrl):
        outF = ROOT.TFile.Open(outUrl,'UPDATE')
        if not outF.cd(self.name):
            outDir = outF.mkdir(self.name)
            outDir.cd()
        for m in self.mc :
            self.mc[m].Write(self.mc[m].GetName(), ROOT.TObject.kOverwrite)
        for m in self.spimpose:
            self.spimpose[m].Write(self.spimpose[m].GetName(), ROOT.TObject.kOverwrite)
        if self.dataH :
            self.dataH.Write(self.dataH.GetName(), ROOT.TObject.kOverwrite)
        if self.data :
            self.data.Write(self.data.GetName(), ROOT.TObject.kOverwrite)
        outF.Close()

    def reset(self):
        for o in self._garbageList:
            try:
                o.Delete()
            except:
                pass

    def show(self, outDir,lumi,noStack=False,saveTeX=False,extraText=None):

        if len(self.mc)<2 and self.dataH is None:
            print '%s has 0 or 1 MC!' % self.name
            return

        if len(self.mc)>0 and self.mc.values()[0].InheritsFrom('TH2') :
            print 'Skipping TH2'
            return

        cwid=1000 if self.wideCanvas else 500
        c = ROOT.TCanvas('c','c',cwid,500)
        c.SetBottomMargin(0.0)
        c.SetLeftMargin(0.0)
        c.SetTopMargin(0)
        c.SetRightMargin(0.00)
        
        #holds the main plot
        c.cd()
        p1 = None
        if self.dataH:
            p1=ROOT.TPad('p1','p1',0.0,0.2,1.0,1.0) if cwid!=1000 else ROOT.TPad('p1','p1',0.0,0.0,1.0,1.0)
            p1.SetRightMargin(0.05)
            p1.SetLeftMargin(0.12)
            p1.SetTopMargin(0.01)
            p1.SetBottomMargin(0.01)
            if self.wideCanvas and len(self.mc)==0 : p1.SetBottomMargin(0.12)
        else:
            p1=ROOT.TPad('p1','p1',0.0,0.0,1.0,1.0)
            p1.SetRightMargin(0.05)
            p1.SetLeftMargin(0.12)
            p1.SetTopMargin(0.1)
            p1.SetBottomMargin(0.12)
        p1.Draw()

        p1.SetGridx(False)
        p1.SetGridy(False) #True)
        self._garbageList.append(p1)
        p1.cd()

        # legend
        iniy=0.9 if self.wideCanvas else 0.85
        dy=0.1 
        ndy= 0.2*max(len(self.mc)-2,1)
        inix,dx =0.45,0.5
        if noStack:
            inix,dx=0.6,0.35
            iniy,dy,ndy=0.95,0.02,len(self.mc)

        leg = ROOT.TLegend(inix, iniy-dy*ndy, inix+dx, iniy+0.06)

        leg.SetBorderSize(0)
        leg.SetFillStyle(0)        
        leg.SetTextFont(42)
        leg.SetTextSize(0.045 if self.wideCanvas else 0.035)
        nlegCols = 0

        if self.dataH is not None:
            if self.data is None: self.finalize()
            leg.AddEntry( self.data, self.data.GetTitle(),'p')
            nlegCols += 1
        for h in self.mc:
            
            #compare
            if noStack:
                refH=self.mc.values()[0]
                if refH!=self.mc[h] and self.doChi2:
                    chi2=refH.Chi2Test( self.mc[h], 'WW CHI2')
                    pval=refH.Chi2Test( self.mc[h], 'WW')     
                    self.mc[h].SetTitle('#splitline{%s}{#chi^{2}=%3.1f (p-val: %3.3f)}'%(self.mc[h].GetTitle(),chi2,pval))
                else:
                    refH.SetLineWidth(2)

            leg.AddEntry(self.mc[h], self.mc[h].GetTitle(), 'f')
            nlegCols += 1
        if nlegCols ==0 :
            print '%s is empty'%self.name
            return

        if not noStack:
            leg.SetNColumns(ROOT.TMath.Min(nlegCols/2,3))

        # Build the stack to plot from all backgrounds
        totalMC = None
        stack = ROOT.THStack('mc','mc')
        for h in self.mc:

            if noStack:
                self.mc[h].SetFillStyle(0)
                self.mc[h].SetLineColor(self.mc[h].GetFillColor())
                
            stack.Add(self.mc[h],'hist')
            
            try:
                totalMC.Add(self.mc[h])
            except:
                totalMC = self.mc[h].Clone('totalmc')
                self._garbageList.append(totalMC)
                totalMC.SetDirectory(0)

        #test for null plots
        if totalMC :
            if totalMC.Integral()==0:
                if self.dataH is None : return
                if self.dataH.Integral()==0: return
        elif self.dataH is None : return
        elif self.dataH.Integral()==0 : return 


        frame = totalMC.Clone('frame') if totalMC is not None else self.dataH.Clone('frame')
        frame.Reset('ICE')
        if noStack:
            maxY=stack.GetStack().At(0).GetMaximum()/1.25
        elif totalMC:
            maxY = totalMC.GetMaximum() 
            if self.dataH:
                if maxY<self.dataH.GetMaximum():
                    maxY=self.dataH.GetMaximum()
        else:
            maxY=self.dataH.GetMaximum()

        frame.GetYaxis().SetRangeUser(self.frameMin,self.frameMax*maxY)

        frame.SetDirectory(0)
        frame.Reset('ICE')
        self._garbageList.append(frame)
        frame.GetYaxis().SetTitleSize(0.045)
        frame.GetYaxis().SetLabelSize(0.04)
        frame.GetYaxis().SetNoExponent()
        frame.GetYaxis().SetTitleOffset(1.3)
        if self.dataH:
            frame.GetXaxis().SetTitleSize(0.0)
            frame.GetXaxis().SetLabelSize(0.0)
        else :
            frame.GetXaxis().SetTitleSize(0.045)
            frame.GetXaxis().SetLabelSize(0.04)
        if self.wideCanvas and totalMC is None : 
            frame.GetXaxis().SetLabelSize(0.03)
            frame.GetXaxis().SetTitleSize(0.035)
        frame.Draw()
        if totalMC is not None   : 
            if noStack: stack.Draw('nostack same')
            else      : stack.Draw('hist same')
        for m in self.spimpose:
            self.spimpose[m].Draw('histsame')
            leg.AddEntry(self.spimpose[m],self.spimpose[m].GetTitle(),'l')
        if self.data is not None : self.data.Draw('p')


        leg.Draw()
        txt=ROOT.TLatex()
        txt.SetNDC(True)
        txt.SetTextFont(42)
        txt.SetTextSize(0.05)
        txt.SetTextAlign(12)
        iniy=0.9 if self.wideCanvas else 0.95
        inix=0.12 if noStack else 0.18
        if lumi<100:
            txt.DrawLatex(inix,iniy,'#bf{CMS} #it{Preliminary} #scale[0.8]{%3.1f pb^{-1} (%s)}' % (lumi,self.com) )
        else:
            txt.DrawLatex(inix,iniy,'#bf{CMS} #it{Preliminary} #scale[0.8]{%3.1f fb^{-1} (%s)}' % (lumi/1000.,self.com) )
        try:
            extraCtr=1
            for extra in extraText.split('\\'):
                txt.DrawLatex(inix,iniy-0.05*extraCtr,'#scale[0.8]{#it{%s}}'%extra)
                extraCtr+=1
        except:
            pass


        #holds the ratio
        c.cd()
        if len(self.mc)>0 and self.dataH:
            p2 = ROOT.TPad('p2','p2',0.0,0.0,1.0,0.2)
            p2.Draw()
            p2.SetBottomMargin(0.4)
            p2.SetRightMargin(0.05)
            p2.SetLeftMargin(0.12)
            p2.SetTopMargin(0.01)
            p2.SetGridx(False)
            p2.SetGridy(True)
            self._garbageList.append(p2)
            p2.cd()
            ratioframe=frame.Clone('ratioframe')
            ratioframe.GetYaxis().SetTitle('Ratio')
            ratioframe.GetYaxis().SetRangeUser(self.ratiorange[0], self.ratiorange[1])
            self._garbageList.append(ratioframe)
            ratioframe.GetYaxis().SetNdivisions(5)
            ratioframe.GetYaxis().SetLabelSize(0.18)        
            ratioframe.GetYaxis().SetTitleSize(0.2)
            ratioframe.GetYaxis().SetTitleOffset(0.25)
            ratioframe.GetXaxis().SetLabelSize(0.15)
            ratioframe.GetXaxis().SetTitleSize(0.2)
            ratioframe.GetXaxis().SetTitleOffset(0.8)
            ratioframe.SetFillStyle(3001)
            ratioframe.SetFillColor(ROOT.kGray+2)
            totalMCnoUnc=totalMC.Clone('totalMCnounc')
            self._garbageList.append(totalMCnoUnc)
            for xbin in xrange(1,totalMC.GetNbinsX()+1):
                ratioframe.SetBinContent(xbin,1)
                val=totalMC.GetBinContent(xbin)
                totalMCnoUnc.SetBinError(xbin,0.)
                if val==0 : continue
                totalUnc=ROOT.TMath.Sqrt((totalMC.GetBinError(xbin)/val)**2+self.mcUnc**2)
                ratioframe.SetBinError(xbin,totalUnc)
            ratioframe.Draw('e2')
            try:
                ratio=self.dataH.Clone('ratio')
                ratio.SetDirectory(0)
                self._garbageList.append(ratio)
                ratio.Divide(totalMCnoUnc)
                gr=ROOT.TGraphAsymmErrors(ratio)
                gr.SetMarkerStyle(self.data.GetMarkerStyle())
                gr.SetMarkerSize(self.data.GetMarkerSize())
                gr.SetMarkerColor(self.data.GetMarkerColor())
                gr.SetLineColor(self.data.GetLineColor())
                gr.SetLineWidth(self.data.GetLineWidth())
                gr.Draw('p')
            except:
                pass

        #all done
        c.cd()
        c.Modified()
        c.Update()

        #save
        for ext in self.plotformats : c.SaveAs(os.path.join(outDir, self.name+'.'+ext))
        if self.savelog:
            p1.cd()
            frame.GetYaxis().SetRangeUser(1,maxY*50)
            p1.SetLogy()
            c.cd()
            c.Modified()
            c.Update()
            for ext in self.plotformats : c.SaveAs(os.path.join(outDir, self.name+'_log.'+ext))

        if saveTeX : self.convertToTeX(outDir=outDir)


    def convertToTeX(self, outDir):
        if len(self.mc)==0:
            print '%s is empty' % self.name
            return

        f = open(outDir+'/'+self.name+'.dat','w')
        f.write('------------------------------------------\n')
        f.write("Process".ljust(20),)
        f.write("Events after each cut\n")
        f.write('------------------------------------------\n')

        tot ={}
        err = {}
        f.write(' '.ljust(20),)
        try:
            for xbin in xrange(1,self.mc.values()[0].GetXaxis().GetNbins()+1):
                pcut=self.mc.values()[0].GetXaxis().GetBinLabel(xbin)
                f.write(pcut.ljust(40),)
                tot[xbin]=0
                err[xbin]=0
        except:
            pass
        f.write('\n')
        f.write('------------------------------------------\n')

        for pname in self.mc:
            h = self.mc[pname]
            f.write(pname.ljust(20),)

            for xbin in xrange(1,h.GetXaxis().GetNbins()+1):
                itot=h.GetBinContent(xbin)
                ierr=h.GetBinError(xbin)
                pval=' & %s'%toLatexRounded(itot,ierr)
                f.write(pval.ljust(40),)
                tot[xbin] = tot[xbin]+itot
                err[xbin] = err[xbin]+ierr*ierr
            f.write('\n')

        f.write('------------------------------------------\n')
        f.write('Total'.ljust(20),)
        for xbin in tot:
            pval=' & %s'%toLatexRounded(tot[xbin],math.sqrt(err[xbin]))
            f.write(pval.ljust(40),)
        f.write('\n')

        if self.dataH is None: return
        f.write('------------------------------------------\n')
        f.write('Data'.ljust(20),)
        for xbin in xrange(1,self.dataH.GetXaxis().GetNbins()+1):
            itot=self.dataH.GetBinContent(xbin)
            pval=' & %d'%itot
            f.write(pval.ljust(40))
        f.write('\n')
        f.write('------------------------------------------\n')
        f.close()



"""
converts a histogram to a graph with Poisson error bars
"""
def convertToPoissonErrorGr(h):
    try:
        htype=h.ClassName()
    except:
        return None
    if htype.find('TH1')<0 : return None

    #check https://twiki.cern.ch/twiki/bin/view/CMS/PoissonErrorBars
    alpha = 1 - 0.6827;
    grpois = ROOT.TGraphAsymmErrors(h);
    for i in xrange(0,grpois.GetN()+1) :
        N = grpois.GetY()[i]
        if N<200 :
            L = 0
            if N>0 : L = ROOT.Math.gamma_quantile(alpha/2,N,1.)
            U = ROOT.Math.gamma_quantile_c(alpha/2,N+1,1)
            grpois.SetPointEYlow(i, N-L)
            grpois.SetPointEYhigh(i, U-N)
        else:
            grpois.SetPointEYlow(i, math.sqrt(N))
            grpois.SetPointEYhigh(i,math.sqrt(N))
    return grpois
