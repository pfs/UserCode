import ROOT
import os
import sys
import optparse
import numpy as np
import pickle
from collections import OrderedDict
from UEAnalysisHandler import VARS,SYSTS

"""
"""
def showSystsSummary(systUpH,systDnH,systSummary,outdir,outname) :
    #show systematics stack
    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetRightMargin(0.02)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.1)
    leg=ROOT.TLegend(0.45,0.95,0.95,0.75)
    leg.SetNColumns(3)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.03)
    nhistos=len(systUpH)
    frame=systUpH[0].Clone('frame')
    frame.Reset('ICE')
    frame.Draw()
    frame.GetYaxis().SetRangeUser(-20.,20.)
    frame.GetYaxis().SetTitle('Cumulative relative uncertainty (%)')
    frame.GetYaxis().SetTitleOffset(1.2)
    frame.GetXaxis().SetTitle('Bin')
    systUpStack=ROOT.THStack()
    systDnStack=ROOT.THStack()
    for i in xrange(0,nhistos):        
        systUpStack.Add(systUpH[nhistos-1-i],'hist') #PFC hist')
        systDnStack.Add(systDnH[nhistos-1-i],'hist') #PFC hist')
        leg.AddEntry(systUpH[nhistos-1-i],systUpH[nhistos-1-i].GetTitle(),"f")
    systUpStack.Draw("nostack same")
    systDnStack.Draw("nostack same")
    leg.Draw()
    l=ROOT.TLine(frame.GetXaxis().GetXmin(),0,frame.GetXaxis().GetXmax(),0)
    l.SetLineColor(1)
    l.Draw()
    tex=ROOT.TLatex()
    tex.SetTextFont(42)
    tex.SetTextSize(0.04)
    tex.SetNDC()
    tex.DrawLatex(0.15,0.9,'#bf{CMS} #it{preliminary}')
    tex.DrawLatex(0.75,0.95,'#scale[0.7]{35.9 fb^{-1} (#sqrt{s}=13 TeV)}')
    c.Modified()
    c.Update()
    for ext in ['png','pdf']:
        c.SaveAs('%s/%s_reluncertainty.%s'%(outdir,outname,ext))

    #save summary to file
    with open('%s/%s_reluncertainty.dat'%(outdir,outname),'w') as fout:
        fout.write('Source & <Rel. unc.> \\\\\n')
        for s in systSummary:
            fout.write('%20s & %3.1f \\\\\n'%(s,systSummary[s]))

"""
"""
def addSystematics(data,fIn,finalSystList,obsAxis,sliceAxis,outdir,outname,isReco=False):

    #normalize to unity
    normalizePerSlice(data,obsAxis,sliceAxis)

    nbinsx=data.GetNbinsX()

    #loop over systematic sources
    systSummary={}
    systUpH,systDnH=[],[]
    for syst in finalSystList:

        #build the up/down envelope of the variations
        systUpH.append( data.Clone(syst+'_up') )
        systUpH[-1].SetTitle(syst)
        ci=len(systUpH)
        systUpH[-1].SetMarkerStyle(1)
        systUpH[-1].SetLineColor(COLORS[ci%len(COLORS)])
        systUpH[-1].SetFillColor(COLORS[ci%len(COLORS)])
        systUpH[-1].SetFillStyle(FILLS[int(ci/len(COLORS))])
        systUpH[-1].Reset('ICE')
        systDnH.append( systUpH[-1].Clone(syst+'_dn') )

        #sub contributions at reco level
        if isReco:
            expSystsH=fIn
            for ybin in xrange(2,expSystsH.GetNbinsY()+1):
                varName=expSystsH.GetYaxis().GetBinLabel(ybin)

                for systKey in finalSystList:
                    for isyst in xrange(0,len(finalSystList[systKey])):
                        isystName=finalSystList[systKey][isyst]
                        if varName!=isystName: continue
                        hv=expSystsH.ProjectionX(finalSystList[systKey][isyst],ybin,ybin)
                        normalizePerSlice(hv,obsAxis,sliceAxis)
                        hv.Add(data,-1)

                        #update envelope
                        for xbin in xrange(1,nbinsx+1):
                            curDeltaUp=systUpH[-1].GetBinContent(xbin)
                            curDeltaDn=systDnH[-1].GetBinContent(xbin)
                            delta=hv.GetBinContent(xbin)
                            if delta>0: systUpH[-1].SetBinContent(xbin,max(delta,curDeltaUp))
                            else:       systDnH[-1].SetBinContent(xbin,min(delta,curDeltaDn))

                        break

        #sub contributions at particle level
        else:
            for v in finalSystList[syst]:

                #compute absolute variations
                hv=fIn.Get('corrected_data%s'%v)
                normalizePerSlice(hv,obsAxis,sliceAxis)
                hv.Add(data,-1)

                #update envelope
                for xbin in xrange(1,nbinsx+1):
                    curDeltaUp=systUpH[-1].GetBinContent(xbin)
                    curDeltaDn=systDnH[-1].GetBinContent(xbin)
                    delta=hv.GetBinContent(xbin)
                    if syst=='FSR' : delta *= 1./ROOT.TMath.Sqrt(2.)

                    if delta>0: systUpH[-1].SetBinContent(xbin,max(delta,curDeltaUp))
                    else:       systDnH[-1].SetBinContent(xbin,min(delta,curDeltaDn))

        #finalize:
        # if single-sided uncertainty, mirror it
        # add uncertainty in quadrature to data
        # determine the median variation per bin
        isSingleSide=True if len(finalSystList[syst])==1 else False
        allVarVals=[]
        for xbin in xrange(1,nbinsx+1):

            cts=data.GetBinContent(xbin)
            deltaUp=systUpH[-1].GetBinContent(xbin)
            deltaDn=systDnH[-1].GetBinContent(xbin)
            if isSingleSide:
                if deltaUp==0 : 
                    deltaUp=abs(deltaDn)
                    systUpH[-1].SetBinContent(xbin,deltaUp)
                if deltaDn==0 : 
                    deltaDn=-abs(deltaUp)
                    systDnH[-1].SetBinContent(xbin,deltaDn)

            #save to determine the median
            if cts>0:
                allVarVals.append( 100*abs(deltaUp)/cts )
                allVarVals.append( 100*abs(deltaDn)/cts )
            else:
                allVarVals.append( 0 )
                allVarVals.append( 0 )

            #increment unertainty in data (use envelope)
            envVar = max(abs(deltaUp),abs(deltaDn))
            data.SetBinError(xbin,
                             ROOT.TMath.Sqrt(data.GetBinError(xbin)**2+envVar**2)
                             )
            
            #accumulate in up/down directions
            cumulUncUp,cumulUncDn=0,0
            if cts>0:
                cumulUncUp,cumulUncDn=(100*deltaUp/cts)**2,(100*deltaDn/cts)**2
            if len(systUpH)>1:
                cumulUncUp += systUpH[-2].GetBinContent(xbin)**2
                cumulUncDn += systDnH[-2].GetBinContent(xbin)**2
            systUpH[-1].SetBinContent(xbin,ROOT.TMath.Sqrt(cumulUncUp))
            systDnH[-1].SetBinContent(xbin,-ROOT.TMath.Sqrt(cumulUncDn))

        #median
        systSummary[syst]=np.median( np.array(allVarVals) )

    #save plot and write file with summary
    showSystsSummary(systUpH,systDnH,systSummary,outdir,outname)
    

"""
"""
def buildPlot(data,signalVars,obsAxis,sliceAxis,opt,outname,isReco=False):

    obs=obsAxis.GetName().split('_')[0]
    sliceVar=sliceAxis.GetName().split('_')[0] if sliceAxis else ''
    nslices=sliceAxis.GetNbins() if sliceAxis else 1
    frame      = ROOT.TH1F('frame',      'frame',     1,obsAxis.GetXmin(),obsAxis.GetXmax())
    frameratio = ROOT.TH1F('frameratio', 'frameratio',1,obsAxis.GetXmin(),obsAxis.GetXmax())

    #prepare signal ratios to display
    nratios=0
    signalVarRatioGr=[]
    for i in xrange(0,len(signalVars)):
        varTitle,signalVarHistoColls=signalVars[i] 
        nratios += 1 + len(signalVarHistoColls)
        ratioColl=[ROOT.TGraphAsymmErrors()]
        ratioColl[-1].SetFillStyle(3001)
        ratioColl[-1].SetFillColor(COLORS[i])
        ratioColl[-1].SetMarkerColor(COLORS[i])
        ratioColl[-1].SetLineColor(COLORS[i])
        ratioColl[-1].SetTitle(varTitle)
        ratioColl[-1].SetName('enveloperatio_%d'%i)
        nsubVars=len(signalVarHistoColls)
        for j in xrange(0,nsubVars):
            ratioColl.append(ROOT.TGraphErrors())
            ratioColl[-1].SetFillStyle(0)
            ratioColl[-1].SetFillColor(0)
            ci = COLORS[i] #1 if nsubVars >1  else COLORS[i]
            ratioColl[-1].SetMarkerColor(ci)
            ratioColl[-1].SetMarkerStyle(MARKERS[j])
            ratioColl[-1].SetLineColor(ci)
            ratioColl[-1].SetLineWidth(2)
            ratioColl[-1].SetTitle(signalVarHistoColls[j][0])
            ratioColl[-1].SetName('nominalratio_%d_%d'%(i,j))
        signalVarRatioGr.append(ratioColl)
    print nratios,'to display'

    for islice in xrange(1,nslices+1):

        idataGr=ROOT.TGraphErrors();
        idataGr.SetMarkerStyle(20)
        idataGr.SetMarkerColor(1)
        idataGr.SetLineColor(1)
        idataGr.SetName('data_%d'%islice)
        title='inc'
        if sliceAxis:
            title='[%d,%d]'%(int(sliceAxis.GetBinLowEdge(islice)),int(sliceAxis.GetBinUpEdge(islice)))
        idataGr.SetTitle(title)

        dataRelUncGr=ROOT.TGraphErrors()
        dataRelUncGr=idataGr.Clone('datarelunc_%d'%islice)
        dataRelUncGr.SetName('datarelunc_%d'%islice)
        dataRelUncGr.SetTitle('Uncertainty')
        dataRelUncGr.SetFillStyle(3004)
        dataRelUncGr.SetFillColor(ROOT.kYellow+3)

        isignalGr=ROOT.TGraphErrors();
        isignalGr.SetFillColor(3005)
        ci=ROOT.kAzure+4
        isignalGr.SetFillColor(ci)
        isignalGr.SetMarkerColor(ci)
        isignalGr.SetLineColor(ci)
        isignalGr.SetMarkerStyle(1)

        signal=None
        for varTitle,signalVarHistoColls in signalVars:
            if varTitle!=MAINMC[0] : continue
            signal=signalVarHistoColls[0][1][0]


        #build final distribution
        nobsBins=obsAxis.GetNbins()
        for xbin in xrange(1,nobsBins+1):
            hcen,hwid=obsAxis.GetBinCenter(xbin),obsAxis.GetBinWidth(xbin)
            
            rawBin=xbin
            if sliceAxis: rawBin+=nobsBins*(islice-1)

            np=idataGr.GetN()

            #signal
            signalCts,signalUnc=signal.GetBinContent(rawBin),signal.GetBinError(rawBin)
            isignalGr.SetPoint     (np,   hcen,      signalCts/(hwid))
            isignalGr.SetPointError(np,   0.5*hwid,  signalUnc/(hwid))

            #data
            dataCts,dataUnc=data.GetBinContent(rawBin),data.GetBinError(rawBin)
            idataGr.SetPoint       (np,   hcen, dataCts/(hwid))
            idataGr.SetPointError  (np,   0,    dataUnc/(hwid))
            if dataCts<=0: continue

            #relative uncertainty from data
            dataRelUnc=dataUnc/dataCts
            dataRelUncGr.SetPoint     (np,   hcen,      1)
            dataRelUncGr.SetPointError(np,   0.5*hwid,  dataRelUnc)

            #signal variations are displayed as ratio only
            deltax=hwid/nratios
            xoffset=hcen-0.5*hwid
            for ivar in xrange(0,len(signalVars)):

                xoffset += deltax

                varTitle,signalVarHistoColls=signalVars[ivar]
                nsubratios=len(signalVarHistoColls)

                centralRatio,centralRatioUnc,uncUp,uncDn=1,0,0,0
                for isubvar in xrange(0,len(signalVarHistoColls)):

                    #build the envelope from all the sub-variatios
                    ienvUp,ienvDn=-999,999
                    contribList=signalVarHistoColls[isubvar][1]
                    for icontrib in xrange(0,len(contribList)):
                        cts,unc=contribList[icontrib].GetBinContent(rawBin),contribList[icontrib].GetBinError(rawBin)
                        icontribRatio=cts/dataCts
                        ienvUp=max(icontribRatio,ienvUp)
                        ienvDn=min(icontribRatio,ienvDn)
                        if isubvar==0 and icontrib==0: 
                            centralRatio=icontribRatio
                            centralRatioUnc=unc/dataCts

                    if icontrib==0 and isubvar==0:
                        signalVarRatioGr[ivar][isubvar+1].SetPoint     (np, xoffset, centralRatio)
                        signalVarRatioGr[ivar][isubvar+1].SetPointError(np, 0,       centralRatioUnc)
                    else:
                        signalVarRatioGr[ivar][isubvar+1].SetPoint     (np, xoffset, 0.5*(ienvUp+ienvDn))
                        signalVarRatioGr[ivar][isubvar+1].SetPointError(np, 0,       0.5*(ienvUp-ienvDn))
                    xoffset += deltax

                    #add to build total uncertainty of the model
                    if ienvUp!=-999 : uncUp += (ienvUp-centralRatio)**2
                    if ienvDn!=999  : uncDn += (ienvDn-centralRatio)**2


                #envelope
                uncUp=ROOT.TMath.Sqrt(uncUp)
                uncDn=ROOT.TMath.Sqrt(uncDn)
                xmax,xmin=xoffset,xoffset-deltax*(len(signalVarHistoColls)+1)
                xcen,xwid=0.5*(xmax+xmin),0.5*(xmax-xmin)
                signalVarRatioGr[ivar][0].SetPoint(np, xcen,centralRatio)
                signalVarRatioGr[ivar][0].SetPointError(np,xwid*0.9,xwid*0.9,uncDn,uncUp) 
            print hcen-0.5*hwid,hcen+0.5*hwid,xoffset,deltax,hwid/deltax

        #show plot
        c=ROOT.TCanvas('c','c',1000,600)
        c.SetTopMargin(0.0)
        c.SetRightMargin(0.0)
        c.SetLeftMargin(0.0)
        c.SetBottomMargin(0.0)

        #main MC
        p1=ROOT.TPad('p1','p1',0.0,0.6,1.0,1.0) 
        p1.SetRightMargin(0.02)
        p1.SetLeftMargin(0.07)
        p1.SetTopMargin(0.1)
        p1.SetBottomMargin(0.01)
        p1.Draw()
        p1.cd()
        p1.SetLogy()
        frame.Draw()
        frame.GetYaxis().SetRangeUser(OBSRANGES[obs][0],OBSRANGES[obs][1])
        frame.GetYaxis().SetTitle('PDF')
        frame.GetYaxis().SetTitleOffset(0.4)
        frame.GetXaxis().SetLabelSize(0)
        frame.GetYaxis().SetTitleSize(0.08)
        frame.GetYaxis().SetLabelSize(0.08)
        isignalGr.Draw('2')        
        idataGr.Draw('ep')
        tex=ROOT.TLatex()
        tex.SetTextFont(42)
        tex.SetTextSize(0.08)
        tex.SetNDC()
        tex.DrawLatex(0.1,0.8,'#bf{CMS} #it{preliminary}')
        if nslices>1: tex.DrawLatex(0.5,0.8,'%s #in %s'%(VARS[sliceVar][0],idataGr.GetTitle()))
        tex.DrawLatex(0.85,0.93,'#scale[0.8]{35.9 fb^{-1} (#sqrt{s}=13 TeV)}')
        leg=ROOT.TLegend(0.7,0.9,0.9,0.65)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.07)
        leg.AddEntry(idataGr,'Data','ep')
        leg.AddEntry(isignalGr,MAINMC[0],'f')
        leg.Draw()
        
        c.cd()
        p2 = ROOT.TPad('p2','p2',0.0,0.0,1.0,0.6)
        p2.Draw()
        p2.SetBottomMargin(0.18)
        p2.SetRightMargin(0.02)
        p2.SetLeftMargin(0.07)
        p2.SetTopMargin(0.01)
        p2.Draw()
        p2.cd()
        frameratio.Draw()
        frameratio.GetYaxis().SetTitleSize(0.06)
        frameratio.GetYaxis().SetLabelSize(0.06)
        frameratio.GetXaxis().SetTitleSize(0.06)
        frameratio.GetXaxis().SetLabelSize(0.06)
        frameratio.GetYaxis().SetTitle('MC/Data')
        frameratio.GetYaxis().SetTitleOffset(0.5)
        frameratio.GetXaxis().SetTitle(VARS[obs][0])
        frameratio.GetYaxis().SetRangeUser(RATIORANGES[obs][0],RATIORANGES[obs][1])
        frameratio.GetYaxis().SetNdivisions(5+100*5)

        for ivar in xrange(0,len(signalVarRatioGr)):
            signalVarRatioGr[ivar][0].Draw('2')

        dataRelUncGr.Draw('2')


        ratiosLeg=[]
        ratiosLeg.append( ROOT.TLegend(0.1,0.95,0.25,0.915) )
        ratiosLeg[-1].SetBorderSize(0)
        ratiosLeg[-1].SetFillStyle(0)
        ratiosLeg[-1].SetTextFont(42)
        ratiosLeg[-1].SetTextSize(0.035)
        ratiosLeg[-1].AddEntry(dataRelUncGr,'Data unc.','f')
        ratiosLeg[-1].Draw()

        for ivar in xrange(0,len(signalVarRatioGr)):
            nsubvars=len(signalVarRatioGr[ivar])
            ratiosLeg.append( ROOT.TLegend(0.2+ivar*0.25,0.95,0.45+ivar*0.25,0.95-nsubvars*0.035) )
            ratiosLeg[-1].SetBorderSize(0)
            ratiosLeg[-1].SetFillStyle(0)
            ratiosLeg[-1].SetTextFont(42)
            ratiosLeg[-1].SetTextSize(0.035)
            ratiosLeg[-1].AddEntry(signalVarRatioGr[ivar][0],signalVarRatioGr[ivar][0].GetTitle(),'f')
            for isubvar in xrange(1,nsubvars):
                signalVarRatioGr[ivar][isubvar].Draw('pz')
                if nsubvars>1:
                    ratiosLeg[-1].AddEntry(signalVarRatioGr[ivar][isubvar],signalVarRatioGr[ivar][isubvar].GetTitle(),'ep')
            ratiosLeg[-1].Draw()

        l=ROOT.TLine(frameratio.GetXaxis().GetXmin(),1,frameratio.GetXaxis().GetXmax(),1)
        l.SetLineColor(ROOT.kBlue)
        l.Draw()
                
        c.cd()
        c.Modified()
        c.Update()
        for ext in ['png','pdf']:
            c.SaveAs('%s/%s_%s_%s.%s'%(opt.output,outname,islice,'smeared' if isReco else 'particle',ext))



def normalizePerSlice(h,obsAxis,sliceAxis):

    nobsBins=obsAxis.GetNbins()
    nslices=sliceAxis.GetNbins() if sliceAxis else 1
    norm=[0]*nslices
    for islice in xrange(1,nslices+1):
        for xbin in xrange(1,nobsBins+1):
            rawBin=xbin
            if sliceAxis: rawBin+=nobsBins*(islice-1)
            norm[islice-1]+=h.GetBinContent(rawBin)
    for islice in xrange(1,nslices+1):
        for xbin in xrange(1,nobsBins+1):
            rawBin=xbin
            if sliceAxis: rawBin+=nobsBins*(islice-1)
            h.SetBinContent( rawBin, h.GetBinContent(rawBin)/norm[islice-1] )
            h.SetBinError( rawBin, h.GetBinError(rawBin)/norm[islice-1] )

"""
"""
def readParticlePlotsFrom(args,opt):

    outdir=args[0].replace('.root','')

    #analysis axes
    analysisaxis=None
    with open(opt.analysisAxis,'r') as cachefile:
        analysisaxis = pickle.load(cachefile)
        
    #list of systematics
    systList={
        'Pileup'         : ['puup','pudn'],
        'Trigger/Sel.'   : ['effup','effdn'],
        'p_{T}(top)  '   : ['toppt'],
        '#mu_{R}/#mu_{F}': ['murup','murdn','mufup','mufdn','qup','qdn'],
        'b-tag'          : ['btagup','btagdn'],
        'JES'            : ['jesup','jesdn'],
        'JER'            : ['jerup','jerdn'],
        'LES'            : ['eesup','eesdn','mesup','mesdn'],
        'Trk. eff.'      : ['tkeff','tkeffeta'],
        'UE'             : ['t#bar{t} UEdn',     't#bar{t} UEup'],
        'FSR'            : ['t#bar{t} fsr dn',   't#bar{t} fsr up'],
        'ME-PS'          : ['t#bar{t} hdamp dn', 't#bar{t} hdamp up'],
        'ISR'            : ['t#bar{t} isr dn',   't#bar{t} isr up'],
        'm_{t}'          : ['t#bar{t} m=169.5',  't#bar{t} m=175.5'],
        #'PY8-HW++'       : ['t#bar{t} Herwig++'],
        'Background'     : ['bckpfakes']
        }
    finalSystList={}
    for syst in systList:
        newVars=[]
        for v in systList[syst]:
            newv=v
            for i in xrange(0,len(SYSTS)):
                if SYSTS[i][0]!=v : continue
                newv='_%d'%i
            newVars.append(newv)
        finalSystList[syst]=newVars

    #build plots
    for obs in OBSERVABLES:

        obsAxis=analysisaxis[(obs,False)]
        
        for s in SLICES:

            outname=obs
            sliceAxis=None if s is None else analysisaxis[(s,False)]
            if sliceAxis: outname += '%s'%s
            if not sliceAxis is None: continue

            #read unfolded data
            nomKey='%s_%s_0_inc'%(obs,s)
            fIn=ROOT.TFile('%s/%s.root'%(args[0],nomKey))
            uePlots['Data']=UEPlot(nomKey,obs,obsAxis)
            uePlots['Data'].addVariation('Data',None, fIn.Get('corrected_data') )
            for syst in finalSystList:
                for v in finalSystList[syst]:
                    uePlots['Data'].addVariation(syst,v,fIn.Get('corrected_data%s'%v))
            fIn.Close()

            #read sets to compare            
            fGen=ROOT.TFile.Open(args[1])
            fSyst=ROOT.TFile.Open(args[2])
            signalVars=[]
            for varTitle,subVars in COMPARISONSETS:
                for x,xvars in subVars:
                    if x!='nominal' : continue

                    key='%s_None_False'%nomKey
                    key='%s/%s_%s'%(key,key,xvars[0])
                    h=None
                    for f in [fGen,fSyst]:
                        try:
                            h=f.Get(key).Clone(ixvar)
                            break
                        except:                       
                            pass
                    if h is None: continue
                    uePlots[compSet]=UEPlot(nomKey+'mc%d'%icomp,compSet,obsAxis)
                    uePlots[compSet].addVariation(x,None,h)

            fGen.Close()
            fSyst.Close()

"""
"""
def readRecoPlotsFrom(args,opt):

    outdir=opt.output

    fIn=ROOT.TFile.Open(args[0])
    fSyst=ROOT.TFile.Open(args[1])

    #analysis axes
    analysisaxis=None
    with open(opt.analysisAxis,'r') as cachefile:
        analysisaxis = pickle.load(cachefile)

    #list of systematics
    systList={
        'Pileup'         : ['puup','pudn'],
        'Trigger/Sel.'   : ['effup','effdn'],
        'b-tag'          : ['btagup','btagdn'],
        'JES'            : ['jesup','jesdn'],
        'JER'            : ['jerup','jerdn'],
        'LES'            : ['eesup','eesdn','mesup','mesdn'],
        'Trk. eff.'      : ['tkeff','tkeffeta']
        }

    #build plots
    for obs in OBSERVABLES:

        obsAxis=analysisaxis[(obs,True)]
        
        for s in SLICES:

            outname=obs
            sliceAxis=None if s is None else analysisaxis[(s,True)]
            if sliceAxis: outname += '%s'%s
            #if not sliceAxis is None: continue
            
            key='%s_%s_inc_None_True'%(obs,s)
        
            #read data, signal and total background
            data,signal,bkg=None,None,None
            t=fIn.Get(key)
            for pkey in t.GetListOfKeys():
                h=t.Get(pkey.GetName())
                if not h.InheritsFrom('TH1') : continue
                if 'Data' in h.GetTitle():
                    data=h.Clone('data')
                elif h.GetTitle() in opt.signal:
                    signal=h.Clone('%s_nominal'%h.GetName())
                else:
                    if bkg is None: bkg=h.Clone('bkg')
                    else : bkg.Add(h)

            data.Add(bkg,-1)
            normalizePerSlice(data,obsAxis,sliceAxis)    

            #experimental systematics
            expSystsKey='%s_%s_inc_syst_True'%(obs,s)
            expSystsH=fIn.Get('{0}/{0}_{1}'.format(expSystsKey,opt.signal))
            addSystematics(signal,expSystsH,systList,obsAxis,sliceAxis,outdir,outname,True)

            #read sets to compare            
            signalVars=[]
            for varTitle,subVars in COMPARISONSETS:

                subVarHistColl=[]
                for x,xvars in subVars:
                    
                    histoColl=[]
                    for ixvar in xvars:
                        if varTitle==MAINMC[0] and ixvar==MAINMC[1]:  
                            histoColl.append(signal.Clone(ixvar))
                        else :
                            histoColl.append( fSyst.Get('{0}/{0}_{1}'.format(key,ixvar)).Clone(ixvar) )
                            normalizePerSlice(histoColl[-1],obsAxis,sliceAxis)
                        histoColl[-1].SetDirectory(0)
                    subVarHistColl.append( (x,histoColl) )
                signalVars.append( (varTitle,subVarHistColl) )

            buildPlot(data,signalVars,obsAxis,sliceAxis,opt,outname,True)

"""
"""
def main():

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(True)

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('--signal',
                      dest='signal',
                      help='signal [%default]',
                      type='string',
                      default='t#bar{t}')
    parser.add_option('--out',
                      dest='output',
                      help='output [%default]',
                      type='string',
                      default='./')
    parser.add_option('--reco',
                      dest='reco',
                      help='reco flag [%default]',
                      default=False,
                      action='store_true')
    parser.add_option('--cfg',
                      dest='analysisAxis',
                      help='cfg with axis definitions [%default]',
                      type='string',
                      default='%s/src/TopLJets2015/TopAnalysis/UEanalysis/analysisaxiscfg.pck'%os.environ['CMSSW_BASE'])
    (opt, args) = parser.parse_args()

    if opt.reco:
        readRecoPlotsFrom(args,opt)
    else:
        readParticlePlotsFrom(args,opt)


"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
