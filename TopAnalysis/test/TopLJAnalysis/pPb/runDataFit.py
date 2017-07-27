#!/usr/bin/env python

lumi=(174.5,8.725)

from TopLJets2015.TopAnalysis.rounding import *
import pickle
import json
import os
import ROOT
import optparse
import os,sys
from roofitTools import *
from parameterizeMCShapes import ALLPDFS
from prepareWorkspace import EVENTCATEGORIES as SELEVENTCATEGORIES

EVENTCATEGORIES=[x for x in SELEVENTCATEGORIES if not '1f' in x]

acceptance={#'e':(0.0508,ROOT.TMath.Sqrt(0.0014**2+(0.031*0.0508)**2)), #eta restricted
            'e' :(0.056,ROOT.TMath.Sqrt(0.0014**2+(0.031*0.056)**2)),
            'mu':(0.060,ROOT.TMath.Sqrt(0.0016**2+(0.031*0.060)**2))}
efficiency={#'e':(0.767,0.0304), #LOOSE ISO
            'e':(0.623,0.0249), #MEDIUM ID+ISO
            #'e':(0.536,0.0206), #TIGHT ID+ISO)
            'mu':(0.915,0.0366)}
ebExp=(0.595,0.0595)
jsf=(1.0,0.034)
WMODEL=None
QCDNORM=None
observables=[('mjj','M(jj)'),('mthad','M(t_{had})'),('mtlep','M(t_{lep})')]

"""
disable RooFit verbosity
"""
def shushRooFit():
    ROOT.RooMsgService.instance().setSilentMode(True);
    ROOT.RooMsgService.instance().getStream(0).removeTopic(ROOT.RooFit.Minimization)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.Minimization)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.ObjectHandling)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.DataHandling)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.Fitting)
    ROOT.RooMsgService.instance().getStream(0).removeTopic(ROOT.RooFit.Plotting)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.Plotting)
    ROOT.RooMsgService.instance().getStream(0).removeTopic(ROOT.RooFit.InputArguments)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.InputArguments)
    ROOT.RooMsgService.instance().getStream(0).removeTopic(ROOT.RooFit.Eval)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.Eval)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.Integration)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.NumIntegration)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.NumIntegration)

"""
"""
def addToFitResults(title,fitResults,varSet):
    fitResults[title]={}
    iter = varSet.createIterator()
    iparam = iter.Next()
    while iparam :
        if not iparam.getAttribute('Constant'):
            fitResults[title][iparam.GetName()]=(iparam.getVal(),iparam.getErrorLo(),iparam.getErrorHi())
        iparam = iter.Next()

"""
"""
def printFitResults(fitResults,impacts,opt):

    #sort out all the parameters
    allVars=set()
    for key in fitResults:
        for pname in fitResults[key]:
            allVars.add(pname)

    #save a neat tex table
    with open('%s/fitResults_wmodel%d_%d.tex'%(opt.output,opt.wModel,opt.fitType),'w') as fOut:
        fOut.write('\\hline\n')
        fOut.write('\multirow{3}{*}{Variable} & \\multicolumn{3}{c}{Fit type} \\\\\n')
        for key in fitResults:
            fOut.write('& %s'%key)
        fOut.write('\\\\\n')
        fOut.write('\\hline\n')
        for pname in allVars:
            fOut.write('%25s'%pname)
            for key in fitResults:
                if not pname in fitResults[key]:
                    fOut.write('& %25s'%'')
                else :
                    val,unclo,unchi=fitResults[key][pname]
                    unc=0.5*(abs(unclo)+abs(unchi))
                    fOut.write('& %25s'%toLatexRounded(val,unc)) #'$%3.2f_{%3.2f}^{%3.2f}$'%fitResults[key][pname])
            fOut.write('\\\\\n')
        fOut.write('\\hline\n')

    #save impacts as a json file
    with open('%s/impacts_wmodel%d_%d.json'%(opt.output,opt.wModel,opt.fitType),'w') as fOut:
        json.dump(impacts,fOut)


"""
"""
def defineCombinedPDFs(w):
    for t in ['2D','3D']:

        simPDF={}
        for ch in ['combined','e','mu']:
            simPDF[ch]=ROOT.RooSimultaneous('model_%s_%s'%(ch,t),'model_%s_%s'%(ch,t),w.cat('sample'))

        for cat in EVENTCATEGORIES:

            ch='e' if cat[0]=='e' else 'mu'

            #multiply PDFs
            sigProd='S_mjj_{0},S_mthad_{0}'.format(cat)
            qcdProd='QCD_mjj_{0},QCD_mthad_{0}'.format(ch)
            wProd  ='W_mjj_{0},W_mthad_{0}'.format(cat)
            if t=='3D':
                sigProd += ',S_mtlep_{0}'.format(cat)
                qcdProd += ',QCD_mtlep_{0}'.format(ch)
                wProd   += ',W_mtlep_{0}'.format(cat)
            w.factory('PROD::S_{1}_{0}({2})'.format(cat,t,sigProd))
            w.factory('PROD::QCD_{0}_{1}({2})'.format(t,ch,qcdProd))
            w.factory('PROD::W_{1}_{0}({2})'.format(cat,t,wProd))

            #sum up so that it can be extended
            w.factory('SUM::model_{1}_{0}(Nsig_{0}*S_{1}_{0},Nqcd_{0}*QCD_{1}_{2},Nw_{0}*W_{1}_{0})'.format(cat,t,ch))

            simPDF['combined'].addPdf(w.pdf('model_{1}_{0}'.format(cat,t)),cat)
            simPDF[ch].addPdf(w.pdf('model_{1}_{0}'.format(cat,t)),cat)

        #import simultaneous pdf to workspace
        for ch in ['combined','e','mu']:
            print 'Importing',simPDF[ch].GetName()
            getattr(w,'import')(simPDF[ch],ROOT.RooFit.RecycleConflictNodes())

"""
"""
def definePDF(w,varName,wModel=0):

    # PDF definition
    simPDF={'combined':ROOT.RooSimultaneous('model_combined_%s'%varName,'model_combined_%s'%varName,w.cat('sample'))}

    #YIELDS
    for ch in ['e','mu']:
        simPDF[ch]=ROOT.RooSimultaneous('model_%s_%s'%(ch,varName),'model_%s_%s'%(ch,varName),w.cat('sample'))
        w.factory("RooFormulaVar::Nsig_{0}('@0*@1*@2*@3',{{acc_{0},eff_{0},lumi,xsec}})".format(ch))
        w.factory("RooFormulaVar::Nsig_{0}1l4j2b('@0*pow(@1,2)',{{Nsig_{0},eb}})".format(ch))
        w.factory("RooFormulaVar::Nsig_{0}1l4j1b1q('@0*2*@1*(1-@1)',{{Nsig_{0},eb}})".format(ch))
        w.factory("RooFormulaVar::Nsig_{0}1l4j2q('@0*pow(1-@1,2)',{{Nsig_{0},eb}})".format(ch))
    w.factory("RooFormulaVar::Nsig_total('@0+@1',{{Nsig_e,Nsig_mu}})")

    #backgrounds yields are indepedent in each category (to be profiled)
    #the number of QCD (W) like is parameterized as fqcd.Nbkg ((1-fqcd).Nw)
    for cat in EVENTCATEGORIES:

        w.factory('Nbkg_{0}[1000,0,100000]'.format(cat))
        if '2b'   in cat : w.factory('fqcd_{0}[0.0,1.0]'.format(cat))
        elif '1b' in cat : w.factory('fqcd_{0}[0.0,1.0]'.format(cat))
        else             : w.factory('fqcd_{0}[0.0,1.0]'.format(cat))
        #if '2b'   in cat : w.factory('fqcd_{0}[1.0]'.format(cat))
        #elif '1b' in cat : w.factory('fqcd_{0}[1.0]'.format(cat))
        #else             : w.factory('fqcd_{0}[1.0]'.format(cat))

        w.factory("RooFormulaVar::Nqcd_{0}('@0*@1',{{Nbkg_{0},fqcd_{0}}})".format(cat))

        #W background shares common parameters for shape in electron and muon samples
        w.factory("RooFormulaVar::Nw_{0}('@0*(1-@1)',{{Nbkg_{0},fqcd_{0}}})".format(cat))
        ch='e' if cat[0]=='e' else 'mu'
        baseCat=cat[1:] if cat[0]=='e' else cat[2:]

        #init QCD normalization constraints
        if QCDNORM:
            val,valUnc=QCDNORM[(ch,baseCat)]
            if val==0 and valUnc==0:
                w.var('fqcd_{0}'.format(cat)).setVal(0)
                w.var('fqcd_{0}'.format(cat)).setConstant(True)
            else:
                w.factory("RooFormulaVar::qcdnorm_{0}_constraint('0.5*pow((@0-@1)/@2,2)',{{Nqcd_{0},NqcdCen_{0}[{1}],NqcdUnc_{0}[{2}]}})".format(cat,val,valUnc))
                w.factory("RooGaussian::qcdnorm_{0}_gconstraint(NqcdCen_{0},Nqcd_{0},NqcdUnc_{0})".format(cat))

        #
        #W+jets component
        #
        mpvName   = 'mpv_w{1}_{0}'.format(baseCat if wModel==0 else cat,varName)
        widthName = 'width_w{1}_{0}'.format(baseCat if wModel==0 else cat,varName)

        #in-situ model
        if wModel<2:            
            minMPV=20  if varName=='mjj' else 150
            maxMPV=120 if varName=='mjj' else 200
            w.factory('{0}[{1},{2},{3}]'.format(mpvName,WMODEL[(varName,baseCat,'MPV')],minMPV,maxMPV))
            w.factory('{0}[{1},5,100]'.format(widthName,WMODEL[(varName,baseCat,'Sigma')]))

            for landauPar,wmodelName,relUnc in [(mpvName,'MPV',0.2), (widthName,'Sigma',0.05)]:
                val=WMODEL[(varName,baseCat,wmodelName)]
                valUnc=relUnc*val
                w.factory("RooFormulaVar::{0}_constraint('0.5*pow((@0-@1)/@2,2)',{{{0},cen_{0}[{1}],unc_{0}[{2}]}})".format(landauPar,val,valUnc))
                w.factory("RooGaussian::{0}_gconstraint({0},cen_{0},unc_{0})".format(landauPar))

        #use pre-defined model
        if wModel==2:
            w.factory('{0}[{1}]'.format(mpvName,WMODEL[(varName,baseCat,'MPV')]))
            w.factory('{0}[{1}]'.format(widthName,WMODEL[(varName,baseCat,'Sigma')]))
            print 'Fixed W model for %s in %s'%(varName,baseCat)
            print '\t mpv=',w.var(mpvName).getVal()
            print '\t width=',w.var(widthName).getVal()

        w.factory('RooLandau::W_{1}_{0}({1},{2},{3})'.format(cat,varName,mpvName,widthName))

        #w.factory('SUM:model_{0}(Nqcd_{0}*QCD_{1},Nw_{0}*W_{0})'.format(cat,varName))
        w.factory('SUM::model_{1}_{0}(Nsig_{0}*S_{1}_{0},Nqcd_{0}*QCD_{1}_{2},Nw_{0}*W_{1}_{0})'.format(cat,varName,ch))
        w.cat('sample').setLabel(cat)
        simPDF['combined'].addPdf(w.pdf('model_{1}_{0}'.format(cat,varName)),cat)
        simPDF[ch].addPdf(w.pdf('model_{1}_{0}'.format(cat,varName)),cat)


    #import simultaneous pdf to workspace
    for ch in simPDF:
        print 'Importing',simPDF[ch].GetName()
        getattr(w,'import')(simPDF[ch],ROOT.RooFit.RecycleConflictNodes())

"""
model QCD using RooKeysPDF and the failing id/iso/btag sideband data
"""
def addQCDModel(w,opt,varName):

    #get sideband data
    for ch in ['e','mu']:
        redData = w.data('data').reduce(ROOT.RooFit.Cut("sample==sample::%s1f4j2q"%ch))
        fullsideBand=redData
        if redData.numEntries()>5000:
            print 'reducing numentries to 5000'
            redData=redData.reduce(ROOT.RooFit.EventRange(0,5000))

        #run RooKeysPdf and add it to the workspace
        keyspdf=ROOT.RooKeysPdf("QCD_%s_%s"%(varName,ch),
                                "QCD_%s_%s"%(varName,ch),
                                w.var(varName),
                                redData,
                                ROOT.RooKeysPdf.NoMirror,
                                1.4)
        getattr(w,'import')(keyspdf,ROOT.RooFit.RecycleConflictNodes())



"""
readout signal models
"""
def addSignalModel(w,opt,varName):

    #get signal workspace from file
    fIn=ROOT.TFile.Open(opt.signal)
    ws=fIn.Get('w')
    fIn.Close()

    #define resolution
    #w.factory('bias_s[0]')
    #w.var('bias_s').SetTitle('#mu_{bias}')
    #w.factory('resol_s[0.1]')
    #w.var('resol_s').SetTitle('#sigma_{resol}')
    #w.factory("Gaussian::S_resolution(%s,bias_s,resol_s)"%opt.varName)

    data=w.data('data')
    for cat in EVENTCATEGORIES:

        if opt.fitType<4:
            for c in ['cor','wro']:
                pdfDef=[x.format(cat) for x in  ALLPDFS[('S_%s'%c,varName)] ]
                for p in pdfDef:
                    if varName=='mjj' and opt.fitType==0: p=p.replace('(mjj','(scaledmjj')
                    pdf=w.factory(p)
                iter = pdf.getParameters(data).createIterator()
                iparam = iter.Next()
                while iparam :
                    try:
                        mcval=ws.var(iparam.GetName()).getVal()
                        iparam.setVal(mcval)
                        iparam.setConstant(True)
                    except:
                        pass
                    iparam = iter.Next()

            resFracName='f_s{1}_{0}'.format(cat,varName)
            w.factory('SUM::S_{1}_{0}( {2}[0,1]*S_cor{1}_{0}, S_wro{1}_{0} )'.format(cat,varName,resFracName))
            w.var(resFracName).setVal(1.0 if opt.fitType==3 else ws.var(resFracName).getVal() )
            w.var(resFracName).setConstant(True)
        else:
            muLimits=[60,90] if varName=='mjj' else [150,180]
            sigmaLimits=[10,30]
            w.factory('mu_scor%s[60,90]'%varName)
            w.factory('sigmaL_scor%s[10,25]'%varName)
            w.factory('sigmaR_scor%s[10,25]'%varName)
            w.factory('RooBifurGauss::S_{1}_{0}({1},mu_scor{1},sigmaL_scor{1},sigmaR_scor{1})'.format(cat,varName))

        #convolve with a common resolution
        #w.factory("FCONV::S_{0}({1}, S_resolution,S_base_{0})".format(cat,opt.varName))

        #frame=w.var(opt.varName).frame()
        #pdf.plotOn(frame)
        #w.pdf('S_{0}'.format(cat)).plotOn(frame)
        #frame.Draw()
        #raw_input()



"""
minimizes the likelihood and profiles the poi
"""
def runFit(pdf,data,poi,obs,w,outDir):

    mc=ROOT.RooStats.ModelConfig('mc_%s'%pdf.GetName(),w)
    mc.SetPdf(pdf)
    mc.SetParametersOfInterest(poi)
    mc.SetObservables(obs)
    getattr(w,'import')(mc)

    pl=ROOT.RooStats.ProfileLikelihoodCalculator(data,mc)
    pl.SetConfidenceLevel(0.683);
    interval=pl.GetInterval()

    ROOT.Math.MinimizerOptions.SetDefaultMinimizer("Minuit2")
    firstPOI = mc.GetParametersOfInterest().first()
    lowerLimit = interval.LowerLimit(firstPOI)
    upperLimit = interval.UpperLimit(firstPOI)
    cenVal     = interval.GetBestFitParameters().find(firstPOI.GetName()).getVal()
    unc=(upperLimit-lowerLimit)*0.5

    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.03)
    c.SetBottomMargin(0.1)
    plot = ROOT.RooStats.LikelihoodIntervalPlot(interval)
    plot.SetRange(lowerLimit-unc*2,upperLimit+unc*2)
    plot.SetNPoints(40)
    plot.Draw("tf1")
    label = ROOT.TLatex()
    label.SetNDC()
    label.SetTextFont(42)
    label.SetTextSize(0.04)
    label.DrawLatex(0.6,0.9,'#bf{CMS} #it{preliminary}')
    label.DrawLatex(0.6,0.85,'%s=%3.1f^{+%3.1f}_{-%3.1f}'%(firstPOI.GetTitle(),cenVal,upperLimit-cenVal,cenVal-lowerLimit))
    c.Modified()
    c.Update()
    for ext in ['png','pdf']:
        c.SaveAs('%s/ll_%s.%s'%(outDir,mc.GetName(),ext))

#simple roofit version
def runSimpleFit(pdf,data,poi,constr=None,pToFix=None):

    #create the log likelihood
    nll=pdf.createNLL(data,
                      ROOT.RooFit.Extended(True),
                      ROOT.RooFit.NumCPU(8))

    #add constraints
    if constr:
        parcels=ROOT.RooArgList(nll)
        iter = constr.createIterator()
        var = iter.Next()
        while var :
            print var.GetName()
            parcels.add(var)
            var = iter.Next()
        nll=ROOT.RooAddition('nllc','nllc',parcels)

    #fix parameters
    if pToFix:
        iter = pdf.getParameters(data).createIterator()
        iparam = iter.Next()
        while iparam :
            try:
                pname=iparam.GetName()
                if pname in pToFix:
                    iparam.setConstant(True)
            except:
                pass
            iparam = iter.Next()

    #maximize the likelihood
    minuit=ROOT.RooMinuit(nll)
    minuit.setStrategy(2)
    minuit.migrad() #minimize with respect to all parameters
    minuit.minos(poi)
    r=minuit.save()
    pll=nll.createProfile(poi)
    return (r,pll)


"""
"""
def performFits(opt):

    #read workspace
    fIn=ROOT.TFile(opt.finalWS)
    w=fIn.Get('w')
    fIn.Close()

    #data to fit
    data=w.data('data')

    #fit results summary
    fitResults={}
    impacts={}

    #parameter of interest
    poi=ROOT.RooArgSet()    
    poi.add(w.var('xsec'))
    paramList=[] #,('eb','#varepsilon_{b}')]
    paramList.append( ('xsec','#sigma(t#bar{t})') )

    #constraints
    constr   = ROOT.RooArgSet()
    allFuncs = w.allFunctions()
    funcIter = allFuncs.createIterator()
    ifunc    = funcIter.Next()
    while ifunc:
        fname=ifunc.GetName()
        if 'constraint' in fname:
            skip=False
            if not opt.fitType in [2,3] and 'mthad' in fname: skip=True
            if not opt.fitType in [3]   and 'mtlep' in fname: skip=True
            if not skip :
                constr.add( w.function(fname) )
                print fname
        ifunc=funcIter.Next()

    #run the 2D/3D fits
    if opt.fitType in [1,2]:

        for ch in ['e','mu','combined']:

            t='2D' if opt.fitType==1 else '3D'
            pdf=w.pdf('model_{0}_{1}'.format(ch,t))
            obsList=ROOT.RooArgSet()
            obsList.add(w.var('mjj'))
            obsList.add(w.var('mthad'))
            if t=='3D': obsList.add(w.var('mtlep'))
            #runFit(pdf,data,poi,obsList,w,opt.output)
            result,pll=runSimpleFit(pdf,data,poi,constr)
            w.saveSnapshot('fitresult_%s'%ch,pdf.getParameters(data))

            addToFitResults('%s_%s'%(ch,t),fitResults,result.floatParsFinal()) #pdf.getParameters(data))

            #pll.plotOn(xsecframe,ROOT.RooFit.Name(ch),ROOT.RooFit.ShiftToZero())

            EVENTCATEGORIES2SHOW=EVENTCATEGORIES if ch=='combined' else [ch+'1l4j'+x for x in ['2q','1b1q','2q']]
            for obs in observables:
                if obs=='mtlep' and t=='2D' : continue
                showFitResult(fitVar=obs[0],
                            data=data,
                            pdf=pdf,
                            categs=EVENTCATEGORIES2SHOW,
                            w=w,
                            showComponents=['S_cor*','S_cor*,S_wro*'],
                            rangeX=(25 if obs[0]=='mjj' else 150,
                                    300 if obs[0]=='mjj' else 400),
                            outDir=opt.output,
                            paramList=paramList,
                            pfix='_%s_%sfit'%(ch,t))

    else:
        w.saveSnapshot('default',w.pdf('model_combined_mjj').getParameters(data))

        for ch in ['combined','mu','e']:
            
            w.loadSnapshot('default')

            #pdf
            pdf=w.pdf('model_%s_mjj'%ch)

            #run the fit to the mjj variable
            #obsList=ROOT.RooArgSet()
            #obsList.add(w.var('mjj'))
            #runFit(pdf,data,poi,obsList,w,opt.output)
            result,pll=runSimpleFit(pdf,data,poi,constr)
            snapshotKey='fitresult_%s'%ch
            w.saveSnapshot(snapshotKey,pdf.getParameters(data))

            #add fit results to the summary, including event yields
            key='%s_1D'%ch
            addToFitResults(key,fitResults,result.floatParsFinal())
            for formName in ['Nsig_e','Nsig_mu','Nsig_total']:
                func=w.function(formName)
                val,unc=func.getVal(),func.getPropagatedError(result)
                fitResults[key][formName]=(val,unc,unc)

            #do the following only for the main fit and if required
            if not opt.impacts: continue

            #fit stat unc. (fix all except xsec and repeat the fit)
            pToFix=[]
            for pname in fitResults[key]:
                if pname!='xsec' and not 'Nsig' in pname: pToFix.append(pname)
            w.loadSnapshot(snapshotKey)
            result,_=runSimpleFit(pdf,data,poi,constr,pToFix)
            fitResults[key]['xsec_statonly']=(
                result.floatParsFinal().find('xsec').getVal(),
                result.floatParsFinal().find('xsec').getErrorLo(),
                result.floatParsFinal().find('xsec').getErrorHi()
                )

            #let all float again
            for pname in pToFix:
                w.var(pname).setConstant(False)
            
            #now do the impacts
            impacts[key]={}
            for pname in ['lumi','eb','eff_e','eff_mu','jsf','acc']:

                impacts[key][pname]=[]

                #get postfit
                val,unclo,unchi=fitResults[key][pname]

                for ivar in xrange(0,2):

                    #set at +/-1 sigma postfit and repeat fit
                    w.loadSnapshot(snapshotKey)
                    w.var(pname).setVal(val+unclo if ivar==0 else val+unchi)
                    w.var(pname).setConstant(True)

                    result,_=runSimpleFit(pdf,data,poi,constr,[pname])

                    #save difference in xsec
                    dR=result.floatParsFinal().find('xsec').getVal()-fitResults['%s_1D'%ch]['xsec'][0]
                    impacts[key][pname].append(dR)

                    #let it float again
                    w.var(pname).setConstant(False)                

    printFitResults(fitResults,impacts,opt)
    origWorkspaceFile=os.path.splitext(os.path.basename(opt.finalWS))[0]
    w.writeToFile('fit_%s_%d.root'%(origWorkspaceFile,opt.fitType))


    #show the likelihoods
    #xsecframe=w.var('xsec').frame(ROOT.RooFit.Bins(10),ROOT.RooFit.Range(0,100))
    #c=ROOT.TCanvas('c','c',500,500)
    #c.SetTopMargin(0.05)
    #c.SetLeftMargin(0.12)
    #c.SetRightMargin(0.03)
    #c.SetBottomMargin(0.1)
    #xsecframe.Draw()
    #xsecframe.GetXaxis().SetTitle('#sigma [nb]')
    #xsecframe.GetYaxis().SetTitle('-2 log #lambda')
    #label = ROOT.TLatex()
    #label.SetNDC()
    #label.SetTextFont(42)
    #label.SetTextSize(0.04)
    #label.DrawLatex(0.6,0.9,'#bf{CMS} #it{preliminary}')
    #c.Modified()
    #c.Update()
    #for ext in ['png','pdf']:
    #    c.SaveAs('%s/pll_%s_%dfit.%s'%(opt.output,ch,opt.fitType))
    #raw_input()

"""
"""
def addPDFToWorkspace(opt):
    #read workspace
    fIn=ROOT.TFile(opt.input)
    w=fIn.Get('w')
    fIn.Close()

    #common to all variables and/or channels
    w.factory('xsec[60,0,300]')    

    w.factory("RooFormulaVar::lumi_constraint('0.5*pow((@0-@1)/@2,2)',{lumiCen[%f],lumi[0,500],lumiUnc[%f]})"%lumi)
    w.factory("RooGaussian::lumi_gconstraint(lumiCen,lumi,lumiUnc)")

    w.factory("RooFormulaVar::accconstraint('0.5*pow((@0-@1)/@2,2)',{accCen[0.0],acc[-5,5],accUnc[1.0]})")
    w.factory("RooGaussian::acc_gconstraint(accCen,acc,accUnc)")

    for ch in ['e','mu']:
        w.factory("RooFormulaVar::acc_%s('%f+%f*@0',{acc})"%(ch,acceptance[ch][0],acceptance[ch][1]))
        w.factory("RooFormulaVar::eff_%s_constraint('0.5*pow((@0-@1)/@2,2)',{effCen_%s[%f],eff_%s[0.,1.0],effUnc_%s[%f]})"%(ch,ch,efficiency[ch][0],ch,ch,efficiency[ch][1]))
        w.factory("RooGaussian::eff_{0}_gconstraint(effCen_{0},eff_{0},effUnc_{0})".format(ch))

    w.factory("RooFormulaVar::ebconstraint('0.5*pow((@0-@1)/@2,2)',{ebCen[%f],eb[0.60,0.0,1.0],ebUnc[%f]})"%ebExp)
    w.factory("RooGaussian::eb_gconstraint(ebCen,eb,ebUnc)")

    w.factory("RooFormulaVar::jsfconstraint('0.5*pow((@0-@1)/@2,2)',{jsfCen[%f],jsf[0.5,1.5],jsfUnc[%f]})"%jsf)
    w.factory("RooGaussian::jsf_gconstraint(jsfCen,jsf,jsfUnc)")

    w.factory("RooFormulaVar::scaledmjj('TMath::Max(0.,@0*@1)',{jsf,mjj})")


    #instantiate PDFs
    for vname,vtit in observables:
        w.var(vname).SetTitle(vtit)

        #model QCD
        addQCDModel(w,opt,vname)

        #add signal model
        addSignalModel(w,opt,vname)

        #define the PDFs
        definePDF(w,vname,opt.wModel)

    #constraints
    constr=ROOT.RooArgList()
    constr.add( w.pdf('eb_gconstraint') )
    constr.add( w.pdf('jsf_gconstraint') )
    constr.add( w.pdf('acc_gconstraint') )
    constr.add( w.pdf('eff_mu_gconstraint') )
    constr.add( w.pdf('eff_e_gconstraint') )
    constr.add( w.pdf('lumi_gconstraint') )
    for cat in EVENTCATEGORIES:
        constr.add( w.pdf('qcdnorm_{0}_gconstraint'.format(cat)) )
    constrProd=ROOT.RooProdPdf('constrProd','constrProd',constr)
    getattr(w,'import')(constrProd,ROOT.RooFit.RecycleConflictNodes())

    defineCombinedPDFs(w)

    w.writeToFile('finalworkspace_wmodel%d.root'%opt.wModel,True)

"""
"""
def main():

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(True) 

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-o', '--output',    dest='output',    default='plots/Data8TeV_pp',          type='string',   help='output directory [%default]')
    parser.add_option('-i', '--input',     dest='input',     default='workspace_Data8TeV_pp.root', type='string',   help='workspace [%default]')
    parser.add_option('-s', '--signal',    dest='signal',    default='pdf_workspace_MC8.16TeV_TTbar_pPb.root', type='string',   help='signal workspace [%default]')
    parser.add_option(      '--wModel',    dest='wModel',    default=0, type=int,   help='W model (0-in-situ; 1-decorrelated in-situ; 2-extrapol. MC) [%default]')
    parser.add_option(      '--fitType',   dest='fitType',   default=0,                  type=int,
        help='0-full signal; 1-full signal 2D; 2-full signal 3D; 3-res from MC; 4-res from CB [%default]')
    parser.add_option(      '--impacts',   dest='impacts',   default=False,                        action='store_true',        help='Run impacts [%default]')
    parser.add_option('-v', '--verbose',   dest='verbose',   default=0,                            type=int,        help='Verbose mode [%default]')
    parser.add_option(      '--finalWorkspace',      dest='finalWS',      default=None,            type='string',   help='final workspace to be used for the fit [%default]')
    (opt, args) = parser.parse_args()

    #keep roofit quite
    if opt.verbose<9 : shushRooFit()

    #load a W model
    global WMODEL
    with open('wmodel.pck','r') as fIn:
        WMODEL=pickle.load(fIn)

    #load the QCD normalization estimation
    global QCDNORM
    with open('qcdnorm.pck','r') as fIn:
        QCDNORM=pickle.load(fIn)

    #create final workspace if not given
    if opt.finalWS is None: addPDFToWorkspace(opt)
    else:                   performFits(opt)


if __name__ == "__main__":
    main()
