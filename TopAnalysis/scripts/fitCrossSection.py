#!/usr/bin/env python

import optparse
import os,sys
import commands
from array import array
import ROOT

POItitles={'r':'#mu=#sigma/#sigma_{th}',
           'BtagEff':'kx#sigma_{varepsilon_{b}}',
           'Mtop':'m_{t} [GeV]'}

"""
common CMS label
"""
def drawCMSlabel(startY=0.97):        
    cmsLabel=ROOT.TLatex()
    cmsLabel.SetTextFont(42)
    cmsLabel.SetTextSize(0.035)
    cmsLabel.SetNDC()
    cmsLabel.DrawLatex(0.12,startY,'#bf{CMS} #it{preliminary}')
    cmsLabel.DrawLatex(0.75,startY,'#scale[0.8]{2.2 fb^{-1} (13 TeV)}')
    cmsLabel.Draw()

"""
Prepares the fitting script
"""
def prepareFitScript(datacard,POIs,unblind=False):
    baseDir=os.path.dirname(datacard)
    fitScriptUrl='%s/runCombine.sh'%baseDir
    fitScript=open(fitScriptUrl,'w')

    fitScript.write('cd %s\n'%baseDir)
    
    fitScript.write('\n# convert datacard to workspace\n')
    fitScript.write('echo \"Converting datacard to workspace\"\n')
    fitScript.write('text2workspace.py %s -m 0 -o workspace.root\n' % os.path.basename(datacard))

    fitScript.write('\n# likelihood scans\n')
    for parameter in POIs:

        minParVal=0.8 if parameter=='r' else -2.0
        maxParVal=1.2 if parameter=='r' else +2.0
        rangeOpt='--setPhysicsModelParameterRanges %s=%f,%f'%(parameter,minParVal,maxParVal)

        poiOpt='' if parameter=='r' else '--redefineSignalPOIs %s'%parameter

        if parameter=='r':
            fitScript.write('\n## max likelihood fit\n')
            fitScript.write('echo \"Running MaxLikelihoodFit for r\"\n')
            fitScript.write('combine workspace.root -M MaxLikelihoodFit -t -1 --expectSignal=1 -m 0 --minimizerTolerance 0.001\n')
            fitScript.write('mv mlfit.root mlfit_exp.root\n')
            
            fitScript.write('\n## impacts\n')
            fitScript.write('echo \"To compute expected impacts re-run runCombine.sh uncommenting the appropriate lines\n')
            fitScript.write('#combineTool.py -M Impacts -d workspace.root -m 0  -t -1 --expectSignal=1 --doInitialFit\n')
            fitScript.write('#combineTool.py -M Impacts -d workspace.root -m 0  -t -1 --expectSignal=1 --doFits\n')
            fitScript.write('#combineTool.py -M Impacts -d workspace.root -m 0 -o impacts.json\n')
            
            fitScript.write('\n## toys\n')
            fitScript.write('echo \n"To run toys  re-run runCombine.sh uncommenting the appropriate lines\n')
            fitScript.write('#rscan=(0.9 1.0 1.1)\n')
            fitScript.write('#for r in ${rscan[@]}; do\n')
            fitScript.write('#\t combine workspace.root -M MaxLikelihoodFit -t 100 --expectSignal=${r} -m ${r} --toysFrequentist --noErrors --minos none --robustFit=1;\n')
            fitScript.write('#done\n\n')

            if unblind:
                fitScript.write('combine workspace.root -M MaxLikelihoodFit -m 0 --setPhysicsModelParameters r=1 --minimizerTolerance 0.001\n')
                fitScript.write('mv mlfit.root mlfit_obs.root\n')
                            
        fitScript.write('\n## function of %s\n'%parameter)
        fitScript.write('echo \"Running likelihood scan for %s\"\n'%parameter)
        fitScript.write('combine workspace.root -M MultiDimFit -P %s -t -1 --expectSignal=1 --algo=grid --points=100 %s %s -m 0\n'%(parameter,rangeOpt,poiOpt))
        fitScript.write('mv higgsCombineTest.MultiDimFit.mH0.root exp_plr_scan_%s.root\n'%parameter)
        fitScript.write('combine workspace.root -M MultiDimFit -P %s -t -1 --expectSignal=1 --algo=grid --points=100 %s %s -m 0 -S 0\n'%(parameter,rangeOpt,poiOpt))
        fitScript.write('mv higgsCombineTest.MultiDimFit.mH0.root exp_plr_scan_stat_%s.root\n'%parameter)
        if unblind:
            fitScript.write('combine workspace.root -M MultiDimFit -P %s --algo=grid --points=100 %s %s -m 0 --saveWorkspace --setPhysicsModelParameters r=1 --minimizerTolerance 0.001\n'%(parameter,rangeOpt,poiOpt))
            fitScript.write('mv higgsCombineTest.MultiDimFit.mH0.root obs_plr_scan_%s.root\n'%parameter)
            fitScript.write('combine obs_plr_scan_%s.root -M MultiDimFit -P %s --algo=grid --points=100 %s %s -m 0 --freezeNuisances all --snapshotName "MultiDimFit"\n'%(parameter,parameter,rangeOpt,poiOpt))
            fitScript.write('mv higgsCombineTest.MultiDimFit.mH0.root obs_plr_scan_stat_%s.root\n'%parameter)

        fitScript.write('echo \"To run pole mass fit uncomment the appropriate lines\n')
        fitScript.write('#text2workspace.py datacard.dat -P TopLJets2015.TopAnalysis.TopPoleMass:topPoleMass -m 172.5 --PO verbose --PO param=${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/python/xsecvsmpole.dat -o workspace_mpole.root\n')
        fitScript.write('#combine workspace_mpole.root -M MultiDimFit -P mpole -t -1 --algo=grid --points=100 --setPhysicsModelParameters mpole=172.5  -m 0 --setPhysicsModelParameterRanges mpole=165,185\n')
        fitScript.write('#mv higgsCombineTest.MultiDimFit.mH0.root exp_plr_scan_mpole.root\n')
        fitScript.write('#combine workspace_mpole.root -M MultiDimFit -P mpole --algo=grid --points=100 --setPhysicsModelParameters mpole=172.5  -m 0 --setPhysicsModelParameterRanges mpole=165,185\n')
        fitScript.write('#mv higgsCombineTest.MultiDimFit.mH0.root obs_plr_scan_mpole.root\n')

    fitScript.write('\n# 2D likelihood scans\n')
    for i in xrange(0,len(POIs)):
        for j in xrange(i+1,len(POIs)):
            minPiVal=0.8 if POIs[i]=='r' else -2.0
            maxPiVal=1.2 if POIs[i]=='r' else +2.0
            minPjVal=0.8 if POIs[j]=='r' else -2.0
            maxPjVal=1.2 if POIs[j]=='r' else +2.0
            rangeOpt='--setPhysicsModelParameterRanges %s=%f,%f:%s=%f,%f'%(POIs[i],minPiVal,maxPiVal,POIs[j],minPjVal,maxPjVal)
            
            fitScript.write('## function of %s,%s\n'%(POIs[i],POIs[j]))
            fitScript.write('echo \"Running 2D likelihood scan for %s vs %s\"\n'%(POIs[i],POIs[j]))
            fitScript.write('combine workspace.root -M MultiDimFit --redefineSignalPOIs %s,%s -P %s -P %s  -t -1 --expectSignal=1 --algo=grid --points=1000 %s -m 0\n'%(POIs[i],POIs[j],POIs[i],POIs[j],rangeOpt)) 
            fitScript.write('mv higgsCombineTest.MultiDimFit.mH0.root  exp_plr_scan_%svs%s.root\n'%(POIs[i],POIs[j]))
            if unblind:
                fitScript.write('combine workspace.root -M MultiDimFit --redefineSignalPOIs %s,%s -P %s -P %s --algo=grid --points=1000  %s -m 0\n'%(POIs[i],POIs[j],POIs[i],POIs[j],rangeOpt)) 
                fitScript.write('mv higgsCombineTest.MultiDimFit.mH0.root obs_plr_scan_%svs%s.root\n'%(POIs[i],POIs[j]))

    fitScript.write('\ncd -\n')
    
    fitScript.close()
    return fitScriptUrl

"""
Opens the output ROOT files from the fits and plots the results for comparison
"""
def show1DLikelihoodScan(resultsSet,parameter='r',output='./'):
   
    #likelihood scans
    nllGrs={}
    colors=[1, ROOT.kOrange,  ROOT.kRed+1, ROOT.kMagenta-9, ROOT.kBlue-7]
    ires=0
    for title,datacard in resultsSet:
        ires+=1
        dir=os.path.dirname(datacard)
        files=[('#splitline{expected}{#scale[0.8]{(stat+syst)} }', 'exp_plr_scan',      1, 1),
               ('#splitline{expected}{#scale[0.8]{(stat)}}',       'exp_plr_scan_stat', 3, 1),
               ('#splitline{observed}{#scale[0.8]{(stat+syst)}}',  'obs_plr_scan',      1, 3),
               ('#splitline{observed}{#scale[0.8]{stat only}}',    'obs_plr_scan_stat', 3, 3)]

        for ftitle,f,lstyle,lwidth in files:

            #if file not available continue
            fIn=ROOT.TFile.Open('%s/%s_%s.root' % (dir,f,parameter) )
            if not fIn : continue
        
            #create new graph for the likelihood scan
            if not ftitle in nllGrs: nllGrs[ftitle]=[]
            nllGrs[ftitle].append(ROOT.TGraph())
            nllGrs[ftitle][-1].SetTitle(title)
            nllGrs[ftitle][-1].SetLineStyle(lstyle)
            nllGrs[ftitle][-1].SetMarkerStyle(1)
            nllGrs[ftitle][-1].SetFillStyle(0)
            nllGrs[ftitle][-1].SetLineWidth(lwidth)
            nllGrs[ftitle][-1].SetLineColor(colors[ires-1])
            nllGrs[ftitle][-1].SetMarkerColor(colors[ires-1])

            #fill graph
            tree=fIn.Get('limit')
            for n in xrange(0,tree.GetEntriesFast()):
                tree.GetEntry(n)
                nll=2*tree.deltaNLL
                if nll>20 : continue
                npoint=nllGrs[ftitle][-1].GetN()
                parVal=getattr(tree,parameter)
                if parameter=='Mtop':
                    parVal=parVal*3+172.5
                nllGrs[ftitle][-1].SetPoint(npoint,parVal,nll)
            nllGrs[ftitle][-1].Sort()
            fIn.Close()

    #show 1D likelihood scan
    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.1)
    c.SetRightMargin(0.05)
    minParVal=0.8 if parameter=='r' else -2.0
    maxParVal=1.2 if parameter=='r' else +2.0
    if parameter=='Mtop' : minParVal,maxParVal=166.5,178.5
    frame=ROOT.TH1F('frame',';%s;-2#DeltalnL'%POItitles[parameter],100,minParVal,maxParVal)
    frame.Draw()
    frame.GetYaxis().SetRangeUser(0,10)
    allLegs=[]
    ileg=0
    for ftitle in nllGrs:
        allLegs.append( ROOT.TLegend(0.15+ileg*0.15,0.92,0.3+ileg*0.15,0.85-0.04*len(nllGrs[ftitle]) ) )
        allLegs[-1].SetTextFont(42)
        allLegs[-1].SetTextSize(0.035)
        allLegs[-1].SetBorderSize(0)
        allLegs[-1].SetFillStyle(1001)
        allLegs[-1].SetFillColor(0)
        allLegs[-1].SetHeader(ftitle)
        for gr in nllGrs[ftitle]:
            gr.Draw('c')
            allLegs[-1].AddEntry(gr,gr.GetTitle(),'l')
        ileg+=1
    for leg in allLegs: leg.Draw()

    cl=ROOT.TLine()
    cl.SetLineStyle(2)
    cl.SetLineColor(ROOT.kGray+2)
    txt=ROOT.TLatex()
    txt.SetTextFont(42)
    txt.SetTextSize(0.025)
    txt.SetTextColor(ROOT.kGray+3)
    for delta,title in [(1.0,'68% CL'),(3.84,'95% CL')]:
        txt.DrawLatex(maxParVal-10*frame.GetBinWidth(1),delta+0.25,title)
        cl.DrawLine(minParVal,delta,maxParVal,delta)

    drawCMSlabel()
   
    c.Modified()
    c.Update()
    for ext in ['png','pdf','C']:
        c.SaveAs('%s/nll1dscan_%s.%s'%(output,parameter,ext))
    
"""
2D likelihood scan
"""
def show2DLikelihoodScan(resultsSet,parameters):

    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.1)
    c.SetRightMargin(0.05)
    
    contours = array('d',[1.0,3.84])


    #likelihood scans
    nllGrs={}
    colors=[1, ROOT.kOrange-1,  ROOT.kRed+1, ROOT.kMagenta-9, ROOT.kBlue-7]
    ires=0
    for title,datacard in resultsSet:
        ires+=1
        dir=os.path.dirname(datacard)
        files=[('#splitline{expected}{#scale[0.8]{(stat+syst)} }', 'exp_plr_scan',      1, 1),
               ('#splitline{expected}{#scale[0.8]{(stat)}}',       'exp_plr_scan_stat', 3, 1),
               ('#splitline{observed}{#scale[0.8]{(stat+syst)}}',  'obs_plr_scan',      1, 3),
               ('#splitline{observed}{#scale[0.8]{stat only}}',    'obs_plr_scan_stat', 3, 3)]

        for ftitle,f,lstyle,lwidth in files:

            if not ftitle in nllGrs: nllGrs[ftitle]=[]

            #if file not available continue
            fIn=ROOT.TFile.Open('%s/%s_%svs%s.root' % (dir,f,parameters[0],parameters[1]) )
            if not fIn : continue
            tree=fIn.Get('limit')

            for ll,ul,tag in [(1-0.99,1.0,'99cl')] : #(1-0.68,1.0,'68cl'),(1-0.95,1.0,'95cl'),(1-0.99,1.0,'99cl')]:
                c.Clear()    
                nllGrs[ftitle].append( ROOT.ll2dContourPlot(tree,parameters[0],parameters[1],ll,ul) )
                nllGrs[ftitle][-1].SetTitle(title)
                nllGrs[ftitle][-1].SetLineStyle(lstyle)
                nllGrs[ftitle][-1].SetMarkerStyle(1)
                nllGrs[ftitle][-1].SetFillStyle(1001)
                nllGrs[ftitle][-1].SetFillColor(colors[ires-1])
                nllGrs[ftitle][-1].SetLineWidth(lwidth)
                nllGrs[ftitle][-1].SetLineColor(colors[ires-1])
                nllGrs[ftitle][-1].SetLineWidth(2)
                nllGrs[ftitle][-1].SetMarkerColor(colors[ires-1])
                            
    #show 2D likelihood scan
    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.1)
    c.SetRightMargin(0.05)
    frame=ROOT.TH2F('frame','frame',10,0.8,1.2,10,-3,3)
    frame.Draw()
    allLegs=[]
    ileg=0
    for ftitle in nllGrs:
        if len(nllGrs[ftitle])==0 : continue
        allLegs.append( ROOT.TLegend(0.15+ileg*0.15,0.92,0.3+ileg*0.15,0.85-0.04*len(nllGrs[ftitle]) ) )
        allLegs[-1].SetTextFont(42)
        allLegs[-1].SetTextSize(0.035)
        allLegs[-1].SetBorderSize(0)
        allLegs[-1].SetFillStyle(1001)
        allLegs[-1].SetFillColor(0)
        allLegs[-1].SetHeader(ftitle)
        for gr in nllGrs[ftitle]:
            gr.Draw('f')
            allLegs[-1].AddEntry(gr,gr.GetTitle(),'f')
        ileg+=1
    for leg in allLegs: leg.Draw()

    drawCMSlabel()
   
    c.Modified()
    c.Update()

"""
compare prefit and postfit nuisances
"""
def compareNuisances(resultsSet,output):
   
    colors=[{'exp':ROOT.kGray,'obs':1}, 
            {'exp':ROOT.kOrange,'obs':ROOT.kOrange-1}]
    postFitNuisGr={}
    nuisCorrelationH={}
    nuisanceList=[]
    dy=1./(len(resultsSet)+2.)
    resCtr=0
    for title,datacard in resultsSet:
        resCtr+=1
        dir=os.path.dirname(datacard)        
        for fit in ['exp','obs']:

            key=(title,fit)
            
            #open file if it exists
            fname='%s/mlfit_%s.root'%(dir,fit)
            inF=ROOT.TFile.Open(fname)
            if inF is None or inF.IsZombie() : continue

            #get (S+B) fit results and store in graph
            fit_s=inF.Get('fit_s')            
            postFitNuisGr[key]=ROOT.TGraphErrors()
            postFitNuisGr[key].SetName('postfitgr_%s'%''.join(key))
            postFitNuisGr[key].SetTitle(title)
            marker=20 if fit=='obs' else 24
            postFitNuisGr[key].SetMarkerStyle(marker)
            postFitNuisGr[key].SetMarkerColor(colors[resCtr-1][fit])
            postFitNuisGr[key].SetLineColor(colors[resCtr-1][fit])
            postFitNuisGr[key].SetLineWidth(2)
            postFitNuisGr[key].SetFillStyle(0)
            npars=fit_s.floatParsFinal().getSize()
            nuisCorrelationH[key]=ROOT.TH1F('nuiscorrelationgr_%s'%''.join(key),';Nuisance;Correlation with #mu=#sigma/#sigma_{th}',npars,0,npars)
            nuisCorrelationH[key].SetLineColor(colors[resCtr-1][fit])
            nuisCorrelationH[key].SetFillColor(colors[resCtr-1][fit])
            fillStyle=1001 if fit=='obs' else 3002
            nuisCorrelationH[key].SetFillStyle(fillStyle)
            nuisCorrelationH[key].SetDirectory(0)
            
            #iterate over nuisances
            for ipar in range(0,npars):
                var=fit_s.floatParsFinal().at(ipar)                
                pname=var.GetName()
                if pname=='r' : continue
                np=postFitNuisGr[key].GetN()
                postFitNuisGr[key].SetPoint(np,var.getVal(),ipar+0.2+resCtr*dy)
                postFitNuisGr[key].SetPointError(np,var.getError(),0)

                nuislabel='#color[%d]{%s}'%((ipar%2)*10+1,pname)
                nuisCorrelationH[key].GetXaxis().SetBinLabel(ipar+1,nuislabel)
                nuisCorrelationH[key].SetBinContent(ipar+1,fit_s.correlation(pname,'r'))
                nuisCorrelationH[key].SetBinError(ipar+1,0)

                #first time around save also the label to display
                if resCtr==1 and fit=='exp':
                    nuisanceList.append( nuislabel )

            #all done with the file
            inF.Close()


    #show correlations
    c=ROOT.TCanvas('c','c',1500,500)
    c.SetLeftMargin(0.1)
    c.SetTopMargin(0.3)
    c.SetRightMargin(0.05)
    c.SetBottomMargin(0.05)
    c.SetGridx(True)
    ires=0
    for title,_ in resultsSet:
        ires+=1
        c.Clear()
        leg=ROOT.TLegend(0.12,0.65,0.3,0.6)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.035)
        leg.AddEntry(nuisCorrelationH[(title,'exp')],nuisCorrelationH[(title,'exp')].GetTitle()+'(exp)','f')
        leg.AddEntry(nuisCorrelationH[(title,'obs')],nuisCorrelationH[(title,'obs')].GetTitle()+'(obs)','f')
        leg.SetNColumns(2)
        nuisCorrelationH[(title,'obs')].Draw('histX+')
        nuisCorrelationH[(title,'obs')].GetXaxis().SetLabelOffset(+0.15)
        nuisCorrelationH[(title,'obs')].GetYaxis().SetRangeUser(-1,1)
        nuisCorrelationH[(title,'exp')].Draw('histX+same')        
        leg.Draw()
        drawCMSlabel(startY=0.65)
        c.RedrawAxis()
        c.Modified()
        c.Update()
        for ext in ['png','pdf','C']:
            c.SaveAs('%s/correlations_%d.%s'%(output,ires,ext))

    #show nuisances
    c=ROOT.TCanvas('c','c',500,1500)
    c.SetTopMargin(0.1)
    c.SetLeftMargin(0.3)
    c.SetBottomMargin(0.1)
    c.SetRightMargin(0.05)
    c.SetGridy(True)
    npars=len(nuisanceList)
    frame=ROOT.TH2F('frame',';N x #sigma_{pre-fit}',1,-3,3,npars,0,npars)
    frame.SetDirectory(0)
    for ipar in xrange(0,npars): frame.GetYaxis().SetBinLabel(ipar+1,nuisanceList[ipar])
    frame.GetXaxis().SetRangeUser(-3,3)
    frame.GetYaxis().SetLabelSize(0.025)
    frame.Draw()

    gr1s=ROOT.TGraph()
    gr1s.SetName('gr1s')
    gr1s.SetMarkerStyle(1)
    gr1s.SetMarkerColor(19)
    gr1s.SetLineColor(19) 
    gr1s.SetFillStyle(1001)
    gr1s.SetFillColor(19) 
    gr1s.SetPoint(0,-1,0)
    gr1s.SetPoint(1,-1,npars)
    gr1s.SetPoint(2,1,npars)
    gr1s.SetPoint(3,1,0)
    gr1s.SetPoint(4,-1,0)
    gr2s=gr1s.Clone('gr2s')
    gr2s.SetMarkerColor(18) 
    gr2s.SetLineColor(18)
    gr2s.SetFillStyle(1001)
    gr2s.SetFillColor(18) 
    gr2s.SetPoint(0,-2,0)
    gr2s.SetPoint(1,-2,npars)
    gr2s.SetPoint(2,2,npars)
    gr2s.SetPoint(3,2,0)
    gr2s.SetPoint(4,-2,0)
    gr2s.Draw('f')
    gr1s.Draw('f')
    
    leg=ROOT.TLegend(0.1,0.93,0.6,0.96)
    leg.SetNColumns(resCtr)
    leg.SetTextFont(42)
    leg.SetTextSize(0.03)
    leg.SetBorderSize(0)
    leg.SetFillStyle(1001)
    leg.SetFillColor(0)
    leg2=ROOT.TLegend(0.1,0.90,0.6,0.93)
    leg2.SetNColumns(resCtr)
    leg2.SetTextFont(42)
    leg2.SetTextSize(0.03)
    leg2.SetBorderSize(0)
    leg2.SetFillStyle(1001)
    leg2.SetFillColor(0)
    for key in postFitNuisGr:
        postFitNuisGr[key].Draw('p')
        title=postFitNuisGr[key].GetTitle()
        if key[1]=='exp':
            leg.AddEntry(postFitNuisGr[key],'%s (exp)'%title,'p')
        else:
            leg2.AddEntry(postFitNuisGr[key],'%s (obs)'%title,'p')
    leg.Draw()
    leg2.Draw()
    
    txt=ROOT.TLatex()
    txt.SetTextFont(42)
    txt.SetTextSize(0.025)
    txt.SetTextColor(ROOT.kGray+3)
    for delta,title in [(1.0,'-1#sigma'),(2,'+2#sigma'),(-1,'-1#sigma'),(-2,'-2#sigma')]:
        txt.DrawLatex(delta-0.2,frame.GetYaxis().GetXmax()+0.5,title)      

    drawCMSlabel()
    c.RedrawAxis()
    c.Modified()
    c.Update()
    for ext in ['png','pdf','C']:
        c.SaveAs('%s/nuisances.%s'%(output,ext))

        
"""
main
"""
def main():
  
    #configuration
    usage = 'usage: %prog category1=datacard1 category2=datacard2 ....'
    parser = optparse.OptionParser(usage)
    parser.add_option('-o', '--output',       dest='output',       help='output directory',       default='./',           type='string')
    parser.add_option(      '--noFit',        dest='noFit',        help='don\'t run the fits',    action='store_true')
    parser.add_option(      '--POIs',         dest='POIs',         help='parameters of interest', default='r',            type='string')
    parser.add_option(      '--unblind',      dest='unblind',      help='unblind',                action='store_true')
    (opt, args) = parser.parse_args()

    print os.path.dirname(os.path.realpath(sys.argv[0]))
    sys.path.insert(0, r'%s/../'% os.path.dirname(os.path.realpath(sys.argv[0])) )

    POIs=opt.POIs.split(',')

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(True) #False)
    
    resultsSet=[]
    for newCat in args:
        cat,datacard=newCat.split('=')
        resultsSet.append( (cat,datacard) )
        if not opt.noFit:
            fitScriptUrl=prepareFitScript(datacard=datacard,POIs=POIs,unblind=opt.unblind)
            print 'Fit script for %s available at %s'%(cat,fitScriptUrl)
            os.system('sh %s'%fitScriptUrl)
            
    compareNuisances(resultsSet=resultsSet,output=opt.output)
    
    for parameter in POIs:
        show1DLikelihoodScan(resultsSet=resultsSet,parameter=parameter,output=opt.output)


    #ROOT.gROOT.LoadMacro('src/RootTools.cc+')

    #for i in xrange(0,len(POIs)):
    #    for j in xrange(i+1,len(POIs)):
    #        show2DLikelihoodScan(resultsSet,parameters=[POIs[i],POIs[j]])
    
            

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
