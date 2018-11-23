#!/usr/bin/env python

import optparse
import os,sys
import commands
from array import array
import ROOT

POItitles={'r':'#mu=#sigma/#sigma_{th}',
           'btagRate':'SF_{b}',
           'Mtop':'m_{t} [GeV]'}

"""
common CMS label
"""
def drawCMSlabel(startY=0.96,label=''):        
    cmsLabel=ROOT.TLatex()
    cmsLabel.SetTextFont(42)
    cmsLabel.SetTextSize(0.035)
    cmsLabel.SetNDC()
    cmsLabel.DrawLatex(0.12,startY,'#bf{CMS} #it{preliminary}')
    cmsLabel.DrawLatex(0.70,startY,'#scale[0.8]{%s}'%label)
    cmsLabel.Draw()

"""
Opens the output ROOT files from the fits and plots the results for comparison
"""
def show1DLikelihoodScan(resultsSet,parameter='r',output='./',label=''):
   
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

            #if file not available continue
            fIn=ROOT.TFile.Open('%s/%s_%s.root' % (dir,f,parameter.replace('Rate','')) )
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
                if parameter=='btagRate' :
                    parVal=parVal*0.1+1.0
                nllGrs[ftitle][-1].SetPoint(npoint,parVal,nll)
            nllGrs[ftitle][-1].Sort()
            fIn.Close()

    #show 1D likelihood scan
    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.1)
    c.SetRightMargin(0.05)
    minParVal=0.6 if parameter=='r' else -2.0
    maxParVal=1.4 if parameter=='r' else +2.0
    if parameter=='btagRate' : minParVal,maxParVal=0.8,1.2
    if parameter=='Mtop'     : minParVal,maxParVal=166.5,178.5
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

    drawCMSlabel(label=label)
    c.Modified()
    c.Update()
    for ext in ['png','pdf','C']:
        c.SaveAs('%s/nll1dscan_%s.%s'%(output,parameter,ext))
    
"""
2D likelihood scan
"""
def show2DLikelihoodScan(resultsSet,parameters,output,label):

    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.1)
    c.SetRightMargin(0.05)
    
    #see Table 33.2 @ http://pdg.lbl.gov/2011/reviews/rpp2011-rev-statistics.pdf
    contours = array('d',[2.30,4.61,5.99,6.18,9.21,11.83])
    #contours = array('d',[1.0,2.71,3.84,4.00,6.63,9.00])

    #likelihood scans
    nllGrs={}
    colors=[1, ROOT.kOrange, ROOT.kRed+1, ROOT.kMagenta-9, ROOT.kBlue-7]
    ires=0
    frame=None
    for title,datacard in resultsSet:
        ires+=1
        dir=os.path.dirname(datacard)
        files=[('#splitline{expected}{#scale[0.8]{(stat+syst)} }', 'exp_plr_scan',      3, 3),
               ('#splitline{observed}{#scale[0.8]{(stat+syst)}}',  'obs_plr_scan',      1, 3)]

        for ftitle,f,lstyle,lwidth in files:

            if not ftitle in nllGrs: nllGrs[ftitle]=[]

            #if file not available continue
            fname='%s/%s_%svs%s.root' % (dir,f,parameters[0],parameters[1].replace('Rate',''))
            fIn=ROOT.TFile.Open(fname)
            if not fIn : 
                print 'Failed to open',fname
                continue
            tree=fIn.Get('limit')


            #fill the 2D likelihood histogram
            xmin,xmax=0,2
            ymin,ymax=0,2
            if parameters[1]=='btagRate' : ymin,ymax=0.8,1.2
            hcont=ROOT.TH2D('hcont','%s;%s;%s'%(ftitle,POItitles[parameters[0]],POItitles[parameters[1]]),
                            50,xmin,xmax,50,ymin,ymax)
            hcont.SetContour(len(contours),contours)
            for ientry in xrange(0,tree.GetEntriesFast()) : 
                tree.GetEntry(ientry)
                x=getattr(tree,parameters[0])
                y=getattr(tree,parameters[1])
                if parameters[1]=='btagRate': y=y*0.1+1.0
                hcont.Fill(x,y,2*tree.deltaNLL)

            #check if there are bins with no entries and interpolate from neighboring bins
            for xbin in xrange(1,hcont.GetNbinsX()+1):
                for ybin in xrange(1,hcont.GetNbinsY()+1):
                    val,unc =hcont.GetBinContent(xbin,ybin), hcont.GetBinError(xbin,ybin)
                    if val!=0 or unc!=0 : continue 
                    avgVal,nUsed=0,0
                    for k in xrange(-1,2):
                        for l in xrange(-1,2):
                            neighborVal=hcont.GetBinContent(xbin+k,ybin+l)
                            if neighborVal==0 : continue
                            avgVal+=neighborVal
                            nUsed+=1
                    if nUsed>0:
                        print 'patching',xbin,ybin
                        hcont.SetBinContent(xbin,ybin,avgVal/nUsed)

            #draw the contours and save them
            c.Clear()
            hcont.Draw('contz list')
            c.Update()
            contGrs = ROOT.gROOT.GetListOfSpecials().FindObject("contours")
            for name,title,idx in [ (f+'68cl','68%CL',0) ] : # (f+'95cl','95%CL',2)]:                
                nllGrs[ftitle].append( contGrs.At(idx).At(0).Clone() )
                nllGrs[ftitle][-1].SetName(name)
                nllGrs[ftitle][-1].SetTitle(title)
                nllGrs[ftitle][-1].SetFillStyle(0)
                nllGrs[ftitle][-1].SetFillColor(0)
                lcolor=colors[len(nllGrs[ftitle])-1]
                nllGrs[ftitle][-1].SetLineColor(lcolor)
                nllGrs[ftitle][-1].SetLineWidth(lwidth)
                nllGrs[ftitle][-1].SetLineStyle(lstyle)
                nllGrs[ftitle][-1].SetMarkerColor(lcolor)
                nllGrs[ftitle][-1].SetMarkerStyle(1)

            if frame is None: 
                frame=hcont.Clone('frame')
                frame.Reset('ICE')
                frame.SetDirectory(0)
            #hcont.Delete()

    #show 2D likelihood scan
    c.Clear()
    frame.Draw()
    frame.GetYaxis().SetTitleSize(0.045)
    frame.GetYaxis().SetTitleOffset(1.1)
    frame.GetXaxis().SetTitleSize(0.045)
    allLegs=[]
    ileg=0
    for ftitle in nllGrs:
        if len(nllGrs[ftitle])==0 : continue
        allLegs.append( ROOT.TLegend(0.15+ileg*0.15,0.93,0.3+ileg*0.15,0.85-0.04*len(nllGrs[ftitle]) ) )
        allLegs[-1].SetTextFont(42)
        allLegs[-1].SetTextSize(0.035)
        allLegs[-1].SetBorderSize(0)
        allLegs[-1].SetFillStyle(1001)
        allLegs[-1].SetFillColor(0)
        allLegs[-1].SetHeader(ftitle)
        for gr in nllGrs[ftitle]:
            gr.Draw('l')
            allLegs[-1].AddEntry(gr,gr.GetTitle(),'l')
        ileg+=1
    for leg in allLegs: leg.Draw()

    drawCMSlabel(label=label)   
    c.Modified()
    c.Update()
    for ext in ['png','pdf','C']:
        print '%s/nll1dscan_%s.%s'%(output,parameter,ext)
        c.SaveAs('%s/nll2d_%svs%s.%s'%(output,parameters[0],parameters[1],ext))
   

"""
compare prefit and postfit nuisances
"""
"""
compare prefit and postfit nuisances
"""
def compareNuisances(resultsSet,output,label):
   
    colors=[{'exp':ROOT.kGray+1,'obs':1}, 
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
            try:
                if inF is None or inF.IsZombie() : continue
            except:
                continue

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
            fillStyle=1001 if fit=='obs' else 3001
            nuisCorrelationH[key].SetFillStyle(fillStyle)
            nuisCorrelationH[key].SetDirectory(0)
            
            #iterate over nuisances
            for ipar in range(0,npars):
                var=fit_s.floatParsFinal().at(ipar)                
                pname=var.GetName()
                if pname=='r' : continue
                np=postFitNuisGr[key].GetN()
                postFitNuisGr[key].SetPoint(np,ipar+0.2+resCtr*dy,var.getVal())
                postFitNuisGr[key].SetPointError(np,0,var.getError())

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
        if (title,'obs') in nuisCorrelationH:
            leg.AddEntry(nuisCorrelationH[(title,'obs')],nuisCorrelationH[(title,'obs')].GetTitle()+'(obs)','f')
            leg.SetNColumns(2)
            nuisCorrelationH[(title,'obs')].Draw('histX+')
            nuisCorrelationH[(title,'obs')].GetXaxis().SetLabelOffset(+0.15)
            nuisCorrelationH[(title,'obs')].GetYaxis().SetRangeUser(-1,1)
            nuisCorrelationH[(title,'exp')].Draw('histX+same')        
        else:
            nuisCorrelationH[(title,'exp')].Draw('histX+')
            nuisCorrelationH[(title,'exp')].GetXaxis().SetLabelOffset(+0.15)
            nuisCorrelationH[(title,'exp')].GetYaxis().SetRangeUser(-1,1)
        leg.Draw()
        drawCMSlabel(startY=0.65,label=label)
        c.RedrawAxis()
        c.Modified()
        c.Update()
        for ext in ['png','pdf','C']:
            c.SaveAs('%s/correlations_%d.%s'%(output,ires,ext))

    #show nuisances
    c.Clear()
    npars=len(nuisanceList)
    frame=ROOT.TH2F('frame',';Nuisance;N x #sigma_{pre-fit}',npars,0,npars,1,-3,3)
    frame.SetDirectory(0)
    for ipar in xrange(0,npars): frame.GetXaxis().SetBinLabel(ipar+1,nuisanceList[ipar])
    frame.GetYaxis().SetRangeUser(-3,3)
    frame.GetXaxis().SetLabelSize(0.025)
    frame.GetXaxis().SetLabelOffset(+0.1)
    frame.Draw('histX+')

    gr1s=ROOT.TGraph()
    gr1s.SetName('gr1s')
    gr1s.SetMarkerStyle(1)
    gr1s.SetMarkerColor(19)
    gr1s.SetLineColor(19) 
    gr1s.SetFillStyle(1001)
    gr1s.SetFillColor(19) 
    gr1s.SetPoint(0,0,-1)
    gr1s.SetPoint(1,npars,-1)
    gr1s.SetPoint(2,npars,1)
    gr1s.SetPoint(3,0,1)
    gr1s.SetPoint(4,-0,-1)
    gr2s=gr1s.Clone('gr2s')
    gr2s.SetMarkerColor(18) 
    gr2s.SetLineColor(18)
    gr2s.SetFillStyle(1001)
    gr2s.SetFillColor(18) 
    gr2s.SetPoint(0,0,-2)
    gr2s.SetPoint(1,npars,-2)
    gr2s.SetPoint(2,npars,2)
    gr2s.SetPoint(3,0,2)
    gr2s.SetPoint(4,0,-2)
    gr2s.Draw('f')
    gr1s.Draw('f')
    
    leg=ROOT.TLegend(0.24,0.69,0.32,0.61)
    leg.SetNColumns(resCtr)
    leg.SetTextFont(42)
    leg.SetTextSize(0.03)
    leg.SetBorderSize(0)
    leg.SetFillStyle(1001)
    leg.SetFillColor(0)
    leg2=ROOT.TLegend(0.32,0.69,0.38,0.61)
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
        txt.DrawLatex(frame.GetXaxis().GetXmax()+0.5,delta-0.2,title)      
    
    drawCMSlabel(startY=0.65,label=label)
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
    parser.add_option('-o', '--output',       dest='output',       help='output directory',       default='./',                    type='string')
    parser.add_option(      '--POIs',         dest='POIs',         help='parameters of interest', default='r',                     type='string')
    parser.add_option(      '--label',        dest='label',        help='plot labels',            default='2.3 fb^{-1} (13 TeV)',  type='string')
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
            
    #for parameter in POIs:
    #    show1DLikelihoodScan(resultsSet=resultsSet,parameter=parameter,label=opt.label,output=opt.output)

    for i in xrange(0,len(POIs)):
        for j in xrange(i+1,len(POIs)):
            show2DLikelihoodScan(resultsSet,parameters=[POIs[i],POIs[j]],output=opt.output,label=opt.label)

    compareNuisances(resultsSet=resultsSet,output=opt.output,label=opt.label)

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
