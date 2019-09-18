import ROOT
from TopLJets2015.TopAnalysis.Plot import fixExtremities, scaleTo, divideByBinWidth

def getPassGenDistributions(pName,gName,plotter,ci,fill=0,marker=20):

    """reads the distributions from the plotter"""

    passH=plotter.Get(pName)
    fixExtremities(passH)
    passH.SetDirectory(0)
    passH.SetFillStyle(fill)
    passH.SetLineColor(ci)
    passH.SetMarkerColor(ci)
    passH.SetLineWidth(2)
    passH.SetMarkerStyle(marker)

    genH=plotter.Get(gName)
    fixExtremities(genH)
    genH.SetDirectory(0)
    genH.SetFillStyle(fill)
    genH.SetLineColor(ci)
    genH.SetMarkerColor(ci)
    genH.SetLineWidth(2)
    genH.SetMarkerStyle(marker)
    
    return passH,genH

def getEfficiencyCurve(passH,genH,marker):
    
    """builds an efficiency curve"""

    effGr=ROOT.TGraphAsymmErrors()
    effGr.SetName(passH.GetName()+'_eff')
    effGr.SetTitle(passH.GetTitle())
    effGr.SetMarkerStyle(marker)
    effGr.SetMarkerColor(passH.GetLineColor())
    effGr.SetLineColor(passH.GetLineColor())
    effGr.SetFillStyle(0)
    effGr.BayesDivide(passH,genH)

    return effGr


def getCurveRatio(gr_den,gr_num):

    """compute the ratio of the two curves"""
    
    compGr=ROOT.TGraphErrors()
    compGr.SetTitle(gr_num.GetTitle())
    compGr.SetName('%s_over_%s'%(gr_num.GetName(),gr_den.GetName()))
    compGr.SetLineColor(gr_num.GetLineColor())
    compGr.SetMarkerColor(gr_num.GetMarkerColor())
    compGr.SetMarkerStyle(gr_num.GetMarkerStyle())

    x_den,x_num,y_den,y_num=ROOT.Double(0),ROOT.Double(0),ROOT.Double(0),ROOT.Double(0)
    for i in xrange(0,gr_den.GetN()):
        gr_den.GetPoint(i,x_den,y_den)
        ey_den=gr_den.GetErrorY(i)

        found=False
        for j in xrange(0,gr_num.GetN()):
            gr_num.GetPoint(j,x_num,y_num)
            ey_num=gr_num.GetErrorY(j)

            if float(x_num)!=float(x_den): continue
            found=True
            break

        if not found:
            print 'Unable to find same point for',float(x_den)
            continue

        if float(y_den)==0 : continue
        newy=float(y_num)/float(y_den)
        newyUnc=ROOT.TMath.Sqrt((float(y_num)*ey_den)**2+(float(y_den)*ey_num)**2)/(float(y_den)**2)

        ip=compGr.GetN()
        compGr.SetPoint(ip,float(x_den),newy)
        compGr.SetPointError(ip,0,newyUnc)

    return compGr


def getCutEfficiencyCurve(h,marker,forward=True):

    """build an efficiency curve for a given cut"""

    eff=ROOT.TGraphErrors()
    eff.SetTitle(h.GetTitle())
    eff.SetMarkerStyle(marker)
    eff.SetMarkerColor(h.GetMarkerColor())
    eff.SetLineColor(h.GetLineColor())
    nbins=h.GetNbinsX()
    err=ROOT.Double(0)        
    total=h.Integral()
    for xbin in xrange(1,nbins+1):        
        n=h.IntegralAndError(xbin,nbins,err)
        eff.SetPoint(xbin-1,h.GetXaxis().GetBinLowEdge(xbin),n/total)
        eff.SetPointError(xbin-1,0,float(err)/total)
    return eff


def showEfficiencyPlot(grColl,xtit,ytit,logx,logy,yran,legPos,outName):

    """ dump efficiency plots"""

    #show plots
    c=ROOT.TCanvas('c','c',500,500)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.02)
    c.SetTopMargin(0.05)
    c.SetBottomMargin(0.1)
    c.SetGridx()
    c.SetGridy()

    c.SetLogx(logx)
    c.SetLogy(logy)
    leg_y0=0.95 if legPos[0]=='t' else 0.45
    leg_x0=0.15 if legPos[1]=='l' else 0.65
    leg=ROOT.TLegend(leg_x0,leg_y0,leg_x0+0.3,leg_y0-0.04*len(grColl))
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    mg=ROOT.TMultiGraph()
    for g in grColl:
        mg.Add(g,'p')
        leg.AddEntry(g,g.GetTitle(),'ep')
    mg.Draw('a')
    if logx:
        mg.GetXaxis().SetMoreLogLabels()
    if xtit:
        mg.GetXaxis().SetTitle(xtit)
    if ytit:
        mg.GetYaxis().SetTitle(ytit)
    if yran:
        mg.GetYaxis().SetRangeUser(*yran)
    leg.Draw()
    tex=ROOT.TLatex()
    tex.SetTextFont(42)
    tex.SetTextSize(0.04)
    tex.SetNDC()
    tex.DrawLatex(0.12,0.96,'#bf{CMS} #it{simulation preliminary}')
    c.Modified()
    c.Update()   
    c.RedrawAxis()
    c.SaveAs('plots/sel/{0}.png'.format(outName))
    c.SaveAs('plots/sel/{0}.pdf'.format(outName))

def showSimpleDistribution(hColl,doLogx,doDivideByBinWidth,outName):

    c=ROOT.TCanvas('c','c',500,500)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.02)
    c.SetTopMargin(0.05)
    c.SetBottomMargin(0.1)
    c.SetGridx()
    c.SetGridy()
    c.SetLogx(doLogx)
    c.SetLogy(True)
    leg=ROOT.TLegend(0.65,0.95,0.95,0.95-0.04*len(hColl))
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    drawOpt='hist'
    for h in hColl:
        scaleTo(h,1.0)
        h.GetYaxis().SetRangeUser(1e-4,1.0)
        h.GetYaxis().SetTitle('PDF')
        if doDivideByBinWidth: 
            divideByBinWidth(h)
            h.GetYaxis().SetTitle('PDF / bin width')
            h.GetYaxis().SetRangeUser(1e-4,0.1)
        h.Draw(drawOpt)
        leg.AddEntry(h,h.GetTitle(),'f')
        drawOpt='histsame'
    leg.Draw()
    tex=ROOT.TLatex()
    tex.SetTextFont(42)
    tex.SetTextSize(0.04)
    tex.SetNDC()
    tex.DrawLatex(0.12,0.96,'#bf{CMS} #it{simulation preliminary}')    
    c.Modified()
    c.Update()
    c.RedrawAxis()
    c.SaveAs('plots/sel/%s.png'%outName)
    c.SaveAs('plots/sel/%s.pdf'%outName)

def runSelectionEfficiencyFor(ch,d):

    bkgROOT=ROOT.TFile.Open('plots/bkg_gen_plotter.root')
    zxROOT=ROOT.TFile.Open('plots/zx_gen_plotter.root')

    procList  =[('DY','#000000',1001,20,bkgROOT),
                ("EWK Zjj",'#f0027f',0,26,zxROOT),
                ("ZH#rightarrowllbb",'#666666',0,20,zxROOT),
                ("qqH(900)#rightarrowZZ#rightarrow2l2#nu",'#8c510a',0,21,zxROOT),
                ("qqH(2000)#rightarrowZZ#rightarrow2l2#nu",'#358CEF',0,23,zxROOT),
                ]
    
    if 'a' in ch:
        procList  =[('#gamma+jets','#000000',1001,20,bkgROOT),
                    ("EWK #gammajj",'#f0027f',0,26,zxROOT)
                    ]

    
    tagsList=['trig','rec'] #'singletrig','trig','rec','2trec']

    genSummary={}
    effSummary={}
    for tag in tagsList:
        genSummary[tag]=[]
        effSummary[tag]=[]
        for p,color,fill,marker,plotter in procList:
            passH,genH = getPassGenDistributions(pName='gen{0}{1}_{2}/gen{0}{1}_{2}_{3}'.format(ch,tag,d,p),
                                                 gName='gen{0}_{1}/gen{0}_{1}_{2}'.format(ch,d,p),
                                                 plotter=plotter,
                                                 ci=ROOT.TColor.GetColor(color),
                                                 fill=fill,
                                                 marker=marker)        
            effGr=getEfficiencyCurve(passH,genH,marker=marker)

            genSummary[tag].append(genH)
            effSummary[tag].append(effGr)

    #relative gains
    #effSummary['trig_gain']=[]
    #effSummary['rec_gain']=[]
    #for i in xrange(0,len(effSummary[tagsList[0]])):
    #    effSummary['trig_gain'].append(getCurveRatio(gr_den=effSummary['singletrig'][i], gr_num=effSummary['trig'][i]))
    #    effSummary['rec_gain'].append(getCurveRatio(gr_den=effSummary['2trec'][i],    gr_num=effSummary['rec'][i]))

    #cut efficiency
    if d=='ptll':
        effSummary['cut']=[]
        for i in xrange(0,len(effSummary[tagsList[0]])):
            effSummary['cut'].append( getCutEfficiencyCurve(h=genSummary['rec'][i],marker=genSummary['rec'][i].GetMarkerStyle()) )

    doLogx=True if d=='ptll' else False
    doDivideByBinWidth=True if d in ['ptll','mll'] else False

    #distribution plots
    showSimpleDistribution(genSummary['rec'],doLogx=False,doDivideByBinWidth=doDivideByBinWidth,outName='gen_%s_%s'%(ch,d))

    #save results 
    fOut=ROOT.TFile.Open('plots/effsummary_%s_%s.root'%(ch,d),'RECREATE')
    for gr in effSummary['rec']: gr.Write()
    fOut.Close()

    #efficiency plots
    for tag,xtit,ytit,logx,logy,yran,legPos in [('trig',      genSummary['rec'][0].GetXaxis().GetTitle(),'Trigger #varepsilon',  doLogx,False,(0.5,1),  'bl'),
                                                #('trig_gain', genSummary['rec'][0].GetXaxis().GetTitle(),'Ratio to single triggers',          doLogx,False,(0.95,1.15),'tr'),
                                                ('rec',       genSummary['rec'][0].GetXaxis().GetTitle(),'Selection #varepsilon',doLogx,False,(0.,1),   'tl' if d=='ptll' else 'tr'),
                                                #('rec_gain',  genSummary['rec'][0].GetXaxis().GetTitle(),'Ratio to tight selection',          doLogx,False,(0.7,2.2),'tr'),
                                                ('cut',       genSummary['rec'][0].GetXaxis().GetTitle(),'Cut efficiency',                    doLogx,True, (1e-3,1), 'bl'),
                                                ]:

        if not tag in effSummary: continue
        if 'pta' in ch and tag=='trig':
            yran=(0,1)
        showEfficiencyPlot(grColl=effSummary[tag],
                           xtit=xtit,
                           ytit=ytit,
                           logx=logx,
                           logy=logy,
                           yran=yran,
                           legPos=legPos,
                           outName='%s_%s_%s'%(tag,ch,d))



ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gROOT.SetBatch(True)

for d in ['ptboson','mll','drll']:
    for ch in ['ee','eez', 'mm','mmz','a']:
        try:
            runSelectionEfficiencyFor(ch,d)
        except Exception as e:
            print e
            pass


