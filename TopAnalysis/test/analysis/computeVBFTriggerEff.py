import optparse
import ROOT
import os,sys
from TopLJets2015.TopAnalysis.Plot import fixExtremities
from computeTransferFactor import getPlotsIn,scaleTo

def parametrizeTurnOn(num,den,turnOn=None):
        
    """parametrizes the turn-on with a given function, and computes the ratio parametrized"""

    if not turnOn: return None

    num.Fit(turnOn,'MR+')
    num.GetListOfFunctions().At(0).SetLineColor(1)
    grnum=ROOT.TGraphErrors()
    for i in xrange(0,num.GetNbinsX()+1):
        grnum.SetPoint(i, num.GetXaxis().GetBinCenter(i+1),0)
    ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(grnum,0.68)
            
    den.Fit(turnOn,'MR+')
    den.GetListOfFunctions().At(0).SetLineColor(ROOT.kBlue)
    grden=ROOT.TGraphErrors()
    for i in xrange(0,den.GetNbinsX()+1):
        grden.SetPoint(i, den.GetXaxis().GetBinCenter(i+1),0)
    ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(grden,0.68)
            
    #compute the parametrized ratio
    data2mcParam=ROOT.TGraphErrors()
    x,ynum,yden=ROOT.Double(0),ROOT.Double(0),ROOT.Double(0)
    for i in xrange(0,grnum.GetN()):
        grden.GetPoint(i,x,yden)
        ydenUnc=grden.GetErrorY(i)
        grnum.GetPoint(i,x,ynum)
        ynumUnc=grnum.GetErrorY(i)
                
        if float(yden)==0 : continue
        ip=data2mcParam.GetN()
        data2mcParam.SetPoint(ip,float(x),float(ynum)/float(yden))
        eratio =  (float(ynum)*ydenUnc)**2
        eratio += (float(yden)*ynumUnc)**2
        eratio = ROOT.TMath.Sqrt(eratio)/(float(yden)**2)
        data2mcParam.SetPointError(ip,0.5*den.GetXaxis().GetBinWidth(i+1),eratio)
    data2mcParam.SetTitle('Data/MC param')
    data2mcParam.SetFillStyle(1001)
    data2mcParam.SetFillColor(ROOT.kGray)
    data2mcParam.SetMarkerColor(ROOT.kGray)
    data2mcParam.SetLineColor(ROOT.kGray)

    return data2mcParam

def computeVBFTriggerEff(f,distsOfInterest):
    
    """computes ratios to evaluate VBF trigger efficiency and fits an erf function"""

    ratios={}
    inF=ROOT.TFile.Open(f)
    for key in inF.GetListOfKeys():
        name=key.GetName()
        if not 'HighPtVBFA' in name : continue
        dist = '_'.join(name.split('_')[1:])
        if not dist in distsOfInterest:continue

        dataNum,mcNum,_=getPlotsIn(inF,name,rebin=False)
        dataDen,mcDen,_=getPlotsIn(inF,name.replace('HighPtVBFA','HighPtOfflineVBFA'),rebin=False)
        for h in [dataNum,mcNum,dataDen,mcDen] : 
            fixExtremities(h)
            scaleTo(h,1.0)

        dataNum.Divide(dataDen)
        dataNum.GetYaxis().SetTitle('High p_{T} VBF #gamma / High p_{T} #gamma')        
        dataNum.SetTitle('Data')
        dataNum.SetMarkerStyle(20)

        mcNum.Divide(mcDen)
        mcNum.GetYaxis().SetTitle('High p_{T} VBF #gamma / High p_{T} #gamma')
        mcNum.SetTitle('MC')
        mcNum.SetMarkerStyle(24)
        mcNum.SetMarkerColor(ROOT.kGray)
        mcNum.SetLineColor(ROOT.kGray)
        mcNum.SetFillStyle(1001)
        mcNum.SetFillColor(ROOT.kGray)

        data2mc=dataNum.Clone('{0}_data2mc'.format(dist))
        data2mc.Divide(mcNum)
        data2mc.SetDirectory(0)
        data2mc.GetYaxis().SetTitle('Data/MC')

        turnOn=None
        if dist in ['mjj','subleadpt'] : 
            turnOn=ROOT.TF1('turnon',
                            '[0]+[1]/(1.+TMath::Exp([2]*(x-[3])))',
                            dataNum.GetXaxis().GetXmin(),
                            dataNum.GetXaxis().GetXmax())
            turnOn.SetParLimits(0,0.,0.8)
            turnOn.SetParLimits(1,0.5,1.2)
            turnOn.SetParLimits(2,-2,2)
            turnOn.SetParLimits(3,500,2000)           
            if dist=='subleadpt': turnOn.SetParLimits(3,20,50)            
        #else:
        #    turnOn=ROOT.TF1('turnon',
        #                    '[0]',
        #                    dataNum.GetXaxis().GetXmin(),
        #                    dataNum.GetXaxis().GetXmax())
        data2mcParam=parametrizeTurnOn(dataNum,mcNum,turnOn)
        ratios[name]=(dataNum,mcNum,data2mc,data2mcParam)

    inF.Close()

    return ratios

def showTriggerEff(triggerEff,outName):

    """shows the trigger efficiency ratio plots"""
    
    data,mc,data2mc,ratioParam=triggerEff


    c=ROOT.TCanvas('c','c',500,500)
    c.SetBottomMargin(0)
    c.SetTopMargin(0)
    c.SetLeftMargin(0)
    c.SetRightMargin(0)
    c.cd()

    p1=ROOT.TPad('p1','p1',0,0.5,1,1.0)
    p1.Draw()
    p1.SetRightMargin(0.03)
    p1.SetLeftMargin(0.12)
    p1.SetTopMargin(0.1)
    p1.SetBottomMargin(0.01)
    p1.SetGridy()
    p1.cd()
    mc.Draw('e2')
    mc.GetYaxis().SetNdivisions(5)
    mc.GetXaxis().SetTitleSize(0)
    mc.GetXaxis().SetLabelSize(0)
    mc.GetYaxis().SetTitleSize(0.08)
    mc.GetYaxis().SetTitleOffset(0.7)
    mc.GetYaxis().SetLabelSize(0.08)
    mc.GetYaxis().SetRangeUser(0.52,1.48)
    data.Draw('e1same')

    leg1=p1.BuildLegend(0.15,0.88,0.6,0.66)
    leg1.SetFillStyle(0)
    leg1.SetBorderSize(0)
    leg1.SetTextFont(42)
    leg1.SetTextSize(0.08)

    txt=ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(42)
    txt.SetTextSize(0.08)
    txt.SetTextAlign(12)
    txt.DrawLatex(0.12,0.95,'#bf{CMS} #it{preliminary}')
    p1.RedrawAxis()

    c.cd()
    p2=ROOT.TPad('p2','p2',0,0,1,0.5)
    p2.SetRightMargin(0.03)
    p2.SetLeftMargin(0.12)
    p2.SetTopMargin(0.01)
    p2.SetBottomMargin(0.18)
    p2.SetGridy()
    p2.Draw()
    p2.cd()
    frame=data2mc.Clone('frame')
    frame.Reset('ICE')
    frame.Draw()
    frame.GetYaxis().SetNdivisions(5)
    frame.GetYaxis().SetRangeUser(0.52,1.48)
    frame.GetXaxis().SetTitleSize(0.08)
    frame.GetXaxis().SetLabelSize(0.08)
    frame.GetYaxis().SetTitleSize(0.08)
    frame.GetYaxis().SetLabelSize(0.08)
    frame.GetYaxis().SetTitleOffset(0.7)
    if ratioParam: ratioParam.Draw('2')
    data2mc.Draw('e1same')

    if ratioParam:
        leg2=ROOT.TLegend(0.15,0.92,0.6,0.7)
        leg2.SetFillStyle(0)
        leg2.SetBorderSize(0)
        leg2.SetTextFont(42)
        leg2.SetTextSize(0.08)
        leg2.AddEntry(ratioParam,ratioParam.GetTitle(),'f')
        leg2.AddEntry(data2mc,data2mc.GetTitle(),'ep')
        leg2.Draw()

    p2.RedrawAxis()

    c.cd()
    c.Modified()
    c.Update()
    for ext in ['png','pdf']:
        c.SaveAs('{0}.{1}'.format(outName,ext))
    

def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-p',  dest='plotter',      help='input plotter [%default]',        default=None,    type='string')
    parser.add_option('-o',  dest='outDir',       help='output directory [%default]',     default=None,    type='string')
    parser.add_option('--titles',  dest='titles', help='titles (CSV list) [%default]',    default=None,    type='string')
    parser.add_option('-t',  dest='triggerBased', help='trigger-based ratio [%default]',  default=False,   action='store_true')
    (opt, args) = parser.parse_args()

    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gROOT.SetBatch(True)

    distsOfInterest=['detajj','dphijj','drvj0','drvj1','forwardeta','centraleta','mjj','nvtx','leadpt','subleadpt','vpt','vy']
    trigEff=computeVBFTriggerEff(opt.plotter,distsOfInterest)
    for name in trigEff:
        fOutName=name
        if opt.outDir: fOutName=os.path.join(opt.outDir,fOutName)
        showTriggerEff(trigEff[name],fOutName+'_trigeff')

if __name__ == "__main__":
    sys.exit(main())
