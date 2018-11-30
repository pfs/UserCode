import optparse
import ROOT
import os
import sys
import numpy as np
from TopLJets2015.TopAnalysis.Plot import fixExtremities

def rebinUnequalBinSize(h,newBins):
    nBin = len(newBins)
    newname = h.GetName()+'_new'
    hnew = h.Rebin(nBin-1, newname, np.array(newBins))
    return hnew

def getPlotsIn(inF,dirname,specList=[],rebin=True,bins=[]):

    """gets the data and mc sums in a given directory"""

    data,mcTotal,h=None,None,None
    mcSpecList=[None]*len(specList)
    for key in inF.Get(dirname).GetListOfKeys():
        name=key.GetName()
        if 'Graph' in name : continue        

        #data histogram
        if name==dirname:
            data_=key.ReadObj()
            if rebin : 
                data_.Rebin(rebin)
            if len(bins)>0: 
                data = rebinUnequalBinSize(data_,bins)
            else:
                data = data_.Clone()
            data.SetDirectory(0)
            fixExtremities(data)

        #predictions
        else:
            h_=key.ReadObj()
            if rebin: 
                h_.Rebin(rebin)
            if len(bins)>0: 
                h = rebinUnequalBinSize(h_,bins)
            else:
                h = h_.Clone();
            proc=name.split('_')[-1]
            
            #check if this specific process should be stored
            try:                
                idx=specList.index(proc)
                mcSpecList[idx]=h.Clone()
                mcSpecList[idx].SetDirectory(0)
                fixExtremities(mcSpecList[idx])
            except:
                pass

            #add to the total
            if not mcTotal:
                mcTotal=h.Clone(dirname+'_mctot')
                mcTotal.SetDirectory(0)
                mcTotal.Reset('ICE')
            mcTotal.Add(h)

    if mcTotal : fixExtremities(mcTotal)

    return data,mcTotal,mcSpecList

def scaleTo(h,ht):
    """ scales h to ht(arget) """
    total_h=h.Integral(1,h.GetNbinsX())
    if type(ht) is float:
        h.Scale(ht/total_h)
    else:
        total_ht=ht.Integral(1,ht.GetNbinsX())
        h.Scale(total_ht/total_h)


def applyTF(h,tf):
    """apply bin-by-bin transfer factor to an histogram"""
    
    htf=h.Clone('{0}_tf'.format(h.GetName()))

    for xbin in xrange(0,h.GetNbinsX()+2):
        val=h.GetBinContent(xbin)
        valUnc=h.GetBinError(xbin)
        tfVal=tf.GetBinContent(xbin)
        tfValUnc=tf.GetBinError(xbin)
        htf.SetBinContent(xbin,val*tfVal)
        htf.SetBinError(xbin,ROOT.TMath.Sqrt((val*tfValUnc)**2+(valUnc*tfVal)**2))

    return htf

def computeTransferFactors(plotter,systPlotter,shapeOnly=True,rebin=True,var=None,binList=None):

    """computes the Z NLO/Z LO, A LO/Z LO, Data/LO MC transfer factors"""
    bins=[]
    if not binList== None:
        bins = [float(i) for i in binList.split(',')]

    tfList={}
    pF=ROOT.TFile.Open(plotter)
    sF=ROOT.TFile.Open(systPlotter)
    for key in pF.GetListOfKeys():
        name=key.GetName()

        if not 'MM' in name: continue
        if '_th' in name or '_exp' in name : continue #skip systs
        if var and not var in name: continue
        out_name=name.replace('MM','')
        tfList[name]={}

        try:

            #DY next-to-leading order is the default
            data,totalMC,nloMC = getPlotsIn(pF,name,['DY'],rebin=rebin,bins=bins)
            dataEE,totalMCEE,nloMCEE = getPlotsIn(pF,name.replace('MM','EE'),['DY'],rebin=rebin,bins=bins)
            data.Add(dataEE)
            totalMC.Add(totalMCEE)
            for i in xrange(0,len(nloMC)): nloMC[i].Add(nloMCEE[i])

            #DY leading order
            _,_,loMC = getPlotsIn(sF,name,['DYlo'],rebin=rebin,bins=bins)
            _,_,loMCEE = getPlotsIn(sF,name.replace('MM','EE'),['DYlo'],rebin=rebin,bins=bins)
            for i in xrange(0,len(loMC)): loMC[i].Add(loMCEE[i])

            #photon leading order is the default
            data_A,totalMClo_A,loMC_A = getPlotsIn(pF,name.replace('MM','A'),['#gamma+jets'],rebin=rebin,bins=bins)

            #transfer factors
            nlo2lo = nloMC[0].Clone('{0}_nlo2lo'.format(out_name))
            if shapeOnly: scaleTo(loMC[0],nloMC[0])
            nlo2lo.Divide(loMC[0])
            nlo2lo.SetMarkerStyle(24)
            nlo2lo.SetMarkerColor(ROOT.kGray)
            nlo2lo.SetLineColor(ROOT.kGray)
            nlo2lo.SetFillStyle(1001)
            nlo2lo.SetFillColor(ROOT.kGray)
            nlo2lo.SetDirectory(0)
            nlo2lo.SetTitle('NLO/LO')
            totalMClo = totalMC.Clone('{0}_Z_lomctot'.format(out_name))
            totalMClo.Add(nloMC[0],-1.0)
            totalMClo.Add(loMC[0],+1.0)            
            data2lo  = data.Clone('{0}_Z_data2lo'.format(out_name))
            if shapeOnly: scaleTo(totalMClo,data)
            data2lo.Divide(totalMClo)
            data2lo.SetMarkerStyle(24)
            data2lo.SetDirectory(0)
            data2lo.SetTitle('Data/LO MC')            
            data2nlo = data.Clone('{0}_Z_data2nlo'.format(out_name))
            if shapeOnly: scaleTo(totalMC,data)
            data2nlo.Divide(totalMC)
            data2nlo.SetMarkerStyle(20)
            data2nlo.SetDirectory(0)
            data2nlo.SetTitle('Data/NLO MC')

            #transfer factors applied to photons
            nloMC_A = applyTF(loMC_A[0],nlo2lo)
            totalMCnlo_A = totalMClo_A.Clone('{0}_A_mctot'.format(out_name))
            totalMCnlo_A.Add(loMC_A[0],-1.0)
            totalMCnlo_A.Add(nloMC_A,+1.0)
            data2lo_A = data_A.Clone('{0}_A_data2lo'.format(out_name))
            if shapeOnly:scaleTo(totalMClo_A,data_A)
            data2lo_A.Divide(totalMClo_A)
            data2lo_A.SetMarkerStyle(24)
            data2lo_A.SetDirectory(0)
            data2lo_A.SetTitle('Data/LO MC')
            data2nlo_A = data_A.Clone('{0}_A_data2nlo'.format(out_name))
            if shapeOnly:scaleTo(totalMCnlo_A,data_A)
            data2nlo_A.Divide(totalMCnlo_A)
            data2nlo_A.SetMarkerStyle(20)
            data2nlo_A.SetDirectory(0)
            data2nlo_A.SetTitle('Data/TF#timesLO MC')

            tfList[name]=(nlo2lo,data2lo,data2nlo,data2lo_A,data2nlo_A)
            
        except Exception as e:
            print e
            pass

    pF.Close()
    sF.Close()

    return tfList

def showTF(tf,outDir):

    """shows the transfer factors"""

    nlo2lo,data2lo,data2nlo,data2lo_A,data2nlo_A=tf

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
    nlo2lo.Draw('e2')
    nlo2lo.GetYaxis().SetTitle('Z ratio')
    nlo2lo.GetYaxis().SetNdivisions(5)
    nlo2lo.GetXaxis().SetTitleSize(0)
    nlo2lo.GetXaxis().SetLabelSize(0)
    nlo2lo.GetYaxis().SetTitleSize(0.08)
    nlo2lo.GetYaxis().SetTitleOffset(0.8)
    nlo2lo.GetYaxis().SetLabelSize(0.08)
    nlo2lo.GetYaxis().SetRangeUser(0.01,1.94)
    data2lo.Draw('e1same')
    data2nlo.Draw('e1same')

    leg1=p1.BuildLegend(0.7,0.88,0.95,0.66)
    leg1.SetFillStyle(0)
    leg1.SetBorderSize(0)
    leg1.SetTextFont(42)
    leg1.SetTextSize(0.06)

    l1=ROOT.TLine()
    l1.SetLineWidth(2)
    l1.SetLineColor(ROOT.kBlue)
    l1.DrawLine(data2lo.GetXaxis().GetXmin(),1,data2lo.GetXaxis().GetXmax(),1)

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
    data2lo_A.Draw('e1')
    data2lo_A.GetYaxis().SetTitle('#gamma ratio')
    data2lo_A.GetYaxis().SetNdivisions(5)
    data2lo_A.GetYaxis().SetRangeUser(0.01,1.94)
    data2lo_A.GetXaxis().SetTitleSize(0.08)
    data2lo_A.GetXaxis().SetLabelSize(0.08)
    data2lo_A.GetYaxis().SetTitleSize(0.08)
    data2lo_A.GetYaxis().SetLabelSize(0.08)
    data2lo_A.GetYaxis().SetTitleOffset(0.8)
    data2nlo_A.Draw('e1same')
    
    leg2=p2.BuildLegend(0.7,0.94,0.95,0.80)
    leg2.SetFillStyle(0)
    leg2.SetBorderSize(0)
    leg2.SetTextFont(42)
    leg2.SetTextSize(0.06)
    
    l2=ROOT.TLine()
    l2.SetLineColor(ROOT.kBlue)
    l2.SetLineWidth(2)
    l2.DrawLine(data2lo_A.GetXaxis().GetXmin(),1,data2lo_A.GetXaxis().GetXmax(),1)

    p2.RedrawAxis()

    c.cd()
    c.Modified()
    c.Update()
    for ext in ['png','pdf']:
        c.SaveAs('{0}.{1}'.format(outDir,ext))


def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-p',  dest='plotter',     help='main plotter [%default]',   default=None,        type='string')
    parser.add_option('-s',  dest='systPlotter', help='syst plotter [%default]',   default=None,        type='string')
    parser.add_option('-o',  dest='outdir',      help='outdir [%default]',         default=None,        type='string')
    parser.add_option('--rebin',  dest='rebin',  help='rebin histograms [%default]', default=None, type=int)
    parser.add_option('--var',  dest='varname',   help='variable name to rebin [%default]', default=None, type='string')
    parser.add_option('--binList',  dest='bins',  help='List of bins with unequal bin sizes [%default]', default=None, type='string')

    (opt, args) = parser.parse_args()

    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gROOT.SetBatch(True)
    
    #get all the transfer factors
    tfList=computeTransferFactors(plotter=opt.plotter,systPlotter=opt.systPlotter,rebin=opt.rebin,var=opt.varname,binList=opt.bins if opt.bins else None)
    
    #save to root and show in neat canvas
    fOutName='tf_%s_plotter.root'%opt.varname
    if opt.outdir : fOutName=os.path.join(opt.outdir,fOutName)
    fOut=ROOT.TFile.Open(fOutName,'RECREATE')
    for name in tfList:
        fOut.cd()
        outd=fOut.mkdir(name)
        outd.cd()
        for h in tfList[name]: h.Write()
        fOutName=name
        if opt.outdir: fOutName=os.path.join(opt.outdir,fOutName)
        if len(tfList[name]) == 0: continue
        showTF(tfList[name],fOutName+'_tfactor')
    fOut.Close()

if __name__ == "__main__":
    sys.exit(main())
