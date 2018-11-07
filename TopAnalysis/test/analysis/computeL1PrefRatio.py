import optparse
import ROOT
import os
import sys
import numpy as np
from itertools import product
from TopLJets2015.TopAnalysis.Plot import fixExtremities

def rebinUnequalBinSize(h,newBins):
    print newBins
    nBin = len(newBins)
    newname = h.GetName()+'_new'
    hnew = h.Rebin(nBin-1, newname, np.array(newBins))
    return hnew

def getPlotsIn(inF,dirname,LPdirname,rebin=True,var=None,bins=[]):

    """gets the data and mc sums in a given directory"""

    data,dataLP=None,None
    if  not inF.Get(LPdirname): return None,None
    LPKeys = inF.Get(LPdirname).GetListOfKeys()
    Keys   = inF.Get(dirname).GetListOfKeys()
    for key, LPkey in product(Keys, LPKeys):
        name=key.GetName()
        nameLP=LPkey.GetName()
        if ('Graph' in name): continue 
        if ('Graph' in nameLP): continue        
        if (not name==dirname) and (not nameLP==LPdirname): continue
        #data histogram
        data_=key.ReadObj()
        dataLP_=LPkey.ReadObj()
        # if data_.InheritsFrom(ROOT.TH2.Class()): continue
        # if dataLP_.InheritsFrom(ROOT.TH2.Class()): continue

        if rebin:
            data_.Rebin()
            dataLP_.Rebin()
        if ( (not bins == [] ) and (var in name)):
            data = rebinUnequalBinSize(data_,bins)
            dataLP = rebinUnequalBinSize(dataLP_,bins)
        else:
            data = data_.Clone()
            data.SetDirectory(0)
            fixExtremities(data)
            dataLP = dataLP_.Clone()
            dataLP.SetDirectory(0)
            fixExtremities(dataLP)

    return data,dataLP


def computeTransferFactors(plotter,shapeOnly=True,rebin=True,var=None,binList=None):

    """computes the L1 prefirig effect"""
    bins = []
    if not binList== None:
        bins = [float(i) for i in binList.split(',')]

    tfList={}
    pF=ROOT.TFile.Open(plotter)
    for key in pF.GetListOfKeys():
        name=key.GetName()
        if not 'MJJ' in name: continue
        if 'LP' in name: continue
        name=key.GetName()
        out_name=name
        nameLP=name.replace('MJJ','MJJLP')

        tfList[name]={}

        try:
            data,dataLP = getPlotsIn(pF,name,nameLP,rebin=rebin,var=var,bins=bins)
            if data==None or dataLP==None: continue
            #transfer factors
            LPtoDef = dataLP.Clone('{0}_LPtoDef'.format(out_name))            
            LPtoDef.Divide(data)
            LPtoDef.SetMarkerStyle(24)
            LPtoDef.SetDirectory(0)
            LPtoDef.SetTitle('')            

            tfList[name]=LPtoDef

            
        except Exception as e:
            print e
            pass
    pF.Close()

    return tfList

def showTF(tf,outDir):

    """shows the transfer factors"""

    LPtoDef=tf

    c=ROOT.TCanvas('c','c',500,500)
    c.SetRightMargin(0.03)
    c.SetLeftMargin(0.12)
    c.SetTopMargin(0.1)
    c.SetBottomMargin(0.11)
    c.SetGridy()
    
    c.cd()

    LPtoDef.Draw('e1')
#    LPtoDef.GetYaxis().SetNdivisions(5)
    LPtoDef.GetXaxis().SetTitleSize(0.05)
    LPtoDef.GetXaxis().SetLabelSize(0.04)
    LPtoDef.GetYaxis().SetTitleSize(0.05)
    LPtoDef.GetYaxis().SetTitleOffset(1)
    LPtoDef.GetYaxis().SetLabelSize(0.04)
    LPtoDef.GetYaxis().SetRangeUser(0.01,2)
    LPtoDef.GetYaxis().SetTitle('L1P Ratio')

    # leg1=c.BuildLegend(0.15,0.88,0.5,0.66)
    # leg1.SetFillStyle(0)
    # leg1.SetBorderSize(0)
    # leg1.SetTextFont(42)
    # leg1.SetTextSize(0.04)

    l1=ROOT.TLine()
    l1.SetLineWidth(2)
    l1.SetLineColor(ROOT.kBlue)
    l1.DrawLine(LPtoDef.GetXaxis().GetXmin(),1,LPtoDef.GetXaxis().GetXmax(),1)

    txt=ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(42)
    txt.SetTextSize(0.04)
    txt.SetTextAlign(12)
    txt.DrawLatex(0.12,0.95,'#bf{CMS} #it{preliminary}')
    c.RedrawAxis()

    

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
    parser.add_option('-o',  dest='outdir',      help='outdir [%default]',         default=None,        type='string')
    parser.add_option('--rebin',  dest='rebin',  help='rebin histograms [%default]', default=False, action='store_true')
    parser.add_option('--var',  dest='varname',   help='variable name to rebin [%default]', default=None, type='string')
    parser.add_option('--binList',  dest='bins',  help='List of bins with unequal bin sizes [%default]', default=None, type='string')

    (opt, args) = parser.parse_args()

    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gROOT.SetBatch(True)
    
    #get all the transfer factors
    tfList=computeTransferFactors(plotter=opt.plotter,rebin=opt.rebin,var=opt.varname,binList=opt.bins)

    #save to root and show in neat canvas
    fOutName='lp_plotter.root'
    if opt.outdir : fOutName=os.path.join(opt.outdir,fOutName)
    fOut=ROOT.TFile.Open(fOutName,'RECREATE')
    for name in tfList:
        fOut.cd()
        outd=fOut.mkdir(name)
        outd.cd()
        if tfList[name] == {}: continue
        tfList[name].Write()
        fOutName=name
        if opt.outdir: fOutName=os.path.join(opt.outdir,fOutName)
        if len(tfList[name]) == 0: continue
        showTF(tfList[name],fOutName+'_lpRatio')
    fOut.Close()

if __name__ == "__main__":
    sys.exit(main())
