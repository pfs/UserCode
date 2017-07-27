import ROOT
ROOT.gROOT.SetBatch(True)
import optparse
import os,sys
import json
import re
from collections import OrderedDict
from math import sqrt

debug = True

"""
steer the script
"""
def main():
    
    cmsLabel='#bf{CMS} #it{preliminary}'
    
    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option(     '--mcUnc',        dest='mcUnc'  ,      help='common MC related uncertainty (e.g. lumi)',        default=0,              type=float)
    parser.add_option(     '--com',          dest='com'  ,        help='center of mass energy',                            default='13 TeV',       type='string')
    parser.add_option('-j', '--json',        dest='json'  ,      help='json with list of files',        default='data/era2016/samples.json',              type='string')
    parser.add_option( '--systJson', dest='systJson', help='json with list of systematics', default='data/era2016/syst_samples.json', type='string')
    parser.add_option('-i', '--inDir',       dest='inDir' ,      help='input directory',                default='unfolding/result',              type='string')
    parser.add_option('-O', '--outDir',      dest='outDir' ,     help='output directory',                default='unfolding/result',              type='string')
    parser.add_option('-o', '--outName',     dest='outName' ,    help='name of the output file',        default='plotter.root',    type='string')
    parser.add_option(      '--silent',      dest='silent' ,     help='only dump to ROOT file',         default=False,             action='store_true')
    parser.add_option(      '--saveTeX',     dest='saveTeX' ,    help='save as tex file as well',       default=False,             action='store_true')
    parser.add_option('-l', '--lumi',        dest='lumi' ,       help='lumi [/pb]',              default=16551.,              type=float)
    parser.add_option('--obs', dest='obs',  default='mult', help='observable [default: %default]')
    parser.add_option('--flavor', dest='flavor',  default='all', help='flavor [default: %default]')
    (opt, args) = parser.parse_args()

    #read lists of samples
    betas = [0.2, 0.5, 1., 2.]
    sbetas = ['02', '05', '10', '20']
    colors = {'all': ROOT.kBlack, 'bottom': ROOT.kRed+1, 'light': ROOT.kBlue+1, 'gluon': ROOT.kGreen+1}
    markers = {'all': 20, 'bottom': 21, 'light': 22, 'gluon': 23}
    fills = {'all': 1001, 'bottom': 3254, 'light': 3245, 'gluon': 3390}
    infiles = {}
    hists = {}
    unchists = {}
    
    dataUnfolded = ROOT.TH1F('dataUnfolded', 'dataUnfolded', 12, 0, 12)
    dataUnfoldedSys = ROOT.TH1F('dataUnfoldedSys', 'dataUnfoldedSys', 12, 0, 12)
    nominalGen = ROOT.TH1F('nominalGen', 'nominalGen', 12, 0, 12)
    FSRUpGen = ROOT.TH1F('FSRUpGen', 'FSRUpGen',       12, 0, 12)
    FSRDownGen = ROOT.TH1F('FSRDownGen', 'FSRDownGen', 12, 0, 12)
    herwigGen = ROOT.TH1F('herwigGen', 'herwigGen',    12, 0, 12)
    
    counter = 0
    for c in range(1,4):
        for sbeta in sbetas:
            counter += 1
            infile = ROOT.TFile.Open('%s/c%i_%s_charged_%s_result.root'%(opt.inDir, c, sbeta, opt.flavor))
            
            inhist = infile.Get('mean')
            dataUnfolded.SetBinContent(counter, inhist.GetBinContent(1))
            dataUnfolded.SetBinError(counter, inhist.GetBinError(1))
            
            inhistErr = infile.Get('meanErr')
            dataUnfoldedSys.SetBinContent(counter, inhistErr.GetBinContent(1))
            dataUnfoldedSys.SetBinError(counter, inhistErr.GetBinError(1))
            dataUnfoldedSys.GetXaxis().SetBinLabel(counter, 'C_{'+str(c)+'}^{('+str(int(sbeta)/10.)+')}')
            
            innominalGen = infile.Get('nominalGen')
            nominalGen.SetBinContent(counter, innominalGen.GetMean())
            nominalGen.SetBinError(counter, innominalGen.GetMeanError())
            
            inFSRUpGen = infile.Get('FSRUpGen')
            FSRUpGen.SetBinContent(counter, inFSRUpGen.GetMean())
            FSRUpGen.SetBinError(counter, inFSRUpGen.GetMeanError())
            
            inFSRDownGen = infile.Get('FSRDownGen')
            FSRDownGen.SetBinContent(counter, inFSRDownGen.GetMean())
            FSRDownGen.SetBinError(counter, inFSRDownGen.GetMeanError())
            
            inHerwigGen = infile.Get('herwigGen')
            herwigGen.SetBinContent(counter, inHerwigGen.GetMean())
            herwigGen.SetBinError(counter, inHerwigGen.GetMeanError())
    
    #plot
    ROOT.gStyle.SetOptStat(0)
    c = ROOT.TCanvas('c','c',500,500)
    c.SetBottomMargin(0.0)
    c.SetLeftMargin(0.0)
    c.SetTopMargin(0)
    c.SetRightMargin(0.00)
    c.cd()
    
    p1=ROOT.TPad('p1','p1',0.0,0.2,1.0,1.0)
    p1.SetRightMargin(0.05)
    p1.SetLeftMargin(0.12)
    p1.SetTopMargin(0.06)
    p1.SetBottomMargin(0.01)
    p1.Draw()
    p1.cd()
    
    dataUnfolded.SetTitle('')
    flavor = opt.flavor
    if flavor in ['light', 'gluon']: flavor += '-enriched'
    dataUnfolded.SetXTitle(flavor+' jets')
    dataUnfolded.GetXaxis().SetTitleSize(0.045)
    dataUnfolded.GetXaxis().SetLabelSize(0.04)
    dataUnfolded.SetYTitle('<C_{N}^{(#beta)}>')
    dataUnfolded.GetYaxis().SetRangeUser(0.0001, 1.5*dataUnfolded.GetMaximum())
    dataUnfolded.GetYaxis().SetTitleSize(0.05)
    dataUnfolded.GetYaxis().SetLabelSize(0.045)
    dataUnfolded.GetYaxis().SetTitleOffset(1.1)
    dataUnfolded.SetLineColor(ROOT.kBlack)
    dataUnfolded.SetMarkerColor(ROOT.kBlack)
    dataUnfolded.SetMarkerStyle(20)
    dataUnfolded.Draw('P X0 E1')
    
    dataUnfoldedSys.SetLineWidth(2)
    dataUnfoldedSys.SetLineColor(ROOT.kBlack)
    dataUnfoldedSys.SetMarkerColor(ROOT.kBlack)
    dataUnfoldedSys.SetMarkerStyle(1)
    dataUnfoldedSys.Draw('SAME P X0 E1')
    
    dataUnfolded.SetLineColor(ROOT.kBlack)
    dataUnfolded.SetMarkerColor(ROOT.kBlack)
    dataUnfolded.SetMarkerStyle(20)
    dataUnfolded.Draw('SAME P X0 E1')
    
    nominalGen.SetLineColor(ROOT.kRed+1)
    nominalGen.SetLineWidth(2)
    nominalGen.SetMarkerColor(ROOT.kRed+1)
    nominalGen.SetMarkerStyle(24)
    nominalGen.Draw('SAME H')
    
    FSRUpGen.SetLineColor(ROOT.kRed+1)
    FSRUpGen.SetMarkerColor(ROOT.kRed+1)
    FSRUpGen.SetMarkerStyle(26)
    FSRUpGen.Draw('SAME P X0 E1')
    
    FSRDownGen.SetLineColor(ROOT.kRed+1)
    FSRDownGen.SetMarkerColor(ROOT.kRed+1)
    FSRDownGen.SetMarkerStyle(32)
    FSRDownGen.Draw('SAME P X0 E1')
    
    herwigGen.SetLineColor(ROOT.kBlue+1)
    herwigGen.SetLineStyle(7)
    herwigGen.SetLineWidth(2)
    herwigGen.SetMarkerColor(ROOT.kBlue+1)
    herwigGen.SetMarkerStyle(25)
    herwigGen.Draw('SAME H')
    
    inix = 0.5
    if (nominalGen.GetMaximumBin() > nominalGen.GetNbinsX()/2.): inix = 0.15
    legend = ROOT.TLegend(inix,0.625,inix+0.35,0.9)
    legend.SetLineWidth(0)
    legend.SetFillStyle(0)
    legend.AddEntry(dataUnfolded, "Data", "ep")
    legend.AddEntry(nominalGen, "Powheg+Pythia 8", "pl")
    legend.AddEntry(FSRUpGen, "#minus FSR up", "p")
    legend.AddEntry(FSRDownGen, "#minus FSR down", "p")
    legend.AddEntry(herwigGen, "Powheg+Herwig++", "pl")
    legend.Draw()
    txt=ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(42)
    txt.SetTextSize(0.05)
    txt.SetTextAlign(12)
    inix = 0.15
    if (nominalGen.GetMaximumBin() > nominalGen.GetNbinsX()/2.): inix = 0.64
    txt.DrawLatex(inix,0.88,cmsLabel)
    txt.DrawLatex(0.7,0.97,'#scale[0.8]{%3.1f fb^{-1} (%s)}' % (opt.lumi/1000.,opt.com) )
    
    c.cd()
    p2 = ROOT.TPad('p2','p2',0.0,0.0,1.0,0.2)
    p2.Draw()
    p2.SetBottomMargin(0.4)
    p2.SetRightMargin(0.05)
    p2.SetLeftMargin(0.12)
    p2.SetTopMargin(0.01)
    p2.cd()
    
    dataUnfoldedRatio=dataUnfolded.Clone('dataUnfoldedRatio')
    dataUnfoldedRatio.Divide(dataUnfolded)
    dataUnfoldedRatio.SetFillStyle(3254)
    dataUnfoldedRatio.SetFillColor(ROOT.kBlack)
    
    dataUnfoldedSysRatio=dataUnfoldedSys.Clone('dataUnfoldedSysRatio')
    dataUnfoldedSysRatio.Divide(dataUnfolded)
    dataUnfoldedSysRatio.SetTitle('')
    dataUnfoldedSysRatio.SetXTitle(dataUnfolded.GetXaxis().GetTitle())
    dataUnfoldedSysRatio.SetYTitle('Ratio ')
    dataUnfoldedSysRatio.SetFillColor(ROOT.kGray)
    dataUnfoldedSysRatio.GetXaxis().SetTitleSize(0.2)
    dataUnfoldedSysRatio.GetXaxis().SetTitleOffset(0.8)
    dataUnfoldedSysRatio.GetXaxis().SetLabelSize(0.18)
    dataUnfoldedSysRatio.GetYaxis().SetTitleSize(0.2)
    dataUnfoldedSysRatio.GetYaxis().SetTitleOffset(0.3)
    dataUnfoldedSysRatio.GetYaxis().SetLabelSize(0.18)
    dataUnfoldedSysRatio.GetYaxis().SetRangeUser(0.9,1.25)
    dataUnfoldedSysRatio.GetYaxis().SetNdivisions(503)
    
    nominalGenRatio=nominalGen.Clone('nominalGenRatio')
    nominalGenRatio.Divide(dataUnfolded)
    FSRUpGenRatio=FSRUpGen.Clone('FSRUpGenRatio')
    FSRUpGenRatio.Divide(dataUnfolded)
    FSRDownGenRatio=FSRDownGen.Clone('FSRDownGenRatio')
    FSRDownGenRatio.Divide(dataUnfolded)
    herwigGenRatio=herwigGen.Clone('herwigGenRatio')
    herwigGenRatio.Divide(dataUnfolded)
    
    dataUnfoldedSysRatio.SetMarkerStyle(0)
    dataUnfoldedSysRatio.Draw('e2')
    dataUnfoldedRatio.SetMarkerStyle(0)
    dataUnfoldedRatio.Draw('e2,same')
    line = dataUnfoldedSysRatio.Clone('line')
    line.SetLineColor(ROOT.kBlack)
    line.SetFillStyle(0)
    for i in range(line.GetNbinsX()+2): line.SetBinContent(i, 1.)
    line.Draw('hist same')
    
    nominalGenRatio.Draw('SAME H')
    FSRUpGenRatio.Draw  ('SAME P X0 E1')
    FSRDownGenRatio.Draw('SAME P X0 E1')
    herwigGenRatio.Draw ('SAME H')
            
    c.Print(opt.outDir+'/meanCvsBeta_'+opt.flavor+'.pdf')
    c.Print(opt.outDir+'/meanCvsBeta_'+opt.flavor+'.png')
        
def normalizeAndDivideByBinWidth(hist):
    hist.Scale(1./hist.Integral())
    for i in range(1, hist.GetNbinsX()+1):
        hist.SetBinContent(i, hist.GetBinContent(i)/hist.GetBinWidth(i))
        hist.SetBinError  (i, hist.GetBinError(i)  /hist.GetBinWidth(i))
    return hist

"""
for execution from another script
"""
if __name__ == "__main__":
    main()
    #sys.exit(main())

