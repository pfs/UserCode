#! /usr/bin/env python
import ROOT
import optparse
import json
import sys
import os
import time
from array import *

"""
Get correlation factor
"""
def getCorrelation(points, obs1, obs2):
    h = ROOT.TH2D('h', '', 100, 0, 1, 100, 0, 1)
    for i in range(len(points[obs1])):
        val1 = points[obs1][i]
        val2 = points[obs2][i]
        
        if obs1 == 'mult': val1 /= 100.
        if obs2 == 'mult': val2 /= 100.
        
        if (val1 < 0.) or (val2 < 0.): continue
        
        h.Fill(val1, val2)
    
    return h.GetCorrelationFactor()

"""
Get data points from tree
"""
def getPointsFromTree(tree, observables):
    points = {}
    for obs in observables:
        points[obs] = []
    for event in tree:
        if event.gen_sel != 1: continue
        for j in range(event.ngj):
            if event.gj_overlap == 1: continue
            for obs in observables:
                val = eval('event.gj_'+obs+'_charged')[j]
                points[obs].append(val)
    
    return points
    
"""
steer
"""
def main():

    cmsLabel='#bf{CMS} #it{preliminary}'
    
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--input',
                            dest='input',   
                            default='analysis.root',
                            help='input file [default: %default]')
    parser.add_option('-o', '--output',
                            dest='output', 
                            default='',
                            help='Output directory [default: %default]')
    parser.add_option('--ro', '--rootoutput',
                            dest='rootoutput',
                            default='correlations.root',
                            help='output root file [default: %default]')
    parser.add_option(     '--com',          dest='com'  ,        help='center of mass energy',                            default='13 TeV',       type='string')
    parser.add_option('-l', '--lumi',        dest='lumi' ,       help='lumi [/pb]',              default=16551.,              type=float)
    (opt, args) = parser.parse_args()
    
    tree = ROOT.TChain('tjsev')
    tree.Add(opt.input)

    rootoutfile = ROOT.TFile(opt.rootoutput, "RECREATE");
    
    observables = ["mult", "width", "ptd", "ptds", "ecc", "tau21", "tau32", "tau43", "zg", "zgxdr", "zgdr", "ga_width", "ga_lha", "ga_thrust", "c1_02", "c1_05", "c1_10", "c1_20", "c2_02", "c2_05", "c2_10", "c2_20", "c3_02", "c3_05", "c3_10", "c3_20"]
    
    nice_observables_root = {"mult": "N", "width": "width", "ptd": "p_{T}D", "ptds": "p_{T}D*", "ecc": "#varepsilon", "tau21": "#tau_{21}", "tau32": "#tau_{32}", "tau43": "#tau_{43}", "zg": "z_{g}", "zgxdr": "z_{g} #times #DeltaR", "zgdr": "z_{g} #DeltaR", "ga_width": "#lambda_{1}^{1}", "ga_lha": "#lambda_{0.5}^{1}", "ga_thrust": "#lambda_{2}^{1}", "c1_02": "C_{1}^{(0.2)}", "c1_05": "C_{1}^{(0.5)}", "c1_10": "C_{1}^{(1.0)}", "c1_20": "C_{1}^{(2.0)}", "c2_02": "C_{2}^{(0.2)}", "c2_05": "C_{2}^{(0.5)}", "c2_10": "C_{2}^{(1.0)}", "c2_20":  "C_{2}^{(2.0)}", "c3_02": "C_{3}^{(0.2)}", "c3_05": "C_{3}^{(0.5)}", "c3_10": "C_{3}^{(1.0)}", "c3_20": "C_{3}^{(2.0)}"}
    
    points = getPointsFromTree(tree, observables)
    
    # 11-class RdBu http://colorbrewer2.org/#type=diverging&scheme=RdBu&n=11
    stops = array('d', [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    red   = array('d')
    green = array('d')
    blue  = array('d')
    colors = [[103,0,31],
              [178,24,43],
              [214,96,77],
              [244,165,130],
              [253,219,199],
              [247,247,247],
              [209,229,240],
              [146,197,222],
              [67,147,195],
              [33,102,172],
              [5,48,97]]
    for color in colors:
        red.append(color[0]/255.)
        green.append(color[1]/255.)
        blue.append(color[2]/255.)
    ROOT.TColor.CreateGradientColorTable(11, stops, red[::-1], green[::-1], blue[::-1], 30)
    
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPaintTextFormat('+2.0f') 
    c = ROOT.TCanvas('c','c',500,500)
    c.SetRightMargin(0.15)
    c.SetLeftMargin(0.12)
    c.SetTopMargin(0.1)
    c.SetTopMargin(0.05)
    c.cd()
    
    h_correlations = ROOT.TH2D('h_correlations', '', len(observables), 0, len(observables), len(observables), 0, len(observables))
    
    for i in range(len(observables)):
        h_correlations.GetXaxis().SetBinLabel(i+1, nice_observables_root[observables[i]])
        h_correlations.GetYaxis().SetBinLabel(i+1, nice_observables_root[observables[i]])
        for j in range(len(observables)):
            correlation = 100.*getCorrelation(points, observables[i], observables[j])
            h_correlations.SetBinContent(i+1, j+1, correlation)
    
    h_correlations.GetXaxis().LabelsOption('v')
    h_correlations.GetZaxis().SetRangeUser(-100., 100.)
    h_correlations.GetZaxis().SetTitle('Correlation [%]     ')
    h_correlations.SetMarkerSize(0.7)
    h_correlations.Draw('colz,text')
    
    ROOT.gPad.Update()
    tl1 = ROOT.TLine(ROOT.gPad.GetUxmin(), ROOT.gPad.GetUymax(),
                     ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax())
    tl1.Draw()
    tl2 = ROOT.TLine(ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymin(),
                     ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax())
    tl2.Draw()
    
    txt=ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(42)
    txt.SetTextSize(0.041)
    txt.SetTextAlign(12)
    txt.DrawLatex(0.16,0.97, cmsLabel)
    txt.DrawLatex(0.63,0.97, '#scale[0.8]{%3.1f fb^{-1} (%s)}' % (opt.lumi/1000.,opt.com) )
    
    c.Print('correlations.pdf')
    c.Print('correlations.png')
    
    observables_low = ["ptds", "ecc", "tau43", "zg", "zgdr"]
    
    h_correlations_low = ROOT.TH2D('h_correlations_low', '', len(observables_low), 0, len(observables_low), len(observables_low), 0, len(observables_low))
    
    for i in range(len(observables_low)):
        h_correlations_low.GetXaxis().SetBinLabel(i+1, nice_observables_root[observables_low[i]])
        h_correlations_low.GetYaxis().SetBinLabel(i+1, nice_observables_root[observables_low[i]])
        for j in range(len(observables_low)):
            correlation = 100.*getCorrelation(points, observables_low[i], observables_low[j])
            h_correlations_low.SetBinContent(i+1, j+1, correlation)
    
    h_correlations_low.GetXaxis().LabelsOption('h')
    h_correlations_low.GetZaxis().SetRangeUser(-100., 100.)
    h_correlations_low.GetZaxis().SetTitle('Correlation [%]     ')
    h_correlations_low.SetMarkerSize(1.)
    h_correlations_low.Draw('colz,text')
    
    ROOT.gPad.Update()
    tl1 = ROOT.TLine(ROOT.gPad.GetUxmin(), ROOT.gPad.GetUymax(),
                     ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax())
    tl1.Draw()
    tl2 = ROOT.TLine(ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymin(),
                     ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax())
    tl2.Draw()
    
    txt=ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(42)
    txt.SetTextSize(0.041)
    txt.SetTextAlign(12)
    txt.DrawLatex(0.16,0.97, cmsLabel)
    txt.DrawLatex(0.63,0.97, '#scale[0.8]{%3.1f fb^{-1} (%s)}' % (opt.lumi/1000.,opt.com) )
    
    c.Print('correlations_low.pdf')
    c.Print('correlations_low.png')

if __name__ == "__main__":
	sys.exit(main())
