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
    flavors = ['all', 'bottom', 'light', 'gluon']
    colors = {'all': ROOT.kBlack, 'bottom': ROOT.kRed+1, 'light': ROOT.kBlue+1, 'gluon': ROOT.kGreen+1}
    markers = {'all': 20, 'bottom': 21, 'light': 22, 'gluon': 23}
    fills = {'all': 1001, 'bottom': 3254, 'light': 3245, 'gluon': 3390}
    infiles = {}
    hists = {}
    unchists = {}
    
    for flavor in flavors:
        infiles[flavor] = ROOT.TFile.Open('%s/%s_charged_%s_result.root'%(opt.inDir, opt.obs, flavor))
        
        hists[flavor] = infiles[flavor].Get('MC13TeV_TTJets_Unfolded').Clone()
        hists[flavor].SetMarkerColor(ROOT.TColor.GetColorDark(colors[flavor]))
        hists[flavor].SetMarkerStyle(markers[flavor])
        hists[flavor].SetLineColor(colors[flavor])
        hists[flavor].SetFillColor(ROOT.TColor.GetColorBright(colors[flavor]))
        hists[flavor].SetFillStyle(fills[flavor])
        
        unchists[flavor] = infiles[flavor].Get('dataUnfoldedSys').Clone()
        unchists[flavor].SetMarkerStyle(0)
        unchists[flavor].SetFillColor(ROOT.TColor.GetColorBright(colors[flavor]))
        unchists[flavor].SetFillStyle(fills[flavor])
        
        if (colors[flavor] == ROOT.kBlack):
            hists[flavor].SetFillColor(ROOT.kGray)
            unchists[flavor].SetFillColor(ROOT.kGray)
    
    #plot
    ROOT.gStyle.SetOptStat(0)
    c = ROOT.TCanvas('c','c',500,500)
    c.cd()
    c.SetRightMargin(0.05)
    c.SetLeftMargin(0.12)
    c.SetTopMargin(0.1)
    c.SetTopMargin(0.05)
    
    unchists['all'].SetTitle('')
    unchists['all'].GetXaxis().SetTitleSize(0.045)
    unchists['all'].GetXaxis().SetTitleOffset(1.)
    unchists['all'].GetXaxis().SetLabelSize(0.035)
    unchists['all'].GetYaxis().SetTitleSize(0.04)
    unchists['all'].GetYaxis().SetTitleOffset(1.4)
    unchists['all'].GetYaxis().SetLabelSize(0.035)
    unchists['all'].GetYaxis().SetRangeUser(0.0001, unchists['all'].GetMaximum()*1.5)
    unchists['all'].GetYaxis().SetTitle(hists['all'].GetYaxis().GetTitle())
    m = re.search('(.*) (\(.+\))', hists['all'].GetXaxis().GetTitle())
    unchists['all'].GetXaxis().SetTitle(m.group(1)[1:])
    unchists['all'].Draw('e2')
    hists['all'].Draw('same')
    unchists['bottom'].Draw('e2,same')
    hists['bottom'].Draw('same')
    unchists['light' ].Draw('e2,same')
    hists['light' ].Draw('same')
    unchists['gluon' ].Draw('e2,same')
    hists['gluon' ].Draw('same')
    
    inix = 0.5
    if (hists['all'].GetMaximumBin() > hists['all'].GetNbinsX()/2.): inix = 0.15
    legend = ROOT.TLegend(inix,0.7,inix+0.35,0.92)
    legend.SetLineWidth(0)
    legend.SetFillStyle(0)
    legend.AddEntry(hists['all'], 'All jets',         'fep')
    legend.AddEntry(hists['bottom'], 'Bottom jets',   'fep')
    legend.AddEntry(hists['light'], 'Light-enriched', 'fep')
    legend.AddEntry(hists['gluon'], 'Gluon-enriched', 'fep')
    legend.Draw()
    txt=ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(42)
    txt.SetTextSize(0.041)
    txt.SetTextAlign(12)
    inix = 0.15
    if (hists['all'].GetMaximumBin() > hists['all'].GetNbinsX()/2.): inix = 0.64
    txt.DrawLatex(inix,0.90,cmsLabel)
    txt.DrawLatex(0.7,0.97,'#scale[0.8]{%3.1f fb^{-1} (%s)}' % (opt.lumi/1000.,opt.com) )
            
    c.Print(opt.outDir+'/'+opt.obs+'_charged_flavors.pdf')
    c.Print(opt.outDir+'/'+opt.obs+'_charged_flavors.png')
        
"""
for execution from another script
"""
if __name__ == "__main__":
    main()
    #sys.exit(main())

