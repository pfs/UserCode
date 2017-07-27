import ROOT
ROOT.gROOT.SetBatch(True)
import optparse
import os,sys
import json
import re
from collections import OrderedDict
from math import sqrt
from array import *
import random
import numpy
import copy

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
    parser.add_option('-i', '--inDir',       dest='inDir' ,      help='input directory',                default='unfolding/fill',              type='string')
    parser.add_option('-O', '--outDir',      dest='outDir' ,     help='output directory',                default='unfolding/result',              type='string')
    parser.add_option('-o', '--outName',     dest='outName' ,    help='name of the output file',        default='plotter.root',    type='string')
    parser.add_option(      '--silent',      dest='silent' ,     help='only dump to ROOT file',         default=False,             action='store_true')
    parser.add_option(      '--saveTeX',     dest='saveTeX' ,    help='save as tex file as well',       default=False,             action='store_true')
    parser.add_option('-l', '--lumi',        dest='lumi' ,       help='lumi [/pb]',              default=16551.,              type=float)
    parser.add_option('--obs', dest='obs',  default='mult', help='observable [default: %default]')
    parser.add_option('--flavor', dest='flavor',  default='all', help='flavor [default: %default]')
    (opt, args) = parser.parse_args()

    #read lists of samples
    samplesList=[]
    jsonList = opt.json.split(',')
    for jsonPath in jsonList:
        jsonFile = open(jsonPath,'r')
        samplesList += json.load(jsonFile, encoding='utf-8', object_pairs_hook=OrderedDict).items()
        jsonFile.close()
    
    #read lists of syst samples
    systSamplesList=None
    if opt.systJson:
        jsonFile = open(opt.systJson,'r')
        systSamplesList=json.load(jsonFile,encoding='utf-8', object_pairs_hook=OrderedDict).items()
        jsonFile.close()

    allSystVars = ['jec_CorrelationGroupMPFInSitu',
                   'jec_RelativeFSR',
                   'jec_CorrelationGroupUncorrelated',
                   'jec_FlavorPureGluon',
                   'jec_FlavorPureQuark',
                   'jec_FlavorPureCharm',
                   'jec_FlavorPureBottom',
                   'jer',
                   'btag_heavy',
                   'btag_light',
                   'csv_heavy',
                   'csv_light',
                   'tracking'
                  ]
    varList = []
    for var in allSystVars:
        varList.append(var+'_up')
        varList.append(var+'_down')
    expSystSamplesList = []
    for var in varList:
        expSystSamplesList.append(['MC13TeV_TTJets_'+var, [832., 0., '', 't#bar{t} '+var]])
    
    data        = None # data histogram
    backgrounds = {}   # background histograms
    nominal     = None # reference unfolding matrix
    systematics = {}   # systematic unfolding matrices
    
    for slist,isSyst in [ (samplesList,False), (systSamplesList,True), (expSystSamplesList,True) ]:
        for tag,sample in slist: 
            if tag in ['MC13TeV_TTJets_cflip']: continue
            if isSyst and not 't#bar{t}' in sample[3] : continue
            isData = sample[1]
            isSig  = not isSyst and sample[3] == 't#bar{t}'
            isBkg  = not (isData or isSig or isSyst)
            xs = sample[0]

            if (not os.path.isfile('%s/%s.root'%(opt.inDir, tag))): continue
            fIn=ROOT.TFile.Open('%s/%s.root'%(opt.inDir, tag))
            if not fIn : continue
            
            if isData:
                #if not 'SingleMuon' in tag: continue
                if not any(x in tag for x in ['2016G', '2016H']): continue
                h=fIn.Get(opt.obs+'_charged_'+opt.flavor+'_responsematrix_py')
                try:
                    data.Add(h)
                except:
                    data=h
                    data.SetDirectory(0)
                if debug: print(tag, data.GetEntries(), data.Integral())
            
            if isBkg:
                #if not tag in ['MC13TeV_SingleTbar_tW', 'MC13TeV_SingleT_tW', 'MC13TeV_SingleTbar_t', 'MC13TeV_SingleT_t']: continue
                h=fIn.Get(opt.obs+'_charged_'+opt.flavor+'_responsematrix_py')
                h.Scale(opt.lumi*xs)
                integral = h.Integral()
                integralError = 0
                integralRelError = 999.
                for i in range(1, h.GetNbinsX()+1):
                    integralError = sqrt(integralError**2 + h.GetBinError(i)**2)
                if integral > 0: integralRelError = integralError/integral
                if debug: print(tag, integral, integralError, integralRelError)
                if integralRelError > 1:
                    if debug: print(tag, 'Relative integral error > 1, skipping sample')
                    continue
                try:
                    backgrounds[tag].Add(h)
                except:
                    backgrounds[tag]=h
                    backgrounds[tag].SetDirectory(0)
                if debug:
                    integral = backgrounds[tag].Integral()
                    integralError = 0
                    integralRelError = 0
                    for i in range(1, backgrounds[tag].GetNbinsX()+1):
                        integralError = sqrt(integralError**2 + backgrounds[tag].GetBinError(i)**2)
                    if integral > 0: integralRelError = integralError/integral
                    print(tag, integral, integralError, integralRelError)
            
            if isSig:
                h=fIn.Get(opt.obs+'_charged_'+opt.flavor+'_responsematrix')
                h.Scale(opt.lumi*xs)
                try:
                    nominal.Add(h)
                except:
                    nominal=h
                    nominal.SetDirectory(0)
                if debug: print(tag, nominal.GetEntries())
                # signal file should contain systematics weights
                for i in range(1, 21):
                    if i in [17, 19]: continue
                    ih=fIn.Get(opt.obs+'_charged_'+opt.flavor+'_wgt'+str(i)+'_responsematrix')
                    itag = tag + '_wgt' + str(i)
                    try:
                        systematics[itag].Add(ih)
                    except:
                        systematics[itag]=ih
                        systematics[itag].SetDirectory(0)
                    if debug: print(itag, systematics[itag].GetEntries())
            
            if isSyst:
                h=fIn.Get(opt.obs+'_charged_'+opt.flavor+'_responsematrix')
                try:
                    systematics[tag].Add(h)
                except:
                    systematics[tag]=h
                    systematics[tag].SetDirectory(0)
                if debug: print(tag, systematics[tag].GetEntries())

    # call unfolding
    
    rootoutfile = ROOT.TFile.Open(opt.outDir+'/'+opt.obs+'_charged_'+opt.flavor+'_result.root', 'RECREATE')
    rootoutfile.cd()
    
    dataUnfolded, dataFoldedBack, dataBkgSub, dataStatCov = unfold('MC13TeV_TTJets', nominal, backgrounds, data)
    
    systematicUnfolded = {}
    for tag,systematic in systematics.iteritems():
        systematicUnfolded[tag] = unfold(tag, systematic, backgrounds, data)[0]
    
    background_systematics = [
        [
            ['singletop'],
            ['MC13TeV_SingleTbar_tW', 0.0537],
            ['MC13TeV_SingleT_tW', 0.0537],
            ['MC13TeV_SingleTbar_t', 0.051],
            ['MC13TeV_SingleT_t', 0.0405]
        ],
        [
            ['wjets'],
            ['MC13TeV_W0Jets', -0.057],
            ['MC13TeV_W1Jets',  0.101],
            ['MC13TeV_W2Jets',  0.328],
        ]
    ]
    
    for bkg_sys in background_systematics:
        for direction in ['up', 'down']:
            backgrounds_var = copy.copy(backgrounds)
            for item in bkg_sys:
                if len(item) == 1: continue
                try:
                  if direction == 'up':   backgrounds_var[item[0]].Scale(1. + item[1])
                  if direction == 'down': backgrounds_var[item[0]].Scale(1. - item[1])
                except: continue
            tag = 'MC13TeV_TTJets_'+bkg_sys[0][0]+'_'+direction
            print(tag)
            systematicUnfolded[tag] = unfold(tag, nominal, backgrounds, data)[0]
    
    systUp=[0.]
    systDown=[0.]
    for i in range(1, dataUnfolded.GetNbinsX()+1):
        systUp.append(dataUnfolded.GetBinError(i))
        systDown.append(dataUnfolded.GetBinError(i))
        for tag,systematic in systematicUnfolded.iteritems():
            diff = systematic.GetBinContent(i) - dataUnfolded.GetBinContent(i)
            if (diff > 0):
                systUp[i] = sqrt(systUp[i]**2 + diff**2)
            else:
                systDown[i] = sqrt(systDown[i]**2 + diff**2)

    dataUnfoldedSys = dataUnfolded.Clone('dataUnfoldedSys')    
    for i in range(1, dataUnfoldedSys.GetNbinsX()+1):
        dataUnfoldedSys.SetBinContent(i, dataUnfolded.GetBinContent(i) + (systUp[i]-systDown[i])/2.)
        dataUnfoldedSys.SetBinError(i, (systUp[i]+systDown[i])/2.)
    
    # Calculate mean value with uncertainties
    meanData = dataUnfolded.GetMean()
    meanDataErr = dataUnfolded.GetMeanError()
    meanUp = meanDataErr
    meanDown = meanDataErr
    for tag,systematic in systematicUnfolded.iteritems():
        diff = systematic.GetMean() - meanData
        if (diff > 0):
            meanUp = sqrt(meanUp**2 + diff**2)
        else:
            meanDown = sqrt(meanDown**2 + diff**2)
    
    print('mean', meanData, meanUp, meanDown)
    
    mean = ROOT.TH1F('mean', 'mean', 1, 0, 1)
    mean.SetBinContent(1, meanData)
    mean.SetBinError(1, meanDataErr)
    meanErr = ROOT.TH1F('meanErr', 'meanErr', 1, 0, 1)
    meanErr.SetBinContent(1, meanData + (meanUp-meanDown)/2.)
    meanErr.SetBinError(1, (meanUp+meanDown)/2.)
    
    # Plot
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
    dataUnfolded.SetXTitle(dataUnfolded.GetXaxis().GetTitle()[9:]+' ('+opt.flavor+')')
    if opt.flavor == 'gluon' or opt.flavor == 'light':
        dataUnfolded.SetXTitle(dataUnfolded.GetXaxis().GetTitle()[:-1]+'-enriched)')
    dataUnfolded.GetXaxis().SetTitleSize(0.045)
    dataUnfolded.GetXaxis().SetLabelSize(0.04)
    m = re.search('(.*) (\(.+\))', dataUnfolded.GetXaxis().GetTitle())
    dataUnfolded.SetYTitle('1/N_{jet} dN_{jet} / d'+m.group(1)[1:])
    dataUnfolded.GetYaxis().SetRangeUser(0.0001, 1.5*dataUnfolded.GetMaximum())
    dataUnfolded.GetYaxis().SetTitleSize(0.05)
    dataUnfolded.GetYaxis().SetLabelSize(0.045)
    dataUnfolded.GetYaxis().SetTitleOffset(1.2)
    dataUnfolded.SetLineColor(ROOT.kBlack)
    dataUnfolded.SetMarkerColor(ROOT.kBlack)
    dataUnfolded.SetMarkerStyle(20)
    dataUnfolded.Draw('P X0 E1')
    
    dataUnfoldedSys.SetLineWidth(2)
    dataUnfoldedSys.SetLineColor(ROOT.kBlack)
    dataUnfoldedSys.SetMarkerColor(ROOT.kBlack)
    dataUnfoldedSys.SetMarkerStyle(1)
    dataUnfoldedSys.Draw('SAME P X0 E1')
    
    nominalGen = normalizeAndDivideByBinWidth(nominal.ProjectionX("nominalGen"))
    nominalGen.SetLineColor(ROOT.kRed+1)
    nominalGen.SetLineWidth(2)
    nominalGen.SetMarkerColor(ROOT.kRed+1)
    nominalGen.SetMarkerStyle(24)
    nominalGen.Draw('SAME H')
    
    nominalGenRatio=nominalGen.Clone('nominalGenRatio')
    nominalGenRatio.Divide(dataUnfolded)
    
    dataUnfoldedRatio=dataUnfolded.Clone('dataUnfoldedRatio')
    dataUnfoldedRatio.Divide(dataUnfolded)
    dataUnfoldedRatio.SetFillStyle(3254)
    #dataUnfoldedRatio.SetFillColor(ROOT.TColor.GetColor('#99d8c9'))
    #dataUnfoldedRatio.SetFillColor(ROOT.kGray+1)
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
    dataUnfoldedSysRatio.GetYaxis().SetRangeUser(0.4,1.6)
    dataUnfoldedSysRatio.GetYaxis().SetNdivisions(503)
    
    FSRUpGen = normalizeAndDivideByBinWidth(systematics['MC13TeV_TTJets_fsrup'].ProjectionX("FSRUpGen"))
    FSRUpGen.SetLineColor(ROOT.kRed+1)
    FSRUpGen.SetMarkerColor(ROOT.kRed+1)
    FSRUpGen.SetMarkerStyle(26)
    FSRUpGen.Draw('SAME P X0 E1')
    FSRUpGenRatio=FSRUpGen.Clone('FSRUpGenRatio')
    FSRUpGenRatio.Divide(dataUnfolded)
    
    FSRDownGen = normalizeAndDivideByBinWidth(systematics['MC13TeV_TTJets_fsrdn'].ProjectionX("FSRDownGen"))
    FSRDownGen.SetLineColor(ROOT.kRed+1)
    FSRDownGen.SetMarkerColor(ROOT.kRed+1)
    FSRDownGen.SetMarkerStyle(32)
    FSRDownGen.Draw('SAME P X0 E1')
    FSRDownGenRatio=FSRDownGen.Clone('FSRDownGenRatio')
    FSRDownGenRatio.Divide(dataUnfolded)
    
    nominalGenFSR = nominalGen.Clone('nominalGenFSR')
    for i in range(1, nominalGenFSR.GetNbinsX()+1):
        up = FSRUpGen.GetBinContent(i)
        down = FSRDownGen.GetBinContent(i)
        nominalGenFSR.SetBinContent(i, (up+down)/2.)
        nominalGenFSR.SetBinError(i, abs(up-down)/2.)
    nominalGenFSR.SetMarkerSize(0)
    nominalGenFSR.SetFillStyle(3245)
    nominalGenFSR.SetFillColor(ROOT.kRed-9)
    #nominalGenFSR.Draw('same,e2')
    nominalGenFSRRatio=nominalGenFSR.Clone('nominalGenFSRRatio')
    nominalGenFSRRatio.Divide(dataUnfolded)
    
    herwigGen = normalizeAndDivideByBinWidth(systematics['MC13TeV_TTJets_herwig'].ProjectionX("herwigGen"))
    herwigGen.SetLineColor(ROOT.kBlue+1)
    herwigGen.SetLineStyle(7)
    herwigGen.SetLineWidth(2)
    herwigGen.SetMarkerColor(ROOT.kBlue+1)
    herwigGen.SetMarkerStyle(25)
    herwigGen.Draw('SAME H')
    herwigGenRatio=herwigGen.Clone('herwigGenRatio')
    herwigGenRatio.Divide(dataUnfolded)
    
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
    
    dataUnfoldedSysRatio.SetMarkerStyle(0)
    dataUnfoldedSysRatio.Draw('e2')
    dataUnfoldedRatio.SetMarkerStyle(0)
    dataUnfoldedRatio.Draw('e2,same')
    line = dataUnfoldedRatio.Clone('line')
    line.SetFillStyle(0)
    line.Draw('hist,same')
    
    nominalGenRatio.Draw('SAME H')
    FSRUpGenRatio.Draw  ('SAME P X0 E1')
    FSRDownGenRatio.Draw('SAME P X0 E1')
    herwigGenRatio.Draw ('SAME H')
    
    c.Print(opt.outDir+'/'+opt.obs+'_charged_'+opt.flavor+'_result.pdf')
    c.Print(opt.outDir+'/'+opt.obs+'_charged_'+opt.flavor+'_result.png')
    #c.Print(opt.outDir+'/'+opt.obs+'_charged_'+opt.flavor+'_result.root')
    rootoutfile.Write()
    
    #Folding test plot
    
    c = ROOT.TCanvas('c','c',500,500)
    c.cd()
    c.SetRightMargin(0.05)
    c.SetLeftMargin(0.12)
    c.SetTopMargin(0.1)
    c.SetTopMargin(0.05)

    dataUnfolded.Draw()
    nominalGen.Draw('same')
    
    nominalReco = normalizeAndDivideByBinWidth(nominal.ProjectionY("nominalReco"))
    nominalReco.SetLineColor(ROOT.kRed+1)
    nominalReco.SetLineStyle(2)
    nominalReco.Draw("SAME,HIST")
    
    normalizeAndDivideByBinWidth(dataBkgSub)
    dataBkgSub.SetLineColor(ROOT.kBlack)
    dataBkgSub.SetLineStyle(0)
    dataBkgSub.Draw("SAME,HIST")
    
    dataFoldedBack.SetMarkerColor(ROOT.kMagenta+1);
    dataFoldedBack.SetMarkerStyle(2);
    dataFoldedBack.Draw("SAME P X0 E1");
    
    inix = 0.5
    if (nominalGen.GetMaximumBin() > nominalGen.GetNbinsX()/2.): inix = 0.15
    
    legend = ROOT.TLegend(inix,0.6,inix+0.35,0.9)
    legend.SetLineWidth(0)
    legend.SetFillStyle(0)
    legend.AddEntry(nominalGen, "MC gen t#bar{t} (Pythia8)", "p")
    legend.AddEntry(nominalReco, "MC reco t#bar{t}", "l")
    legend.AddEntry(dataBkgSub, "data (bg-sub)", "l")
    legend.AddEntry(dataUnfolded, "data unfolded", "ep")
    legend.AddEntry(dataFoldedBack, "data folded back", "p")
    legend.Draw()
    
    c.Print(opt.outDir+'/'+opt.obs+'_charged_'+opt.flavor+'_test.pdf')
    c.Print(opt.outDir+'/'+opt.obs+'_charged_'+opt.flavor+'_test.png')
    
    #Uncertainty plot
    c = ROOT.TCanvas('c','c',500,500)
    c.cd()
    c.SetRightMargin(0.05)
    c.SetLeftMargin(0.12)
    c.SetTopMargin(0.1)
    c.SetTopMargin(0.05)
    
    legend = ROOT.TLegend(0.15,0.65,0.85,0.875)
    legend.SetNColumns(2)
    legend.SetLineWidth(0)
    legend.SetFillStyle(0)
    dataUnfoldedSysRatio.GetYaxis().SetNdivisions()
    dataUnfoldedSysRatio.GetXaxis().SetTitleSize(0.045)
    dataUnfoldedSysRatio.GetXaxis().SetTitleOffset(1.)
    dataUnfoldedSysRatio.GetXaxis().SetLabelSize(0.04)
    dataUnfoldedSysRatio.GetYaxis().SetTitleSize(0.045)
    dataUnfoldedSysRatio.GetYaxis().SetTitleOffset(1.2)
    dataUnfoldedSysRatio.GetYaxis().SetLabelSize(0.04)
    dataUnfoldedSysRatio.GetYaxis().SetRangeUser(0.8,1.5)
    dataUnfoldedSysRatio.Draw('e2')
    legend.AddEntry(dataUnfoldedSysRatio, 'Total uncertainty', 'pf')
    dataUnfoldedRatio.SetFillColor(ROOT.kGray+1)
    dataUnfoldedRatio.Draw('e2,same')
    legend.AddEntry(dataUnfoldedRatio, 'Statistical uncertainty', 'pf')
    colors = [ROOT.kRed+1, ROOT.kGreen+1, ROOT.kBlue+1, ROOT.kCyan+1, ROOT.kMagenta+1]
    color = 0
    uncmap = {'fsrup':'FSR up', 'fsrdn':'FSR down', 'herwig':'Herwig++', 'tracking_up':'Tracking up', 'tracking_down':'Tracking down'}
    for tag,systematic in systematicUnfolded.iteritems():
        if not tag in ['MC13TeV_TTJets_fsrup', 'MC13TeV_TTJets_fsrdn', 'MC13TeV_TTJets_herwig', 'MC13TeV_TTJets_tracking_up', 'MC13TeV_TTJets_tracking_down']: continue
        #if any(x in tag for x in ['wgt', 'jec', 'btag', 'csv', 'jer']): continue
        #if 'weight' in tag: continue
        systematic.Divide(dataUnfolded)
        systematic.SetLineColor(colors[color])
        systematic.SetMarkerColor(colors[color])
        systematic.SetMarkerStyle(21+color)
        color += 1
        legend.AddEntry(systematic, uncmap.get(tag.replace('MC13TeV_TTJets_', ''), tag), 'pl')
        systematic.Draw('same')
    legend.Draw()
    txt=ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(42)
    txt.SetTextSize(0.041)
    txt.SetTextAlign(12)
    txt.DrawLatex(0.16,0.91,cmsLabel)
    txt.DrawLatex(0.7,0.97,'#scale[0.8]{%3.1f fb^{-1} (%s)}' % (opt.lumi/1000.,opt.com) )
            
    c.Print(opt.outDir+'/'+opt.obs+'_charged_'+opt.flavor+'_syst.pdf')
    c.Print(opt.outDir+'/'+opt.obs+'_charged_'+opt.flavor+'_syst.png')
    c.Print(opt.outDir+'/'+opt.obs+'_charged_'+opt.flavor+'_syst.root')
    
    # Covariance from TUnfold
    c.SetRightMargin(0.15)
    c.SetLeftMargin(0.12)
    c.SetTopMargin(0.1)
    c.SetTopMargin(0.05)
    
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
    
    dataStatCov.SetTitle('')
    dataStatCov.GetZaxis().SetRangeUser(-max(abs(dataStatCov.GetMaximum()), abs(dataStatCov.GetMinimum())),
                                         max(abs(dataStatCov.GetMaximum()), abs(dataStatCov.GetMinimum())))
    dataStatCov.Draw('colz')
    
    matrix = []
    for i in range(1, dataStatCov.GetNbinsX()+1):
        column = []
        for j in range(1, dataStatCov.GetNbinsY()+1):
            column.append(dataStatCov.GetBinContent(i, j))
        matrix.append(column)
    X = numpy.array(matrix)
    print(X)
    print(numpy.linalg.det(X))
    X_inv = numpy.linalg.inv(X)
    print(X_inv)
    
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
    txt.DrawLatex(0.16,0.91,cmsLabel)
    txt.DrawLatex(0.63,0.97,'#scale[0.8]{%3.1f fb^{-1} (%s)}' % (opt.lumi/1000.,opt.com) )
    
    c.Print(opt.outDir+'/'+opt.obs+'_charged_'+opt.flavor+'_cov_tunfold.pdf')
    c.Print(opt.outDir+'/'+opt.obs+'_charged_'+opt.flavor+'_cov_tunfold.png')
    

def unfold(Mtag, Morig, backgrounds, data):
    # Background subtraction
    dataBkgSub = data.Clone(Mtag+'_dataBkgSub')
    
    M = Morig.Clone('M')

    # TODO: background uncertainties
    for tag,background in backgrounds.iteritems():
        dataBkgSub.Add(background, -1.)
    #if debug:
    #    print('dataBkgSub', dataBkgSub.GetEntries(), dataBkgSub.Integral())
    #    for i in range(dataBkgSub.GetNbinsX()+2):
    #        print('dataBkgSub bin content', dataBkgSub.GetBinContent(i), 'bin error', dataBkgSub.GetBinError(i))
    
    sigRecoNonGen = M.ProjectionY("sigRecoNonGen", 0, 0);
    sigReco       = M.ProjectionY("sigReco");
    
    for i in range(dataBkgSub.GetNbinsX()+2):
        sf = 1.
        if sigReco.GetBinContent(i) > 0:
          sf = sigRecoNonGen.GetBinContent(i) / sigReco.GetBinContent(i)
        dataBkgSub.SetBinContent(i, (1.-sf)*dataBkgSub.GetBinContent(i))
        dataBkgSub.SetBinError  (i, (1.-sf)*dataBkgSub.GetBinError(i))
        M.SetBinContent(0, i, 0.)
    #if debug:
    #    print('dataBkgSub', dataBkgSub.GetEntries(), dataBkgSub.Integral())
    #    for i in range(dataBkgSub.GetNbinsX()+2):
    #        print('dataBkgSub bin content', dataBkgSub.GetBinContent(i), 'bin error', dataBkgSub.GetBinError(i))
    
    # Do unfolding
    unfold = ROOT.TUnfoldDensity(M, ROOT.TUnfold.kHistMapOutputHoriz, ROOT.TUnfold.kRegModeCurvature, ROOT.TUnfold.kEConstraintArea, ROOT.TUnfoldDensity.kDensityModeBinWidth)
    
    if (unfold.SetInput(dataBkgSub) >= 10000):
        print('SetInput bad return value >= 10000')
        return
    
    unfold.DoUnfold(0.);
    
    # Retrieve results
    dataUnfolded = normalizeAndDivideByBinWidth(unfold.GetOutput(Mtag+"_Unfolded"))
    
    #if debug:
    #    for i in range(1, dataUnfolded.GetNbinsX()+1):
    #        print('dataUnfolded', i, dataUnfolded.GetBinContent(i), dataUnfolded.GetBinError(i))
    dataFoldedBack = normalizeAndDivideByBinWidth(unfold.GetFoldedOutput(Mtag+"_FoldedBack"))
    
    # get error matrix (input distribution [stat] errors only)
    dataStatCov = unfold.GetEmatrixTotal(Mtag+"_StatCov");
    dataUnfoldedOrig = unfold.GetOutput(Mtag+"_UnfoldedOrig")
    dataStatCov.Scale(1./dataUnfoldedOrig.Integral()**2)
    
    return dataUnfolded, dataFoldedBack, dataBkgSub, dataStatCov



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

