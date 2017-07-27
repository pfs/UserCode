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
    parser.add_option('-i', '--inDir',       dest='inDir' ,      help='input directory',                default='unfolding/result',              type='string')
    parser.add_option('', '--inDirToys',       dest='inDirToys' ,      help='input toy directory',                default='unfolding/toys',              type='string')
    parser.add_option('-O', '--outDir',      dest='outDir' ,     help='output directory',                default='unfolding/covariance',              type='string')
    parser.add_option('-o', '--outName',     dest='outName' ,    help='name of the output file',        default='plotter.root',    type='string')
    parser.add_option(      '--silent',      dest='silent' ,     help='only dump to ROOT file',         default=False,             action='store_true')
    parser.add_option(      '--saveTeX',     dest='saveTeX' ,    help='save as tex file as well',       default=False,             action='store_true')
    parser.add_option('-l', '--lumi',        dest='lumi' ,       help='lumi [/pb]',              default=16551.,              type=float)
    parser.add_option('--obs', dest='obs',  default='mult', help='observable [default: %default]')
    parser.add_option('--flavor', dest='flavor',  default='all', help='flavor [default: %default]')
    (opt, args) = parser.parse_args()
    
    observables = ["mult", "width", "ptd", "ptds", "ecc", "tau21", "tau32", "tau43", "zg", "zgxdr", "zgdr", "ga_width", "ga_lha", "ga_thrust", "c1_02", "c1_05", "c1_10", "c1_20", "c2_02", "c2_05", "c2_10", "c2_20", "c3_02", "c3_05", "c3_10", "c3_20"]
    
    nice_observables_tex = {"mult": "N", "width": "width", "ptd": "$p_{T}D$", "ptds": "$p_{T}D^{s}$", "ecc": "$\\varepsilon$", "tau21": "$\\tau_{21}$", "tau32": "$\\tau_{32}$", "tau43": "$\\tau_{43}$", "zg": "$z_{g}$", "zgxdr": "$z_{g} \\times \\Delta R$", "zgdr": "$z_{g} \\Delta R$", "ga_width": "$\\lambda_{1}^{1}$ (width)", "ga_lha": "$\\lambda_{0.5}^{1}$ (LHA)", "ga_thrust": "$\\lambda_{2}^{1}$ (thrust)", "c1_02": "$C_{1}^{(0.2)}$", "c1_05": "$C_{1}^{(0.5)}$", "c1_10": "$C_{1}^{(1.0)}$", "c1_20": "$C_{1}^{(2.0)}$", "c2_02": "$C_{2}^{(0.2)}$", "c2_05": "$C_{2}^{(0.5)}$", "c2_10": "$C_{2}^{(1.0)}$", "c2_20":  "$C_{2}^{(2.0)}$", "c3_02": "$C_{3}^{(0.2)}$", "c3_05": "$C_{3}^{(0.5)}$", "c3_10": "$C_{3}^{(1.0)}$", "c3_20": "$C_{3}^{(2.0)}$"}
    
    observables_low = ["ptds", "ecc", "tau43", "zg", "zgdr"]
    
    flavors = ['all', 'bottom', 'light', 'gluon']

    sumNominal = 0.
    sumFSRUp = 0.
    sumFSRDown = 0.
    sumHerwig = 0.
    
    sumLowNominal = {}
    sumLowFSRUp   = {}
    sumLowFSRDown = {}
    sumLowHerwig  = {}
    sumLowNominal['total'] = 0.
    sumLowFSRUp  ['total'] = 0.
    sumLowFSRDown['total'] = 0.
    sumLowHerwig ['total'] = 0.
    for flavor in flavors:
        sumLowNominal[flavor] = 0.
        sumLowFSRUp  [flavor] = 0.
        sumLowFSRDown[flavor] = 0.
        sumLowHerwig [flavor] = 0.
    
    with open('%s/table.tex'%(opt.outDir), 'w') as tex:
        for obs in observables:
            for flavor in flavors:
        
                # statistical covariance
                
                toyfile = '%s/%s_charged_%s_toys.root'%(opt.inDirToys, obs, flavor)
                fInToy = ROOT.TFile.Open(toyfile)
                if not fInToy: continue
                
                pseudoresults = []
                counter = 0
                while True:
                    pseudoresult = []
                    h = fInToy.Get('Unfolded_' + str(counter))
                    if not h: break
                    for i in range(1, h.GetNbinsX()+1):
                        pseudoresult.append(h.GetBinContent(i)/h.GetBinWidth(i))
                    integral = sum(pseudoresult)
                    for i in range(len(pseudoresult)):
                        pseudoresult[i] = pseudoresult[i]/integral
                    pseudoresults.append(pseudoresult)
                    counter += 1
                
                print('Imported', counter, 'toy experiments')
                
                x = numpy.array(pseudoresults).T
                #print(x)
                statcov = numpy.cov(x)
                #print(statcov)
                statcov_reduced = numpy.delete(statcov, 3, 0)
                statcov_reduced = numpy.delete(statcov_reduced, 3, 1)
                #print(statcov_reduced)
                #print(numpy.linalg.det(statcov_reduced))
                #print(numpy.linalg.inv(statcov_reduced))
                
                #rootoutfile = ROOT.TFile.Open(opt.outDir+'/'+obs+'_charged_'+flavor+'_cov.root', 'RECREATE')
                #rootoutfile.cd()
                
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
                c = ROOT.TCanvas('c','c',500,500)
                c.SetRightMargin(0.15)
                c.SetLeftMargin(0.12)
                c.SetTopMargin(0.1)
                c.SetTopMargin(0.05)
                c.cd()
                
                h = fInToy.Get('Unfolded_0')
                dataStatCovNorm = ROOT.TH2D('dataStatCovNorm', '', h.GetNbinsX(), h.GetXaxis().GetXbins().GetArray(), h.GetNbinsX(), h.GetXaxis().GetXbins().GetArray())
                axistitle = h.GetXaxis().GetTitle().replace('generated ', '')
                dataStatCovNorm.SetXTitle(axistitle + ' (%s)'%(flavor))
                dataStatCovNorm.SetYTitle(axistitle + ' (%s)'%(flavor))
                
                for i in range(1, dataStatCovNorm.GetNbinsX()+1):
                    for j in range(1, dataStatCovNorm.GetNbinsY()+1):
                        dataStatCovNorm.SetBinContent(i, j, statcov[i-1][j-1])
                
                dataStatCovNorm.GetZaxis().SetRangeUser(-max(abs(statcov.min()), abs(statcov.max())),
                                                         max(abs(statcov.min()), abs(statcov.max())))
                dataStatCovNorm.Draw('colz')
                
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
                txt.DrawLatex(0.16,0.91, cmsLabel)
                txt.DrawLatex(0.63,0.97, '#scale[0.8]{%3.1f fb^{-1} (%s)}' % (opt.lumi/1000.,opt.com) )
                txt.DrawLatex(0.16,0.85, 'Statistical covariance')
                
                c.Print(opt.outDir+'/'+obs+'_charged_'+flavor+'_cov_stat.pdf')
                c.Print(opt.outDir+'/'+obs+'_charged_'+flavor+'_cov_stat.png')
                
                
                # Read lists of syst samples
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
                               'tracking',
                               'singletop',
                               'wjets'
                              ]
                varList = []
                for var in allSystVars:
                    varList.append([var+'_up', var+'_down'])
                varList += [['evtgen'],
                            ['m171v5', 'm173v5'],
                            ['herwig'],
                            ['isrup', 'isrdn'],
                            ['fsrup', 'fsrdn'],
                            ['hdampup', 'hdampdn'],
                            ['ueup', 'uedn'],
                            ['erdON'],
                            ['qcdBased'],
                            ['wgt1', 'wgt2'],
                            ['wgt3', 'wgt4'],
                            ['wgt5', 'wgt6'],
                            ['wgt7', 'wgt8'],
                            ['wgt9'],
                            ['wgt10', 'wgt11'],
                            ['wgt12'],
                            ['wgt13', 'wgt14'],
                            ['wgt15', 'wgt18'],
                            ['wgt16', 'wgt20'],
                           ]
                
                resultfile = '%s/%s_charged_%s_result.root'%(opt.inDir, obs, flavor)
                fIn=ROOT.TFile.Open(resultfile)
                
                # reference
                hnominal    = fIn.Get('MC13TeV_TTJets_Unfolded')
                nominal     = []
                for i in range(1, hnominal.GetNbinsX()+1):
                    nominal.append(hnominal.GetBinContent(i))
                
                systcov = copy.copy(statcov)
                systcov -= systcov
                
                for var in varList:
                    if len(var) == 1:
                        hsyst = fIn.Get('MC13TeV_TTJets_'+var[0]+'_Unfolded')
                        if not hsyst:
                            print(var, 'not loaded')  
                            continue
                        
                        syst = []
                        for i in range(1, hsyst.GetNbinsX()+1):
                            syst.append(hsyst.GetBinContent(i))
                        
                        x = numpy.array([nominal, syst]).T
                        cov = numpy.cov(x)
                        systcov += cov
                    if len(var) == 2:
                        hsyst_up = fIn.Get('MC13TeV_TTJets_'+var[0]+'_Unfolded')
                        hsyst_dn = fIn.Get('MC13TeV_TTJets_'+var[1]+'_Unfolded')
                        if not hsyst_up or not hsyst_dn:
                            print(var, 'not loaded')  
                            continue
                        
                        up = []
                        dn = []
                        maxdelta = []
                        for i in range(1, hsyst_up.GetNbinsX()+1):
                            up.append(hsyst_up.GetBinContent(i))
                            dn.append(hsyst_dn.GetBinContent(i))
                            maxdelta.append(max(abs(hsyst_up.GetBinContent(i) - nominal[i-1]), abs(hsyst_dn.GetBinContent(i) - nominal[i-1])))
                        
                        x = numpy.array([up, dn]).T
                        cov = numpy.cov(x)
                        
                        for i in range(len(maxdelta)):
                            for j in range(len(maxdelta)):
                                cov[i][j] = maxdelta[i] * maxdelta[j] * numpy.sign(cov[i][j])
                        
                        systcov += cov
                
                dataSystCovNorm = dataStatCovNorm.Clone('dataSystCovNorm')
                dataSystCovNorm.Reset()
                
                for i in range(1, dataSystCovNorm.GetNbinsX()+1):
                    for j in range(1, dataSystCovNorm.GetNbinsY()+1):
                        dataSystCovNorm.SetBinContent(i, j, systcov[i-1][j-1])
                
                dataSystCovNorm.GetZaxis().SetRangeUser(-max(abs(systcov.min()), abs(systcov.max())),
                                                         max(abs(systcov.min()), abs(systcov.max())))
                dataSystCovNorm.Draw('colz')
                
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
                txt.DrawLatex(0.16,0.91, cmsLabel)
                txt.DrawLatex(0.63,0.97, '#scale[0.8]{%3.1f fb^{-1} (%s)}' % (opt.lumi/1000.,opt.com) )
                txt.DrawLatex(0.16,0.85, 'Systematic covariance')
                
                c.Print(opt.outDir+'/'+obs+'_charged_'+flavor+'_cov_syst.pdf')
                c.Print(opt.outDir+'/'+obs+'_charged_'+flavor+'_cov_syst.png')
                
                dataCovNorm = dataStatCovNorm.Clone('dataCovNorm')
                dataCovNorm.Reset()
                
                cov = statcov + systcov
                
                for i in range(1, dataCovNorm.GetNbinsX()+1):
                    for j in range(1, dataCovNorm.GetNbinsY()+1):
                        dataCovNorm.SetBinContent(i, j, cov[i-1][j-1])
                
                dataCovNorm.GetZaxis().SetRangeUser(-max(abs(cov.min()), abs(cov.max())),
                                                     max(abs(cov.min()), abs(cov.max())))
                dataCovNorm.Draw('colz')
                
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
                txt.DrawLatex(0.16,0.91, cmsLabel)
                txt.DrawLatex(0.63,0.97, '#scale[0.8]{%3.1f fb^{-1} (%s)}' % (opt.lumi/1000.,opt.com) )
                txt.DrawLatex(0.16,0.85, 'Total covariance')
                
                c.Print(opt.outDir+'/'+obs+'_charged_'+flavor+'_cov.pdf')
                c.Print(opt.outDir+'/'+obs+'_charged_'+flavor+'_cov.png')
                
                #print(numpy.linalg.det(cov))
                cov_reduced = numpy.delete(cov, 0, 0)
                cov_reduced = numpy.delete(cov_reduced, 0, 1)
                #print(cov_reduced)
                #print(numpy.linalg.det(cov_reduced))
                #print(numpy.linalg.inv(cov_reduced))
                
                chi2Nominal = returnChi2(fIn, cov_reduced, nominal, 'nominalGen')
                chi2FSRUp   = returnChi2(fIn, cov_reduced, nominal, 'FSRUpGen')
                chi2FSRDown = returnChi2(fIn, cov_reduced, nominal, 'FSRDownGen')
                chi2Herwig  = returnChi2(fIn, cov_reduced, nominal, 'herwigGen')
                
                sumNominal += chi2Nominal
                sumFSRUp   += chi2FSRUp
                sumFSRDown += chi2FSRDown
                sumHerwig  += chi2Herwig
                
                if obs in observables_low:
                    sumLowNominal['total']  += chi2Nominal
                    sumLowFSRUp  ['total']  += chi2FSRUp
                    sumLowFSRDown['total']  += chi2FSRDown
                    sumLowHerwig ['total']  += chi2Herwig
                    sumLowNominal[flavor] += chi2Nominal
                    sumLowFSRUp  [flavor] += chi2FSRUp
                    sumLowFSRDown[flavor] += chi2FSRDown
                    sumLowHerwig [flavor] += chi2Herwig
                
                tex.write('%s & %s & %.1f & %.1f & %.1f & %.1f \\\\\n'%(nice_observables_tex[obs], flavor, chi2Nominal, chi2FSRUp, chi2FSRDown, chi2Herwig))
        tex.write('\\hline\nTotal &  & %.1f & %.1f & %.1f & %.1f \\\\\n'%(sumNominal, sumFSRUp, sumFSRDown, sumHerwig))
        tex.write('\\hline\nTotal low corr. &  & %.1f & %.1f & %.1f & %.1f \\\\\n'%(sumLowNominal['total'], sumLowFSRUp['total'], sumLowFSRDown['total'], sumLowHerwig['total']))
        for flavor in flavors:
            tex.write('-- %s &  & %.1f & %.1f & %.1f & %.1f \\\\\n'%(flavor, sumLowNominal[flavor], sumLowFSRUp[flavor], sumLowFSRDown[flavor], sumLowHerwig[flavor]))
            probLowNominal = ROOT.TMath.Prob(sumLowNominal[flavor], len(observables_low))
            probLowFSRUp   = ROOT.TMath.Prob(sumLowFSRUp[flavor], len(observables_low))
            probLowFSRDown = ROOT.TMath.Prob(sumLowFSRDown[flavor], len(observables_low))
            probLowHerwig  = ROOT.TMath.Prob(sumLowHerwig[flavor], len(observables_low))
            tex.write('$P(\\chi^{2})$ &  & %.3f & %.3f & %.3f & %.3f \\\\\n'%(probLowNominal, probLowFSRUp, probLowFSRDown, probLowHerwig))
    
def returnChi2(fIn, cov_reduced, nominal, prediction):
    hpred = fIn.Get(prediction)
    pred  = []
    for i in range(1, hpred.GetNbinsX()+1):
        pred.append(hpred.GetBinContent(i))
    diff = []
    for i in range(len(nominal)):
        diff.append(pred[i] - nominal[i])
    
    chi2 = numpy.array(diff[1:]).T.dot(numpy.linalg.inv(cov_reduced).dot(numpy.array(diff[1:])))
    ndf  = hpred.GetNbinsX()-1
    prob = ROOT.TMath.Prob(chi2, ndf)
    return chi2/ndf

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

