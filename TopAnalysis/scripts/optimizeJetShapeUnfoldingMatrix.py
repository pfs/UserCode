#! /usr/bin/env python
import ROOT
import optparse
import json
import sys
import os
import time
from array import *
from math import sqrt

"""
Analysis loop
"""
def optimize(inputfile, output, obs, reco, ptcut, rootoutput):
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadTopMargin(0.07);
    ROOT.gStyle.SetPadBottomMargin(0.13);
    ROOT.gStyle.SetPadLeftMargin(0.16)
    ROOT.gStyle.SetPadRightMargin(0.15)
    ROOT.gStyle.SetTitleSize(0.03, "t")
    ROOT.gStyle.SetTitleXOffset(1.5);
    ROOT.gStyle.SetTitleYOffset(1.75);
    ROOT.gStyle.SetPaintTextFormat("4.2f")
    
    # Viridis palette reversed + white
    stops = array('d', [0.0, 0.05, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0000])
    red   = array('d', [26./255., 51./255.,  43./255.,  33./255.,  28./255.,  35./255.,  74./255., 144./255., 246./255., 1., 1.])
    green = array('d', [9./255., 24./255.,  55./255.,  87./255., 118./255., 150./255., 180./255., 200./255., 222./255., 1., 1.])
    blue  = array('d', [30./255., 96./255., 112./255., 114./255., 112./255., 101./255.,  72./255.,  35./255.,   0./255., 1., 1.])
    ROOT.TColor.CreateGradientColorTable(11, stops, red[::-1], green[::-1], blue[::-1], 255)
    #ROOT.gStyle.SetPalette(53)
    
    tree = ROOT.TChain('tjsev')
    tree.AddFile(inputfile)
    totalEntries = tree.GetEntries()
    
    print("\nObservable: " + obs)
    print("Reco method: " + reco)
    
    #initial histogram
    # jet pt
    nbins   = 50
    lowbin  = 0
    highbin = 1
    sigmaFactor = 0.5
    
    if (obs == "mult"):
        label = "N"
        nbins   = 30
        highbin = 30
    if (obs == "width"):
        label = "width"
        highbin = 0.25
    if (obs == "ptd"):
        label = "p_{T}D"
    if (obs == "ecc"):
        label = "eccentricity"
        sigmaFactor = 0.25
    if (obs == "tau21"):
        label = "#tau_{21}"
        #nbins   = 60
        #lowbin  = 0.2
        #highbin = 0.8
    if (obs == "tau32"):
        label = "#tau_{32}"
    if (obs == "tau43"):
        label = "#tau_{43}"
    if (obs == "zg"):
        label = "z_{g}"
        lowbin = 0.1
        highbin = 0.5
        nbins = 40
        sigmaFactor = 0.25
    if (obs == "zgdr"):
        label = "z_{g} #DeltaR"
        highbin = 0.5
        sigmaFactor = 0.25
    if (obs == "zgxdr"):
        label = "z_{g} #times #DeltaR"
        highbin = 0.25
    if (obs == "ga_width"):
        label = "#lambda_{ 1}^{1} (width)"
    if (obs == "ga_thrust"):
        label = "#lambda_{ 2}^{1} (thrust)"
        highbin = 0.5
    if (obs == "ga_lha"):
        label = "#lambda_{ 0.5}^{1} (LHA)"
    if (obs == "c1_02"):
        label = "C_{ 1}^{ (0.2)}"
        highbin = 0.5
    if (obs == "c1_05"):
        label = "C_{ 1}^{ (0.5)}"
        highbin = 0.3
        nbins = 60
    if (obs == "c1_10"):
        label = "C_{ 1}^{ (1.0)}"
        highbin = 0.2
        nbins = 40
    if (obs == "c1_20"):
        label = "C_{ 1}^{ (2.0)}"
        highbin = 0.1
    if (obs == "c2_02"):
        label = "C_{ 2}^{ (0.2)}"
        highbin = 0.7
        nbins = 35
    if (obs == "c2_05"):
        label = "C_{ 2}^{ (0.5)}"
        highbin = 0.4
        nbins = 40
    if (obs == "c2_10"):
        label = "C_{ 2}^{ (1.0)}"
        highbin = 0.25
        nbins = 50
    if (obs == "c2_20"):
        label = "C_{ 2}^{ (2.0)}"
        highbin = 0.15
        nbins = 30
    if (obs == "c3_02"):
        label = "C_{ 3}^{ (0.2)}"
        highbin = 0.7
        nbins = 35
    if (obs == "c3_05"):
        label = "C_{ 3}^{ (0.5)}"
        highbin = 0.4
        nbins = 40
    if (obs == "c3_10"):
        label = "C_{ 3}^{ (1.0)}"
        highbin = 0.25
    if (obs == "c3_20"):
        label = "C_{ 3}^{ (2.0)}"
        highbin = 0.15
        nbins = 30
    
    label1 = reco + " particles, p_{T}>" + ptcut + " GeV;" + label
    labels = reco + " particles, p_{T}>" + ptcut + " GeV;generated  " + label + ";reconstructed  " + label
    
    print("Starting with bin width " + str((highbin-lowbin)/float(nbins)))
    
    f = ROOT.TFile(output + "/" + obs + "_" + reco +"_" + ptcut +".root", "RECREATE");
    
    h = ROOT.TH2F("h", labels, nbins, lowbin, highbin, nbins, lowbin, highbin)
    
    fillHist(h, tree, obs, reco, ptcut)
    
    c = ROOT.TCanvas('c', 'c', 500, 450)
    c.cd()
    
    h.SetMinimum(-1e-10)
    h.Draw("colz")
    #h.Fit("pol1")
    c.Print(output + "/" + obs + "_" + reco +"_" + ptcut +"_h.eps")
    
    h_gen  = h.ProjectionX()
    h_reco = h.ProjectionY()
    #print(h_gen.GetMaximum(), h_reco.GetMaximum())
    h_gen.GetYaxis().SetRangeUser(0., 1.5*h_reco.GetMaximum())
    h_gen.SetTitle(label1)
    h_gen.SetLineColor(ROOT.kRed+1)
    h_gen.Draw()
    h_reco.SetLineStyle(7)
    h_reco.Draw("same")
    leg = ROOT.TLegend(0.5,0.7,0.85,0.9)
    leg.SetLineWidth(0)
    leg.SetFillStyle(0)
    leg.AddEntry(h_gen, "generated", "l")
    leg.AddEntry(h_reco, "reconstructed", "l")
    leg.Draw()
    c.Print(output + "/" + obs + "_" + reco +"_" + ptcut +"_h1.eps")
    del h_gen
    del h_reco
    
    indices  = []
    diagonal = []
    gensums  = []
    recosums = []
    
    for g in range(1, h.GetNbinsX()+1):
        gensum = 0
        indices.append(g)
        diagonal.append(h.GetBinContent(g, g))
        for r in range(1, h.GetNbinsY()+1):
            gensum += h.GetBinContent(g, r)
        gensums.append(gensum)
    
    for r in range(1, h.GetNbinsY()+1):
        recosum = 0
        for g in range(0, h.GetNbinsX()+1):
            recosum += h.GetBinContent(g, r)
        recosums.append(recosum)

    #print diagonal
    #print gensums
    #print recosums
    
    purities    = []
    stabilities = []
    
    for i in range(len(gensums)):
        purity    = -1
        stability = -1
        if recosums[i] > 0:
            purity = diagonal[i]/recosums[i]
        purities.append(purity)
        if gensums[i] > 0:
            stability = diagonal[i]/gensums[i]
        stabilities.append(stability)
        
    #print purities
    #print stabilities
    
    #splitForEqualPurity(h, gensums, recosums, indices)
    #bins = splitForMinPurity(h, gensums, recosums, indices)
    
    bins = splitForMinSigma(h, output, obs, reco, ptcut, sigmaFactor)
    
    c = ROOT.TCanvas('c', 'c', 500, 450)
    c.cd()
    
    print(bins)
    print("Number of bins = " + str(len(bins)-1))
    
    hnorm = ROOT.TH2F("hnorm", labels, nbins, lowbin, highbin, nbins, lowbin, highbin)
    for g in range(1, h.GetNbinsX()+1):
        for r in range(1, h.GetNbinsY()+1):
            if recosums[r-1] > 0:
                hnorm.SetBinContent(g, r, h.GetBinContent(g, r)/recosums[r-1])
            else: hnorm.SetBinContent(g, r, 0)
    
    hnorm.SetMinimum(-1e-10)
    #hnorm.SetMaximum(1)
    hnorm.Draw("colz")
    #hnorm.Fit("pol1")
    c.Print(output + "/" + obs + "_" + reco +"_" + ptcut +"_hnorm.eps")
    
    # reco bin splitting
    divisor = 2
    minWidth = 0.
    if (obs == "mult"): minWidth = 1.0
    bins2 = []
    for i in range(len(bins)-1):
        if (abs(bins[i]-bins[i+1]) >= minWidth*divisor):
          for j in range(divisor):
              bins2.append(bins[i] + abs(bins[i]-bins[i+1])/divisor*j)
    bins2.append(bins[-1])
    
    print(bins2)
    
    # original reco bins
    bins3 = []
    for i in range(nbins+1):
        bins3.append(lowbin + (highbin-lowbin)/nbins*i)
    
    #print bins3
    
    binArray = array('d', bins)
    bin2Array = array('d', bins2) # put bins2 for reco bin split
    
    hopt = ROOT.TH2F("hopt", labels, len(binArray)-1, binArray, len(binArray)-1, binArray)
    fillHist(hopt, tree, obs, reco, ptcut)
    
    hopt.SetMinimum(-1e-10)
    hopt.Draw("colz")
    c.Print(output + "/" + obs + "_" + reco +"_" + ptcut +"_hopt.eps")
    
    h_gen  = hopt.ProjectionX()
    h_reco = hopt.ProjectionY()
    for i in range(1, h_gen.GetNbinsX()+1):
        h_gen.SetBinContent(i, h_gen.GetBinContent(i)/h_gen.GetBinWidth(i))
        h_reco.SetBinContent(i, h_reco.GetBinContent(i)/h_reco.GetBinWidth(i))
    h_gen.GetYaxis().SetRangeUser(0., 1.5*h_reco.GetMaximum())
    h_gen.GetXaxis().SetTitle(label)
    h_gen.SetLineColor(ROOT.kRed+1)
    h_gen.Draw()
    h_reco.SetLineStyle(7)
    h_reco.Draw("same")
    leg = ROOT.TLegend(0.5,0.7,0.85,0.9)
    leg.SetLineWidth(0)
    leg.SetFillStyle(0)
    leg.AddEntry(h_gen, "generated", "l")
    leg.AddEntry(h_reco, "reconstructed", "l")
    leg.Draw()
    c.Print(output + "/" + obs + "_" + reco +"_" + ptcut +"_hopt1.eps")
    del h_gen
    del h_reco
    
    optgensums = []
    for g in range(1, hopt.GetNbinsX()+1):
        gensum = 0
        for r in range(1, hopt.GetNbinsY()+1):
            gensum += hopt.GetBinContent(g, r)
        optgensums.append(gensum)
    
    optrecosums = []
    for r in range(1, hopt.GetNbinsY()+1):
        recosum = 0
        for g in range(1, hopt.GetNbinsX()+1):
            recosum += hopt.GetBinContent(g, r)
        optrecosums.append(recosum)
    
    hoptpur  = ROOT.TH1F("hoptpur", label1, len(binArray)-1, binArray)
    hoptsta  = ROOT.TH1F("hoptsta", label1, len(binArray)-1, binArray)
    hoptnorm = ROOT.TH2F("hoptnorm", labels, len(binArray)-1, binArray, len(binArray)-1, binArray)
    for g in range(1, hopt.GetNbinsX()+1):
        for r in range(1, hopt.GetNbinsY()+1):
            purity    = 0
            if optrecosums[r-1] > 0:
              purity    = hopt.GetBinContent(g, r)/optrecosums[r-1]
            hoptnorm.SetBinContent(g, r, purity)
            if g == r:
              hoptpur.SetBinContent(g, purity)
              stability = 0
              if optgensums[g-1] > 0:
                stability = hopt.GetBinContent(g, r)/optgensums[g-1]
              hoptsta.SetBinContent(g, stability)
    
    hoptnorm.SetMinimum(-1e-10)
    #hoptnorm.SetMaximum(1)
    hoptnorm.SetMarkerColor(ROOT.kWhite)
    hoptnorm.SetMarkerSize(1.5)
    hoptnorm.Draw("colz,text")
    c.Print(output + "/" + obs + "_" + reco +"_" + ptcut +"_hoptnorm.eps")
    
    hoptpur.GetYaxis().SetRangeUser(0., 1.)
    hoptpur.SetLineColor(ROOT.kRed+1)
    hoptpur.Draw()
    hoptsta.SetLineStyle(7)
    hoptsta.Draw("same")
    leg = ROOT.TLegend(0.5,0.75,0.85,0.9)
    leg.SetLineWidth(0)
    leg.SetFillStyle(0)
    leg.AddEntry(hoptpur, "purity", "l")
    leg.AddEntry(hoptsta, "stability", "l")
    leg.Draw()
    c.Print(output + "/" + obs + "_" + reco +"_" + ptcut +"_hoptpursta.eps")
    
    responsematrix = ROOT.TH2F("responsematrix", labels, len(binArray)-1, binArray, len(bin2Array)-1, bin2Array)
    fillHist(responsematrix, tree, obs, reco, ptcut)
    responsematrix.ProjectionX()
    responsematrix.ProjectionY()
    
    f.Write();
    
    del h
    del hopt


def fillHist(h, tree, obs, reco, ptcut):
    for event in tree:
        if event.gen_sel*event.reco_sel == -1: continue
        for j in range(event.nj):
            if event.j_gj[j] >= 0:
                #TODO jet criteria: no overlap, eta<2
                
                #j_p4 = ROOT.TLorentzVector()
                #j_p4.SetPtEtaPhiM(event.j_pt[j], event.j_eta[j], event.j_phi[j], event.j_m[j])
                #
                #g = event.j_gj[j]
                #gj_p4 = ROOT.TLorentzVector()
                #gj_p4.SetPtEtaPhiM(event.gj_pt[g], event.gj_eta[g], event.gj_phi[g], event.gj_m[g])
                #
                ##h.Fill(event.gj_pt[event.j_gj[j]], event.j_pt[j]) # pt
                ##h.Fill(event.gj_m[event.j_gj[j]], event.j_m[j]) # m
                ##h.Fill(gj_p4.M(), j_p4.M()) # m
                #h.Fill(gj_p4.M()/gj_p4.E(), j_p4.M()/j_p4.E()) # m/E
                
                #ga = []
                #for i in range(len(event.gj_ga)): ga.append(event.gj_ga[i])
                #print ga
                
                #o = 0 #mult
                #if (obs == "width"): o = 36
                #if (obs == "ptd"):   o = 18
                #if (obs == "mass"):  o = 63
                #
                #c = 0 #charged
                #if (reco == "all")   : c = 1
                #if (reco == "puppi") : c = 2
                #
                #p = 0 #pt>0.5
                #if (ptcut == "1.0") : p = 3
                #if (ptcut == "1.5") : p = 6
                #
                #i = event.j_gj[j]
                #if (obs == "ptd"): h.Fill(sqrt(event.gj_ga[i*81+o+c+p]), sqrt(event.j_ga[j*81+o+c+p]))
                #else :             h.Fill(event.gj_ga[i*81+o+c+p], event.j_ga[j*81+o+c+p])
                
                i = event.j_gj[j]
                valReco = eval('event.j_'+obs+'_'+reco)[j]
                valGen  = eval('event.gj_'+obs+'_'+reco)[i]
                h.Fill(valGen, valReco)



def splitForMinPurity(h, gensums, recosums, indices):
    for i in range(len(gensums)):
        if len(gensums) <= 1: break
        if sum(gensums[:i]) == 0: continue
        if sum(gensums[i:]) == 0: continue
        if len(recosums) <= 1: break
        if sum(recosums[:i]) == 0: continue
        if sum(recosums[i:]) == 0: continue
        sumdiagonal1 = 0
        for j in range(i):
            for k in range(i):
                sumdiagonal1 += h.GetBinContent(indices[j], indices[k])
        sumpurity1 = sumdiagonal1/sum(recosums[:i])
        sumstability1 = sumdiagonal1/sum(gensums[:i])
        
        threshold = 0.5
        
        if sumpurity1 > threshold and sumstability1 > threshold:
            
            print(h.GetXaxis().GetBinLowEdge(indices[0]), h.GetXaxis().GetBinUpEdge(indices[i-1]), sumpurity1, sumstability1)
            
            return [h.GetXaxis().GetBinLowEdge(indices[0])] + splitForMinPurity(h, gensums[i:], recosums[i:], indices[i:])
            
    return [h.GetXaxis().GetBinUpEdge(indices[-1])]
            
            #return True

def splitForMinSigma(h, output, obs, reco, ptcut, factor = 0.5):
    f1 = ROOT.TF1("f1", "gaus", h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax());
    f1.SetParameters(100, (h.GetXaxis().GetXmax()-h.GetXaxis().GetXmin())/2., 0.01);
    slices = ROOT.TObjArray()
    h.FitSlicesX(f1, 1, h.GetNbinsX(), 0, "QNRLM", slices)
    
    c2 = ROOT.TCanvas('c2', 'c2', 500, 450)
    c2.cd()
    
    slices[0].Draw()
    c2.Print(output + "/" + obs + "_" + reco +"_" + ptcut +"_slices_N.eps")
    slices[1].GetYaxis().SetRangeUser(h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax())
    slices[1].Draw()
    c2.Print(output + "/" + obs + "_" + reco +"_" + ptcut +"_slices_mean.eps")
    slices[2].GetYaxis().SetRangeUser(0., h.GetXaxis().GetXmax())
    slices[2].Draw()
    c2.Print(output + "/" + obs + "_" + reco +"_" + ptcut +"_slices_sigma.eps")
    
    bins  = [h.GetXaxis().GetXmin()]
    exact = [h.GetXaxis().GetXmin()]
    
    for i in range(1, h.GetNbinsX()+1):
      mean  = slices[1].GetBinContent(i)
      sigma = slices[2].GetBinContent(i) * factor
      #print(bins)
      #print(exact)
      #print(mean)
      #print(sigma)
      if (mean - sigma) > exact[-1]:
        #if (mean + sigma > h.GetXaxis().GetXmax()):
        #  bins.append(h.GetXaxis().GetXmax())
        #  break
        if (mean + sigma < h.GetXaxis().GetXmax()):
          exact.append(mean+sigma)
          newbinedge = h.GetXaxis().GetBinUpEdge((h.GetXaxis().FindBin(mean+sigma)))
          if (newbinedge > bins[-1]): bins.append(newbinedge)
    
    bins[-1] = h.GetXaxis().GetXmax()
    
    return bins


def splitForEqualPurity(h, gensums, recosums, indices):
    for i in range(len(gensums)):
        if len(gensums) <= 1: break
        if sum(gensums[:i]) == 0: continue
        if sum(gensums[i:]) == 0: continue
        if len(recosums) <= 1: break
        if sum(recosums[:i]) == 0: continue
        if sum(recosums[i:]) == 0: continue
        sumdiagonal1 = 0
        for j in range(i):
            for k in range(i):
                sumdiagonal1 += h.GetBinContent(indices[j], indices[k])
        sumpurity1 = sumdiagonal1/sum(gensums[:i])
        sumstability1 = sumdiagonal1/sum(recosums[:i])
        sumdiagonal2 = 0
        for j in range(i,len(gensums)):
            for k in range(i,len(gensums)):
                sumdiagonal2 += h.GetBinContent(indices[j], indices[k])
        sumpurity2 = sumdiagonal2/sum(gensums[i:])
        sumstability2 = sumdiagonal2/sum(recosums[i:])
        #print sumpurity1
        #print sumpurity2
        
        threshold = 0.4
        
        if sumpurity1 > sumpurity2:
            if sumpurity1 < threshold or sumpurity2 < threshold or sumstability1 < threshold or sumstability2 < threshold: return False
            
            success1 = splitForPurity(h, gensums[:i], recosums[:i], indices[:i]) 
            success2 = splitForPurity(h, gensums[i:], recosums[i:], indices[i:])
            
            if not success1:
                print(h.GetXaxis().GetBinLowEdge(indices[0]), h.GetXaxis().GetBinUpEdge(indices[i-1]), sumpurity1, sumstability1)
            if not success2:
                print(h.GetXaxis().GetBinLowEdge(indices[i]), h.GetXaxis().GetBinUpEdge(indices[-1]), sumpurity2, sumstability2)
            
            return True
            
            #return indices[i], sumpurity1, sumstability1, sumpurity2, sumstability2
        
    
"""
steer
"""
def main():
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--input',
                            dest='input',   
                            default='analysis.root',
                            help='input file [default: %default]')
    parser.add_option('--obs',
                            dest='obs',   
                            default='mult',
                            help='observable [default: %default]')
    parser.add_option('-o', '--output',
                            dest='output', 
                            default='unfolding',
                            help='Output directory [default: %default]')
    parser.add_option('-r', '--reco',
                            dest='reco',
                            default='charged',
                            help='Use charged/puppi/all particles [default: %default]')
    parser.add_option('-p', '--ptcut',
                            dest='ptcut',
                            default='1.0',
                            help='Use particles with pt>0.5/1.0/1.5 GeV [default: %default]')
    parser.add_option('-a', '--all',
                            dest='all',
                            action="store_true",default=False,
                            help='Run all plots [default: %default]')
    parser.add_option('-O', '--rootoutput',
                            dest='rootoutput',
                            action="store_true",default=False,
                            help='Run all plots [default: %default]')
    (opt, args) = parser.parse_args()

    os.system('mkdir -p %s' % opt.output)

    if (not opt.all): optimize(opt.input, opt.output, opt.obs, opt.reco, opt.ptcut, opt.rootoutput)
    else:
        observables = ["mult", "width", "ptd", "ecc", "tau21", "tau32", "tau43", "zg", "zgxdr", "zgdr", "ga_width", "ga_lha", "ga_thrust", "c1_02", "c1_05", "c1_10", "c1_20", "c2_02", "c2_05", "c2_10", "c2_20", "c3_02", "c3_05", "c3_10", "c3_20"]
        reco        = ["charged", "puppi", "all"]
        #ptcuts      = {"0.5", "1.0", "1.5"}
        
        for o in observables:
            for r in reco:
                #for p in ptcuts:
                optimize(opt.input, opt.output, o, r, opt.ptcut, opt.rootoutput)
        

if __name__ == "__main__":
	sys.exit(main())
