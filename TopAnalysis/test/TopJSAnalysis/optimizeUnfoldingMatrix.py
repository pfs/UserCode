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
    if inputfile == 'eos':
        tree.Add('/eos/user/m/mseidel/analysis/TopJetShapes/b312177/Chunks/MC13TeV_TTJets_0.root')
        tree.Add('/eos/user/m/mseidel/analysis/TopJetShapes/b312177/Chunks/MC13TeV_TTJets_1.root')
        #tree.Add('/eos/user/m/mseidel/analysis/TopJetShapes/b312177/Chunks/MC13TeV_TTJets_2.root')
        #tree.Add('/eos/user/m/mseidel/analysis/TopJetShapes/b312177/Chunks/MC13TeV_TTJets_3.root')
        #tree.Add('/eos/user/m/mseidel/analysis/TopJetShapes/b312177/Chunks/MC13TeV_TTJets_4.root')
        #tree.Add('/eos/user/m/mseidel/analysis/TopJetShapes/b312177/Chunks/MC13TeV_TTJets_5.root')
        #tree.Add('/eos/user/m/mseidel/analysis/TopJetShapes/b312177/Chunks/MC13TeV_TTJets_6.root')
        #tree.Add('/eos/user/m/mseidel/analysis/TopJetShapes/b312177/Chunks/MC13TeV_TTJets_7.root')
        #tree.Add('/eos/user/m/mseidel/analysis/TopJetShapes/b312177/Chunks/MC13TeV_TTJets_8.root')
        #tree.Add('/eos/user/m/mseidel/analysis/TopJetShapes/b312177/Chunks/MC13TeV_TTJets_9.root')
    else: tree.Add(inputfile)
    totalEntries = tree.GetEntries()
    
    print("Total entries: " + str(totalEntries))
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
    if (obs == "ptds"):
        label = "scaled p_{T}D"
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
    
    basename = obs + "_" + reco + "_"
    
    h = ROOT.TH2F(basename+"h", labels, nbins, lowbin, highbin, nbins, lowbin, highbin)
    
    points = getPointsFromTree(tree, obs, reco, ptcut)
    #fillHist(h, tree, obs, reco, ptcut)
    fillHistFromPoints(h, points)
    
    c = ROOT.TCanvas('c', 'c', 500, 450)
    c.cd()
    
    h.SetMinimum(-1e-10)
    h.Draw("colz")
    #h.Fit("pol1")
    
    plotformats = ['.png', '.pdf']
    plotbasename = output + "/" + basename
    for p in plotformats: c.Print(plotbasename + "h" + p)
    
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
    for p in plotformats: c.Print(plotbasename + "h1" + p)
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
    print("Number of gen bins = " + str(len(bins)-1))
    
    hnorm = ROOT.TH2F(basename+"hnorm", labels, nbins, lowbin, highbin, nbins, lowbin, highbin)
    for g in range(1, h.GetNbinsX()+1):
        for r in range(1, h.GetNbinsY()+1):
            if recosums[r-1] > 0:
                hnorm.SetBinContent(g, r, h.GetBinContent(g, r)/recosums[r-1])
            else: hnorm.SetBinContent(g, r, 0)
    
    hnorm.SetMinimum(-1e-10)
    #hnorm.SetMaximum(1)
    hnorm.Draw("colz")
    #hnorm.Fit("pol1")
    for p in plotformats: c.Print(plotbasename + "hnorm" + p)
    
    # reco bin splitting
    divisor = 2
    #if (obs == "mult"): minWidth = 1.0
    bins2 = []
    for i in range(len(bins)-1):
        for j in range(divisor):
            step = abs(bins[i]-bins[i+1])/divisor*j
            if (j != 0 and obs == 'mult' and step < 1): continue
            bins2.append(bins[i] + step)
    bins2.append(bins[-1])
    
    print(bins2)
    print("Number of reco bins = " + str(len(bins2)-1))
    
    # original reco bins
    bins3 = []
    for i in range(nbins+1):
        bins3.append(lowbin + (highbin-lowbin)/nbins*i)
    
    #print bins3
    
    binArray = array('d', bins)
    bin2Array = array('d', bins2) # put bins2 for reco bin split
    
    hopt = ROOT.TH2F(basename+"hopt", labels, len(binArray)-1, binArray, len(binArray)-1, binArray)
    #fillHist(hopt, tree, obs, reco, ptcut)
    fillHistFromPoints(hopt, points)
    
    hopt.SetMinimum(-1e-10)
    hopt.Draw("colz")
    for p in plotformats: c.Print(plotbasename + "hopt" + p)
    
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
    for p in plotformats: c.Print(plotbasename + "hopt1" + p)
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
    
    hoptpur  = ROOT.TH1F(basename+"hoptpur", label1, len(binArray)-1, binArray)
    hoptsta  = ROOT.TH1F(basename+"hoptsta", label1, len(binArray)-1, binArray)
    hoptnorm = ROOT.TH2F(basename+"hoptnorm", labels, len(binArray)-1, binArray, len(binArray)-1, binArray)
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
    for p in plotformats: c.Print(plotbasename + "hoptnorm" + p)
    
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
    for p in plotformats: c.Print(plotbasename + "hoptpursta" + p)
    
    responsematrix = ROOT.TH2F(basename+"responsematrix", labels, len(binArray)-1, binArray, len(bin2Array)-1, bin2Array)
    #fillHist(responsematrix, tree, obs, reco, ptcut)
    fillHistFromPoints(responsematrix, points)
    responsematrix.SetMinimum(-1e-10)
    responsematrix.Draw("colz")
    for p in plotformats: c.Print(plotbasename + "responsematrix" + p)
    
    responsematrix.ProjectionX()
    responsematrix.ProjectionY()
    
    rootoutput.Write();
    
    del h
    del hopt


def getPointsFromTree(tree, obs, reco, ptcut):
    points = []
    for event in tree:
        if event.gen_sel*event.reco_sel != 1: continue
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
                
                '''
                if not 'ptds' in obs:
                    valReco = eval('event.j_'+obs+'_'+reco)[j]
                    valGen  = eval('event.gj_'+obs+'_'+reco)[i]
                
                # TODO: remove special recipes with next ntuples
                # additional multiplicity cuts
                multReco = eval('event.j_mult_'+reco)[j]
                multGen  = eval('event.gj_mult_'+reco)[i]
                minMult = 0
                if 'c1' in obs: minMult = 3
                if 'c2' in obs: minMult = 4
                if 'c3' in obs: minMult = 5
                if 'tau' in obs: minMult = int(obs[3]) + 2
                if multReco < minMult: valReco = -1
                if multGen  < minMult: valGen  = -1
                if 'ptd' in obs and not 'ptds' in obs:
                    valReco = valReco**2
                    valGen  = valGen**2
                if 'ptds' in obs:
                    if multReco > 1:
                        valReco = eval('event.j_ptd_'+reco)[j]
                        valReco = (valReco**2 - 1./multReco) * multReco/(multReco-1)
                        if valReco > 0: valReco = sqrt(valReco)
                    else: valReco = -1
                    if multGen > 1:
                        valGen  = eval('event.gj_ptd_'+reco)[i]
                        valGen = (valGen**2 - 1./multGen) * multGen/(multGen-1)
                        if valGen > 0: valGen = sqrt(valGen)
                    else: valGen = -1
                if 'zgdr' in obs or 'zgxdr' in obs:
                    zcut = 0.1
                    if (eval('event.j_zg_'+reco)[j] < zcut): valReco = -1
                    if (eval('event.gj_zg_'+reco)[i] < zcut): valGen = -1
                '''
                #h.Fill(valGen, valReco)
                points.append([valGen, valReco])
    return points

def fillHistFromPoints(h, points):
    for point in points:
        h.Fill(point[0], point[1])

def fillHist(h, tree, obs, reco, ptcut):
    for event in tree:
        if event.gen_sel*event.reco_sel != 1: continue
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
                
                if not 'ptds' in obs:
                    valReco = eval('event.j_'+obs+'_'+reco)[j]
                    valGen  = eval('event.gj_'+obs+'_'+reco)[i]
                
                # TODO: remove special recipes with next ntuples
                # additional multiplicity cuts
                multReco = eval('event.j_mult_'+reco)[j]
                multGen  = eval('event.gj_mult_'+reco)[i]
                minMult = 0
                if 'c1' in obs: minMult = 3
                if 'c2' in obs: minMult = 4
                if 'c3' in obs: minMult = 5
                if 'tau' in obs: minMult = int(obs[3]) + 2
                if multReco < minMult: valReco = -1
                if multGen  < minMult: valGen  = -1
                if 'ptd' in obs and not 'ptds' in obs:
                    valReco = valReco**2
                    valGen  = valGen**2
                if 'ptds' in obs:
                    if multReco > 1:
                        valReco = eval('event.j_ptd_'+reco)[j]
                        valReco = (valReco**2 - 1./multReco) * multReco/(multReco-1)
                        if valReco > 0: valReco = sqrt(valReco)
                    else: valReco = -1
                    if multGen > 1:
                        valGen  = eval('event.gj_ptd_'+reco)[i]
                        valGen = (valGen**2 - 1./multGen) * multGen/(multGen-1)
                        if valGen > 0: valGen = sqrt(valGen)
                    else: valGen = -1
                if 'zgdr' in obs or 'zgxdr' in obs:
                    zcut = 0.1
                    if (eval('event.j_zg_'+reco)[j] < zcut): valReco = -1
                    if (eval('event.gj_zg_'+reco)[i] < zcut): valGen = -1
                
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
    plotformats = ['.png', '.pdf']
    plotbasename = output + "/" + obs + "_" + reco + "_"
    for p in plotformats: c2.Print(plotbasename + "slices_N" + p)
    slices[1].GetYaxis().SetRangeUser(h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax())
    slices[1].Draw()
    for p in plotformats: c2.Print(plotbasename + "slices_mean" + p)
    slices[2].GetYaxis().SetRangeUser(0., h.GetXaxis().GetXmax())
    slices[2].Draw()
    for p in plotformats: c2.Print(plotbasename + "slices_sigma" + p)
    
    bins  = [h.GetXaxis().GetXmin()]
    exact = [h.GetXaxis().GetXmin()]
    
    integral = h.Integral()
    print('integral', integral)
    binfraction = 0.
    for i in range(1, h.GetNbinsX()+1):
      binfraction += h.ProjectionX('px', i, i).Integral()/integral
      mean  = slices[1].GetBinContent(i)
      meanError = slices[1].GetBinError(i)
      if (mean == 0 or meanError/mean > 0.1): continue
      if (obs == 'tau21' and h.GetXaxis().GetBinCenter(i) < 0.1): continue
      if (obs == 'tau32' and h.GetXaxis().GetBinCenter(i) < 0.2): continue
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
          if (newbinedge > bins[-1] and binfraction > 0.01):
              bins.append(newbinedge)
              print('accepted bin', binfraction)
              binfraction = 0.
    
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
                            default='unfolding/optimize',
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
    parser.add_option('--ro', '--rootoutput',
                            dest='rootoutput',
                            default='output.root',
                            help='output root file [default: %default]')
    (opt, args) = parser.parse_args()

    os.system('mkdir -p %s' % opt.output)
    rootoutfile = ROOT.TFile(opt.output + "/" + opt.rootoutput, "UPDATE");
    
    if opt.obs == 'all': observables = ["mult", "width", "ptd", "ptds", "ecc", "tau21", "tau32", "tau43", "zg", "zgxdr", "zgdr", "ga_width", "ga_lha", "ga_thrust", "c1_02", "c1_05", "c1_10", "c1_20", "c2_02", "c2_05", "c2_10", "c2_20", "c3_02", "c3_05", "c3_10", "c3_20"]
    else: observables = opt.obs.split(',')

    if len(observables) == 1: optimize(opt.input, opt.output, opt.obs, opt.reco, opt.ptcut, rootoutfile)
    else:
        reco        = ["charged"]
        #reco        = ["charged", "puppi", "all"]
        #ptcuts      = {"0.5", "1.0", "1.5"}
        
        for o in observables:
            for r in reco:
                #for p in ptcuts:
                optimize(opt.input, opt.output, o, r, opt.ptcut, rootoutfile)
        

if __name__ == "__main__":
	sys.exit(main())
