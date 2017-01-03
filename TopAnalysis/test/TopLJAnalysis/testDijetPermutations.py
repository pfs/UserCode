#!/usr/bin/env python

import sys
import ROOT
import numpy as np

"""
associate different metrics to the dijet systems in an event
"""
def dijetsWithMetric(jets):

    #sort by decreasing pT
    sorted(jets, key=lambda p4 : p4.Pt(), reverse=True)
    
    #build all possible dijet pairs with different metrics for ranking
    dijets=[]
    for i in xrange(0,len(jets)):
        for j in xrange(i+1,len(jets)):
            alpha = [ 60./(jets[i].Pt()+jets[j].Pt()), 60./(jets[i]+jets[j]).Pt() ]
            for p in [-2, -1, 0, 1 ,2 ]:
                for metric in ['dr','deta','dphi']:
                    kti2p=ROOT.TMath.Power(jets[i].Pt(),2*p)
                    ktj2p=ROOT.TMath.Power(jets[j].Pt(),2*p)
                    deltaij=jets[i].DeltaR(jets[j]) if metric=='dr' else ROOT.TMath.Abs(jets[i].Eta()-jets[j].Eta())
                    if metric=='dphi' : deltaij=ROOT.TVector2.Phi_mpi_pi(jets[i].Phi()-jets[j].Phi())
                    dij=ROOT.TMath.Min(kti2p,ktj2p)*ROOT.TMath.Power(deltaij,2)
                    diB=kti2p
                    alpha.append(dij/diB)
            dijets.append( (jets[i]+jets[j], alpha) )

    return dijets

"""
summarize resolutions
"""
def summarizeRankingPerformance(h,window,relativeresol):

    resolH,offPeakFracH=[],[]
    for i in xrange(0,len(h)):
        title=h[i][0]
        nbinsx=h[i][1].GetNbinsX()
        resolH.append( ROOT.TH1F('resol%d'%i,'%s;Ranking mode;Resolution'%title,nbinsx,0,nbinsx) )
        if relativeresol : resolH[i].GetYaxis().SetTitle('Relative resolution')
        resolH[i].SetMarkerStyle(20+i)
        offPeakFracH.append( ROOT.TH1F('offpeakfrac%d'%i,'%s;Ranking mode;N(off-peak)/N(on-peak)'%title,nbinsx,0,nbinsx) )
        offPeakFracH[i].SetMarkerStyle(20+i)
        
        for xbin in xrange(0,nbinsx):
            py=h[i][1].ProjectionY('py',xbin+1,xbin+1)
            resolH[i].GetXaxis().SetBinLabel(xbin+1,h[i][1].GetXaxis().GetBinLabel(xbin+1))
            offPeakFracH[i].GetXaxis().SetBinLabel(xbin+1,h[i][1].GetXaxis().GetBinLabel(xbin+1))

            #determine maximum
            ycenBin=py.GetMaximumBin()
            ycen=py.GetXaxis().GetBinCenter(ycenBin)

            #count off-peak to on-peak fraction
            offPeak,onPeak=0.,0.
            offPeakUnc,onPeakUnc=0.,0.
            for ybin in xrange(0,py.GetNbinsX()):
                yval=py.GetXaxis().GetBinCenter(ybin+1)
                cts,ctsUnc=py.GetBinContent(ybin+1),py.GetBinError(ybin+1)
                if ROOT.TMath.Abs(yval-ycen)>window : 
                    offPeak += cts
                    offPeakUnc += ctsUnc**2
                else : 
                    onPeak += cts
                    onPeakUnc += ctsUnc
            frac=offPeak/onPeak
            fracUnc=ROOT.TMath.Sqrt(onPeakUnc*(offPeak**2)+offPeakUnc*(onPeak**2))/onPeak
            offPeakFracH[i].SetBinContent(xbin+1,frac)
            offPeakFracH[i].SetBinError(xbin+1,fracUnc)
            
            #resolution
            resol=py.GetRMS()
            if relativeresol : resol=resol/ycen
            resolUnc=resol*py.GetRMSError()/py.GetRMS()
            resolH[i].SetBinContent(xbin+1,resol)
            resolH[i].SetBinError(xbin+1,resolUnc)

            py.Delete()

    #
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    c=ROOT.TCanvas('c','c',800,500)
    c.SetLeftMargin(0)
    c.SetRightMargin(0)
    c.SetTopMargin(0)
    c.SetBottomMargin(0)
    p1 = ROOT.TPad('p1','p1',0.0,0.5,1.0,0.0)
    p1.SetRightMargin(0.05)
    p1.SetLeftMargin(0.12)
    p1.SetTopMargin(0.008)
    p1.SetBottomMargin(0.12)
    p1.Draw()
    p1.cd()
    for i in xrange(0,len(resolH)):
        drawOpt='e1' if i==0 else 'e1same'
        resolH[i].Draw(drawOpt)
        resolH[i].GetXaxis().SetLabelSize(0.06)
        resolH[i].GetXaxis().SetTitleSize(0.06)
        resolH[i].GetYaxis().SetLabelSize(0.06)
        resolH[i].GetYaxis().SetTitleSize(0.06)
    resolH[0].GetYaxis().SetRangeUser(resolH[0].GetMinimum()*0.9,
                                      resolH[-1].GetMaximum()*1.1)
    c.cd()
    p2 = ROOT.TPad('p2','p2',0.0,0.5,1.0,1.0)
    p2.SetBottomMargin(0.005)
    p2.SetRightMargin(0.05)
    p2.SetLeftMargin(0.12)
    p2.SetTopMargin(0.05)
    p2.Draw()
    p2.cd()
    leg=ROOT.TLegend(0.65,0.8,0.9,0.6)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.06)
    for i in xrange(0,len(resolH)):
        drawOpt='e1' if i==0 else 'e1same'
        offPeakFracH[i].Draw(drawOpt)
        offPeakFracH[i].GetXaxis().SetLabelSize(0.00)
        offPeakFracH[i].GetXaxis().SetTitleSize(0.00)
        offPeakFracH[i].GetYaxis().SetLabelSize(0.06)
        offPeakFracH[i].GetYaxis().SetTitleSize(0.06)
        leg.AddEntry(offPeakFracH[i],offPeakFracH[i].GetTitle(),'p')
    offPeakFracH[0].GetYaxis().SetRangeUser(offPeakFracH[0].GetMinimum()*0.9,
                                            offPeakFracH[-1].GetMaximum()*1.1)
    label = ROOT.TLatex()
    label.SetNDC()
    label.SetTextFont(42)
    label.SetTextSize(0.08)
    label.DrawLatex(0.65,0.85,'#bf{CMS} #it{simulation preliminary}')
    leg.Draw()
    c.cd()
    c.Modified()
    c.Update()
    for ext in ['png','pdf']:
        c.SaveAs('%s_resol.%s'%(h[-1][1].GetName(),ext))


"""
"""
def main():

    rankQuantiles=[100,80,70,60]
    histos={}
    for r in rankQuantiles:
        pf=""
        if r!=100 : pf='q%d'%r
        histos['ranked_mjj'+pf]      = ROOT.TH2F('ranked_mjj'+pf,      ';Ranking mode;Dijet invariant mass [GeV]', 17,0,17,50,0,200)
        histos['ranked_mw'+pf]       = ROOT.TH2F('ranked_mw'+pf,       ';Ranking mode;Closest W mass [GeV]',       17,0,17,50,0,200)
        histos['ranked_mjjresol'+pf] = ROOT.TH2F('ranked_mjjresol'+pf, ';Ranking mode;M(jj)-M(W) [GeV]',           17,0,17,50,-50,50)
        histos['ranked_drresol'+pf]  = ROOT.TH2F('ranked_drresol'+pf,  ';Ranking mode;#DeltaR(jj,W)',              17,0,17,50,0,1)
    
    for xbin,label in [(1,'#Sigmap_{T}'),   (2,'p_{T}(jj)'),
                       (3,'#DeltaR,p=-2'), (4,'#Delta#eta,p=-2'), (5,'#Delta#phi,p=-2'),
                       (6,'#DeltaR,p=-1'), (7,'#Delta#eta,p=1'),  (8,'#Delta#phi,p=1'),
                       (9,'#DeltaR,p=0'),  (10,'#Delta#eta,p=0'), (11,'#Delta#phi,p=0'),
                       (12,'#DeltaR,p=+1'), (13,'#Delta#eta,p=+1'), (14,'#Delta#phi,p=+1'),
                       (15,'#DeltaR,p=+2'),(16,'#Delta#eta,p=+2'), (17,'#Delta#phi,p=+2')]:
        for r in rankQuantiles:
            pf=""
            if r!=100 : pf='q%d'%r
            histos['ranked_mjj'+pf].GetXaxis().SetBinLabel(xbin,label)
            histos['ranked_mw'+pf].GetXaxis().SetBinLabel(xbin,label)
            histos['ranked_mjjresol'+pf].GetXaxis().SetBinLabel(xbin,label)
            histos['ranked_drresol'+pf].GetXaxis().SetBinLabel(xbin,label)
    for h in histos:
        histos[h].Sumw2()
        histos[h].SetDirectory(0)
        
    rankSummary={}
    for mode in xrange(0,17):
        rankSummary[mode]=([],[])

    #loop over the events
    data=ROOT.TChain('data')
    for url in sys.argv[1].split(','):
        data.AddFile(url)
    nentries=data.GetEntries()
    for i in xrange(0,nentries):
        data.GetEntry(i)

        #get reco and gen jets which are not b-tagged
        jets,genJets=[],[]
        matchedGenJetsIdx=[]
        for j in xrange(0,data.nj):
            if data.j_btag[j]: continue        
            jets.append(ROOT.TLorentzVector())
            jets[-1].SetPtEtaPhiM( data.j_pt[j], data.j_eta[j], data.j_phi[j], data.j_m[j])

            try:
                for k in xrange(0,data.ngj):
                    if k in matchedGenJetsIdx : continue
                    gjp4=ROOT.TLorentzVector()
                    gjp4.SetPtEtaPhiM( data.gj_pt[k], data.gj_eta[k], data.gj_phi[k], data.gj_m[k])
                    if gjp4.DeltaR(jets[-1])>0.4: continue
                    matchedGenJetsIdx.append(k)
                    genJets.append( ROOT.TLorentzVector(gjp4) )
                    break
            except:
                pass

        #get Ws
        Wbosons=[]
        try:
            for k in xrange(0,data.ngp):
                Wbosons.append(ROOT.TLorentzVector())
                Wbosons[-1].SetPtEtaPhiM( data.gp_pt[k], data.gp_eta[k], data.gp_phi[k], data.gp_m[k])
        except:
            pass

        #two jets are reco level
        if len(jets)<2 : continue
                    
        #fill mass and direction resolutions
        weight=data.w
        if len(genJets)>=2 and len(Wbosons)==2:
            genDijets=dijetsWithMetric(genJets)
            for mode in xrange(0,17):
                jj,rank=sorted(genDijets, key=lambda jjranks : jjranks[1][mode])[0]
                mjj=jj.M()

                dRjjw1,dRjjw2=jj.DeltaR(Wbosons[0]),jj.DeltaR(Wbosons[1])
                widx=0 if dRjjw1<dRjjw2 else 1
                mw=Wbosons[widx].M()
                drjj2w=jj.DeltaR(Wbosons[widx])
                histos['ranked_mjj'].Fill(mode,mjj,weight)
                histos['ranked_mw'].Fill(mode,mw,weight)
                histos['ranked_mjjresol'].Fill(mode,mjj-mw,weight)
                histos['ranked_drresol'].Fill(mode,drjj2w,weight)

                rankSummary[mode][0].append(rank[mode])
                rankSummary[mode][1].append((mjj,mw,drjj2w))

    #compute ranking 0.80,0.70 percentile and cut all above that value
    print rankQuantiles,'working points are listed below:'
    for mode in rankSummary:
        cutVal=np.percentile(rankSummary[mode][0],rankQuantiles, axis=0)
        print mode,histos['ranked_mjj'].GetXaxis().GetBinLabel(mode+1),cutVal
        for i in xrange(0,len(rankSummary[mode][0])):
            val=rankSummary[mode][0][i]
            mjj,mw,drjj2w = rankSummary[mode][1][i]
            for k in xrange(1,len(cutVal)):
                if val>cutVal[k]: continue
                pf='q%d'%rankQuantiles[k]
                histos['ranked_mjj'+pf].Fill(mode,mjj,weight)
                histos['ranked_mw'+pf].Fill(mode,mw,weight)
                histos['ranked_mjjresol'+pf].Fill(mode,mjj-mw,weight)
                histos['ranked_drresol'+pf].Fill(mode,drjj2w,weight)

    #build summary plots
    for name,window,relResol in [('mjj',15,True),('mjjresol',5,True),('drresol',0.1,False)]:
        summarizeRankingPerformance([
                ('60wp',      histos['ranked_%sq60'%name]),                
                ('70wp',      histos['ranked_%sq70'%name]),
                ('80wp',      histos['ranked_%sq80'%name]),
                ('inclusive', histos['ranked_%s'%name])
                ],
                                    window,
                                    relResol)

    #save to file
    fOut=ROOT.TFile('workspace.root','RECREATE')
    for h in histos:
        histos[h].SetDirectory(fOut)
        histos[h].Write()
    fOut.Close()
            

if __name__ == "__main__":
    main()
