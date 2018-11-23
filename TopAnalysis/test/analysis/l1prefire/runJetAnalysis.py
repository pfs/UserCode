import ROOT
import os
import sys
from array import array

MAXMUFRAC=0.5

def runJetAnalysis(t):
    """fills histograms to determine prefiring probability"""

    jetPt=[30.0, 36.0, 43.0, 52.0, 63.0, 75.0, 91.0, 109.0, 131.0, 158.0, 190.0, 228.0, 274.0, 397.0, 574.0, 830.0, 1200.0]
    jetEta=[0, 1., 1.5, 2., 2.25, 2.5, 2.75, 3., 3.5]
    histos={
        'njets'     : ROOT.TH1D('njets',';Jet multiplicity;Events',10,0,10),
        'jetpt'     : ROOT.TH1D('jetpt',';Transverse momentum [GeV];Jets',len(jetPt)-1,array('d',jetPt)),
        'jeteta'    : ROOT.TH1D('jeteta',';Pseudo-rapidity;Jets',len(jetEta)-1,array('d',jetEta)),
        'jetetapt'  : ROOT.TH2D('jetpteta',';Pseudo-rapidity;Transverse momentum [GeV];Events',len(jetEta)-1,array('d',jetEta),len(jetPt)-1,array('d',jetPt)),
        'jetetaphi' : ROOT.TH2D('jetetaphi',';Pseudo-rapidity;#phi [rad];Events',57, -3, 3, 72, -3.1415, 3.1415),
        'jetemf'    : ROOT.TH1D('jetemf',';e.m. fraction;Events',10,0,1),
        }
    for key in histos:
        histos[key].SetDirectory(0)
        histos[key].Sumw2()

    def fillHisto(val,name,cat=None):
        key='{0}_{1}'.format(cat,name) if cat else name
        if not key in histos:
            histos[key]=histos[name].Clone(key)
            histos[key].Reset('ICE')
        if isinstance(val,float) or isinstance(val,int):
            histos[key].Fill(val)
        else:
            histos[key].Fill(*val)

    for ientry in xrange(0,1000): #t.GetEntries()):
        t.GetEntry(ientry)
        #print ientry,'/',t.GetEntries()
        
        #select events
        hasHotSpotJet=False
        hasVeryFwdJet=True

        selJets=[]
        for j in xrange(0,t.jet_p4.size()):

            pdgid=t.jet_pdgid[j]
            if pdgid!=0: continue

            emfrac=t.jet_neutralEmFrac[j]
            mufrac=t.jet_muFrac[j]
            if mufrac>MAXMUFRAC : continue

            eta,phi=t.jet_p4[j].Eta(),t.jet_p4[j].Phi()            
            #if t.run>=280919 and t.run<=284044:
            #    if abs(eta+2.81)<0.2 and abs(phi-2.07)<0.2:
            #        hasHotSpotJet=True

            p4=ROOT.TLorentzVector(t.jet_p4[j].Px(),t.jet_p4[j].Py(),t.jet_p4[j].Pz(),t.jet_p4[j].E())
            fillHisto(p4.Pt(),                 'jetpt')
            fillHisto(abs(p4.Eta()),           'jeteta')
            fillHisto(p4.Pt(),                 'jetpt')
            fillHisto([p4.Eta(),p4.Pt()],      'jetetapt')
            fillHisto([abs(p4.Eta()),p4.Phi()],'jetetaphi')
            fillHisto(emfrac,                  'jetemf')
            if abs(eta)>=2:
                if abs(eta)<3.25 : selJets.append(p4)
                else             : hasVeryFwdJet=True

        njets=len(selJets)
        fillHisto(njets,'njets')
        if njets==0 : continue     

    for key in histos:
        histos[key].Draw()
        raw_input(key)

"""
        #match to L1EG
        for ieg in xrange(0,t.L1EG_p4.size()):

            if t.L1EG_p4[ieg].Pt()<30: continue
            isIso=(t.L1EG_iso[ieg]&0x1)
            if not isIso: continue

            p4=ROOT.TLorentzVector(t.L1EG_p4[ieg].Px(),t.L1EG_p4[ieg].Py(),t.L1EG_p4[ieg].Pz(),t.L1EG_p4[ieg].E())
            for j in selJets:
                if p4.DeltaR(j)>0.4 : continue
                histos["jetmatchl1eg"].fill(p4.Eta(),p4.Phi())



      if ( L1EG_p4[iEG].Pt() > 30. and (L1EG_iso[iEG] & 0x1) and DeltaR(L1EG_p4[iEG], jet) < 0.4 ) {
        match_bx |= bit(L1EG_bx[iEG]+2);
        if ( jet.Pt() > 30. and L1EG_bx[iEG] == -1 ) {
          hJet30EGEtaPhi_->Fill(L1EG_p4[iEG].Eta(), L1EG_p4[iEG].Phi(), jweight);
        }
      }
      if ( etaBinCut(jet) and L1EG_bx[iEG] == -1 and (L1EG_iso[iEG] & 0x1) ) {
        if ( DeltaR(L1EG_p4[iEG], jet) < 0.4 ) {
          if ( L1EG_p4[iEG].Pt() > 20. ) match_thr |= bit(0);
          if ( L1EG_p4[iEG].Pt() > 30. ) match_thr |= bit(1);
          if ( L1EG_p4[iEG].Pt() > 40. ) match_thr |= bit(2);
        }
      }
      if ( etaBinCut(jet) and (L1EG_iso[iEG] & 0x1) and L1EG_p4[iEG].Pt() > 20. ) {
        if ( L1EG_bx[iEG] == -1 ) hJetEGdeltaR_bxm1_->Fill(DeltaR(L1EG_p4[iEG], jet), jweight);
        if ( L1EG_bx[iEG] ==  0 ) hJetEGdeltaR_bx0_->Fill(DeltaR(L1EG_p4[iEG], jet), jweight);
        if ( L1EG_bx[iEG] ==  1 ) hJetEGdeltaR_bx1_->Fill(DeltaR(L1EG_p4[iEG], jet), jweight);
      }
    }

    hJetPtEtaEGeffDenom_->Fill(std::abs(jet.Eta()), jet.Pt(), jweight);
    if ( match_bx & bit(-2+2) ) hJetPtEtaEGeffNum_bxm2_->Fill(std::abs(jet.Eta()), jet.Pt(), jweight);
    if ( match_bx & bit(-1+2) ) hJetPtEtaEGeffNum_bxm1_->Fill(std::abs(jet.Eta()), jet.Pt(), jweight);
    if ( match_bx & bit( 0+2) ) hJetPtEtaEGeffNum_bx0_ ->Fill(std::abs(jet.Eta()), jet.Pt(), jweight);
    if ( match_bx & bit( 1+2) ) hJetPtEtaEGeffNum_bx1_ ->Fill(std::abs(jet.Eta()), jet.Pt(), jweight);
    if ( match_bx & bit( 2+2) ) hJetPtEtaEGeffNum_bx2_ ->Fill(std::abs(jet.Eta()), jet.Pt(), jweight);

    if ( etaBinCut(jet) ) {
      hJetEGm1thrDenom_->Fill(jet.Pt(), jweight);
      if ( match_thr & bit(0) ) hJetEGm1thrNum_EGlow_->Fill(jet.Pt(), jweight);
      if ( match_thr & bit(1) ) hJetEGm1thrNum_EGmed_->Fill(jet.Pt(), jweight);
      if ( match_thr & bit(2) ) hJetEGm1thrNum_EGhigh_->Fill(jet.Pt(), jweight);
    }
  }

  if ( forwardJets.size() == 1 and not vetoEvent ) {
    LorentzVector jet = forwardJets[0];
    if ( hjetKinReweight_ != nullptr ) {
      weight *= hjetKinReweight_->GetBinContent(hjetKinReweight_->FindBin(std::abs(jet.Eta()), jet.Pt()));
    }
    hJetPtEtaFinOReffDenom_->Fill(std::abs(jet.Eta()), jet.Pt(), weight);
    if ( etaBinCut(jet) ) hJetL1ADenom_->Fill(jet.Pt(), weight);
    if ( L1GtBx[0] ) {
      hJetPtEtaFinOReffNum_bxm2_->Fill(std::abs(jet.Eta()), jet.Pt(), weight);
      if ( etaBinCut(jet) ) hJetL1ANum_bxm2_->Fill(jet.Pt(), weight);
    }
    if ( L1GtBx[1] ) {
      hJetPtEtaFinOReffNum_bxm1_->Fill(std::abs(jet.Eta()), jet.Pt(), weight);
      if ( etaBinCut(jet) ) hJetL1ANum_bxm1_->Fill(jet.Pt(), weight);
    }
    if ( L1GtBx[2] ) {
      hJetPtEtaFinOReffNum_bx0_->Fill(std::abs(jet.Eta()), jet.Pt(), weight);
      if ( etaBinCut(jet) ) hJetL1ANum_bx0_->Fill(jet.Pt(), weight);
    }
    if ( L1GtBx[3] ) {
      hJetPtEtaFinOReffNum_bx1_->Fill(std::abs(jet.Eta()), jet.Pt(), weight);
      if ( etaBinCut(jet) ) hJetL1ANum_bx1_->Fill(jet.Pt(), weight);
    }
    if ( L1GtBx[4] ) {
      hJetPtEtaFinOReffNum_bx2_->Fill(std::abs(jet.Eta()), jet.Pt(), weight);
      if ( etaBinCut(jet) ) hJetL1ANum_bx2_->Fill(jet.Pt(), weight);
    }
  }

  return kTRUE;
}
"""

def main():
    dirList=['/eos/cms/store/cmst3/group/top/RunIIReReco/l1prefire/Data13TeV_2017C_JetHT']

    t=ROOT.TChain('prefiringVBFAna/l1prefire')
    for d in dirList:
        for f in os.listdir(d):
            t.AddFile(os.path.join(d,f))

    runJetAnalysis(t)


if __name__ == "__main__":
    sys.exit(main())
