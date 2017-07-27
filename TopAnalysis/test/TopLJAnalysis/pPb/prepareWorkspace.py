#!/usr/bin/env python2.7

import ROOT
import optparse
import os,sys
import random
from collections import OrderedDict
from kinTools import *
from datasets import *

EVENTCATEGORIES=[ch+cat for ch in ['e','mu'] for cat in ['1l4j2b','1l4j1b1q','1l4j2q','1f4j2q']]


"""
"""
def getMC(tree,doPseudoTop=False):
    partColl=[]
    baseID = 1000     if doPseudoTop else 1
    nmc    = tree.npt if doPseudoTop else tree.ngp
    tag    = 'pt'     if doPseudoTop else 'gp'
    for j in xrange(0,nmc):
        pdgid=getattr(tree,'%s_pdgid'%tag)[j]
        if abs(pdgid)<=5*baseID or abs(pdgid)==11*baseID or abs(pdgid)==13*baseID :
            if abs(pdgid)<5*baseID :
                if baseID==1 and j<2 : continue #sometimes incoming partons are stored
                pdgid=1
            partColl.append( Particle(pdgid,
                                      getattr(tree,'%s_pt'%tag)[j],
                                      getattr(tree,'%s_eta'%tag)[j],
                                      getattr(tree,'%s_phi'%tag)[j],
                                      getattr(tree,'%s_m'%tag)[j]) )
    return partColl

"""
"""
def getJets(tree,minPt=25.,maxEta=2.4,mcTruth=None,shiftJES=0,jerProf=None):
    partColl=[]
    matchedGen=[]
    for j in xrange(0,tree.nj):
        abseta=abs(tree.j_eta[j])
        if tree.j_pt[j]<minPt or abseta>maxEta : continue
        pdgid = 1 if tree.j_btag[j]==0 else 5
        partColl.append( Particle(pdgid,
                                  tree.j_pt[j],
                                  tree.j_eta[j],
                                  tree.j_phi[j],
                                  tree.j_m[j]) )
        partColl[-1].setRank(tree.j_btag[j])
        
        #residual JES
        if shiftJES!=0: scaleP4(partColl[-1].p4,shiftJES)
        
        #JER smearing
        if   abseta<0.5 : ptSF,ptSF_err=1.109,0.008
        elif abseta<0.8 : ptSF,ptSF_err=1.138,0.013
        elif abseta<1.1 : ptSF,ptSF_err=1.114,0.013
        elif abseta<1.3 : ptSF,ptSF_err=1.123,0.024
        elif abseta<1.7 : ptSF,ptSF_err=1.084,0.011
        elif abseta<1.9 : ptSF,ptSF_err=1.082,0.035
        elif abseta<2.1 : ptSF,ptSF_err=1.140,0.047
        elif abseta<2.3 : ptSF,ptSF_err=1.067,0.053
        elif abseta<2.5 : ptSF,ptSF_err=1.117,0.041

        #inflate due to difference between high and low HF energy
        ptSF_err =ROOT.TMath.Sqrt(ptSF_err**2+0.098**2) 

        cjer=1
        try:

            resol=-1
            for gj in xrange(0,tree.ngj):
                if tree.gj_pt[gj]<1 : continue
                gjet=Particle(0,tree.j_pt[gj],tree.j_eta[gj], tree.j_phi[gj],tree.gj_m[gj])
                if gjet.p4.DeltaR( partColl[-1].p4 ) >0.2 : continue
                resol=(tree.j_pt[j]/tree.gj_pt[gj]-1)
                cjer   = ROOT.TMath.Max(0,tree.gj_pt[gj]+ptSF*(tree.j_pt[j]-tree.gj_pt[gj]))/tree.j_pt[j]
                break

            #fill jer profile
            if jerProf.InheritsFrom('TH1') :          
                if resol>0 : jerProf.Fill(abs(partColl[-1].p4.Eta()),resol)
                #testP4=ROOT.TLorentzVector(partColl[-1].p4)
                #scaleP4(testP4,cjer)
                #jerProf.Fill(abs(partColl[-1].p4.Eta()),testP4.Pt()/partColl[-1].p4.Pt())

        except:
            cjerUp=jerProf.Eval(abs(partColl[-1].p4.Eta()))
            cjerDn=1./cjerUp
            pass

        #add scale uncertainties
        jesUnc=0.028
        if partColl[-1].p4.Pt()<40 : jesUnc=ROOT.TMath.Sqrt(0.02**2+jesUnc**2)
        partColl[-1].setScaleUnc('jesup',(1.0+jesUnc))
        partColl[-1].setScaleUnc('jesdn',(1.0-jesUnc))
        partColl[-1].setScaleUnc('jerup',1 if cjer==0 else max(cjer,1./cjer))
        partColl[-1].setScaleUnc('jerdn',1 if cjer==0 else min(cjer,1./cjer))

        #match parton level, if available
        minDR=999
        imatch=-1
        try:
            for i in xrange(0,len(mcTruth)):
                p=mcTruth[i]
                dR=p.DeltaR(partColl[-1])
                if dR>minDR : continue
                minDR=dR
                imatch=i
            if imatch<0 or imatch in matchedGen: continue
            if minDR>0.4 : continue
            matchedGen.append(imatch)
            partColl[-1].setMCtruth(mcTruth[imatch])
        except:
            pass

    return partColl


"""
starts the workspace
"""
def startWorkspace(opt):

    #start the workspace
    ws  = ROOT.RooWorkspace('w')

    #list of variables of interest
    args=ROOT.RooArgSet()
    varList=['mjj[0,500]',  
             'dmjj_jes[0,2]',
             'dmjj_jer[0,2]',
             'mthad[0,500]', 
             'mtlep[0,500]',
             'corwjj[0,1]', 
             'corthad[0,1]', 
             'cortlep[0,1]' ]
    for var in varList :
        args.add( ws.factory(var) )

    #categories for the data
    sample=ROOT.RooCategory('sample','sample')
    for name in EVENTCATEGORIES : sample.defineType(name)
    getattr(ws,'import')(sample)
    args.add(ws.cat('sample'))

    #define the dataset
    data=ROOT.RooDataSet('data','data',args)
    return ws,data,args

"""
books control distributions and loops over the events to fill the dataset for the fit
"""
def fillDataAndPlots(data,ws,args,tree,leptonSel,mcTruth=False,wjjOrder='drjj',thadOrder='dm2tlep',shiftJES=0,jerProf=None,etaRestr=False):

    p4=ROOT.TLorentzVector(0,0,0,0)

    baseCat='1l4j'
    for lid in leptonSel:
        if lid>1000 :
            baseCat='1f4j'
            break

    plots={}
    plots['recacc']=ROOT.TH1F('recacc',';Variation;Acceptance',3,0,3)
    plots['jerprofile']=ROOT.TProfile('jerprofile',';Pseudo-rapidity;#delta p_{T}/p_{T}',5,0,2.5)
    for cat in EVENTCATEGORIES:

        plots['jeteff_%s'%cat] = ROOT.TH1F('jeteff_%s'%cat,';Matched jet multiplicity;Events',11,0,11)
        xbin=1
        plots['jeteff_%s'%cat].GetXaxis().SetBinLabel(xbin,'inc')
        for step in ['inc','wjj']:
            for label in ['1b','2b','1q','2q']:
                xbin+=1
                plots['jeteff_%s'%cat].GetXaxis().SetBinLabel(xbin,step+label)
        plots['jeteff_%s'%cat].GetXaxis().SetBinLabel(10,'thad')
        plots['jeteff_%s'%cat].GetXaxis().SetBinLabel(11,'tlep')

        for t in ['','_cor','_wro']:
            plots['drjj_%s%s'%(cat,t)]  = ROOT.TH1F('drjj_%s%s'%(cat,t), ';#DeltaR; Events',  6,0, 6)
            plots['mjj_%s%s'%(cat,t)]   = ROOT.TH1F('mjj_%s%s'%(cat,t),  ';W_{jj} mass [GeV]; Events / 10 GeV',  20,0, 200)            
            plots['mthad_%s%s'%(cat,t)] = ROOT.TH1F('mthad_%s%s'%(cat,t),';t_{had} mass [GeV]; Events / 10 GeV', 30,50,350)
            plots['mtlep_%s%s'%(cat,t)] = ROOT.TH1F('mtlep_%s%s'%(cat,t),';t_{lep} mass [GeV]; Events / 10 GeV', 30,50,350)
            if t=='':
                plots['minmlb_%s%s'%(cat,t)]   = ROOT.TH1F('minmlb_%s%s'%(cat,t),  ';min M_{lb} [GeV]; Events / 10 GeV',  20,0, 200)
                plots['maxmlb_%s%s'%(cat,t)]   = ROOT.TH1F('maxmlb_%s%s'%(cat,t),  ';max M_{lb} [GeV]; Events / 10 GeV',  20,0, 200)
                plots['mtw_%s'%cat] = ROOT.TH1F('mtw_%s'%cat,';M_{T}(l,MET) [GeV]; Events / 10 GeV', 20,0,200)
                plots['met_%s'%cat] = ROOT.TH1F('met_%s'%cat,';MET [GeV]; Events / 10 GeV', 20,0,200)
                plots['ntracks_%s'%cat] = ROOT.TH1F('ntracks_%s'%cat,';Track multiplicity; Events', 20,0,200)
                plots['ntrackshp_%s'%cat] = ROOT.TH1F('ntrackshp_%s'%cat,';UE track multiplicity; Events', 10,0,100)

        #check this point fwd
        plots['pt_l_%s'%cat]  = ROOT.TH1F('pt_l_%s'%cat,';Lepton transverse momentum [GeV];Events',20,0,100)
        plots['y_l_%s'%cat]   = ROOT.TH1F('y_l_%s'%cat,';Lepton rapidity;Events',20,-2.1,2.1)

        plots['pt_b_%s'%cat]  = ROOT.TH1F('pt_b_%s'%cat,';b transverse momentum [GeV];b-jets',20,0,100)
        plots['y_b_%s'%cat]   = ROOT.TH1F('y_b_%s'%cat,';b rapidity;b-jets',20,-2.1,2.1)

        for t in ['jj','l']:
            plots['y_w%s_%s'%(t,cat)]   = ROOT.TH1F('y_w%s_%s'%(t,cat),  '; W_{%s} rapidity;Events'%t,20,-2.1,2.1)
            plots['pt_w%s_%s'%(t,cat)]  = ROOT.TH1F('pt_w%s_%s'%(t,cat), '; W_{%s} transverse momentum [GeV]; Events'%t,50,0,500)



    #loop over all the events
    totalEntries=tree.GetEntries()
    for i in xrange(0,totalEntries):

        tree.GetEntry(i)
        if i%100==0 : sys.stdout.write('\r [ %d/%d ] done' %(i,totalEntries))

        #lepton selection
        if len(leptonSel)!=0 and not tree.l_id in leptonSel : continue
        leptonCat='e' if (abs(tree.l_id)==11 or abs(tree.l_id)==1100) else 'mu'
        if tree.l_pt<30 or abs(tree.l_eta)>2.1: continue
        if etaRestr and leptonCat=='e' and abs(tree.l_eta)>1.5: continue

        #apply filtering, if available (older or pp versions) don't have this
        try:
            if tree.pPAprimaryVertexFilter<1 or tree.HBHENoiseFilterResultRun2Loose<1 : continue
        except:
            pass

        #readout the event
        mcTruth=getMC(tree=tree) if mcTruth else None
        lepton=Particle(tree.l_id,tree.l_pt,tree.l_eta,tree.l_phi,tree.l_m)
        try:
            for p in mcTruth:
                if p.DeltaR(lepton)>0.4 : continue
                lepton.setMCtruth(p)
                break
        except:
            pass
        MET=Particle(0,tree.met_pt,0.,tree.met_phi,0.)
        jets=getJets(tree=tree,minPt=25,maxEta=2.4,mcTruth=mcTruth,shiftJES=shiftJES,jerProf=jerProf if jerProf else plots['jerprofile'])
        jets.sort(key=lambda x: x.rankVal, reverse=True)

        #for the acceptance systematics
        njets,njetsUp,njetsDn=0,0,0
        for j in jets:
            if j.p4.Pt()>25 : njets+=1
            if j.p4.Pt()*j.scaleUnc['jesup'] > 25 : njetsUp+=1
            if j.p4.Pt()*j.scaleUnc['jesdn'] > 25 : njetsDn+=1
        if njets>=4   : plots['recacc'].Fill(1)
        if njetsUp>=4 : plots['recacc'].Fill(2)
        if njetsDn>=4 : plots['recacc'].Fill(3)

        if len(jets)<4 : continue

        #separate the jets according to the b-tagging decision
        bJets=jets[0:2]
        lightJets=[]
        for j in xrange(2,len(jets)):
            #if jets[j].rankVal!=0 : continue
            lightJets.append(jets[j])
        if len(lightJets)<2 : continue

        cat=leptonCat+baseCat
        if bJets[0].rankVal*bJets[1].rankVal==9 : cat += '2b'
        elif bJets[0].rankVal==3 or bJets[1].rankVal==3 : cat += '1b1q'
        else : cat += '2q'

        if baseCat=='1f4j' and not '1f4j2q' in cat : continue

        #W->qq' system
        WjjCandidates=buildWjj(lightJets=lightJets,orderBy=wjjOrder)
        Wjj=WjjCandidates[0]

        #W->lnu
        Wlnu=buildWlnu(lepton,MET)

        #ttbar hypothesis
        ttCandidates=buildTTbar(bJets=bJets,Wjj=Wjj,Wlnu=Wlnu)

        minmlb=min((lepton.p4+bJets[0].p4).M(),(lepton.p4+bJets[1].p4).M())
        maxmlb=max((lepton.p4+bJets[0].p4).M(),(lepton.p4+bJets[1].p4).M())

        thad,tlep = ttCandidates[0]
        if thadOrder=='dm2tlep':
            chi=[ abs(x.p4.M()-y.p4.M()) for x,y in ttCandidates ]
            if chi[0]>chi[1] : thad,tlep = ttCandidates[1]
        if thadOrder=='dm2pdg':
            chi=[ abs(x.p4.M()-172.5) for x,_ in ttCandidates ]
            if chi[0]>chi[1] : thad,tlep = ttCandidates[1]
        if thadOrder=='dr':
            chi=[ x.daughters[0].DeltaR(x.daughters[1]) for x,_ in ttCandidates]
            if chi[0]>chi[1] : thad,tlep = ttCandidates[1]

        plots['drjj_%s'%cat].Fill( Wjj.daughters[0].DeltaR(Wjj.daughters[1]) )
        plots['mjj_%s'%cat].Fill(Wjj.p4.M())
        plots['mtlep_%s'%cat].Fill(tlep.p4.M())
        plots['mthad_%s'%cat].Fill(thad.p4.M())
        plots['minmlb_%s'%cat].Fill(minmlb)
        plots['maxmlb_%s'%cat].Fill(maxmlb)

        #check matching efficiency
        goodWjj,goodTlep,goodThad=True,True,True
        if mcTruth:
            nbmatch,nqmatch=0,0
            for b in bJets     : nbmatch += 1 if b.mcTruth and abs(b.mcTruth.id)==5 else 0
            for j in lightJets : nqmatch += 1 if j.mcTruth and abs(j.mcTruth.id)==1 else 0

            nqmatchAfterWjj=0
            for j in Wjj.daughters : nqmatchAfterWjj += 1 if j.mcTruth and abs(j.mcTruth.id)==1 else 0
            plots['mjj_%s_%s'%(cat,'cor' if nqmatchAfterWjj==2 else 'wro')].Fill(Wjj.p4.M())
            plots['drjj_%s_%s'%(cat,'cor' if nqmatchAfterWjj==2 else 'wro')].Fill( Wjj.daughters[0].DeltaR(Wjj.daughters[1]) )
            goodWjj=True if nqmatchAfterWjj==2 else False

            if tlep.daughters[0].mcTruth is None or tlep.daughters[1].daughters[0].mcTruth is None:
                goodTlep=False
            else:
                bmatchid=tlep.daughters[0].mcTruth.id
                blepid=tlep.daughters[1].daughters[0].mcTruth.id
                if not bmatchid*blepid in [-55,-65]: goodTlep=False
            plots['mtlep_%s_%s'%(cat,'cor' if goodTlep else 'wro')].Fill(tlep.p4.M())

            goodThad=True if goodWjj else False
            if thad.daughters[0].mcTruth is None:
                goodThad=False
            elif abs(thad.daughters[0].mcTruth.id)!=5 :
                goodThad=False
            plots['mthad_%s_%s'%(cat,'cor' if goodThad else 'wro')].Fill(thad.p4.M())

            plots['jeteff_%s'%cat].Fill(0)
            if nbmatch>0 :
                plots['jeteff_%s'%cat].Fill(nbmatch)
                plots['jeteff_%s'%cat].Fill(nbmatch+4)
            if nqmatch>0 :
                plots['jeteff_%s'%cat].Fill(min(nqmatch,2)+2)
            if nqmatchAfterWjj>0:
                plots['jeteff_%s'%cat].Fill(nqmatchAfterWjj+6)
            if goodThad: plots['jeteff_%s'%cat].Fill(9)
            if goodTlep: plots['jeteff_%s'%cat].Fill(10)


        #add entry to dataset
        ws.cat('sample').setLabel(cat)
        args.setCatLabel(cat)
        argValList=[('mjj',    Wjj.p4.M()),
                    ('dmjj_jes',0.5*(abs(1-Wjj.scaleUnc['jesup'])+abs(1-Wjj.scaleUnc['jesdn'])) ),
                    ('dmjj_jer',0.5*(abs(1-Wjj.scaleUnc['jerup'])+abs(1-Wjj.scaleUnc['jerdn'])) ),
                    ('mthad',  thad.p4.M()),
                    ('mtlep',  tlep.p4.M()),
                    ('corwjj', 1 if goodWjj  else 0),
                    ('corthad',1 if goodThad else 0),
                    ('cortlep',1 if goodTlep else 0)]
        for argName,argVal in argValList:
            args.find(argName).setVal(argVal)
        data.add(args)

        #ue tracks
        plots['ntracks_%s'%cat].Fill(tree.ntracks)
        plots['ntrackshp_%s'%cat].Fill(tree.ntracks-tree.ntracks_hp)
        plots['met_%s'%cat].Fill(MET.p4.Pt())
        mtw=ROOT.TMath.Sqrt(2*lepton.p4.Pt()*MET.p4.Pt()*(1-ROOT.TMath.Cos(lepton.p4.DeltaPhi(MET.p4))))
        plots['mtw_%s'%cat].Fill(mtw)
        plots['y_l_%s'%cat].Fill(lepton.p4.Rapidity())
        plots['pt_l_%s'%cat].Fill(lepton.p4.Pt())
        for b in bJets:
            plots['y_b_%s'%cat].Fill(b.p4.Rapidity())
            plots['pt_b_%s'%cat].Fill(b.p4.Pt())
        plots['y_wjj_%s'%cat].Fill(Wjj.p4.Rapidity())
        plots['pt_wjj_%s'%cat].Fill(Wjj.p4.Pt())
        plots['y_wl_%s'%cat].Fill(Wlnu.p4.Rapidity())
        plots['pt_wl_%s'%cat].Fill(Wlnu.p4.Pt())


        #add entry to dataset
        #ws.cat('sample').setLabel(cat)
        #args.setCatLabel(cat)
        #argValList=[('y_l',lepton.Rapidity()),
        #            ('y_bhad',bhad.Rapidity()),
        #            ('y_blep',blep.Rapidity() if blep else -999.),
        #            ('m_wjj',wjj.M()),
        #            ('dr_jj',drjj),
        #            ('pt_wjj',wjj.Pt()),
        #            ('y_wjj',wjj.Rapidity()),
        #            ('pt_wl',wl.Pt()),
        #            ('y_wl',wl.Rapidity()),
        #            ('m_bwjj',bwjj.M()),
        #            ('y_bwjj',bwjj.Rapidity()),
        #            ('m_bwl',bwl.M() if bwl else -999),
        #            ('y_bwl',bwl.Rapidity() if bwl else -999),
        #            ('dy_bwjjbwl',abs(bwl.Rapidity()-bwjj.Rapidity()) if bwl else -999),
        #            ('sy_bwjjbwl',abs(bwl.Rapidity()+bwjj.Rapidity()) if bwl else -999),
        #           ]
        #for argName,argVal in argValList:
        #    args.find(argName).setVal(argVal)
        #data.add(args)

    return plots


"""
display the results of the fit
"""
def showControlPlots(allPlots,fOut=None,outDir='plots/'):

    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.04)
    c.SetBottomMargin(0.12)

    for p in allPlots[0]:
        c.Clear()
        totalH,stackH=None,None
        for coll in allPlots:
            if not p in coll: continue
            if totalH is None:
                totalH=coll[p].Clone('%s_total'%p)
                totalH.SetDirectory(0)
                totalH.SetFillStyle(0)
                totalH.SetFillColor(0)
                totalH.SetLineWidth(2)
                totalH.SetLineColor(1)
                totalH.SetTitle('total')
                totalH.SetMarkerStyle(20)
                totalH.Reset('ICE')
            if coll[p].Integral()==0 : continue
            totalH.Add( coll[p] )
            if stackH is None:
                stackH=ROOT.THStack('%s_stack'%p,'')
            stackH.Add( coll[p], 'hist' )

        if totalH.Integral()==0 : continue

        is2D=True if totalH.InheritsFrom('TH2') else False

        leg = ROOT.TLegend(0.14,0.88,0.5,0.75)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.04)
        leg.AddEntry(totalH,'total','p')

        if is2D:
            totalH.Draw('coltext')
        else:
            frame=totalH.Clone('frame')
            frame.Reset('ICE')
            frame.Draw()
            frame.GetYaxis().SetTitleOffset(1.4)
            frame.GetYaxis().SetRangeUser(0,totalH.GetMaximum()*1.5)
            stackH.Draw('histsame')
            totalH.Draw('same')

            for i in xrange(0,stackH.GetStack().GetEntriesFast()):
                h=stackH.GetStack().At(i)
                leg.AddEntry(h,h.GetTitle(),'f')

        leg.Draw()

        label = ROOT.TLatex()
        label.SetNDC()
        label.SetTextFont(42)
        label.SetTextSize(0.04)
        label.DrawLatex(0.15,0.9,'#bf{CMS} #it{preliminary}')
        if '_pp' in outDir:
            label.DrawLatex(0.58,0.96,'#scale[0.8]{881 pb^{-1} (pp at #sqrt{s}=8.0 TeV)}')
        else:
            label.DrawLatex(0.58,0.96,'#scale[0.8]{180 nb^{-1} (pPb at #sqrt{s}=8.16 TeV)}')
        tag='#mu+jets'
        if '1f'   in p : tag ='non-iso '+tag
        if '4j'   in p : tag='#geq4j,'
        if '2q'   in p : tag+='=0b'
        if '1b1q' in p : tag+='=1b'
        if '2b'   in p : tag+='#geq2b'
        label.DrawLatex(0.68,0.9,'#scale[0.8]{#it{%s}}'%tag)

        c.cd()
        c.Modified()
        c.Update()
        c.SaveAs('%s/%s.png'%(outDir,p))

        if fOut:
            fOut.cd()
            totalH.SetDirectory(fOut)
            totalH.Write(p)



"""
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-d', '--data',       dest='data',      default='MC8.16TeV_TTbar_pPb_Pohweg',   type='string',    help='dataset to use [%default]')
    parser.add_option('-v', '--verbose',    dest='verbose',   default=0,                              type=int,         help='Verbose mode [%default]')
    parser.add_option(      '--shiftJES',   dest='shiftJES',  default=0.97,                           type=float,       help='shift jes [%default]')
    parser.add_option(      '--wjjOrder',   dest='wjjOrder',  default='drjj',                         type='string',    help='wjj ordering (drjj,mjj,sumpt) [%default]')
    parser.add_option(      '--thadOrder',  dest='thadOrder', default='dm2tlep',                      type='string',    help='thad ordering (dr,dm2tlep,dm2pdg) [%default]')
    parser.add_option(      '--jerProf',    dest='jerProf',   default=None,                      type='string',    help='plotter with JER profile [%default]')
    parser.add_option(      '--etaRestr',   dest='etaRestr',   default=False,   action="store_true", help='restrict electron selection to barrel [%default]')
    (opt, args) = parser.parse_args()

    jerProf=None
    if opt.jerProf:
        print 'Reading JER profile from',opt.jerProf
        fIn=ROOT.TFile.Open(opt.jerProf)
        jerProf=ROOT.TGraph(fIn.Get('jerprofile'))        
        fIn.Close()

    #init the workspace
    ws,data,args=startWorkspace(opt)

    #prepare control plots and fill dataset
    allPlots=[]
    dataset=defineDataset(data=opt.data)
    for coll in dataset:
        urlList,color,leptonSel,mcTruth=dataset[coll]
        tree=ROOT.TChain('data')
        for url in urlList: 
            tree.AddFile(url)
            print url
        allPlots.append( fillDataAndPlots(data,ws,args,tree,leptonSel,mcTruth,opt.wjjOrder,opt.thadOrder,opt.shiftJES,jerProf,opt.etaRestr) )

        color=ROOT.TColor.GetColor(color)
        for p in allPlots[-1]:
            if 'Data' in opt.data:
                allPlots[-1][p].SetBinErrorOption(ROOT.TH1.kPoisson)
            allPlots[-1][p].SetFillStyle(1001)
            allPlots[-1][p].SetFillColor(color)
            allPlots[-1][p].SetLineColor(color)
            allPlots[-1][p].SetTitle(coll)
            allPlots[-1][p].Sumw2()
            allPlots[-1][p].SetDirectory(0)

    pfix='_etarestr' if opt.etaRestr else ''

    #import data and save workspace
    getattr(ws,'import')(data)
    ws.writeToFile('workspace_%s%s.root'%(opt.data,pfix),True)

    #show results
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(True)
    outDir='plots/%s%s'%(opt.data,pfix)
    os.system('mkdir -p %s'%outDir)
    fOut=ROOT.TFile.Open('%s/controlplots.root'%outDir,'RECREATE')
    showControlPlots(allPlots,fOut,outDir=outDir)
    fOut.Close()


if __name__ == "__main__":
    main()
