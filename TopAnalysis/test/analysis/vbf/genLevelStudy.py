import ROOT
import sys
import os
from TopLJets2015.TopAnalysis.HistoTool import *
from TopLJets2015.TopAnalysis.myProgressBar import *

def getChain(url,maxFiles):

    """returns the chain of events"""
    sumw=0.
    t=ROOT.TChain('analysis/data')
    ifile=0
    for f in [os.path.join(url,x) for x in os.listdir(url)]:
        t.AddFile(f)
        inF=ROOT.TFile.Open(f)
        sumw += inF.Get('analysis/fidcounter').GetBinContent(1,1)
        inF.Close()
        ifile+=1
        if ifile>=maxFiles: break
    return t,sumw

def fillP4(pt,eta,phi,m):

    """fills a tlorentzvector"""

    p4=ROOT.Math.PtEtaPhiMVector(pt,eta,phi,m)
    return p4

def fillGenLevelPlots(t,sumw=1.0,leadPt=50,subLeadPt=40,minMjj=500,minCenJetPt=15,out_name='histos'):

    """makes a simple generator-level selection and fills control histograms"""

    ht=HistoTool()
    ht.add(ROOT.TH1F("sel", ";Selection step;Events", 5,0,5))
    for xbin,xtit in enumerate(['gen','1#gamma','#geq2j','lpT','hpT']):
        ht.get('sel').GetXaxis().SetBinLabel(xbin+1,xtit)
    ht.add(ROOT.TH1F("pt",         ";Boson p_{T}[GeV];Events",            25,50,550))  
    ht.add(ROOT.TH1F("eta",        ";Boson #eta;Events",                  25,-2.5,2.5))
    ht.add(ROOT.TH1F("mindraj",    ";min#DeltaR(#gamma,j);Events",        25, 0.4,5))  
    ht.add(ROOT.TH1F("leadpt",     ";Leading jet p_{T} [GeV];Events",     25,40,800))
    ht.add(ROOT.TH1F("subleadpt",  ";Sub-leading jet p_{T} [GeV];Events", 25,40,500))  
    ht.add(ROOT.TH1F("leadeta",    ";Leading jet |#eta|;Events",          40,-5,5))  
    ht.add(ROOT.TH1F("subleadeta", ";Sub-leading jet |#eta|;Events",      40,-5,5))  
    ht.add(ROOT.TH1F("mjj",        ";Dijet invariant mass [GeV];Events",  40,minMjj,4000))  
    ht.add(ROOT.TH1F("detajj" ,    ";#Delta#eta(jj);Events",              20,0,8)) 
    ht.add(ROOT.TH1F("balance",    ";System p_{T} balance [GeV];Events",  20,0,300)) 
    ht.add(ROOT.TH1F("ystar",      ";#eta^{*}=#eta-(1/2)[#eta(1)+#eta(2)];Events",  25,-3,3)) 
    ht.add(ROOT.TH1F("ncentj",     ";Number of central jets;Events",                3,0,3))
    ht.add(ROOT.TH1F("centjystar", ";#eta^{*}=#eta-(1/2)[#eta(1)+#eta(2)];Events",  25,-3,3))

    maxEvts=t.GetEntries()
    print 'Processing',maxEvts,'events'
    ht.fill( (0,sumw),      'sel',  ['inc'],'')
    for i in range(maxEvts):

        if i%1000==0 : drawProgressBar(float(i)/maxEvts)

        t.GetEntry(i)
        wgt=t.g_w[0]
        
        #select photon+jets
        photon_particle=None
        alljets=[]
        for ip in range(t.ng):
            pid=t.g_id[ip]
            pt,eta,phi,m=t.g_pt[ip],t.g_eta[ip],t.g_phi[ip],t.g_m[ip]
            if pid==22 and abs(eta)<2.5 and pt>75 : 
                photon_particle=fillP4(pt,eta,phi,0)
            if abs(pid)<=5 and pt>minCenJetPt and abs(eta)<4.7: 
                alljets.append( fillP4(pt,eta,phi,m) )

        #check if photon is trigger-able
        if not photon_particle: continue
        ptPhoton=photon_particle.Pt()
        etaPhoton=photon_particle.Eta()
        isBarrel=True if abs(etaPhoton)<1.5 else False        
        if ptPhoton<200 and not isBarrel:
            continue
        ht.fill( (1,wgt),      'sel',  ['inc'],'')

        #select jets
        alljets=[x for x in alljets if ROOT.Math.VectorUtil.DeltaR(x,photon_particle)>0.4]
        jets=[x for x in alljets if x.Pt()>subLeadPt]
        if len(jets)<2 : continue
        if jets[0].Pt()<leadPt : continue
        mjj=(jets[0]+jets[1]).M()
        if mjj<minMjj : continue
        ht.fill( (2,wgt),      'sel',  ['inc'],'')

        #assign final category
        cats=[]
        if ptPhoton>75 and isBarrel: 
            ht.fill( (3,wgt),      'sel',  ['inc'],'')
            cats.append(0)
        if ptPhoton>200: 
            ht.fill( (4,wgt),      'sel',  ['inc'],'')
            cats.append(1)
        pcats=['inc'] + ['cat%d'%c for c in cats]

        #plot some kinematics
        mindraj=min( [ROOT.Math.VectorUtil.DeltaR(jets[ij],photon_particle) for ij in range(2)] )
        maxY=max(jets[0].Eta(),jets[1].Eta())
        minY=min(jets[0].Eta(),jets[1].Eta())
        syjj=(maxY+minY)
        dyjj=abs(maxY-minY)
        ptBal=(jets[0]+jets[1]+photon_particle).Pt()
        ystar=photon_particle.Eta()-0.5*syjj
        cen_jets=[x for x in alljets if x.Pt()>20 and  x.Rapidity()>minY+0.2 and x.Rapidity()<maxY-0.2]
        ncenj=len(cen_jets)
        cenjystar=None
        if ncenj>0:
            cenjystar=cen_jets[0].Eta()-0.5*syjj

        ht.fill( (ptPhoton,wgt),      'pt',         pcats,'')
        ht.fill( (etaPhoton,wgt),     'eta',        pcats,'')
        ht.fill( (mindraj,wgt),       'mindraj',    pcats,'')
        ht.fill( (jets[0].Pt(),wgt),  'leadpt',     pcats,'')
        ht.fill( (jets[1].Pt(),wgt),  'subleadpt',  pcats,'')
        ht.fill( (jets[0].Eta(),wgt), 'leadeta',    pcats,'')
        ht.fill( (jets[1].Eta(),wgt), 'subleadeta', pcats,'')
        ht.fill( (mjj,wgt),           'mjj',        pcats,'')
        ht.fill( (dyjj,wgt),          'detajj',     pcats,'')
        ht.fill( (ptBal,wgt),         'balance',    pcats,'')
        ht.fill( (ystar,wgt),         'ystar',      pcats,'')
        ht.fill( (ncenj,wgt),         'ncentj',     pcats,'')
        ht.fill( (cenjystar,wgt),     'centjystar', pcats,'')

    ht.postNormalize(1./sumw)
    ht.writeToFile(out_name+'.root')


def fillGenLevelPlotsPacked(args):

    """ unpacks the information needed to fill generator level plots """

    url,out,maxFiles=args
    os.system('mkdir -p ' + out)
    t,sumw=getChain(url,maxFiles=maxFiles)
    out_name=os.path.join(out,url.split('/')[-1])
    print 'Starting fillGenLevelPlotsPacked with',url,' #events=',t.GetEntries(),' output @',out_name
    fillGenLevelPlots(t=t,sumw=sumw,out_name=out_name)


def main():

    url=sys.argv[1]
    out=sys.argv[2]
    maxFiles=1
    if len(sys.argv)>3: maxFiles=int(sys.argv[3])

    task_list=[]
    for f in os.listdir(url):
        if not 'EWKAJJ' in f and not 'IntAJJ' in f: continue
        task_list.append( (os.path.join(url,f),out,maxFiles) )

    import multiprocessing as MP
    pool = MP.Pool(8)
    pool.map(fillGenLevelPlotsPacked,task_list)


if __name__ == "__main__":
    sys.exit(main())
