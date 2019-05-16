import ROOT
import os
from runExclusiveAnalysis import getTracksPerRomanPot,buildDiproton
from TopLJets2015.TopAnalysis.HistoTool import *
from TopLJets2015.TopAnalysis.myProgressBar import *
from TopLJets2015.TopAnalysis.Plot import *


def buildControlPlots(tag,excSel):

    """loop over data events and make control plots"""

    baseDir='/eos/cms/store/cmst3/user/psilva/ExclusiveAna/ab05162/Chunks/'
    
    if len(tag)==1:
        tree=ROOT.TChain('tree')
        for f in os.listdir(baseDir):
            if not 'Data13TeV_2017%s_DoubleMuon'%tag in f : continue
            tree.AddFile(os.path.join(baseDir,f))
    elif 'MC13TeV_ppxz_m' in tag:
        tree=ROOT.TChain('analysis/data')
        for x in [120,130,140,150]:
            tree.AddFile(os.path.join(baseDir,'%s_x%d_0.root'%(tag,x)))

    nentries=tree.GetEntries()
    print 'Analysing',nentries,'events for tag',tag
    if nentries==0:
        return

    #book histograms
    ht=HistoTool()
    ht.add(ROOT.TH1F('n','RP;Proton multiplicity;PDF', 5,0,5))
    ht.add(ROOT.TH1F('csi',';#xi;PDF',                 50,0,0.3))
    ht.add(ROOT.TH2F('csi2d',';#xi(1);#xi(2);PDF',     50,0,0.3,50,0,0.3))
    ht.add(ROOT.TH1F('mpp',';m_{pp} [GeV];PDF',        100,0,2500))
    ht.add(ROOT.TH2F('mpp2d',';m_{pp}(far) [GeV];m_{pp}(near) PDF',                 100,0,2500,100,0,2500))
    ht.add(ROOT.TH1F('dmpp','Near-Far;m_{pp}(near)-m_{pp}(far) [GeV];PDF',          100,-200,200))
    ht.add(ROOT.TH1F('mmass',';Missing mass [GeV];Events',                          100,-1000,3000))
    ht.add(ROOT.TH2F('mmass2d',';Missing mass (near) [GeV];Missing mass (far);PDF', 100,-1000,3000,100,-1000,3000))
    ht.add(ROOT.TH1F('dmmass','Near-Far;#Delta missing mass [GeV];Events',          100,-200,200))
    ht.add(ROOT.TH1F('nvtx',';Vertex multiplicity;Events',                          50,0,50))
    ht.add(ROOT.TH1F('xangle',';LHC crossing angle [#murad];Events',                4,120,160))
    ht.add(ROOT.TH1F('phiboson',';#phi(boson);Events',                              50,-3.15,3.15))
    ht.add(ROOT.TH1F('yboson',';Boson rapidity;Events',                             50,-2,2))
    ht.add(ROOT.TH1F('dyboson',';y(boson)-y(central);Events',                       50,-2,2))
    

    for i in range(nentries):
        tree.GetEntry(i)
        if not tree.isZ : continue
        if tree.bosonpt>10 : continue
        if i%1000==0 : drawProgressBar(float(i)/float(nentries))

        #get basic information from the event
        boson=ROOT.TLorentzVector(0,0,0,0)
        boson.SetPtEtaPhiM(tree.bosonpt,tree.bosoneta,tree.bosonphi,tree.mboson)
        rp023,rp123=getTracksPerRomanPot(tree)
        rp003,rp103=getTracksPerRomanPot(tree,False,False)
        
        #build up category flags
        passPix=False
        passStrip=False
        passMatch=False        
        csi_rp023, csi_rp123 = -1,-1
        csi_rp003, csi_rp103 = -1,1
        if excSel:
            if len(rp023)==1 and len(rp123)==1:
                passPix=True
                csi_rp023=rp023[0]
                csi_rp123=rp123[0]
            if len(rp003)==1 and len(rp103)==1 : 
                passStrip=True
                csi_rp003=rp003[0]
                csi_rp103=rp103[0]

        else:

            passPix   = True if len(rp023)>0 and len(rp123)>0 else False
            passStrip = True if len(rp003)==1 and len(rp103)==1 else False

            #take the max. csi for pixels as starting point
            csi_rp023=max(rp023) if len(rp023)>0 else -1
            csi_rp123=max(rp123) if len(rp123)>0 else -1

            #if matching to strips is found use that instead
            if len(rp003)==1:
                csi_rp003=rp003[0]
                for x in rp023:
                    if abs(csi_rp003-x)>0.02: continue
                    csi_rp023=x
                    break
            if len(rp103)==1:
                csi_rp103=rp103[0]
                for x in rp123:
                    if abs(csi_rp103-x)>0.02: continue
                    csi_rp123=x
                    break

        if passPix and passStrip:
            if abs(csi_rp003-csi_rp023) <0.02:
                if abs(csi_rp103-csi_rp123)<0.02:
                    passMatch=True
            
        #diproton kinematics
        pp_far=buildDiproton([[csi_rp023],[csi_rp123]]) if passPix else None
        mmass_far=(pp_far-boson).M() if pp_far else None
        pp_near=buildDiproton([[csi_rp003],[csi_rp103]]) if passStrip else None
        mmass_near=(pp_near-boson).M() if pp_near else None


        #categories
        cats=['inc']
        if passPix                             : cats += ['passPix']
        if             passStrip               : cats += ['passStrip']
        if passPix and passStrip               : cats += ['passPixandpassStrip']
        if passPix and passStrip and passMatch : cats += ['passPixandpassStripmatched']
        if             passStrip and passMatch : cats += ['passStripmatched']
        if passPix               and passMatch : cats += ['passPixmatched']
        
        #global variables
        ht.fill((tree.nvtx,1),        'nvtx',     cats)
        ht.fill((tree.beamXangle,1),  'xangle',   cats)
        ht.fill((boson.Rapidity(),1), 'yboson',   cats)
        ht.fill((boson.Phi(),1),      'phiboson', cats)
        ht.fill((boson.Pt(),1),       'ptboson',  cats)

        #individual RP plots
        for rpinfo,rp in [(rp003,'RP003'),(rp023,'RP023'),(rp103,'RP103'),(rp123,'RP123')]:
            ht.fill((len(rpinfo),1),'n',cats,rp)
            for x in rpinfo:
                ht.fill((x,1),'csi',cats,rp)

        #RP correlations
        for csi1,csi2,rpcor in [(csi_rp003,csi_rp023,'pos'),(csi_rp103,csi_rp123,'neg')]:
            if csi1<0 or csi2<0 : continue
            ht.fill((csi1,csi2,1),'csi2d',cats,rpcor)

        #diproton kinematics
        for pp,mmass,csi0,csi1,pptag in [(pp_near,mmass_near,csi_rp003,csi_rp103,'near'),
                                         (pp_far, mmass_far, csi_rp023,csi_rp123,'far')]:
            if pp is None: continue
            ht.fill( (pp.M(),1), 'mpp',   cats, pptag )
            ht.fill( (mmass,1),  'mmass', cats, pptag )

            ycen=0.5*ROOT.TMath.Log(csi0/csi1)
            ht.fill((boson.Rapidity()-ycen,1), 'dyboson',   cats)


        #diproton correlations
        if pp_near and pp_far:
            ht.fill( (pp_near.M()-pp_far.M(),1), 'dmpp',    cats)
            ht.fill( (pp_near.M(),pp_far.M(),1), 'mpp2d',   cats)
            ht.fill( (mmass_near-mmass_far,1),   'dmmass',  cats)
            ht.fill( (mmass_near,mmass_far,1),   'mmass2d', cats)

    recTag='exc' if excSel else 'inc'
    outURL='RPcontrol_%s_era%s.root'%(recTag,tag)
    ht.writeToFile(outURL)
    print 'Results can be found in',outURL

def drawPlots(data,mc,period):

    """draw final plots comparing distributions in data and in mc"""

    files={'data':data,'mc':mc}

    base='passStrip'
    baseTitle='=2 strip'
    base='passPix'
    baseTitle='=2 px'
    for d,pfix in [
        ('xangle',''),    
        #('nvtx',''),      ('yboson',''),  ('dyboson',''),
        #('n','RP003'),    ('n','RP023'),    ('n','RP103'),   ('n','RP123'),
        #('csi','RP003'),  ('csi','RP023'),  ('csi','RP103'), ('csi','RP123'),
        ('mpp','far'),    ('mpp','near'),   ('mmass','far'),  ('mmass','near'),
        #('dmpp',''),      ('dmmass',''),
        ]:

        hpercat={}

        for tag in ['','matched']:
            for x in files:

                for cat in [base,'passPixandpassStrip']:
                    try:
                        hpercat[ (x,tag,cat) ]=files[x].Get('%s_%s%s%s'%(d,cat,tag,pfix))
                        hpercat[ (x,tag,cat) ].SetDirectory(0)
                        hpercat[ (x,tag,cat) ].Sumw2()
                        hpercat[ (x,tag,cat) ].GetYaxis().SetTitle('PDF')
                    except:
                        pass
                
                #remove overlap
                #hpercat[(x,tag,base)].Add(hpercat[(x,tag,'passPixandpassStrip')],-1)
                
                #build the total
                hpercat[(x,tag,'tot')]=hpercat[(x,tag,base)].Clone('tot_%s_%s'%(d,tag))
                hpercat[(x,tag,'tot')].Add(hpercat[(x,tag,'passPixandpassStrip')])                
                hpercat[(x,tag,'tot')].SetDirectory(0)
                ntot=hpercat[(x,tag,'tot')].Integral()

                continue
                for cat in [base,'passPixandpassStrip']:
                    hpercat[ (x,tag,cat) ].Scale(1./ntot)

            if tag=='':
                print base,hpercat[ ('data','',base) ].Integral()
                print base+'passPixandpassStrip',hpercat[('data','','passPixandpassStrip')].Integral()
                print '\t',hpercat[ ('data','',base) ].Integral()/hpercat[('data','','passPixandpassStrip')].Integral()
            continue


            if tag=='matched' : continue


            p=Plot('%s%s_%s_era%s'%(d,pfix,tag,period),com='13 TeV')
            p.savelog=True
            p.range=[1e-3,1]
            p.doPoissonErrorBars=False
            p.add(hpercat[('data',tag,base)],          title=baseTitle,     color=633,       isData=False, spImpose=False, isSyst=False)
            p.add(hpercat[('data',tag,'passPixandpassStrip')], title='=2 px+strip', color="#fdc086", isData=False, spImpose=False, isSyst=False)
            try:
                hpercat[('mc',tag,'tot')].Scale(1./hpercat[('mc',tag,'tot')].Integral())
                p.add(hpercat[('mc',tag,'tot')],               title='signal',   color=1,         isData=False, spImpose=True,  isSyst=False)
            except:
                pass
            p.show(outDir='./', lumi=37500, noStack=False)

        continue
        #efficiency
        effgr={}
        for x,ms in [('data',20),('mc',24)]:
            effgr[x]=ROOT.TGraphAsymmErrors()
            effgr[x].Divide(hpercat[(x,'matched','tot')],hpercat[(x,'','tot')])
            effgr[x].SetMarkerStyle(ms)
            effgr[x].SetTitle(x)
            effgr[x].SetName(x)
        p=Plot('%s%s_eff_era%s'%(d,pfix,period),com='13 TeV')
        p.ytit='Efficiency'
        p.xtit=hpercat[(x,'matched','tot')].GetXaxis().GetTitle()
        p.savelog=False
        p.range=[1e-3,1]
        #p.add(effgr['mc'],   title='MC',   color=ROOT.kCyan+1,  isData=False, spImpose=True,  isSyst=False)        
        p.add(effgr['data'], title='Data', color=1,             isData=False, spImpose=True,  isSyst=False)        
        p.show(outDir='./', lumi=37500, noStack=True)

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetPaintTextFormat("4.2f");
if len(sys.argv)==1:
    for period in 'BCDEF':
        buildControlPlots(period,excSel=True)
        #buildControlPlots(period,excSel=False)
elif len(sys.argv)==2:
    arg=sys.argv[1]
    buildControlPlots(arg,excSel=True)
    #buildControlPlots(arg,excSel=False)
else:
    data,mc=sys.argv[1:3]
    for period in 'BDCEF':
        drawPlots(ROOT.TFile.Open('%s%s.root'%(data,period)),ROOT.TFile.Open(mc),period)
