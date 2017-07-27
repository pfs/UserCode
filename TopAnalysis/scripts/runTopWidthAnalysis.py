#! /usr/bin/env python
import ROOT
import optparse
import json
import sys
import os
from TopLJets2015.TopAnalysis.storeTools import getEOSlslist
import MT2Calculator

"""
Take the ratio of two Breit-Wigner functions at fixed mass as a reweighting factor
"""
def weightTopWidth(tmassList,bwigner,targetWidth,origWidth=1.324):
    bwigner.FixParameter(2,origWidth)
    origNorm=bwigner.Integral(0,300)

    bwigner.FixParameter(2,targetWidth)
    targetNorm=bwigner.Integral(0,300)

    wgt=1.0
    for m in tmassList:
        bwigner.FixParameter(2,origWidth)
        origVal=bwigner.Eval(m)
        bwigner.FixParameter(2,targetWidth)
        targetVal=bwigner.Eval(m)
        wgt *= (targetVal/targetNorm) / (origVal/origNorm)
    return wgt

"""
Analysis loop
"""
def runTopWidthAnalysis(fileName,
                        outFileName,
                        widthList=[0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.5,4.0],
                        systs=['',
                               'puup','pudn',
                               'btagup','btagdn',
                               'ltagup','ltagdn',
                               'jerup','jerdn',
                               'eesup','eesdn',
                               'mesup','mesdn',
                               'trigup','trigdn',
                               'eselup','eseldn',
                               'mselup','mseldn',
                               'bfragup','bfragdn','petersfrag',
                               'semilepup','semilepdn'],
                        thsysts=[]):

    print '....analysing',fileName,'with output @',outFileName

    #open file
    isData=False if 'MC13TeV' in fileName else True
    puNormSF=[1.0,1.0,1.0]
    if not isData:
        fIn=ROOT.TFile.Open(fileName)
        puCorrH=fIn.Get('puwgtctr')
        nonWgt=puCorrH.GetBinContent(1)
        for xbin in xrange(2,5):
            wgt=puCorrH.GetBinContent(xbin)
            if wgt>0 : puNormSF[xbin-2]=nonWgt/wgt
        fIn.Close()

    #instantiate the tree
    tree=ROOT.TChain('twev')
    tree.AddFile(fileName)

    #check if this is data beforehand
    if isData:
        widthList=[1.0]
        systs=['']
    else:
        
        #add jet energy scale
        tree.GetEntry(0)
        for i in xrange(0,tree.j_jes[0].size()):
            systs += ['jesup_%d'%i,'jesdn_%d'%i]

        #generator level systematics for ttbar
        if 'TTJets' in fileName or 'SingleT_tW' in fileName:
            thsysts += ['topptup','topptdn','nloproddec']
            ngenWgts=tree.nw-15
            print 'Adding',ngenWgts,'gen-level weights for systs'
            for i in xrange(1,ngenWgts):  thsysts += ['gen%d'%i]
        else:
            widthList=[1.0]

    #powheg samples have running width 2%/1GeV
    #see  https://github.com/jfernan2/genproductions/tree/8b309da3427fb5fdcc2dadc1860774f20adb517c/bin/Powheg/production
    smWidth=1.31
    smMass=172.5
    if '166v5' in fileName : smMass,smWidth=166.5,1.16
    if '169v5' in fileName : smMass,smWidth=169.5,1.23
    if '171v5' in fileName : smMass,smWidth=171.5,1.28    
    if '173v5' in fileName : smMass,smWidth=173.5,1.34
    if '175v5' in fileName : smMass,smWidth=175.5,1.39
    if '178v5' in fileName : smMass,smWidth=178.5,1.48

    #define the relativistic Breit-Wigner function
    bwigner=ROOT.TF1('bwigner',
                     '[0]*([1]*[2]*sqrt([1]*[1]*([1]*[1]+[2]*[2]))/sqrt([1]*[1]+sqrt([1]*[1]*([1]*[1]+[2]*[2]))))/(TMath::Power(x*x-[1]*[1],2)+TMath::Power([1]*[2],2))',
                     0,300)
    bwigner.SetParName(0,"N")
    bwigner.SetParameter(0,1.0)
    bwigner.SetParName(1,"m_{t}")
    bwigner.FixParameter(1,smMass)
    bwigner.SetParName(2,"#Gamma_{t}")
    bwigner.FixParameter(2,smWidth)

    #read the MCFM NLO(prod+dec)/NLO(prod) weights
    mcfmIn=ROOT.TFile.Open('${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/MCFM_todk2tota.root')
    todk2totaGr={'wro':ROOT.TGraph( mcfmIn.Get('wrong_combinations/dist_NLO prod+decay_ratio') ),
                 'cor':ROOT.TGraph( mcfmIn.Get('correct_combinations/dist_NLO prod+decay_ratio') ) }
    mcfmIn.Close()

    #book histograms
    observablesH={}

    #MC truth control histograms
    for w in widthList:
        var='tmass_w%d'%int(100*w)
        observablesH[var]=ROOT.TH1F(var,';Top quark mass [GeV];Top quarks',100,150,200)
        for assig in ['cor','wro']:
            var='%sgenmlbvsmtop_w%d'%(assig,int(100*w))
            observablesH[var]=ROOT.TH2F(var,';Mass(lepton,jet) [GeV];Top mass;l+j pairs',30,0,300,100,150,200)

    #RECO level histograms
    for j in ['EE','MM','EM']:
        var=j+'_evcount'
        observablesH[var]=ROOT.TH1F(var,';Category;Events',2,0,2)
        observablesH[var].GetXaxis().SetBinLabel(1,'=1b')
        observablesH[var].GetXaxis().SetBinLabel(2,'#geq2b')
        for b in ['1b','2b']:
            var=j+b+'_ptlb'
            observablesH[var]=ROOT.TH1F(var,';p_{T}(lepton,jet) [GeV];l+j pairs',30,0,300)
            var=j+b+'_drlb'
            observablesH[var]=ROOT.TH1F(var,';#Delta R(lepton,jet);l+j pairs',30,0,6)
            var=j+b+'_dphilb'
            observablesH[var]=ROOT.TH1F(var,';#Delta #phi(lepton,jet)[rad];l+j pairs',30,-3.15,3.15)
            var=j+b+'_met'
            observablesH[var]=ROOT.TH1F(var,';Missing transverse energy [GeV];Events',30,0,300)
            var=j+b+'_mll'
            observablesH[var]=ROOT.TH1F(var,';Dilepton invariant mass [GeV];Events',30,0,300)
            var=j+b+'_njets'
            observablesH[var]=ROOT.TH1F(var,';Jet multiplicity;Events',6,2,8)
            var=j+b+'_njetsnonboost'
            observablesH[var]=ROOT.TH1F(var,';Jet multiplicity;Events',6,2,8)
            var=j+b+'_njetsboost'
            observablesH[var]=ROOT.TH1F(var,';Jet multiplicity;Events',6,2,8)

            for i in ['lowpt','highpt']:

                for w in widthList:
                    var=j+b+i+'_incmlb_w%d'%int(100*w)
                    observablesH[var]=ROOT.TH1F(var,';Mass(lepton,jet) (Inclusive) [GeV];l+j pairs',30,0,300)

                    var=j+b+i+'_incmlb_w%d_exp'%int(100*w)
                    nsysts=len(systs)
                    observablesH[var]=ROOT.TH2F(var,';Mass(lepton,jet) (Inclusive) [GeV];Systematics;l+j pairs',30,0,300,nsysts+1,0,nsysts+1)
                    for isyst in xrange(0,nsysts):
                        observablesH[var].GetYaxis().SetBinLabel(isyst+1,systs[isyst])

                    var=j+b+i+'_incmlb_w%d_gen'%int(100*w)
                    nsysts=len(thsysts)
                    observablesH[var]=ROOT.TH2F(var,';Mass(lepton,jet) (Inclusive) [GeV];Systematics;l+j pairs',30,0,300,nsysts+1,0,nsysts+1)
                    for isyst in xrange(0,nsysts):
                        observablesH[var].GetYaxis().SetBinLabel(isyst+1,thsysts[isyst])

                    if w!=1.0 : continue
                    var=j+b+i+'_pairing'
                    observablesH[var]=ROOT.TH1F(var,';Pairing;l+j pairs',2,0,2)
                    observablesH[var].GetXaxis().SetBinLabel(1,'correct')
                    observablesH[var].GetXaxis().SetBinLabel(2,'wrong')

    for var in observablesH:
        observablesH[var].SetDirectory(0)
        observablesH[var].Sumw2()

    #loop over events in the tree and fill histos
    totalEntries=tree.GetEntries() 
    for i in xrange(0,totalEntries):

        tree.GetEntry(i)

        if i%100==0 : sys.stdout.write('\r [ %d/100 ] done' %(int(float(100.*i)/float(totalEntries))) )

        #focus on dilepton analysis for the moment
        evcat=''
        if abs(tree.cat)==11*11 : evcat='EE'
        if abs(tree.cat)==11*13 : evcat='EM'
        if abs(tree.cat)==13*13 : evcat='MM'
        if evcat==''            : continue

        #base event weight
        baseEvWeight=puNormSF[0]*tree.weight[0]

        #determine weighting factors for the width
        tops={}
        tmassList=[]
        bpartonList=[]
        for it in xrange(0,tree.nt):
            tid=tree.t_id[it]
            if abs(tid)==5:
                bpartonList.append( ROOT.TLorentzVector() )
                bpartonList[-1].SetPtEtaPhiM(tree.t_pt[it],tree.t_eta[it],tree.t_phi[it],tree.t_m[it])
            if it<2 and abs(tid)==6:
                tops[ tid ] = ROOT.TLorentzVector()
                tops[ tid ].SetPtEtaPhiM(tree.t_pt[it],tree.t_eta[it],tree.t_phi[it],tree.t_m[it])
                tmassList.append( tops[tid].M() )
        widthWeight={}
        for w in widthList:
            widthWeight[w]=weightTopWidth(tmassList,bwigner,w*smWidth,smWidth)

            #paranoid check for the reweighted mass based on the Breit-Wigner
            var='tmass_w%d'%int(100*w)
            for mtop in tmassList:
                observablesH[var].Fill(mtop,baseEvWeight*widthWeight[w])


        #preselect the b-jets (central, b-tag up, b-tag dn, l-tag up, l-tag dn, jer up, jer dn, jes_1 up, jes_1 dn, ...)
        bjets     = [ [], [], [], [], [], [], [] ]
        otherjets = [ [], [], [], [], [], [], [] ]
        for ijs in xrange(0,tree.j_jes[0].size()):
            bjets     += [ [],[] ]
            otherjets += [ [],[] ]

        #define b-tag variations (bit to use, do heavy flavour shift)
        #if heavy flavour shift = None the standard b-tag decision is used
        btagVars=[ (0,None), (1,True), (2,True), (1,False), (2,False) ]
        for ij in xrange(0,tree.nj):

            jp4=ROOT.TLorentzVector()
            jp4.SetPtEtaPhiM(tree.j_pt[ij],tree.j_eta[ij],tree.j_phi[ij],tree.j_m[ij])

            for ibvar in xrange(0,len(btagVars)):

                #reset b-tag bit to 0 if flavour is not the one to vary
                ibit, shiftHeavyFlav = btagVars[ibvar]
                if shiftHeavyFlav is not None:
                    hadFlav=abs(tree.gj_hadflav[ij])
                    if shiftHeavyFlav :
                        if hadFlav!=4 and hadFlav!=5 : ibit=0
                    else:
                        if hadFlav==4 or hadFlav==5 : ibit=0

                #get b-tag decision
                btagVal=((tree.j_btag[ij] >> ibit) & 0x1)

                if btagVal > 0:
                    
                    bjets[ibvar].append( (ij,jp4) )

                    #if standard b-tagging decision is being used, store JES/JER variations
                    if shiftHeavyFlav is not None : continue

                    jres=ROOT.TMath.Abs(1-tree.j_jer[ij])
                    bjets[5].append( (ij,ROOT.TLorentzVector(jp4)*(1+jres)) )
                    bjets[6].append( (ij,ROOT.TLorentzVector(jp4)*(1-jres)) )

                    for ijs in xrange(0,tree.j_jes[ij].size()):
                        jscale=abs(tree.j_jes[ij][ijs])
                        bjets[7+ijs*2].append( (ij,ROOT.TLorentzVector(jp4)*(1+jscale)) )
                        bjets[8+ijs*2].append( (ij,ROOT.TLorentzVector(jp4)*(1-jscale)) )

                elif btagVal == 0 :
                    otherjets[ibvar].append( (ij,jp4) )
                    
                    #if standard b-tagging decision is being used, store JES/JER variations
                    if shiftHeavyFlav is not None : continue

                    jres=ROOT.TMath.Abs(1-tree.j_jer[ij])
                    otherjets[5].append( (ij,ROOT.TLorentzVector(jp4)*(1+jres)) )
                    otherjets[6].append( (ij,ROOT.TLorentzVector(jp4)*(1-jres)) )

                    for ijs in xrange(0,tree.j_jes[ij].size()):
                        jscale=abs(tree.j_jes[ij][ijs])
                        otherjets[7+ijs*2].append( (ij,ROOT.TLorentzVector(jp4)*(1+jscale)) )
                        otherjets[8+ijs*2].append( (ij,ROOT.TLorentzVector(jp4)*(1-jscale)) )

        #build the dilepton
        dilepton=ROOT.TLorentzVector()
        for il in xrange(0,2):
            stdlp4=ROOT.TLorentzVector()
            stdlp4.SetPtEtaPhiM(tree.l_pt[il],tree.l_eta[il],tree.l_phi[il],tree.l_m[il])
            dilepton+=stdlp4

        #quarkonia/low mass dilepton veto
        if dilepton.M()<15: continue

        #b-tagging category
        nbtags=len(bjets[0])
        if nbtags>2 : nbtags=2
        btagcat='%db'%nbtags

        #full Mll variable is used to re-scale DY (Rin/Rout)
        if nbtags>0:
            var=evcat+btagcat+'_mll'
            observablesH[var].Fill(dilepton.M(),baseEvWeight)

        #remove Z/quarkonia candidates from this point forward
        if abs(tree.cat)==11*11 or abs(tree.cat)==13*13:
            if ROOT.TMath.Abs(dilepton.M()-91)<15 : continue

        #control variables
        if nbtags>0:
            var=evcat+"_evcount"
            observablesH[var].Fill(nbtags-1,baseEvWeight)
            var=evcat+btagcat+'_met'
            observablesH[var].Fill(tree.met_pt,baseEvWeight)
            var=evcat+btagcat+'_njets'
            observablesH[var].Fill(tree.nj,baseEvWeight)

            #check if event has a boosted pair
            hasBoostedPair=False
            for il in xrange(0,2):
                lp4=ROOT.TLorentzVector()
                lp4.SetPtEtaPhiM(tree.l_pt[il],tree.l_eta[il],tree.l_phi[il],tree.l_m[il])
                for ib in xrange(0,nbtags):
                    ij,jp4 = bjets[0][ib]
                    ptlb=(lp4+jp4).Pt()
                    if ptlb>100 : hasBoostedPair=True
            var=evcat+btagcat+'_njetsboost' if hasBoostedPair else evcat+btagcat+'_njetsnonboost'
            observablesH[var].Fill(tree.nj,baseEvWeight)

        #pair with the leptons and for all possible variations
        for il in xrange(0,2):

            stdlp4=ROOT.TLorentzVector()
            stdlp4.SetPtEtaPhiM(tree.l_pt[il],tree.l_eta[il],tree.l_phi[il],tree.l_m[il])
            lscale=tree.l_les[il]

            for isyst in xrange(0,len(systs)):

                s=systs[isyst]

                #event weight
                evWeight=puNormSF[0]*tree.weight[0]

                #base lepton kinematics
                lp4=ROOT.TLorentzVector(stdlp4)

                #experimental uncertainties
                ijhyp=0
                if s=='btagup' : ijhyp=1
                if s=='btagdn' : ijhyp=2
                if s=='ltagup' : ijhyp=3
                if s=='ltagdn' : ijhyp=4
                if s=='jerup'  : ijhyp=5
                if s=='jerdn'  : ijhyp=6                
                if 'jes' in s:
                    ijs=int(s.split('_')[1])
                    idx=7 if 'dn' in s else 8
                    ijhyp=idx+2*ijs
                if s=='eesup' and abs(tree.l_id[il])==11:  lp4 *= (1.0+lscale)
                if s=='eesdn' and abs(tree.l_id[il])==11 : lp4 *= (1.0-lscale)
                if s=='mesup' and abs(tree.l_id[il])==13:  lp4 *= (1.0+lscale)
                if s=='mesdn' and abs(tree.l_id[il])==13 : lp4 *= (1.0-lscale)
                if s=='puup'   : evWeight=puNormSF[1]*tree.weight[1]
                if s=='pudn'   : evWeight=puNormSF[2]*tree.weight[2]
                if s=='trigup':  evWeight=puNormSF[0]*tree.weight[3]
                if s=='trigdn':  evWeight=puNormSF[0]*tree.weight[4]
                if s=='eselup':  evWeight=puNormSF[0]*tree.weight[5]
                if s=='eseldn':  evWeight=puNormSF[0]*tree.weight[6]
                if s=='mselup':  evWeight=puNormSF[0]*tree.weight[7]
                if s=='mseldn':  evWeight=puNormSF[0]*tree.weight[8]
                if s=='semilepdn' : evWeight=puNormSF[0]*tree.weight[10]
                if s=='semilepup' : evWeight=puNormSF[0]*tree.weight[11]
                if s=='bfragup'   : evWeight=puNormSF[0]*tree.weight[12]
                if s=='bfragdn'   : evWeight=puNormSF[0]*tree.weight[13]
                if s=='petersfrag': evWeight=puNormSF[0]*tree.weight[14]
                
                #require two jets
                njets=len(bjets[ijhyp])+len(otherjets[ijhyp])
                if njets<2 : continue

                #btag hypothesis
                nbtags=len(bjets[ijhyp])
                if nbtags<1 : continue
                if nbtags>2 : nbtags=2
                btagcat='%db'%nbtags

                for ib in xrange(0,nbtags):
                    
                    ij,jp4 = bjets[ijhyp][ib]

                    #RECO kinematics of the l,b system
                    mlb=(lp4+jp4).M()
                    ptlb=(lp4+jp4).Pt()
                    ptCat='lowpt' if ptlb<100 else 'highpt'
                    drlb=lp4.DeltaR(jp4)
                    dphilb=lp4.DeltaPhi(jp4)
                    
                    #fill the nominal histos
                    if isyst==0:
                        var=evcat+btagcat
                        observablesH[var+'_ptlb'].Fill(ptlb,evWeight)
                        observablesH[var+'_drlb'].Fill(drlb,evWeight)
                        observablesH[var+'_dphilb'].Fill(dphilb,evWeight)

                    #observable and associated systematics
                    for w in widthList:

                        if isyst==0:
                            var=evcat+btagcat+ptCat+'_incmlb_w%d'%int(100*w)
                            observablesH[var].Fill(mlb,evWeight*widthWeight[w])

                        var=evcat+btagcat+ptCat+'_incmlb_w%d_exp'%int(100*w)
                        observablesH[var].Fill(mlb,isyst,evWeight*widthWeight[w])
                        
                    #for the nominal variation do MC truth and the theory systematics, if available
                    if isyst!=0: continue

                    pairFullyMatchedAtGen = (tree.gl_id[il]!=0 and abs(tree.gj_flav[ij])==5)
                    assignmentType,tmass,genmlb,genmlb_parton=1,0.0,0.0,-1.0
                    if pairFullyMatchedAtGen and tree.nt>0:

                        #MC truth  kinematkcs
                        glp4=ROOT.TLorentzVector()
                        glp4.SetPtEtaPhiM(tree.gl_pt[il],tree.gl_eta[il],tree.gl_phi[il],tree.gl_m[il])
                        gjp4=ROOT.TLorentzVector()
                        gjp4.SetPtEtaPhiM(tree.gj_pt[ij],tree.gj_eta[ij],tree.gj_phi[ij],tree.gj_m[ij])
                        genmlb=(glp4+gjp4).M()
                        
                        for ibp in xrange(0,len(bpartonList)):
                            dR=gjp4.DeltaR(bpartonList[ibp])
                            if dR>0.5 : continue
                            genmlb_parton=(glp4+bpartonList[ibp]).M()

                        #correctness of the assignment can be checked by op. charge
                        if tree.gl_id[il]*tree.gj_flav[ij]<0 : assignmentType=0

                        #top mass (parton level)
                        try:
                            if tree.gl_id[il]<0 : tmass=tops[6].M()
                            else                : tmass=tops[-6].M()
                        except:
                            pass

                    #save MC truth distribution
                    if pairFullyMatchedAtGen:
                        for w in widthList:
                            var='cor' if assignmentType==0 else 'wro'
                            var+='genmlbvsmtop_w%d'%int(100*w)
                            observablesH[var].Fill(genmlb,tmass,evWeight*widthWeight[w])

                    #emulate reweighting to NLO prod+dec based on MCFM calculations
                    pairWeightAtNLO=1.0
                    if pairFullyMatchedAtGen and genmlb_parton>0 :
                        if assignmentType==0 : pairWeightAtNLO=todk2totaGr['cor'].Eval(ROOT.TMath.Min(300.,genmlb_parton))
                        else :                 pairWeightAtNLO=todk2totaGr['wro'].Eval(ROOT.TMath.Min(300.,genmlb_parton))
                            
                    #theory systematics
                    for ith in xrange(0,len(thsysts)):                            
                        s=thsysts[ith]
                        thEvWeight=evWeight
                        if s=='topptup': thEvWeight=puNormSF[0]*tree.weight[9]
                        if s=='topptdn': thEvWeight=puNormSF[0]*tree.weight[9]
                        if 'gen' in s  : thEvWeight=puNormSF[0]*tree.weight[10+int(s[3:])]
                        if s=='nloproddec'  : thEvWeight = evWeight*pairWeightAtNLO
                        for w in widthList:
                            var=evcat+btagcat+ptCat+'_incmlb_w%d_gen'%int(100*w)
                            observablesH[var].Fill(mlb,ith,thEvWeight*widthWeight[w])

    #save results
    fOut=ROOT.TFile.Open(outFileName,'RECREATE')
    for var in observablesH: observablesH[var].Write()
    fOut.Close()


"""
Wrapper for when the analysis is run in parallel
"""
def runTopWidthAnalysisPacked(args):
    try:
        fileNames,outFileName=args
        runTopWidthAnalysis(fileNames,outFileName)
    except : # ReferenceError:
        print 50*'<'
        print "  Problem with", name, "continuing without"
        print 50*'<'
        return False

"""
Create analysis tasks
"""
def createAnalysisTasks(opt):

    onlyList=opt.only.split(',')

    ## Local directory
    file_list=[]
    if os.path.isdir(opt.input):
        for file_path in os.listdir(opt.input):
            if file_path.endswith('.root'):
                file_list.append(os.path.join(opt.input,file_path))
    elif opt.input.startswith('/store/'):
        file_list = getEOSlslist(opt.input)
    elif '.root' in opt.input:
        file_list.append(opt.input)

    #list of files to analyse
    tasklist=[]
    for filename in file_list:
        baseFileName=os.path.basename(filename)
        tag,ext=os.path.splitext(baseFileName)
        if len(onlyList)>0:
            processThis=False
            for filtTag in onlyList:
                if filtTag in tag:
                    processThis=True
            if not processThis : continue
        tasklist.append((filename,'%s/%s'%(opt.output,baseFileName)))

    #loop over tasks
    if opt.queue=='local':
        if opt.jobs>1:
            print ' Submitting jobs in %d threads' % opt.jobs
            import multiprocessing as MP
            pool = MP.Pool(opt.jobs)
            pool.map(runTopWidthAnalysisPacked,tasklist)
        else:
            for fileName,outFileName in tasklist:
                runTopWidthAnalysis(fileName,outFileName)
    else:
        cmsswBase=os.environ['CMSSW_BASE']
        FarmDirectory = '%s/%s'%(cmsswBase,opt.farm)
        os.system('mkdir -p %s'%FarmDirectory)
        
        with open ('%s/condor.sub'%FarmDirectory,'w') as condor:
            condor.write('executable = {0}/$(jobName).sh\n'.format(FarmDirectory))
            condor.write('output     = {0}/output_$(jobName).out\n'.format(FarmDirectory))
            condor.write('error      = {0}/output_$(jobName).err\n'.format(FarmDirectory))
            condor.write('+JobFlavour = workday\n')
            condor.write('RequestCpus = 4\n')

            for fileName,_ in tasklist:
                jobName='%s'%(os.path.splitext(os.path.basename(fileName))[0])
                jobScript='%s/%s.sh'%(FarmDirectory,jobName)
                with open(jobScript,'w') as job:
                    job.write('#!/bin/bash\n')
                    job.write('WORKDIR=`pwd`\n')
                    job.write('echo "Working directory is ${WORKDIR}"\n')
                    job.write('cd %s\n'%cmsswBase)
                    job.write('eval `scram r -sh`\n')
                    job.write('cd ${WORKDIR}\n')
                    job.write('python {0}/src/TopLJets2015/TopAnalysis/scripts/runTopWidthAnalysis.py -o {1} -q local -i {2}\n'.format(cmsswBase,opt.output,fileName))
                    job.write('echo "All done"\n')

                os.system('chmod u+x %s'%jobScript)
                condor.write('jobName=%s\n'%jobName)
                condor.write('queue 1\n')

            os.system('condor_submit %s/condor.sub'%FarmDirectory)


"""
steer
"""
def main():
	usage = 'usage: %prog [options]'
	parser = optparse.OptionParser(usage)
	parser.add_option('-i', '--input',
                          dest='input',
                          default='/afs/cern.ch/user/p/psilva/work/TopWidth',
                          help='input directory with the files [default: %default]')
	parser.add_option('--jobs',
                          dest='jobs',
                          default=1,
                          type=int,
                          help='# of jobs to process in parallel the trees [default: %default]')
	parser.add_option('--only',
                          dest='only',
                          default='',
                          type='string',
                          help='csv list of tags to process')
	parser.add_option('-o', '--output',
                          dest='output',
                          default='analysis',
                          help='Output directory [default: %default]')
	parser.add_option('--farm',
                          dest='farm',
                          default='TOP17010ANA',
                          help='farm directory name [default: %default]')
	parser.add_option('-q', '--queue',
                          dest='queue',
                          default='local',
                          help='Batch queue to use [default: %default]')
	(opt, args) = parser.parse_args()

        ROOT.FWLiteEnabler.enable()
	os.system('mkdir -p %s' % opt.output)

        createAnalysisTasks(opt)


if __name__ == "__main__":
	sys.exit(main())
