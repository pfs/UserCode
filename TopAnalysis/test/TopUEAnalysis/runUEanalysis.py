#!/usr/bin/env python

import glob
import sys
import os
import optparse
import ROOT
import numpy as np
import array as array
import pickle

from UEEventCounter import *
from UEAnalysisHandler import *
from TopLJets2015.TopAnalysis.storeTools import getEOSlslist

#GLOBAL VARIABLES TO DEFINE THE ANALYSIS

BASEQUANTILES=[ 100*x/50 for x in xrange(0,51) ]

"""
determines the resolutions needed to define the migration matrices
"""
def determineSliceResolutions(opt):

    #build the chain
    t=ROOT.TChain('tue')
    for f in opt.input.split(','): 
        if len(f):
            t.AddFile(f)

    #loop the available events and fill resolution arrays for events passing the selection cuts
    varVals={}
    for var in VARS.keys() : varVals[var]=[[],[]]
    ue=UEEventCounter(ptthreshold=[float(x) for x in opt.ptThr.split(',')])
    totalEntries=t.GetEntries()
    for i in xrange(0,totalEntries):
        t.GetEntry(i)
        if i%1000==0 :
            sys.stdout.write('\r [ %d/100 ] done' %(int(float(100.*i)/float(totalEntries))) )
            sys.stdout.flush()

        #count particles in the event
        ue.count(t=t,debug=False,isMC=True)

        #require a pure event selected at reco and gen levels
        if not ue.rec_passSel[0] or not ue.gen_passSel: continue

        #fill resolution arrays
        for var in VARS.keys():
            try:
                val=getattr(t,'gen_%s'%var)
                deltaVal=getattr(t,var)[0]-val
            except:
                val=getattr(ue,'gen_%s'%var) 
                deltaVal=getattr(ue,'rec_%s'%var)[0] - val
            varVals[var][0].append( val )
            varVals[var][1].append( deltaVal )

    #prepare canvas for the resolution plots
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetPalette(ROOT.kTemperatureMap)
    c=ROOT.TCanvas('c','c',550,500)
    c.SetTopMargin(0.05)
    c.SetRightMargin(0.12)
    c.SetBottomMargin(0.1)

    #determine resolution for 10 quantiles
    varResolutions={}
    for var in varVals:

        #use quantiles and extremes to define the 2D resolution histogram
        genvarQ=[]
        if var=='nj': 
            genvarQ = [0,1,2,5]
        else : 
            genvarQ = np.percentile( np.array(varVals[var][0]), BASEQUANTILES )
        dvarQ = np.percentile( np.array(varVals[var][1]), [2.5,97.5])
        h2d=ROOT.TH2F(var,
                      ';%s (gen. level);#Delta(reco-gen);'%VARS[var][0],
                      len(genvarQ)-1,
                      array.array('d',genvarQ),50,dvarQ[0],dvarQ[1])
        h2d.GetZaxis().SetNdivisions(5)
        for i in xrange(0,len(varVals[var][0])):
            h2d.Fill(varVals[var][0][i],varVals[var][1][i])
        
        #fit a gausian to determine resolution per quantile
        resGr=ROOT.TGraphErrors()
        resGr.SetName(var+'_resol')
        resGr.SetMarkerStyle(20)
        for xbin in xrange(1,h2d.GetNbinsX()+1):
            
            tmp=h2d.ProjectionY('tmp',xbin,xbin)

            total=tmp.Integral(0,tmp.GetNbinsX()+1)
            if total==0 : continue

            #require minimum stats for resolution fit
            if total>10 :
                
                mean,sigma=tmp.GetMean(),tmp.GetRMS()
                if not var in ['C','D','aplanarity','sphericity']:
                    tmp.Fit('gaus','MQ+')
                    gaus=tmp.GetFunction('gaus')
                    mean=gaus.GetParameter(1)
                    sigma=gaus.GetParameter(2)
                npts=resGr.GetN()
                resGr.SetPoint(npts,h2d.GetXaxis().GetBinCenter(xbin),mean)
                resGr.SetPointError(npts,0,sigma)

            #normalize entries in quantile
            for ybin in xrange(1,h2d.GetNbinsY()+1):
                val=h2d.GetBinContent(xbin,ybin)
                err=h2d.GetBinError(xbin,ybin)
                h2d.SetBinContent(xbin,ybin,100.*val/total)
                h2d.SetBinError(xbin,ybin,100.*err/total)
            tmp.Delete()

        #show plot
        c.Clear()
        c.SetRightMargin(0.12)
        h2d.Draw('colz')
        h2d.GetYaxis().SetTitleOffset(1.1)
        resGr.Draw('p')
        tex=ROOT.TLatex()        
        tex.SetTextFont(42)
        tex.SetTextSize(0.035)
        tex.SetNDC()
        tex.DrawLatex(0.1,0.96,'#bf{CMS} #it{simulation preliminary} %s'%VARS[var][0])
        tex.DrawLatex(0.8,0.96,'#sqrt{s}=13 TeV')
        logx=True if ('chavg' in var or 'chflux' in var or 'ptttbar' in var or 'sumpt' in var or 'ptpos' in var or 'ptll' in var) else False
        c.SetLogx(logx)
        c.Modified()
        c.Update()
        for ext in ['png','pdf']: c.SaveAs('%s/%s_resol.%s'%(opt.out,var,ext))

        #inclusive resolution
        c.Clear()
        c.SetRightMargin(0.02)
        tmp=h2d.ProjectionY('tmp',1,h2d.GetNbinsX())
        tmp.Scale(1./tmp.Integral())
        tmp.SetMarkerStyle(20)
        tmp.Draw()
        tmp.GetYaxis().SetTitle('PDF')
        tmp.GetYaxis().SetTitleOffset(1.2)
        mean,sigma=tmp.GetMean(),tmp.GetRMS()
        tmp.Fit('gaus','MRQ+','same',mean-sigma,mean+sigma)
        gaus=tmp.GetFunction('gaus')
        mean=gaus.GetParameter(1)
        sigma=gaus.GetParameter(2)
        tex=ROOT.TLatex()        
        tex.SetTextFont(42)
        tex.SetTextSize(0.035)
        tex.SetNDC()
        tex.DrawLatex(0.12,0.9,'#bf{CMS} #it{simulation preliminary} %s'%VARS[var][0])
        tex.DrawLatex(0.84,0.96,'#sqrt{s}=13 TeV')
        tex.DrawLatex(0.12,0.85,'#mu=%3.2f #sigma=%3.2f'%(mean,sigma))
        c.SetLogx(False)
        c.Modified()
        c.Update()
        for ext in ['png','pdf']: c.SaveAs('%s/%s_incresol.%s'%(opt.out,var,ext))


        varResolutions[var]=(h2d,resGr)


    #compute correlations at generator level
    ROOT.gStyle.SetPaintTextFormat("4.0f")
    hcorr=ROOT.TH2F('slicecorr',';Variable; Variable; Correlation (%)',len(varVals),0,len(varVals),len(varVals),0,len(varVals))
    hcorr.GetZaxis().SetNdivisions(5)
    xbin=0
    for var1 in varVals:
        xbin+=1
        hcorr.GetXaxis().SetBinLabel(xbin,VARS[var1][0])
        ybin=0
        for var2 in varVals:
            rho=np.corrcoef(varVals[var1][0],varVals[var2][0])[0][1]
            ybin+=1
            hcorr.GetYaxis().SetBinLabel(ybin,VARS[var2][0])
            hcorr.SetBinContent(xbin,ybin,rho*100)
    c.Clear()
    hcorr.Draw('colztext')
    hcorr.GetYaxis().SetTitleOffset(1.25)
    tex=ROOT.TLatex()        
    tex.SetTextFont(42)
    tex.SetTextSize(0.035)
    tex.SetNDC()
    tex.DrawLatex(0.1,0.96,'#bf{CMS} #it{simulation preliminary}')
    tex.DrawLatex(0.8,0.96,'#sqrt{s}=13 TeV')
    c.Modified()
    c.Update()
    for ext in ['png','pdf']: c.SaveAs('%s/slicecorr.%s'%(opt.out,ext))

    #all done, save to pickle file
    with open(os.path.join(opt.out,'sliceresolutions.pck'), 'w') as cachefile:
        pickle.dump(varResolutions, cachefile, pickle.HIGHEST_PROTOCOL)
        pickle.dump(hcorr,          cachefile, pickle.HIGHEST_PROTOCOL)


"""
defines the variables to slice the UE measurement and saves their reco quantiles
"""
def defineAnalysisBinning(opt):

    #read the resolutions
    varResolutions=None
    with open(os.path.join(opt.out,'sliceresolutions.pck'),'r') as cachefile:
        varResolutions = pickle.load(cachefile)

    #
    # VARIABLES
    # readout resolution curves and determine the bins for the observables
    #
    print 'Defining axes for the slice variables'
    varAxes={}
    for var in varResolutions:

        #special case for jet multiplicity, ptll and ptttbar
        if var=='nj':
            genBin=[0,1,2,3]
            recBin=[0,1,2,3]
        elif var=='ptll' : 
            genBin=[0,50,100,150,300]
            recBin=[0,50,100,150,300]
        elif var in ['ptttbar'] : 
            genBin=[0,40,80,120,300]
            recBin=[0,40,80,120,300]
        else:
             
            #get resolution map and quantiles       
            genBin=[0]
            resolGr=varResolutions[var][1]
            nSigmaForBins=2.5
            if var in ['C','D','aplanarity','sphericity'] : nSigmaForBins=1.5
            #if var in ['chavgpt','chavgpz'] : nSigmaForBins=3
            if VARS[var][1] and var!='chmult' : nSigmaForBins=5

            lastAcceptResol=resolGr.GetErrorY(0)
            for n in xrange(1,resolGr.GetN()-1):

                #center value, bias and resolution (don't correct bias)
                xgen_i,delta_i=ROOT.Double(0),ROOT.Double(0)        
                resolGr.GetPoint(n,xgen_i,delta_i)
                dx=xgen_i-genBin[-1]
                if dx<0 : continue
                if dx<nSigmaForBins*lastAcceptResol : continue
                genBin.append( xgen_i )
                lastAcceptResol=resolGr.GetErrorY(n-1)
            genBin.append( genBin[-1]+lastAcceptResol*2 )

            
            recBin=[0]
            for i in xrange(1,len(genBin)):
                dBin=genBin[i]-genBin[i-1]
                if i>1 : recBin.append( genBin[i-1]+0.25*dBin )
                recBin.append( genBin[i-1]+0.75*dBin )
            recBin.append( genBin[-1]+0.25*(genBin[-1]+genBin[-2]) )

        if var in ['chavgpt','chavgpz']  : 
            genBin[0],recBin[0]=0.9,0.9
            del genBin[1]
            del recBin[1]
        if var in ['chflux','chfluxz']:
            genBin[0],recBin[0]=0.9,0.9
        if var in ['C','D','sphericity'] : 
            genBin[-1],recBin[-1]=1,1
        if var in ['aplanarity']         : 
            genBin[-1],recBin[-1]=0.5,0.5

        #save binning in histos
        varAxes[(var,False)] = ROOT.TAxis(len(genBin)-1,array.array('d',genBin))
        varAxes[(var,False)].SetName('%s_genSlices'%var)
        varAxes[(var,True)]  = ROOT.TAxis(len(recBin)-1,array.array('d',recBin))
        varAxes[(var,True)].SetName('%s_recSlices'%var)
        print '%30s'%('Bin definition for '+var)
        for level in [False,True]:
            nbins=varAxes[(var,level)].GetNbins()
            print '%30s'%'','[%3.1f,%3.1f] / %d @ %s'%(varAxes[(var,level)].GetBinLowEdge(1),
                                                       varAxes[(var,level)].GetBinUpEdge(nbins),
                                                       nbins,
                                                       'rec' if level else 'gen')
            print '%30s'%'',(recBin if level else genBin)
        print ''

    #
    # DEFINE MIGRATION MATRICES, GEN/REC LEVEL HISTOS
    #
    histos={}    
    print 'Saving gen/rec level histos and migration matrix templates'
    print 'Observables:',OBSVARS
    print 'Event axes:',EVAXES
    print 'Slice vars:',SLICEVARS

    #histograms to fill
    for obs in OBSVARS:

        for sliceVar in SLICEVARS+[None]:

            nsliceBins=1
            if not sliceVar is None: nsliceBins=varAxes[(sliceVar,False)].GetNbins()+1

            #check if distributions can also be obtained in different regions
            axes=['inc']
            if sliceVar in EVAXES : 
                if VARS[obs][3]:
                    axes+=[0,1,2]
                
            #create as many histograms as bins in the gen level axis for the slicing variable
            for islice in xrange(0,nsliceBins):
                
                #iterate over the axis
                for axis in axes:

                    nbins={}
                    for level in [False,True]:
                        nbins[level]=varAxes[(obs,level)].GetNbins()

                        #distribution
                        key=(obs,sliceVar,islice,axis,None,level)                
                        name='_'.join(map(str,key))
                        histos[ key ] = ROOT.TH1F(name,name,nbins[level],0,nbins[level])
                        histos[ key ].SetDirectory(0)
                        histos[ key ].Sumw2()

                        if not level: continue

                        #fakes (only for reco)
                        key=(obs,sliceVar,islice,axis,'fakes',level)                
                        name='_'.join(map(str,key))
                        histos[ key ] = ROOT.TH1F(name,name,nbins[level],0,nbins[level])
                        histos[ key ].SetDirectory(0)
                        histos[ key ].Sumw2()

                        #experimental systs (only for reco)
                        key=(obs,sliceVar,islice,axis,'syst',level)
                        name='_'.join(map(str,key))
                        histos[ key ] = ROOT.TH2F(name,name,nbins[level],0,nbins[level],len(SYSTS),-0.5,len(SYSTS)-0.5)
                        histos[ key ].SetDirectory(0)
                        histos[ key ].Sumw2()
                        for ybin in xrange(0,len(SYSTS)):
                            histos[ key ].GetYaxis().SetBinLabel(ybin+1,SYSTS[ybin][0])

                    #migration matrix (1 per syst variation)
                    for i in xrange(0,len(SYSTS)):
                        key=(obs,sliceVar,islice,axis,i,'mig')
                        name='_'.join(map(str,key))
                        histos[ key ] = ROOT.TH2F(name,name,nbins[False],0,nbins[False],nbins[True],0,nbins[True])
                        histos[ key ].SetTitle( SYSTS[i][0] )
                        histos[ key ].SetDirectory(0)
                        histos[ key ].Sumw2()

    #all done, save to pickle file
    with open(os.path.join(opt.out,'analysisaxiscfg.pck'), 'w') as cachefile:
        pickle.dump(varAxes,     cachefile, pickle.HIGHEST_PROTOCOL)
    with open(os.path.join(opt.out,'analysiscfg.pck'), 'w') as cachefile:
        pickle.dump(varAxes,     cachefile, pickle.HIGHEST_PROTOCOL)
        pickle.dump(histos,      cachefile, pickle.HIGHEST_PROTOCOL)
    print 'Analysis cfg saved in',os.path.join(opt.out,'analysiscfg.pck')
    print 'Only axis definitions also available in',os.path.join(opt.out,'analysisaxiscfg.pck')


"""
loops over a set of files with common name to fill the migration matrices
"""        
def runUEAnalysis(inF,outF,cfgDir):

    print '[runAnalysis] %s -> %s'%(inF,outF)
    
    #configure from pickle file
    ueHandler=UEAnalysisHandler(os.path.join(cfgDir,'analysiscfg.pck'))

    varList=SYSTS
    if not 'MC13TeV_TTJets' in inF : varList=[('',   0,0,False)]
    isMC=True if 'MC13TeV' in inF else False
    ALLSLICEVARS=[None]+SLICEVARS

    #loop over the tree to fill histos
    ue=UEEventCounter(EVAXES,varList=varList)
    t=ROOT.TChain('tue')
    t.AddFile(inF)
    totalEntries=t.GetEntries()
    for i in xrange(0,totalEntries):

        t.GetEntry(i)

        if i%100==0 :
            sys.stdout.write('\r [ %d/100 ] done' %(int(float(100.*i)/float(totalEntries))) )
            sys.stdout.flush()

        #count particles
        ue.count(t=t,isMC=isMC)

        #
        for ivar in xrange(0,len(varList)):
            _,_,varIdx,_ = varList[ivar]

            #fill histos (if sliceVar is None, non-sliced histos are filled)
            for sliceVar in ALLSLICEVARS:
            
                #check if sliceVar is non null and configure values to use
                sliceVarVals=None
                try:
                    sliceVarVals = (sliceVar, getattr(t,'gen_%s'%sliceVar), getattr(t,sliceVar)[varIdx] )
                except:
                    try:
                        sliceVarVals = (sliceVar, getattr(ue,'gen_%s'%sliceVar), getattr(ue,'rec_%s'%sliceVar)[ivar] )
                    except:
                        pass
                
                #loop over UE observables
                for obs in OBSVARS:
                    ueHandler.fillHistos(sliceVarVals=sliceVarVals, obs=obs, ue=ue, ivar=ivar)


    #save histos to ROOT file
    fOut=ROOT.TFile.Open(outF,'RECREATE')
    for h in ueHandler.histos:        
        ueHandler.histos[h].SetDirectory(fOut)
        ueHandler.histos[h].Write()
    fOut.Close()

"""
Wrapper for when the analysis is run in parallel
"""
def runUEAnalysisPacked(args):
    try:
        fileName,outDir,varIdx,wgtIdx,cfgDir=args
        runUEAnalysis(fileName,outDir,varIdx,wgtIdx,cfgDir)
    except : 
        print 50*'<'
        print "  Problem with", name, "continuing without"
        print 50*'<'
        return False


def main():

    ROOT.gROOT.SetBatch(True)

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--in',    dest='input',  help='input',                       default='MC13TeV_TTJets_dilpowheg_0.root',   type='string')
    parser.add_option(      '--ptThr', dest='ptThr',  help='ptThreshold gen,reco',        default='1.0,0.9',   type='string')
    parser.add_option('-s', '--step',  dest='step',   help='step',                        default=1,   type=int)
    parser.add_option('-j', '--jobs',  dest='jobs',   help='jobs to run in parallel',     default=1,   type=int)
    parser.add_option('-o', '--out',   dest='out',    help='output',                      default='./UEanalysis',   type='string')
    parser.add_option(      '--only',  dest='only',   help='csv list of tags to process', default='',  type='string')
    parser.add_option('-q', '--queue', dest='queue',  help='Batch queue to use [default: %default]', default='local')
    (opt, args) = parser.parse_args()

    onlyList=opt.only.split('v')
    os.system('mkdir -p %s'%opt.out)

    if opt.step==0:
        determineSliceResolutions(opt)

    if opt.step==1:
        defineAnalysisBinning(opt)

    if opt.step==2:

        #prepare output
        outDir=opt.out+'/analysis/Chunks'
        os.system('mkdir -p %s'%outDir)

        #create the tasklist
        file_list=[]
        if os.path.isdir(opt.input):
            for file_path in os.listdir(opt.input):
                if file_path.endswith('.root'):
                    file_list.append(os.path.join(opt.input,file_path))
        elif opt.input.startswith('/store/'):
            file_list = getEOSlslist(opt.input)
        elif '.root' in opt.input:
            file_list.append(opt.input)


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
            tasklist.append((filename,'%s/%s'%(outDir,baseFileName),opt.out))

        #run jobs locally
        if opt.queue=='local':
            if opt.jobs>1:
                print ' Running %d jobs in %d threads' % (len(tasklist), opt.jobs)
                import multiprocessing as MP
                pool = MP.Pool(opt.jobs)
                pool.map(runUEAnalysisPacked,tasklist)
            else:
                for fileName,outfile,cfgDir in tasklist:
                    runUEAnalysis(fileName,outfile,cfgDir)

        #submit jobs
        else:
            print 'Running %d jobs to the batch'%len(tasklist)
            cmsswBase=os.environ['CMSSW_BASE']
            FarmDirectory='%s/FARM-UEANA'%cmsswBase
            os.system('mkdir -p %s'%FarmDirectory)
            print 'Scripts and logs will be available @ %s'%FarmDirectory

            condorScript='%s/condor.sub'%FarmDirectory
            with open (condorScript,'w') as condor:

                condor.write('executable = {0}/$(jobName).sh\n'.format(FarmDirectory))
                condor.write('output     = {0}/output_$(jobName).out\n'.format(FarmDirectory))
                condor.write('error      = {0}/output_$(jobName).err\n'.format(FarmDirectory))
                condor.write('+JobFlavour = workday\n')

                for fileName,_,cfgDir in tasklist:
                    
                    jobName='%s'%(os.path.splitext(os.path.basename(fileName))[0])

                    jobScript='%s/%s.sh'%(FarmDirectory,jobName)
                    with open(jobScript,'w') as job:
                        job.write('#!/bin/bash\n')
                        job.write('WORKDIR=`pwd`\n')
                        job.write('echo "Working directory is ${WORKDIR}"\n')
                        job.write('cd %s\n'%cmsswBase)
                        job.write('eval `scram r -sh`\n')
                        job.write('cd ${WORKDIR}\n')
                        job.write('python {0}/src/TopLJets2015/TopAnalysis/test/TopUEAnalysis/runUEanalysis.py -o {1} -q local -s 2 -i {2}\n'.format(cmsswBase,cfgDir,fileName))
                        job.write('echo "All done"\n')

                    os.system('chmod u+x %s'%jobScript)
                    condor.write('jobName=%s\n'%jobName)
                    condor.write('queue 1\n')

                #localRun='python %s/src/TopLJets2015/TopAnalysis/test/TopUEAnalysis/runUEanalysis.py -i %s -o %s -q local -s 2'%(cmsswBase,fileName,cfgDir)
                #cmd='bsub -q %s %s/src/TopLJets2015/TopAnalysis/scripts/wrapLocalAnalysisRun.sh \"%s\"' % (opt.queue,cmsswBase,localRun)
                #os.system(cmd)
            os.system('condor_submit %s'%condorScript)


"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
