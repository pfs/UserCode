#! /usr/bin/env python
import ROOT
import optparse
import json
import sys
import os
import time
from array import *
from math import sqrt
from math import isnan

flavorMap = {'all': -1, 'light': 1, 'bottom': 5, 'gluon': 0}

def fillHist(hmap, tree, nweights):
    for event in tree:
        weights = []
        for i in range(nweights+1):
            if not isnan(event.weight[i]):
                weights.append(event.weight[i])
            else: weights.append(0.)
        for obs,hmap_reco in hmap.iteritems():
            for reco,hmap_flavor in hmap_reco.iteritems():
                for flavor,hmap_weight in hmap_flavor.iteritems():
                    iflavor = flavorMap[flavor]
                    for j in range(event.nj):
                        if (event.reco_sel != 1
                            or (iflavor > -1 and event.j_flavor[j] != iflavor)):
                            valReco = -1
                        else:
                            valReco = eval('event.j_'+obs+'_'+reco)[j]
                        g = event.j_gj[j]
                        if g >= 0: #match
                            if (event.gen_sel != 1
                                or (iflavor > -1 and event.gj_flavor[g] != iflavor)):
                                valGen = -1
                            else:
                                valGen = eval('event.gj_'+obs+'_'+reco)[g]
                        else: valGen = -1
                        for weightindex,h in hmap_weight.iteritems():
                            weight = weights[0]
                            if (weightindex > 0): weight *= weights[weightindex]
                            h.Fill(valGen, valReco, weight)
                    for g in range(event.ngj):
                        if event.gj_j[g] < 0: #no match
                            if (event.gen_sel != 1
                                or (iflavor > -1 and event.gj_flavor[g] != iflavor)):
                                valGen = -1
                            else:
                                valGen = eval('event.gj_'+obs+'_'+reco)[g]
                            for weightindex,h in hmap_weight.iteritems():
                                weight = weights[0]
                                if (weightindex > 0): weight *= weights[weightindex]
                                h.Fill(valGen, -1, weight)


"""
steer
"""
def main():
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-t', '--task',
                            dest='task',   
                            default='fill',
                            help='task: optimize,fill [default: %default]')
    parser.add_option('-i', '--input',
                            dest='input',   
                            default='analysis.root',
                            help='input file, if directory the script will run in batch mode [default: %default]')
    parser.add_option('--obs',
                            dest='obs',   
                            default='mult',
                            help='observable [default: %default]')
    parser.add_option('-o', '--output',
                            dest='output', 
                            default='unfolding/fill',
                            help='Output directory [default: %default]')
    parser.add_option('--ri', '--rootinput',
                            dest='rootinput',
                            default='unfolding/optimize/output.root',
                            help='output root file [default: %default]')
    parser.add_option('--skipexisting', 
                            dest='skipexisting',
                            help='skip jobs with existing output files  [%default]',
                            default=False, action='store_true')
    parser.add_option('-q', '--queue', dest='queue',  help='Batch queue to use [default: %default]', default='longlunch')
    parser.add_option(      '--only',  dest='only',   help='csv list of samples to process (exact name, applies only for batch mode) [%default]', default=None, type='string')
    parser.add_option(      '--skip',  dest='skip',   help='csv list of samples to exclude (exact name, applies only for batch mode) [%default]', default=None, type='string')
    parser.add_option('--nw', '--nweights', dest='nweights', help='number of weights to run [%default]', default=0, type='int')
    (opt, args) = parser.parse_args()

    if opt.input.endswith('.root'):
        rootoutfilepath = opt.output+'/'+os.path.basename(opt.input)
        if (opt.skipexisting and os.path.isfile(rootoutfilepath)):
            print("skip existing file " + rootoutfilepath)
            return
        
        print("Filling response matrix from " + opt.input)
        rootinfile = ROOT.TFile(opt.rootinput, "READ");
        
        keys = []
        for tkey in rootinfile.GetListOfKeys():
            keys.append(tkey.GetName())
        #observables = ["mult", "width"]
        observables = ["mult", "width", "ptd", "ptds", "ecc", "tau21", "tau32", "tau43", "zg", "zgxdr", "zgdr", "ga_width", "ga_lha", "ga_thrust", "c1_02", "c1_05", "c1_10", "c1_20", "c2_02", "c2_05", "c2_10", "c2_20", "c3_02", "c3_05", "c3_10", "c3_20"]
        
        reco = ["charged"]
        
        flavors = ['all', 'light', 'bottom', 'gluon']
        
        hmap = {}
        for o in observables:
            hmap[o] = {}
            for r in reco:
                if o+'_'+r+'_responsematrix' in keys:
                    hmap[o][r] = {}
                    obj = rootinfile.Get(o+'_'+r+'_responsematrix').Clone()
                    obj.Reset()
                    for f in flavors:
                        hmap[o][r][f] = {}
                        hmap[o][r][f][0] = obj.Clone(o+'_'+r+'_'+f+'_responsematrix')
                        for w in range(opt.nweights):
                            hmap[o][r][f][w+1] = obj.Clone(o+'_'+r+'_'+f+'_wgt'+str(w+1)+'_responsematrix')
        
        #rootinfile.Close()
        os.system('mkdir -p %s' % opt.output)
        rootoutfile = ROOT.TFile(rootoutfilepath, "RECREATE");    
        rootoutfile.cd()
        
        tree = ROOT.TChain('tjsev')
        if opt.input == 'eos':
            opt.input = '/eos/user/m/mseidel/analysis/TopJetShapes/b312177/Chunks/MC13TeV_TTJets_0.root'
        if opt.input == 'eosdata':
            opt.input = '/eos/user/m/mseidel/analysis/TopJetShapes/b312177/Chunks/Data13TeV_SingleMuon_2016G_27.root'
        if (tree.Add(opt.input) == 0): return
        
        fillHist(hmap, tree, opt.nweights)

        for obs,hmap_reco in hmap.iteritems():
            for reco,hmap_flavor in hmap_reco.iteritems():
                for flavor,hmap_weight in hmap_flavor.iteritems():
                    for weight,h in hmap_weight.iteritems():
                        h.Write()
                        h.ProjectionX().Write()
                        h.ProjectionY().Write()
    else: # got input directory, go to batch
        cmsswBase = os.environ['CMSSW_BASE']
        workdir   = cmsswBase+'/src/TopLJets2015/TopAnalysis/'
        #parse selection lists
        onlyList=[]
        try:
            onlyList=opt.only.split(',')
        except:
            pass
        skipList=[]
        try:
            skipList=opt.skip.split(',')
        except:
            pass
        #create the tasklist
        inputlist=[]
        if os.path.isdir(opt.input):
            for file_path in os.listdir(opt.input):
                if file_path.endswith('.root'):
                    rootoutfilepath = workdir+opt.output+'/Chunks/'+os.path.basename(file_path)
                    if (opt.skipexisting and os.path.isfile(rootoutfilepath)): continue
                    #filter tags
                    if len(onlyList)>0:
                        if (not os.path.basename(file_path).rsplit('_',1)[0] in onlyList): continue
                    if len(skipList)>0:
                        if (os.path.basename(file_path).rsplit('_',1)[0] in skipList): continue
                    inputlist.append(os.path.join(opt.input,file_path))
        #print(inputlist)
        #FIXME old
        #print 'Running %d jobs to %s'%(len(inputlist),opt.queue)
        #njob = 1
        #for inputfile in inputlist:
        #    outputdir=workdir+opt.output
        #    localRun='python %s/src/TopLJets2015/TopAnalysis/test/TopJSAnalysis/fillUnfoldingMatrix.py --input %s --output %s --nweights %d'%(cmsswBase,inputfile,outputdir,opt.nweights)
        #    if (opt.skipexisting): localRun += " --skipexisting"
        #    cmd='bsub -q %s %s/src/TopLJets2015/TopAnalysis/scripts/wrapLocalAnalysisRun.sh \"%s\"' % (opt.queue,cmsswBase,localRun)
        #    print 'Submitting job %d/%d: %s'%(njob,len(inputlist),os.path.basename(inputfile))
        #    njob+=1
        #    os.system(cmd)
            
        #FIXME new
        FarmDirectory = '%s/FARM%s'%(cmsswBase,os.path.basename(opt.output))
        os.system('mkdir -p %s'%FarmDirectory)
        
        print 'Preparing %d tasks to submit to the batch'%len(inputlist)
        print 'Executables and condor wrapper are stored in %s'%FarmDirectory

        with open ('%s/condor.sub'%FarmDirectory,'w') as condor:

            condor.write('executable = {0}/$(cfgFile).sh\n'.format(FarmDirectory))
            condor.write('output     = {0}/output_$(cfgFile).out\n'.format(FarmDirectory))
            condor.write('error      = {0}/output_$(cfgFile).err\n'.format(FarmDirectory))
            condor.write('+JobFlavour = "{0}"\n'.format(opt.queue))

            jobNb=0
            for inputfile in inputlist:

                jobNb+=1
                outF=workdir+opt.output+'/Chunks/'+os.path.basename(inputfile)
                cfgFile='%s'%(os.path.splitext(os.path.basename(outF))[0])

                condor.write('cfgFile=%s\n'%cfgFile)
                condor.write('queue 1\n')
                
                with open('%s/%s.sh'%(FarmDirectory,cfgFile),'w') as cfg:

                    cfg.write('#!/bin/bash\n')
                    cfg.write('WORKDIR=`pwd`\n')
                    cfg.write('echo "Working directory is ${WORKDIR}"\n')
                    cfg.write('cd %s\n'%cmsswBase)
                    cfg.write('eval `scram r -sh`\n')
                    cfg.write('cd ${WORKDIR}\n')
                    localOutF=os.path.basename(outF)
                    nweights=opt.nweights
                    if localOutF.rsplit('_',1)[0] == 'MC13TeV_TTJets': nweights = 20
                    binningfile = '%s/src/TopLJets2015/TopAnalysis/unfolding/optimize/output.root'%(cmsswBase)
                    runOpts='--input %s --rootinput %s --output . --nweights %d'\
                        %(inputfile, binningfile, nweights)
                    cfg.write('python %s/src/TopLJets2015/TopAnalysis/test/TopJSAnalysis/fillUnfoldingMatrix.py %s\n'%(cmsswBase,runOpts))
                    if outF!=localOutF:
                        cfg.write('mv -v ${WORKDIR}/%s %s\n'%(localOutF,outF))

                os.system('chmod u+x %s/%s.sh'%(FarmDirectory,cfgFile))

        print 'Submitting jobs to condor, flavour "%s"'%(opt.queue)
        os.system('condor_submit %s/condor.sub'%FarmDirectory)

if __name__ == "__main__":
	sys.exit(main())
