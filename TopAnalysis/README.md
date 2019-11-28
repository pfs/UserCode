# TopLJets2015

## Installation instructions

These installation instructions correspond to the 2017 data/MC production.
To install execute the following in your work area.
Notice: if you are not creating the ntuples, you can skip the part of the instructions 
marked with the `##OPTIONAL/##END OPTIONAL` markers.
If compilation fails for some reason repeat the scram b...

```
cmsrel CMSSW_10_6_0
cd CMSSW_10_6_0/src
cmsenv

#EGM id
git cms-merge-topic cms-egamma:EgammaPostRecoTools
scram b -j 8

#B-fragmentation analyzer
mkdir TopQuarkAnalysis 
cd TopQuarkAnalysis
git clone -b 94x https://gitlab.cern.ch/psilva/BFragmentationAnalyzer.git
scram b -j 8
cd -

#This package
cd $CMSSW_BASE/src
git clone https://github.com/pfs/TopLJets2015.git -b 106_protonreco
cd TopLJets2015

scram b -j 8
```

## Running ntuple creation and checking the selection

The ntuplizer is steered with test/runMiniAnalyzer_cfg.py.
It takes several options from command line (see cfg for details).
To run locally the ntuplizer, for testing purposes do something like:

```
cmsRun test/runMiniAnalyzer_cfg.py runOnData=False era=era2017 outFilename=MC13TeV_TTJets.root noParticleLevel=True
cmsRun test/runMiniAnalyzer_cfg.py runOnData=True  era=era2017 outFilename=Data13TeV_DoubleEG_UL.root
cmsRun test/runL1PrefireAna_cfg.py runOnData=True  era=era2017 outFilename=Data13TeV_SinglePhoton_l1prefire.root
```

To submit the ntuplizer to the grid start by setting the environment for crab3.
More details can be found in [CRAB3CheatSheet](https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3CheatSheet#Environment_setup)

```
source /cvmfs/cms.cern.ch/crab3/crab.sh
```
The following script helps submitting a list of files described in a json file.
Partial submission can be made adding "-o csv_list" as an option.
Adding "-s" will trigger the submission to the grid (otherwise the script only writes down the crab cfg files)

```
python scripts/submitToGrid.py -j data/era2017/samples.json -c ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/runMiniAnalyzer_cfg.py --only MC -s
python scripts/submitToGrid.py -j data/era2017/samples.json -c ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/runMiniAnalyzer_cfg.py --only 2017B_Single,2017C_Single,2017D_Single,2017E_Single,2017F_Single --addParents -s
python scripts/submitToGrid.py -j data/era2017/samples.json -c ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/runMiniAnalyzer_cfg.py --only 2017H --addParents --rawParents --lumi /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Final/Cert_306896-307082_13TeV_PromptReco_Collisions17_JSON_LowPU.txt  -s
python scripts/submitToGrid.py -j data/era2017/samples.json -c ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/runL1PrefireAna_cfg.py --addParents --only JetHT,SinglePhoton,SingleElectron,MuonEG,DoubleEG --lfn /store/group/cmst3/group/top/psilva/l1prefire -w grid_prefire -s
python scripts/submitToGrid.py -j data/era2016/samples.json -c ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/runMiniAnalyzer_cfg.py -w grid_2016 --lumi /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Final/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt --era era2016
python scripts/submitToGrid.py -j data/era2018/samples.json -c ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/runMiniAnalyzer_cfg.py -w grid_2018 --lumi FIXME --era era2018
```

As soon as ntuple production starts to finish, to move from crab output directories to a simpler directory structure which can be easily parsed by the local analysis runThe merging can be run locally if needed by using the checkProductionIntegrity.py script

```
python scripts/mergeGridOutputs.py -i /store/cmst3/group/top/psilva/ab05162/ -o /store/cmst3/group/top/RunIIReReco/ab05162/
python scripts/mergeGridOutputs.py -i /store/cmst3/group/top/grid_2016/113427a -o /store/cmst3/group/top/RunIIReReco/113427a_2016
```

## Luminosity

After ntuples are processed, you can create the list of runs/lumi sections processed using crab as:
```
a=(`find grid/ -maxdepth 1 | grep crab_Data `)
for i in ${a[@]}; do
    crab kill ${i};
    crab status ${i};
    crab report ${i}; 
done
``` 
In case of failed jobs the missing lumis can be processed with the following script to wrap the tedious process of 
updating the cfg with a finer grain luminosity per job and the missing lumis json
```
for i in ${a[@]}; do
    python scripts/runMissingLumiSecs.py ${i}
done
```
You can then run the brilcalc tool to get the integrated luminosity in total and per run 
(see http://cms-service-lumi.web.cern.ch/cms-service-lumi/brilwsdoc.html for more details).

```
export PATH=$HOME/.local/bin:/cvmfs/cms-bril.cern.ch/brilconda/bin:$PATH
brilcalc lumi -b "STABLE BEAMS" -u /pb -i processedLumis.json --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json
```

The following script runs brilcalc inclusively and per trigger path, and stores the results in a ROOT file with the total integrated lumi per run.
It takes a bit to run, depending on the number of triggers configured to use in the analysis

```
python scripts/convertLumiTable.py --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json 
python scripts/convertLumiTable.py -o data/era2016/ -y 2016 --lumi /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Final/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt
```

## Running ntuplizer on condor

You can use this option also, although there won't be a lumi report in the end to run brilcalc on.
Check the options with
```
python scripts/submitLocalNtuplizer.py -h

```

After jobs run, check the integrity by doing
```
python scripts/checkLocalNtuplizerInteg.py condor_file
```


## Preparing the analysis 

Correction and uncertainty files are stored under data by era directories (e.g data/era2017) in order no to mix different periods.

* Pileup weighting. To update the pileup distributions run the script below. It will store the data pileup distributions for different min.bias cross section in data/pileupWgts.root
```
python scripts/runPileupEstimation.py --out data/era2017/pileupWgts.root
python scripts/runPileupEstimation.py --out data/era2016/pileupWgts.root \
       --json /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Final/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt \
       --puJson /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt
```
* B-tagging. To apply corrections to the simulation one needs the expected efficiencies stored somwewhere. The script below will project the jet pT spectrum from the TTbar sample before and after applying b-tagging, to compute the expecte efficiencies. The result will be stored in data/expTageff.root
```
python scripts/saveExpectedBtagEff.py -i /store/cmst3/group/top/RunIIReReco/ab05162/MC13TeV_2017_TTJets      -o data/era2017/expectedBtagEff.root;
python scripts/saveExpectedBtagEff.py -i /store/cmst3/group/top/RunIIReReco/2016/0c522df/MC13TeV_2016_TTJets -o data/era2016/expectedBtagEff.root;
```
* MC normalization. This will loop over all the samples available in EOS and produce a normalization cache (weights to normalize MC). The file will be available in data/genweights.pck
```
python scripts/produceNormalizationCache.py -i /store/cmst3/group/top/RunIIReReco/ab05162      -o data/era2017/genweights_ab05162.root
python scripts/produceNormalizationCache.py -i /store/cmst3/group/top/RunIIReReco/2016/0c522df -o data/era2016/genweights_0c522df.root
```
The lepton/photon trigger/id/iso efficiencies should also be placed under data/era2017. 
The src/EfficiencyScaleFactorsWrapper.cc  should then be updated to handle the reading of the ROOT files and the application of the scale factors
event by event.

## Updating the code

Commit your changes regularly with
```
git commit -a -m'comment on the changes made'
```
Push to your forked repository
```
git push git@github.com:MYGITHUBLOGIN/TopLJets2015.git
```
From the github area of the repository cleak on the green button "Compare,review and create a pull request" to create the PR to merge with your colleagues.

# Local analyses

The ROOT trees created by MiniAnalyzer.cc can be analyzed with a simple executable.
See some examples under `src`.
The new executable should be included in `bin/analysisWrapper.cc` so that it can be used with the runLocalAnalysis.py script
which allows to run over single files or full directories.
See examples under test/ in ```steer*Analysis.sh```.
To plot the output of the local analysis you can run the following:
```
python scripts/plotter.py -i analysis/   -j data/era2017/samples.json  -l 12870
```
# BDT training
This part currently works only for the VBF analysis. The signal and background trees must have been produced in the previous session by enabling "--mvatree" with "SEL" option.
The training is done per category with the following command:
```
python scripts/trainVbfBDT --vbf nt=50:mns=5:md=3:abb=0.6:nc=1 --ext LowVPtHighMJJ --sig signal.root --bkg backgrounds.root --cat A:VBF --card LowVPtHighMJJCard
```
The example above is for "LowVPtHighMJJ" category. Once happy with the training, 
```
cp vbf/weights/LowVPtHighMJJ_BDT_VBF0LowVPtHighMJJ* test/analysis/VBF_weights
```
Update src/VBFVectorBoson.cc accordingly and compile.

## Transformed BDT
To make the background BDT distribution flat, following steps must be followed:
   * Produce plots with the compiled version of the code including new BDT weights (see Local analysis) and
     ```
     cd test/analysis/VBF_weights
     python getInverseCDFFromPlotter.py PATH_TO_plotter.root
     ```
   * Reproduce the plots to have the proper flat BDT distribution of background

# Preparation of the data cards and workspaces

This part currently works only for the VBF analysis. It assumes that there are root files in the working directory which includes the plots created in the previous section. For the example below, the root file is called plotter_LowVPtHighMJJ2016.root. Running the script below will make two data cards and a single corresponding workspace for signal region (here gamma+jets) and control region (here Z+jets) of the "LowVPtHighMJJ" the category:
```
python scripts/makeWorkspace.py --Chan LowVPtHighMJJ --Hist vbfmva --nBin 20 --year 2016 --shapeOnly --doSignalPH
```
The input histogram will be rebinned to have five bins. If you remove "--doSignalPH", the signal process in the signal and control region will NOT be connected via the transfer function (TF). Note that the TF part is not developed fully since it is statistics limited. "--shapeOnly" is to disentangle the shape and rate systematics for non-TF version.
# Photon fake rate estimation and distributions
The method is explained in AN-18-046. The FR is measured in bins of mjj separately for photons in the barrel and endcap. It is measured in a control region (CR) which has the same selection as the signal region except the 2nd jet that is required to ```fail``` the loose pileup identification.
## Input data for FR measurement
   * Photon data (HighVPt categories): add ```--CR``` to the SEL option in steerVBFVectorBoson.sh
   * Jet data (LowVPt category): run steerVBFVectorBoson.sh with SELJETHT in which ```--CR``` is a default input
## Input templates for FR measurement
   * Tight template: the output of running default SEL option in steerVBFVectorBoson.sh on GJet samples
   * Fake template: the output of running SELJETHT option with ```-q QCDTemp``` argument in steerVBFVectorBoson.sh
## Fake rate estimation
   * Run the command below on the output of the aforementioned steps:
```
python scripts/createFakeRatio.py --fGdata PHOTON-DATA --fJdata JETHT-DATA --fJQCD JETHT-DATA-QCDTEMP --fGMC GJET-MC --cats HighVPtHighMJJA:HighVPtLowMJJA:LowVPtHighMJJA
```
   * Copy the output ```fakeRatios.root``` to data/eraYEAR
## Fake rate application
Here the distributions of not-tight photons are scaled by the estimated FR. The distributions can be used then in plotting. Need to activate ```--SRfake``` for the SEL option in steerVBFVectorBoson.sh and run on SinglePhoton data

