# TopLJets2015

## Installation instructions

These installation instructions correspond to the 2016 data/MC Moriond17 re-reco.
To install execute the following in your work area.

```
cmsrel CMSSW_8_0_28
cd CMSSW_8_0_28/src 
cmsenv

#EGM electron regression+smearer
git cms-init
git cms-merge-topic cms-egamma:EGM_gain_v1
cd EgammaAnalysis/ElectronTools/data
git clone -b Moriond17_gainSwitch_unc https://github.com/ECALELFS/ScalesSmearings.git

#MET
git cms-merge-topic cms-met:METRecipe_8020 -u
git cms-merge-topic cms-met:METRecipe_80X_part2 -u

#pseudo-top rivet based and b-frag utils
git cms-merge-topic -u intrepid42:pseudotoprivet_80x
cd TopQuarkAnalysis
git clone ssh://git@gitlab.cern.ch:7999/CMS-TOPPAG/BFragmentationAnalyzer.git
cd -

#muon calibration
git clone -o Analysis https://github.com/bachtis/analysis.git -b KaMuCa_V4 KaMuCa

#ntuplizer
git clone -b 80x_rereco_rev git@github.com:pfs/TopLJets2015.git

#compile
scram b -j 8
```

## Running ntuple creation and checking the selection

The ntuplizer is steered with test/runMiniAnalyzer_cfg.py.
It takes several options from command line (see cfg for details).
To run locally the ntuplizer, for testing purposes do something like:

```
cmsRun test/runMiniAnalyzer_cfg.py runOnData=False era=era2016 outFilename=MC13TeV_TTJets.root
cmsRun test/runMiniAnalyzer_cfg.py runOnData=True  era=era2016 outFilename=Data13TeV_MuonEG.root
```

The default files point to the ones used in the TOP synchronization exercise
(see details in )
The output files can be analyzed with a simple executable (also a skeleton for the implementation of new analysis) as detailed below.
The output file will contain cutflow histograms which can be used to check the selection against the synchronization twiki.
Check the implementation of this analysis in src/TOPSynchExercise.cc. 
Other analysis should also be implemented there.

```
analysisWrapper --in MC13TeV_TTJets.root   --out mc_synch.root   --method TOPSynchExercise::RunTOPSynchExercise
analysisWrapper --in Data13TeV_MuonEG.root --out data_synch.root --method TOPSynchExercise::RunTOPSynchExercise
```

A wrapper script is also provided to run over directories with trees or single files and allowing to schedule the execution to the job.
See examples under scripts/ in ```steer*Analysis.sh```.


## Submitting ntuple creation through the grid

To submit the ntuplizer to the grid start by setting the environment for crab3.
More details can be found in [CRAB3CheatSheet](https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3CheatSheet#Environment_setup)

```
source /cvmfs/cms.cern.ch/crab3/crab.sh
```
The following script helps submitting a list of files described in a json file.
Partial submission can be made adding "-o csv_list" as an option.
Adding "-s" will trigger the submission to the grid (otherwise the script only writes down the crab cfg files)

```
python scripts/submitToGrid.py -j data/era2016/samples.json -c ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/runMiniAnalyzer_cfg.py 
```

As soon as ntuple production starts to finish, to move from crab output directories to a simpler directory structure which can be easily parsed by the local analysis runThe merging can be run locally if needed by using the checkProductionIntegrity.py script

```
python scripts/submitCheckProductionIntegrity.py -i /store/cmst3/group/top/psilva/b312177 -o /store/cmst3/group/top/ReReco2016/b312177
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
The following script runs brilcalc inclusively and per trigger path, and stores the results in a ROOT file with the total integrated lumi per run.
It takes a bit to run, depending on the number of triggers configured to use in the analysis
```
export PATH=$HOME/.local/bin:/afs/cern.ch/cms/lumi/brilconda/bin:$PATH
python scripts/convertLumiTable.py -o data/era2016/
```

## Preparing the analysis 

Correction and uncertainty files are stored under data by era directories (e.g. data/era2015, data/era2016) in order no to mix different periods.

* Pileup weighting. To update the pileup distributions run the script below. It will store the data pileup distributions for different min.bias cross section in data/pileupWgts.root
```
python scripts/runPileupEstimation.py --json /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt --out data/era2016/pileupWgts.root
```
* B-tagging. To apply corrections to the simulation one needs the expected efficiencies stored somwewhere. The script below will project the jet pT spectrum from the TTbar sample before and after applying b-tagging, to compute the expecte efficiencies. The result will be stored in data/expTageff.root
```
python scripts/saveExpectedBtagEff.py -i /store/cmst3/group/top/ReReco2016/b32c02e/MC13TeV_TTJets -o data/era2016/expTageff.root;
```
* MC normalization. This will loop over all the samples available in EOS and produce a normalization cache (weights to normalize MC). The file will be available in data/genweights.pck
```
python scripts/produceNormalizationCache.py -i /store/cmst3/group/top/ReReco2016/b312177 -o data/era2016/genweights.root
```
The lepton trigger/id/iso efficiencies should also be placed under data/era2016. 
The src/LeptonEfficiencyWrapper.cc  should then be updated to handle the reading of the ROOT files and the application of the scale factors
event by event.

## Plotting
To plot the output of the local analysis you can run the following:
```
python scripts/plotter.py -i analysis_muplus/   -j data/era2016/samples.json  -l 12870
```

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
