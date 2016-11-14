# TopLJets2015

## Analysis twiki
Please keep 
https://twiki.cern.ch/twiki/bin/view/Main/TopLJ2015Analysis
up to date with the on-going tasks and results

## Installation instructions
To execute in your lxplus work area.
```
cmsrel CMSSW_8_0_11 
cd CMSSW_8_0_11/src 
cmsenv
git cms-init
#EGM smearer
git remote add -f -t ecal_smear_fix_80X emanueledimarco https://github.com/emanueledimarco/cmssw.git
git cms-addpkg EgammaAnalysis/ElectronTools
git checkout -b from-52f192a 52f192a
cd EgammaAnalysis/ElectronTools/data
git clone -b ICHEP2016_v2 https://github.com/ECALELFS/ScalesSmearings.git
cd -
#Pseudo-top producer with Markus fix
git remote add -f -t bugfix-pseudotop-bjets-81X intrepid42 https://github.com/intrepid42/cmssw.git
git cms-addpkg TopQuarkAnalysis/TopEventProducers
git fetch intrepid42
git checkout intrepid42/bugfix-pseudotop-bjets-81X
#analysis code
git clone git@github.com:pfs/TopLJets2015.git
cd TopLJets2015/TopAnalysis
git checkout 80x_dev
scram b -j 8
```

## Running ntuple creation
First time create a symbolic link to the jet energy corrections files
```
ln -s data/era2016/Spring16_25nsV6_DATA.db
ln -s data/era2016/Spring16_25nsV6_MC.db
ln -s data/era2016/RoccoR_13tev.txt 
```
To run locally the ntuplizer, for testing purposes
```
cmsRun test/runMiniAnalyzer_cfg.py runOnData=False/True outFilename=MiniEvents.root
```
To submit a list of samples, described in a json file to the grid you can use the following script.
```
python scripts/submitToGrid.py -j data/era2016/samples.json -c ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/runMiniAnalyzer_cfg.py --lfn my_output_directory_in_eos -s
```
Partial submission can be made adding "-o csv_list" as an option
Don't forget to init the environment for crab3 (e.g. https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3CheatSheet#Environment_setup)
```
source /cvmfs/cms.cern.ch/crab3/crab.sh
```
As soon as ntuple production starts to finish, to move from crab output directories to a simpler directory structure which can be easily parsed by the local analysis runThe merging can be run locally if needed by using the checkProductionIntegrity.py script
```
python scripts/submitCheckProductionIntegrity.py -i /store/group/phys_top/psilva/8db9ad6 -o /store/cmst3/user/psilva/LJets2016/8db9ad6
```

## Preparing the analysis 

Correction and uncertainty files are stored under data by era directories (e.g. data/era2015, data/era2016) in order no to mix different periods.
After ntuples are processed start by creating the json files with the list of runs/luminosity sections processed, e.g. as:
```
crab report grid/crab_Data13TeV_DoubleMuon_2016B
``` 
Then you can merge the json files for the same dataset to get the full list of run/lumi sections to analyse
```
mergeJSON.py grid/crab_Data13TeV_DoubleMuon_2016B/results/processedLumis.json grid/crab_Data13TeV_DoubleMuon_2015C/results/processedLumis.json grid/crab_Data13TeV_DoubleMuon_2015D/results/processedLumis.json --output data/era2016/Data13TeV_DoubleMuon_lumis.json
```
You can then run the brilcalc tool to get the integrated luminosity in total and per run (see https://twiki.cern.ch/twiki/bin/view/CMS/2015LumiNormtag for more details).
```
export PATH=$HOME/.local/bin:/afs/cern.ch/cms/lumi/brilconda-1.0.3/bin:$PATH
brilcalc lumi --normtag /afs/cern.ch/user/l/lumipro/public/normtag_file/normtag_DATACERT.json -i data/era2016/Data13TeV_DoubleMuon_lumis.json
```
Use the table which is printed out to update the "lumiPerRun" method in ReadTree.cc.
That will be used to monitor the event yields per run in order to identify outlier runs.
* Pileup weighting. To update the pileup distributions run the script below. It will store the data pileup distributions for different min.bias cross section in data/pileupWgts.root
```
python scripts/runPileupEstimation.py --json data/era2016/Data13TeV_DoubleMuon_lumis.json --out data/era2016/pileupWgts.root
```
* B-tagging. To apply corrections to the simulation one needs the expected efficiencies stored somwewhere. The script below will project the jet pT spectrum from the TTbar sample before and after applying b-tagging, to compute the expecte efficiencies. The result will be stored in data/expTageff.root
```
for i in "_powheg" "_herwig" "_scaledown" "_scaleup"; do
    python scripts/saveExpectedBtagEff.py -i /store/cmst3/user/psilva/LJets2016/8db9ad6/MC13TeV_TTJets${i} -o data/era2016/expTageff${i}.root;
done
mv data/era2016/expTageff_powheg.root data/era2016/expTageff.root
```
* MC normalization. This will loop over all the samples available in EOS and produce a normalization cache (weights to normalize MC). The file will be available in data/genweights.pck
```
python scripts/produceNormalizationCache.py -i /store/cmst3/user/psilva/LJets2016/8db9ad6 -o data/era2016/genweights.root
```
You're now ready to start locally the analysis.


## Running locally the analysis for testing

The analysis (histogram filling, final selection) is in src/ReadTree.cc.
Recompile (scram b) everytime you change it so that you can test the new features.
To test the code on a single file to produce plots.
```
python scripts/runLocalAnalysis.py -i MiniEvents.root
```
To run the code on a set of samples stored in EOS you can run it as shown below.
If "-q queue_name" is appended the jobs are submitted to the batch system instead of running locally. 
To check the status of your jobs run "bjobs" and then "bpeek job_number" if you want to inspect how the job is running in the cluster.
If "-n n_jobs" is passed the script runs locally using "n_jobs" parallel threads.
```
python scripts/runLocalAnalysis.py -i /store/cmst3/user/psilva/LJets2015/7e62835 -n 8 --runSysts -o analysis_muplus   --ch 13   --charge 1
```
If you want to suppress the mails sent automatically after job completion please do
```
export LSB_JOB_REPORT_MAIL=N
```
before submitting the jobs to the batch. After the jobs have run you can merge the outputs with
```
./scripts/mergeOutputs.py analysis_muplus
```
To plot the output of the local analysis you can run the following:
```
python scripts/plotter.py -i analysis_muplus/   -j data/era2016/samples.json  -l 12870
```
After the plotters are created one can run the QCD estimation normalization, by fitting the MET distribution.
The script will also produce the QCD templates using the data from the sideband region. It runs as
```
python scripts/runQCDEstimation.py --iso analysis_muplus/plots/plotter.root --noniso analysis_munoniso/plots/plotter.root --out analysis_muplus/
```
The output is a ROOT file called Data_QCDMultijets.root which can now be used in addition to the predictions of all the other backgrounds.
To include it in the final plots you can run the plotter script again (see instructions above).

## Submitting the full analysis to the batch system

A script wraps up the above procedure for all the signal and control regions used in the analyis.
To use it you can use the following script
```
sh scripts/steerAnalysis.sh <DISTS/MERGE/PLOT/BKG>
```

## Cross section fitting

We use the Higgs combination tool to perform the fit of the production cross section.
(cf. https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideHiggsAnalysisCombinedLimit for details of the release to use).
It currently has to be run from a CMSSW_7_1_5 release. To create the datacard you can run the following script
```
python scripts/createDataCard.py -i analysis_muplus/plots/plotter.root -o  analysis_muplus/datacard  -q analysis_muplus/.qcdscalefactors.pck -d nbtags
```
The script can be used to create the datacard from any histogram stored in plotter.root.
For the systematic variations it expects a 2D histogram named as HISTONAMEshapes_{exp,gen} filled with alternative variations of the shape,
being exp/gen used for experimental/generator-level systematics.
Additional systematics from alternative samples can also be used to build the datacards using the --systInput option.
Other options are available to choose the categories to use.
The datacards can be further combined using the standard combineCards.py script provided by the Higgs Combination package.

To run the fits and show the results you can use the following script.
```
python scripts/fitCrossSection.py "#mu^{+}"=analysis_muplus/datacard/datacard.dat -o analysis_muplus/datacard &
```
If --noFit is passed it displays the results of the last fit. The script is a wrapper used to run combine 
to perform the fit with and without systematics, produce the post-fit nuisance parameters summary
and the likelihood scans.
For the standard analysis one can re-use the steerAnalysis.sh script with two options CinC/SHAPE
will run the Cut-in-Categories/Shape analyses.
```
sh scripts/steerAnalysis.sh CinC/SHAPE
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
From the github area of the repository cleak on the green button "Compare,review and create a pull request"
to create the PR to merge with your colleagues.
