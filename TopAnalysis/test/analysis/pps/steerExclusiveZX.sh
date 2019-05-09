#!/bin/bash

WHAT=$1; 
if [ "$#" -ne 1 ]; then 
    echo "steerExclusiveZX.sh <SEL/MERGE/PLOT/WWW>";
    echo "   SEL           - launches selection jobs to the batch, output will contain summary trees and control plots"; 
    echo "   MERGESEL      - merge output"
    echo "   PLOTSEL        - make plots"
    echo "   WWWSEL        - move plots to web-based are"
    echo "   TRAINPUDISCR  - train pileup discriminator"
    echo "   RUNPRED       - run pileup discriminator prediction"
    echo "   PREPAREANA    - prepare bank of events for event mixing from the summary trees"
    echo "   COLLECTMIX    - collects all the mixing events found in PREPAREANA"
    echo "   ANA           - run analysis on the summary trees"
    echo "   PLOTANA       - plot analysis results"
    echo "   OPTIMSTATANA  - optimize the statistical analysis"
    echo "   DEFINESTATANA - define the final datacards based on the results of the optimization"
    exit 1; 
fi

#to run locally use local as queue + can add "--njobs 8" to use 8 parallel jobs
queue=tomorrow
githash=ab05162
eosdir=/store/cmst3/group/top/RunIIReReco/${githash}
outdir=/store/cmst3/user/psilva/ExclusiveAna/final/${githash}
signal_json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/signal_samples.json
plot_signal_json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/plot_signal_samples.json
samples_json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/samples.json
jetht_samples_json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/jetht_samples.json
zx_samples_json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/zx_samples.json
zbias_samples_json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/zbias_samples.json
RPout_json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/golden_noRP.json
wwwdir=~/www/ExclusiveAna
inputfileTag=MC13TeV_2017_GGH2000toZZ2L2Nu
inputfileTag=MC13TeV_2017_GGToEE_lpair
inputfileTag=Data13TeV_2017F_MuonEG
inputfileTag=MC13TeV_2017_GJets_HT400to600
#inputfileTag=Data13TeV_2017B_DoubleMuon
inputfileTag=Data13TeV_2017B_ZeroBias
inputfileTESTSEL=${eosdir}/${inputfileTag}/Chunk_0_ext0.root
lumi=41833
ppsLumi=37500
lumiUnc=0.025


RED='\e[31m'
NC='\e[0m'
case $WHAT in

    TESTSYNCH )
        
        python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/runLocalAnalysis.py \
            -i /store/cmst3/user/psilva/ExclusiveAna/synch/Data13TeV_2017_DoubleMuon_synch.root \
            -o testsynch.root --genWeights genweights_${githash}.root \
            --njobs 1 -q local --debug \
            --era era2017 -m ExclusiveZX::RunExclusiveZX --ch 0 --runSysts;

        ;;

    TESTSEL )
        python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/runLocalAnalysis.py \
            -i ${inputfileTESTSEL} --tag ${inputfileTag} \
            -o testsel.root --genWeights genweights_${githash}.root \
            --njobs 1 -q local --debug \
            --era era2017 -m ExclusiveZX::RunExclusiveZX --ch 0 --runSysts;
        ;;

    SEL )
        baseOpt="-i ${eosdir} --genWeights genweights_${githash}.root"
        baseOpt="${baseOpt} -o ${outdir} -q ${queue} --era era2017 -m ExclusiveZX::RunExclusiveZX --ch 0 --runSysts"
        #baseOpt="${baseOpt} --exactonly"        
	python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/runLocalAnalysis.py ${baseOpt} \
            --only ${samples_json},${zx_samples_json}; #,${zbias_samples_json};
	;;

    CHECKSELINTEG )
        python scripts/checkLocalAnalysisInteg.py ../../../FARM${githash}/ analysis/data tree
        ;;


    MERGESEL )
	mergeOutputs.py /eos/cms/${outdir} True;
	;;

    PLOTSEL )
        lumiSpecs="--lumiSpecs a:2642"
        kFactorList="--procSF #gamma+jets:1.33"
	commonOpts="-i /eos/cms/${outdir} --puNormSF puwgtctr -l ${lumi} --mcUnc ${lumiUnc} ${kFactorList} ${lumiSpecs}"
	python scripts/plotter.py ${commonOpts} -j ${samples_json}    -O plots/sel --only mboson,mtboson,pt,eta,met,jets,nvtx,ratevsrun --saveLog; 
        python scripts/plotter.py ${commonOpts} -j ${samples_json}    --rawYields --silent --only gen -O plots/ -o bkg_plotter.root ; 
	python scripts/plotter.py ${commonOpts} -j ${zx_samples_json} --rawYields --silent --only gen -O plots/;
        python test/analysis/pps/computeDileptonSelEfficiency.py 
        mv *.{png,pdf} plots/sel/
        mv effsummary* plots/
	;;

    WWWSEL )
	mkdir -p ${wwwdir}/presel
	cp plots/sel/*.{png,pdf} ${wwwdir}/presel
	cp $CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/index.php ${wwwdir}/presel
	;;

    TRAINPUDISCR )
        echo "Please remember to use a >10_3_X release for this step"

        commonOpts="--trainFrac 0.3  --RPout test/analysis/pps/golden_noRP.json -o /eos/cms/${outdir}/train_results"
        python test/analysis/pps/trainPUdiscriminators.py ${commonOpts} -s "isZ && evcat==13*13 && bosonpt<10 && trainCat>=0" 
        #python test/analysis/pps/trainPUdiscriminators.py ${commonOpts} -s "hasZBTrigger && trainCat>=0" --zeroBiasTrain

        ;;

    RUNPRED )
        predin=/eos/cms/${outdir}/Chunks
        predout=/eos/cms/${outdir}/Chunks/pudiscr 
        condor_prep=runpred_condor.sub
        echo "executable  = ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/pps/wrapPUDiscrTrain.sh" > $condor_prep
        echo "output      = ${condor_prep}.out" >> $condor_prep
        echo "error       = ${condor_prep}.err" >> $condor_prep
        echo "log         = ${condor_prep}.log" >> $condor_prep
        echo "requirements = (OpSysAndVer =?= \"SLCern6\")"  >> $condor_prep
        echo "arguments   = ${predout} \$(chunk)" >> $condor_prep
        echo "queue chunk matching (${predin}/*.root)" >> $condor_prep
        condor_submit $condor_prep
        ;;

    PREPAREANA )       
        step=0
        predin=/eos/cms/${outdir}/Chunks
        predout=/eos/cms/${outdir}/mixing
        condor_prep=runana_condor.sub
        echo "executable  = ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/pps/wrapAnalysis.sh" > $condor_prep
        echo "output      = ${condor_prep}.out" >> $condor_prep
        echo "error       = ${condor_prep}.err" >> $condor_prep
        echo "log         = ${condor_prep}.log" >> $condor_prep
        echo "arguments   = ${step} ${predout} ${predin} \$(chunk)" >> $condor_prep
        echo "queue chunk matching (${predin}/Data*.root)" >> $condor_prep
        condor_submit $condor_prep
        
        #python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/runExclusiveAnalysis.py --step 0 --jobs 8 \
        #    --json ${samples_json} --RPout ${RPout_json} -o ${predout} -i ${predin};

        ;;

    COLLECTMIX )
        python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/collectEventsForMixing.py /eos/cms/${outdir}
        ;;

    ANA )
        step=1
        predin=/eos/cms/${outdir}/Chunks
        predout=/eos/cms/${outdir}/analysis
        condor_prep=runana_condor.sub
        mix_file=/eos/cms/${outdir}/mixing/mixbank.pck
        echo "executable  = ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/pps/wrapAnalysis.sh" > $condor_prep
        echo "output      = ${condor_prep}.out" >> $condor_prep
        echo "error       = ${condor_prep}.err" >> $condor_prep
        echo "log         = ${condor_prep}.log" >> $condor_prep
        echo "arguments   = ${step} ${predout} ${predin} \$(chunk)" ${mix_file} >> $condor_prep
        echo "queue chunk matching (${predin}/*.root)" >> $condor_prep
        condor_submit $condor_prep

        #run locally
        #python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/runExclusiveAnalysis.py --step 1 --jobs 8 \
        #    --json ${samples_json} --RPout ${RPout_json} -o ${predout} --mix ${mix_file} -i ${predin};
        
        ;;

    ANASIG )
        
        python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/runExclusiveAnalysis.py --step 1 --jobs 8 \
            --json ${signal_json} --RPout ${RPout_json} --mix /eos/cms/${outdir}/mixing/mixbank.pck \
            -i /eos/cms/${outdir}/Chunks -o anasig/;

         cp -v anasig/Chunks/*.root /eos/cms/${outdir}/analysis/

        ;;

    MERGEANA )
        mergeOutputs.py /eos/cms/${outdir}/analysis;
        ;;

    PLOTANA )

        lptalumi=2642
        lumiSpecs="lpta:${lptalumi},lptaneg:${lptalumi},lptapos:${lptalumi},lptahpur:${lptalumi},lptahpur120:${lptalumi},lptahpur130:${lptalumi},lptahpur140:${lptalumi},lptahpur150:${lptalumi}"
        lumiSpecs="${lumiSpecs},lpta120neg:${lptalumi},lpta130neg:${lptalumi},lpta140neg:${lptalumi},lpta150neg:${lptalumi},lpta120pos:${lptalumi},lpta130pos:${lptalumi},lpta140pos:${lptalumi},lpta150pos:${lptalumi}"
        lumiSpecs="${lumiSpecs},lptahpur120neg:${lptalumi},lptahpur130neg:${lptalumi},lptahpur140neg:${lptalumi},lptahpur150neg:${lptalumi},lptahpur120pos:${lptalumi},lptahpur130pos:${lptalumi},lptahpur140pos:${lptalumi},lptahpur150pos:${lptalumi}"
	baseOpts="-i /eos/cms/${outdir}/analysis --lumiSpecs ${lumiSpecs} --procSF #gamma+jets:1.33 -l ${ppsLumi} --mcUnc ${lumiUnc} ${lumiSpecs} ${kFactorList}"
        
        plots=xangle_eeZhpur,xangle_mmZhpur,xangle_emhpur,xangle_lptahpur
        commonOpts="${baseOpts} -j ${samples_json} --signalJson ${plot_signal_json} -O plots/analysis/zx_yields" 
        python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/plotter.py ${commonOpts} --only ${plots} --strictOnly --saveTeX --rebin 4;

        plots=""
        for c in eeZ mmZ em lpta hpta eeZhpur mmZhpur emhpur lptahpur hptahpur; do         
            for d in acopl ptll yll l1pt l2pt l1eta l2eta; do
                plots="${plots},${d}_${c}"                
            done            
            for d in xangle ntk; do
                for s in pos neg; do                    
                    plots="${plots},${d}_${c}${s}"
                done
            done
            for d in csi; do            
                for x in 120 130 140 150; do
                    for s in pos neg; do
                        plots="${plots},${d}_${c}${x}${s}"
                    done
                done        
            done
            
            for x in 120 130 140 150; do
                for d in rho mpp ypp mmass nextramu ptll yll; do
                    plots="${plots},${d}_${c}${x}"
                done
                for r in HF HE EB EE; do
                    for d in PFMult PFHt PFPz; do
                        plots="${plots},${d}${r}_${c}${x}"
                    done
                done
                for d in csi csi2d mpp2d ypp2d; do
                    for s in pos neg; do
                        plots="${plots},${d}_${c}${x}${s}"
                    done
                done
            done
        done
        #python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/plotter.py ${commonOpts} --only ${plots} --strictOnly --signalJson ${plot_signal_json} --saveLog; # --normToData;
        ;;

    WWWANA )
        mkdir -p ${wwwdir}/ana
        cp /eos/cms/${outdir}/analysis/plots/*.{png,pdf,dat} ${wwwdir}/ana
        cp $CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/index.php ${wwwdir}/ana
        ;;

    OPTIMSTATANA )

        #run combine on condor
        for m in 800 1000 1200 1400 1600; do

            echo "Will first prepare the datacards/shapes files for m=${m}"
            python test/analysis/pps/prepareOptimScanCards.py ${m}

            echo "Will launch combine runs to condor for m=${m}"
            condor_prep="condor_optim_m${m}.sub"
            echo "executable  = ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/pps/wrapOptim.sh" > $condor_prep
            echo "output      = ${condor_prep}.out" >> $condor_prep
            echo "error       = ${condor_prep}.err" >> $condor_prep
            echo "log         = ${condor_prep}.log" >> $condor_prep
            echo "+JobFlavour =\"workday\""> $condor_prep
            for x in 120 130 140 150; do
                echo "arguments   = ${m} ${x} ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/analysis/stat_m${m}" >> $condor_prep
                echo "queue 1" >> $condor_prep
            done
            condor_submit $condor_prep
        done
        ;;

    DEFINESTATANA)
        python test/analysis/pps/compareOptimResults.py analysis/
        #python test/analysis/pps/plotLimits.py 1000=optimresults_1000.pck,1200=optimresults_1200.pck,1400=optimresults_1400.pck
        #python test//analysis/pps/compareOptimResults.py
        ;;

esac
