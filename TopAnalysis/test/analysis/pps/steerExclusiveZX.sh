#!/bin/bash

WHAT=$1; 
if [ "$#" -ne 1 ]; then 
    echo "steerExclusiveZX.sh <SEL/MERGE/PLOT/WWW>";
    echo "   SEL           - launches selection jobs to the batch, output will contain summary trees and control plots"; 
    echo "   MERGESEL      - merge output"
    echo "   PLOTSEL        - make plots"
    echo "   TRAINPUDISCR  - train pileup discriminator"
    echo "   RUNPRED       - run pileup discriminator prediction"
    echo "   PREPAREANA    - prepare bank of events for event mixing from the summary trees"
    echo "   COLLECTMIX    - collects all the mixing events found in PREPAREANA"
    echo "   ANA/ANASIG    - run analysis on the summary trees"
    echo "   CHECKANA      - check analysis integrity and re-run locally jobs which have failed"
    echo "   PLOTANA       - plot analysis results"
    echo "   OPTIMSTATANA  - optimize the statistical analysis"
    echo "   DEFINESTATANA - define the final datacards based on the results of the optimization"
    echo "   WWW           - move plots to web-based area"
    exit 1; 
fi

#to run locally use local as queue + can add "--njobs 8" to use 8 parallel jobs
queue=tomorrow
githash=ab05162
eosdir=/store/cmst3/group/top/RunIIReReco/${githash}
#githash=b949800_newppscalib
#eosdir=/store/cmst3/group/top/RunIIReReco/2017/${githash}
signal_dir=/store/cmst3/group/top/RunIIReReco/2017/vxsimulations_19Aug
outdir=/store/cmst3/user/psilva/ExclusiveAna/final/${githash}
signal_json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/signal_samples.json
plot_signal_json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/plot_signal_samples.json
samples_json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/samples.json
jetht_samples_json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/jetht_samples.json
zx_samples_json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/zx_samples.json
zbias_samples_json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/zbias_samples.json
RPout_json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/golden_noRP.json
wwwdir=/eos/user/p/psilva/www/ExclusiveAna
inputfileTag=MC13TeV_2017_GGH2000toZZ2L2Nu
inputfileTag=MC13TeV_2017_GGToEE_lpair
inputfileTag=Data13TeV_2017F_MuonEG
inputfileTag=MC13TeV_2017_GJets_HT400to600
#inputfileTag=Data13TeV_2017B_DoubleMuon
inputfileTag=Data13TeV_2017B_ZeroBias
inputfileTESTSEL=${eosdir}/${inputfileTag}/Chunk_0_ext0.root
lumi=41529
ppsLumi=37193
lptalumi=2642
lptappslumi=2288
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
        lumiSpecs="--lumiSpecs a:${lptalumi}"
        kFactorList="--procSF #gamma+jets:1.33"
	commonOpts="-i /eos/cms/${outdir} --puNormSF puwgtctr -l ${lumi} --mcUnc ${lumiUnc} ${kFactorList} ${lumiSpecs}"
	python scripts/plotter.py ${commonOpts} -j ${samples_json}    -O /eos/cms/${outdir}/plots/sel -o plotter.root --only mboson,mtboson,pt,eta,met,jets,nvtx,ratevsrun --saveLog; 
        python scripts/plotter.py ${commonOpts} -j ${samples_json}    --rawYields --silent --only gen -O /eos/cms/${outdir}/plots/ -o plots/bkg_gen_plotter.root; 
	python scripts/plotter.py ${commonOpts} -j ${zx_samples_json} --rawYields --silent --only gen -O /eos/cms/${outdir}/plots/ -o plots/zx_gen_plotter.root;
        python test/analysis/pps/computeDileptonSelEfficiency.py /eos/cms/${outdir}/plots/       
	;;

    TRAINPUDISCR )
        echo "Please remember to use a >10_3_X release for this step"

        commonOpts="--trainFrac 0.3  --RPout ${RPout_json} -o /eos/cms/${outdir}/train_results"
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

    CHECKMIX )
        python test/analysis/pps/checkFinalNtupleInteg.py /eos/cms/${outdir} 0 0
        ;;

    COLLECTMIX )
        python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/collectEventsForMixing.py /eos/cms/${outdir}
        ;;


    TESTANA )

        predin=/eos/cms/${outdir}/Chunks
        #file=Data13TeV_2017D_DoubleEG_2.root
        #file=MC13TeV_2017_DY50toInf_fxfx_0.root
        file=Data13TeV_2017C_DoubleEG_0.root

        #predin=/eos/cms/${signal_dir}
        #file=Z_m_X_1440_xangle_130_2017_preTS2.root

        predout=/eos/cms/${outdir}/analysis
        mix_file=/eos/cms/${outdir}/mixing/mixbank.pck
        
        #run locally
        python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/runExclusiveAnalysis.py --step 1 --jobs 1 \
            --json ${samples_json},${signal_json} --RPout ${RPout_json} -o ${predout} --mix ${mix_file} -i ${predin} --only ${file} --maxEvents 10000;
        
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
        echo "arguments   = ${step} ${predout} ${predin} \$(chunk) ${mix_file}" >> $condor_prep
        echo "queue chunk matching (${predin}/*.root)" >> $condor_prep
        condor_submit $condor_prep

        #run locally
        #python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/runExclusiveAnalysis.py --step 1 --jobs 8 \
        #    --json ${samples_json} --RPout ${RPout_json} -o ${predout} --mix ${mix_file} -i ${predin};
        
        ;;

    ANASIG )
        step=1
        predin=/eos/cms/${signal_dir}
        predout=/eos/cms/${outdir}/analysis
        condor_prep=runanasig_condor.sub
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
        #    --json ${signal_json} --RPout ${RPout_json} --mix /eos/cms/${outdir}/mixing/mixbank.pck \
        #    -i /eos/cms/${signal_dir} -o anasig/;
        # cp -v anasig/Chunks/*.root /eos/cms/${outdir}/analysis/

        ;;

    CHECKANA )
        #0-just check
        #1-run locally
        #2-submit to condor
        python test/analysis/pps/checkFinalNtupleInteg.py /eos/cms/${outdir} 1 1
        ;;

    MERGEANA )
        mergeOutputs.py /eos/cms/${outdir}/analysis;
        ;;

    PLOTANAPERERA )

        plots=""
        for evcat in a ee mm em; do
            for p in ptll mll; do 
                plots="${p}_${evcat},${plots}"
            done
        done        
        for era in B C D E F; do
            alumi=${lptalumi}
            eralumi=${lumi}
            if [ "${era}" = "B" ]; then
                alumi=`echo ${alumi}*0.115 | bc`
                eralumi=`echo ${eralumi}*0.115 | bc`
            elif [ "${era}" = "C" ]; then
                alumi=`echo ${alumi}*0.233 | bc`
                eralumi=`echo ${eralumi}*0.233 | bc`
            elif [ "${era}" = "D" ]; then
                alumi=`echo ${alumi}*0.103 | bc`
                eralumi=`echo ${eralumi}*0.103 | bc`
            elif [ "${era}" = "E" ]; then
                alumi=`echo ${alumi}*0.22 | bc`
                eralumi=`echo ${eralumi}*0.22 | bc`
            elif [ "${era}" = "F" ]; then
                alumi=`echo ${alumi}*0.329 | bc`
                eralumi=`echo ${eralumi}*0.329 | bc`
            fi

	    baseOpts="-i /eos/cms/${outdir}/analysis --lumiSpecs a:${alumi} --procSF #gamma+jets:1.33 -l ${eralumi} --mcUnc ${lumiUnc}"
            era_json=test/analysis/pps/test_samples_${era}.json;
            commonOpts="${baseOpts} -j ${era_json} --signalJson ${plot_signal_json} -O /eos/cms/${outdir}/analysis/plots${era}"
            python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/plotter.py ${commonOpts} --only ${plots} --strictOnly;
        done
        ;;

    PLOTSENS )

        for xangle in 120 130 140 150; do
            python test/analysis/pps/estimateLocalSensitivity.py \
                -i /eos/cms/${outdir}/analysis \
                -o /eos/cms/${outdir}/analysis/plots \
                --xangle ${xangle};
        done
        ;;

    PLOTANA )

        lumiSpecs="mmrpin:${ppsLumi},eerpin:${ppsLumi},emrpin:${ppsLumi},a:${lptalumi},arpin:${lptappslumi}"        
        for c in neg pos hpur hpurpos hpurneg hpur120xangle hpur130xangle hpur140xangle hpur150xangle; do
            lumiSpecs="${lumiSpecs},mmrpin${c}:${ppsLumi},eerpin${c}:${ppsLumi},emrpin${c}:${ppsLumi},arpin${c}:${lptappslumi}";
        done
	baseOpts="-i /eos/cms/${outdir}/analysis --lumiSpecs ${lumiSpecs} --procSF #gamma+jets:1.33 -l ${lumi} --mcUnc ${lumiUnc} ${lumiSpecs} ${kFactorList}"
        commonOpts="${baseOpts} -j ${samples_json} --signalJson ${plot_signal_json} -O /eos/cms/${outdir}/analysis/plots --saveLog"

        plots=xangle_eerpinhpur,xangle_mmrpinhpur,xangle_emrpinhpur,xangle_arpinhpur
        python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/plotter.py ${commonOpts} --only ${plots} --strictOnly --saveTeX --rebin 4;

        cats=(
            "" 
            "rpin" 
            "rpinpos" "rpinneg"
            "rpinhpur" 
            "rpinhpurneg" "rpinhpurpos" 
            "rpinhpur120xangle" "rpinhpur130xangle" "rpinhpur140xangle" "rpinhpur150xangle" 
        )
        for ch in mm ee a em; do
           plots=""

            for c in "${cats[@]}"; do            
                evcat=${ch}${c}
                
                for p in ptll mll nvtx rho xangle mll yll etall ptll l1eta l1pt l2eta l2pt acopl costhetacs; do
                    plots="${p}_${evcat},${plots}"
                done

                if [[ $evcat == *"rpin"* ]]; then
                    for p in nextramu extramupt extramueta mmass ppcount ypp2d ypp mpp mpp2d mmass_full; do
                        plots="${p}_${evcat},${plots}"
                    done
                fi
                
                if [[ $evcat == *"pos"* || $evcat == *"neg"* ]]; then
                    for p in ntk csi csi2d ntk; do
                        plots="${p}_${evcat},${plots}"
                    done
                fi            
            done
            python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/plotter.py ${commonOpts} --only ${plots} --strictOnly  -o plots/plotter_${ch}.root &
        done
        
        ;;

    WWW )

	mkdir -p ${wwwdir}/presel
	cp /eos/cms/${outdir}/plots/*.{png,pdf,dat} ${wwwdir}/presel
	cp $CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/index.php ${wwwdir}/presel

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
            echo "+JobFlavour =\"workday\"" >> $condor_prep
            for x in 120 130 140 150; do
                echo "arguments   = ${CMSSW_BASE} ${m} ${x} ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/analysis/stat_m${m}" >> $condor_prep
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
