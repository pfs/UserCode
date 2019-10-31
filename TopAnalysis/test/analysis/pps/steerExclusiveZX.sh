#!/bin/bash

WHAT=$1; 
if [ "$#" -ne 1 ]; then 
    echo "steerExclusiveZX.sh <SEL/MERGE/PLOT/WWW>";
    echo "   SEL           - launches selection jobs to the batch, output will contain summary trees and control plots"; 
    echo "   MERGESEL        - merge output"
    echo "   PLOTSEL         - make plots"
    echo "   TRAINPUDISCR    - train pileup discriminator"
    echo "   RUNPRED         - run pileup discriminator prediction"
    echo "   PREPAREMIX      - prepare bank of events for event mixing from the summary trees"
    echo "   COLLECTMIX      - collects all the mixing events found in PREPAREANA"
    echo "   ANA/ANASIG      - run analysis on the summary trees"
    echo "   CHECKANA        - check analysis integrity and re-run locally jobs which have failed"
    echo "   PLOTANA         - plot analysis results"
    echo "   OPTIMSTATANA    - optimize the statistical analysis"
    echo "   FINALIZESTATANA - define the final datacards based on the results of the optimization"
    echo "   WWW             - move plots to web-based area"
    exit 1; 
fi

#to run locally use local as queue + can add "--njobs 8" to use 8 parallel jobs
queue=tomorrow
githash=ab05162
eosdir=/store/cmst3/group/top/RunIIReReco/${githash}
#githash=b949800_newppscalib
#eosdir=/store/cmst3/group/top/RunIIReReco/2017/${githash}
signal_dir=/store/cmst3/group/top/RunIIReReco/2017/vxsimulations_30Sep
outdir=/store/cmst3/user/psilva/ExclusiveAna/final/${githash}
#anadir=${outdir}/analysis_0p05
anadir=${outdir}/analysis_0p04
signal_json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/signal_samples.json
plot_signal_json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/plot_signal_samples.json
plot_signal_ext_json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/plot_signal_samples_ext.json
samples_json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/samples.json
jetht_samples_json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/jetht_samples.json
zx_samples_json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/zx_samples.json
zbias_samples_json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/zbias_samples.json
RPout_json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/golden_noRP.json
wwwdir=/eos/user/p/psilva/www/ExclusiveAna_0p04
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
        kFactorList="--procSF #gamma+jets:1.4"
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

    TESTPREPAREMIX )

        predin=/eos/cms/${outdir}/Chunks
        file=Data13TeV_2017D_DoubleMuon_2.root
        predout=/eos/cms/${anadir}/mixing
        mkdir -p $predout
        
        #run locally
        python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/runExclusiveAnalysis.py --step 0 --jobs 1 \
            --json ${samples_json},${signal_json} --RPout ${RPout_json} -o ${predout} -i ${predin} --only ${file} --maxEvents 10000;
        
        ;;

    PREPAREMIX )       
        step=0
        predin=/eos/cms/${outdir}/Chunks
        predout=/eos/cms/${anadir}/mixing
        condor_prep=runana_condor.sub
        echo "executable  = ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/pps/wrapAnalysis.sh" > $condor_prep
        echo "output      = ${condor_prep}.out" >> $condor_prep
        echo "error       = ${condor_prep}.err" >> $condor_prep
        echo "log         = ${condor_prep}.log" >> $condor_prep
        echo "arguments   = ${CMSSW_BASE} ${step} ${predout} ${predin} \$(chunk)" >> $condor_prep
        echo "queue chunk matching (${predin}/Data*.root)" >> $condor_prep
        condor_submit $condor_prep
        
        #python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/runExclusiveAnalysis.py --step 0 --jobs 8 \
        #    --json ${samples_json} --RPout ${RPout_json} -o ${predout} -i ${predin};

        ;;


    CHECKMIX )
        python test/analysis/pps/checkFinalNtupleInteg.py /eos/cms/${outdir}/Chunks /eos/cms/${anadir}/mixing/Chunks 1 0
        ;;

    COLLECTMIX )
        python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/collectEventsForMixing.py /eos/cms/${anadir}
        ;;


    TESTANA )

        predin=/eos/cms/${outdir}/Chunks
        #file=Data13TeV_2017D_DoubleEG_2.root
        #file=MC13TeV_2017_DY50toInf_fxfx_0.root
        file=Data13TeV_2017C_DoubleEG_0.root

        predin=/eos/cms/${signal_dir}
        file=Z_m_X_1440_xangle_130_2017_preTS2.root
        file=gamma_m_X_1440_xangle_130_2017_preTS2.root
        
        predout=/eos/cms/${anadir}
        #mix_file=/eos/cms/${anadir}/mixing/mixbank.pck
        mix_file=/eos/cms//store/cmst3/user/psilva/ExclusiveAna/final/ab05162/mixing/mixbank.pck
        #run locally
        python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/runExclusiveAnalysis.py --step 1 --jobs 1 \
            --json ${samples_json},${signal_json} --RPout ${RPout_json} -o ${predout} --mix ${mix_file} -i ${predin} --only ${file} --maxEvents 10000;
        
        ;;


    ANA )
        step=1
        predin=/eos/cms/${outdir}/Chunks
        predout=/eos/cms/${anadir}
        condor_prep=runana_condor.sub
        mix_file=/eos/cms/${anadir}/mixing/mixbank.pck
        echo "executable  = ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/pps/wrapAnalysis.sh" > $condor_prep
        echo "output       = ${condor_prep}.out" >> $condor_prep
        echo "error        = ${condor_prep}.err" >> $condor_prep
        echo "log          = ${condor_prep}.log" >> $condor_prep
        echo "+JobFlavour = \"tomorrow\"">> $condor_prep
        echo "request_cpus = 4" >> $condor_prep
        echo "arguments    = ${CMSSW_BASE} ${step} ${predout} ${predin} \$(chunk) ${mix_file}" >> $condor_prep
        echo "queue chunk matching (${predin}/*.root)" >> $condor_prep
        condor_submit $condor_prep

        #run locally
        #python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/runExclusiveAnalysis.py --step 1 --jobs 8 \
        #    --json ${samples_json} --RPout ${RPout_json} -o ${predout} --mix ${mix_file} -i ${predin};
        
        ;;

    ANASIG )
        step=1
        predin=/eos/cms/${signal_dir}
        predout=/eos/cms/${anadir}
        condor_prep=runanasig_condor.sub
        mix_file=/eos/cms/${anadir}/mixing/mixbank.pck
        echo "executable  = ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/pps/wrapAnalysis.sh" > $condor_prep
        echo "output      = ${condor_prep}.out" >> $condor_prep
        echo "error       = ${condor_prep}.err" >> $condor_prep
        echo "log         = ${condor_prep}.log" >> $condor_prep
        echo "+JobFlavour = \"tomorrow\"">> $condor_prep
        echo "request_cpus = 4" >> $condor_prep
        echo "arguments   = ${CMSSW_BASE} ${step} ${predout} ${predin} \$(chunk)" ${mix_file} >> $condor_prep
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
        python test/analysis/pps/checkFinalNtupleInteg.py /eos/cms/${outdir}/Chunks /eos/cms/${anadir}/Chunks 2 1 /eos/cms/${anadir}/mixing/mixbank.pck
        ;;
   
    CHECKANASIG )
        python test/analysis/pps/checkFinalNtupleInteg.py /eos/cms/${signal_dir} /eos/cms/${anadir}/Chunks 1 1 /eos/cms/${anadir}/mixing/mixbank.pck
        ;;

    MERGEANA )
        mergeOutputs.py /eos/cms/${anadir};
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

	    baseOpts="-i /eos/cms/${anadir} --lumiSpecs a:${alumi} --procSF #gamma+jets:1.4 -l ${eralumi} --mcUnc ${lumiUnc}"
            era_json=test/analysis/pps/test_samples_${era}.json;
            commonOpts="${baseOpts} -j ${era_json} --signalJson ${plot_signal_json} -O /eos/cms/${anadir}/plots${era}"
            python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/plotter.py ${commonOpts} --only ${plots} --strictOnly;
        done
        ;;

    PLOTSENS )

        for xangle in 120 130 140 150; do
            python test/analysis/pps/estimateLocalSensitivity.py \
                -i /eos/cms/${anadir} \
                -o /eos/cms/${anadir}/plots \
                --xangle ${xangle};
        done
        ;;

    PLOTANA )

        lumiSpecs="mmrpin:${ppsLumi},eerpin:${ppsLumi},emrpin:${ppsLumi},a:${lptalumi},arpin:${lptappslumi},arpinhpur:${lptappslumi}"   
        for c in neg pos hpur hpurpos hpurneg hpur120xangle hpur130xangle hpur140xangle hpur150xangle; do
            lumiSpecs="${lumiSpecs},mmrpin${c}:${ppsLumi},eerpin${c}:${ppsLumi},emrpin${c}:${ppsLumi},arpin${c}:${lptappslumi}";
        done
	baseOpts="-i /eos/cms/${anadir} --lumiSpecs ${lumiSpecs} --procSF #gamma+jets:1.4 -l ${lumi} --mcUnc ${lumiUnc} ${lumiSpecs} ${kFactorList}"
        commonOpts="${baseOpts} -j ${samples_json} --signalJson ${plot_signal_ext_json} -O /eos/cms/${anadir}/plots"

        plots=xangle_arpinhpur #xangle_eerpinhpur,xangle_mmrpinhpur,xangle_emrpinhpur,xangle_arpinhpur
        python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/plotter.py ${commonOpts} --only ${plots} --strictOnly --saveTeX --rebin 4;

        commonOpts="${baseOpts} -j ${samples_json} --signalJson ${plot_signal_json} -O /eos/cms/${anadir}/plots"
        cats=(
            "" 
            "rpin"
            "rpinpos" "rpinneg"
            "rpinhpur" 
            "rpinhpurneg" "rpinhpurpos" 
            "rpinhpur120xangle" "rpinhpur130xangle" "rpinhpur140xangle" "rpinhpur150xangle" 
        )
        for ch in a; do #mm ee a em; do
           plots=""

            for c in "${cats[@]}"; do            
                evcat=${ch}${c}
                
                for p in ptll mll nvtx rho xangle mll yll etall ptll l1eta l1pt l2eta l2pt acopl costhetacs met njets mpf zjb zj2b nch; do                  
                    plots="${p}_${evcat},${plots}"
                done

                if [[ $evcat == *"rpin"* ]]; then
                    for p in nextramu extramupt extramueta mmass ppcount ypp2d ypp mpp pzpp mpp2d mmass_full; do
                        plots="${p}_${evcat},${plots}"
                    done
                fi
                
                if [[ $evcat == *"pos"* || $evcat == *"neg"* ]]; then
                    for p in ntk csi csi2d; do
                        plots="${p}_${evcat},${plots}"
                    done
                fi            
            done
            python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/plotter.py ${commonOpts} --only ${plots} --strictOnly  -o plots/plotter_${ch}.root --saveLog &
        done
        
        ;;

    BKGVALIDATION )
        python test/analysis/pps/doBackgroundValidation.py --doPerAngle -i /eos/cms/${anadir} -o /eos/cms/${outdir}/bkg &
        python test/analysis/pps/doBackgroundValidation.py --doPerAngle -i /eos/cms/${anadir} -o /eos/cms/${outdir}/bkg_ptll0 --selCuts "bosonpt>=0. && l1pt>30 && l2pt>20" &
        python test/analysis/pps/doBackgroundValidation.py --doPerAngle -i /eos/cms/${anadir} -o /eos/cms/${outdir}/bkg_ptll60 --selCuts "bosonpt>60 && l1pt>30 && l2pt>20" &
        ;;

    WWW )

	cp $CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/index.php ${wwwdir}
	mkdir -p ${wwwdir}/presel
	cp /eos/cms/${outdir}/plots/*.{png,pdf,dat} ${wwwdir}/presel
	cp /eos/cms/${outdir}/plots/sel/*.{png,pdf,dat} ${wwwdir}/presel
	cp $CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/index.php ${wwwdir}/presel

        mkdir -p ${wwwdir}/ana
        cp /eos/cms/${anadir}/plots/*.{png,pdf,dat} ${wwwdir}/ana
        cp $CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/index.php ${wwwdir}/ana

        mkdir -p ${wwwdir}/bkg
        cp /eos/cms/${anadir}/plots/*.{png,pdf,dat} ${wwwdir}/bkg
        cp $CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/index.php ${wwwdir}/bkg

        ;;

    TESTSTATANA)

        echo "Generating a datacard takes a bit as it'll project the shapes for a given set of cuts"
        echo "You can run locally with python test/analysis/pps/generatedBinnedWorkspace.py and your preferred set of cuts"
        echo "Running with the default values for ${anadir} and output @ ppvx_analysis/test"
        python test/analysis/pps/generatedBinnedWorkspace.py -i /eos/cms/${anadir} -o ppvx_analysis/test

        ;;

    OPTIMSTATANA )

        #afs needs to be used here...
        python test/analysis/pps/prepareOptimScanCards.py -o ppvx_analysis -i /eos/cms/${anadir}

        ;;

    FINALIZESTATANA )
        #python test/analysis/pps/compareOptimResults.py ppvx_analysis/
        python test/analysis/pps/compareOptimResults.py ppvx_analysis/ 45
        ;;

    
    INJECTSIGNAL )

        #afs needs to be used here...
        for m in 800 1200 1400; do
            python test/analysis/pps/prepareOptimScanCards.py -o ppvx_analysis_${m}  -i /eos/cms/${anadir} --injectMass ${m} --just 2,11,37;
        done

        ;;


esac
