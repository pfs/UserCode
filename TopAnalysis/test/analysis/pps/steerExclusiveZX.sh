#!/bin/bash

WHAT=${1}; 
ALLOWPIX=${2}

if [ "$#" -ne 2 ]; then 
    echo "steerExclusiveZX.sh <SEL/MERGE/PLOT/WWW> <ALLOWPIX>";
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
    echo "   {PREPARE,CHECK,RUN,SUMMARIZE}OPTIMSTATANA - optimize the statistical analysis"
    echo "   WWW             - move plots to web-based area"
    exit 1; 
fi

#to run locally use local as queue + can add "--njobs 8" to use 8 parallel jobs
queue=tomorrow

#max. pixels allowed per arm
if [[ $ALLOWPIX == *"1,2"* ]]; then
    pfix=""
else
    pfix="_1exc"
fi
if [[ $ALLOWPIX == *"-"* ]]; then
    pfix="${pfix}single"
fi


githash=2017_unblind_multi

datadir=/store/cmst3/group/top/RunIIUL/2017/6bfa3f2e/
datajson=${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/pps/datasamples.json
RPout_json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/golden_noRP.json

mcdir=/store/cmst3/group/top/RunIIReReco/ab05162
mcjson=${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/pps/mcsamples.json
cleanmcjson=${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/pps/mcsamples_nowqcd.json
zxjson=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/zx_samples.json
genweights=genweights_ab05162.root

sddir=/store/cmst3/group/top/RunIIReReco/2017/sdz
sdjson=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/pps_sd_samples.json
genweights_sd=genweights_sdz.root

signaldir=/store/cmst3/group/top/RunIIReReco/2017/vxsimulations_7Sep2020 
signaljson=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/signal_samples.json
signalpostts2json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/signal_samples_postTS2.json
fullsimsignaljson=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/signal_samples_fullsim.json

outdir=/store/cmst3/user/psilva/ExclusiveAna/final/${githash}

anadir=${outdir}/analysis
wwwdir=/eos/user/p/psilva/www/EXO-19-009/analysis${pfix}

plot_signal_json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/plot_signal_samples.json
plot_signal_ext_json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/plot_signal_samples_ext.json
plot_signal_sdz_json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/plot_signal_sdz_samples.json

zbias_samples_json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/zbias_samples.json

lumi=41529
ppsLumi=37193
lptalumi=2642
lptappslumi=2288
lumiUnc=0.025

RED='\e[31m'
NC='\e[0m'
case $WHAT in

    TESTSEL )
        
        inputfileTag=MC13TeV_2017_DY50toInf_fxfx #MC13TeV_2017_GJets_HT400to600
        inputfileTESTSEL=${mcdir}/${inputfileTag}/Chunk_0_ext0.root
        
        #inputfileTag=Data13TeV_2017B_SinglePhoton
        #inputfileTESTSEL=${datadir}/${inputfileTag}/Chunk_0_ext0.root

        python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/runLocalAnalysis.py \
            -i ${inputfileTESTSEL} --tag ${inputfileTag} \
            -o testsel.root --genWeights ${genweights} \
            --njobs 1 -q local --debug \
            --era era2017 -m ExclusiveZX::RunExclusiveZX --ch 0 --runSysts;
        ;;

    FULLSIMSIG )

        for ch in ee; do # ee mm; do
            inputfileTag=MC13TeV_Z${ch}_m_X_950_xangle_120_2017_postTS2_fullsim
            inputfileTESTSEL=${datadir}/${inputfileTag}/Chunk_0_ext0.root
            selout=/eos/cms/${outdir}/Chunks/${inputfileTag}_0.root

            python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/runLocalAnalysis.py \
                -i ${inputfileTESTSEL} --tag ${inputfileTag} \
                -o ${selout} --genWeights  genweights_f439f08_ul.root \
                --njobs 1 -q local \
                --era era2017 -m ExclusiveZX::RunExclusiveZX --ch 0 --runSysts;

            mix_file=/eos/cms/${anadir}/mixing/
            anaout=./mixtest${pfix}
            addOpt="--effDir test/analysis/pps"
            python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/runExclusiveAnalysis.py --step 1 --jobs 1 \
                --json ${fullsimsignaljson} --RPout ${RPout_json} -o ${anaout} \
                --mix ${mix_file} -i `dirname ${selout}` --only ${inputfileTag}_0.root ${addOpt} \
                --allowPix ${ALLOWPIX};
        done

        ;;


    FULLSIMSD )
        
        for tag in `ls /eos/cms/${sddir}`; do
            echo ${tag};

            inputfile=${sddir}/${tag}/Chunk_0_ext0.root
            selout=/eos/cms/${outdir}/Chunks/${tag}_0.root
            python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/runLocalAnalysis.py \
                -i ${inputfile} --tag ${tag} \
                -o ${selout} --genWeights  ../../test/analysis/pps/genweights_sdz.root \
                --njobs 1 -q local \
                --era era2017 -m ExclusiveZX::RunExclusiveZX --ch 0 --runSysts;
                    
        done

        step=1
        predin=/eos/cms/${outdir}/Chunks/
        predout=/eos/cms/${anadir}${pfix}
        condor_prep=runanasdsig${pfix}_condor.sub
        mix_file=/eos/cms/${anadir}/mixing/
        echo "executable  = ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/pps/wrapAnalysis.sh" > $condor_prep
        echo "output      = ${condor_prep}.out" >> $condor_prep
        echo "error       = ${condor_prep}.err" >> $condor_prep
        echo "log         = ${condor_prep}.log" >> $condor_prep
        # echo "requirements = (OpSysAndVer =?= \"SLCern6\")" >> $condor_prep
        echo "+AccountingGroup = \"group_u_CMST3.all\"" >> $condor_prep
        echo "+JobFlavour = \"tomorrow\"">> $condor_prep
        echo "request_cpus = 4" >> $condor_prep
        echo "arguments   = ${CMSSW_BASE} ${step} ${predout} ${predin} \$(chunk) ${mix_file} ${ALLOWPIX}" >> $condor_prep
        echo "queue chunk matching (${predin}/MC13TeV_SDZ*.root)" >> $condor_prep
        condor_submit $condor_prep        

        ;;

    
    SEL )
        baseCmd="$CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/runLocalAnalysis.py"
        baseCmd="${baseCmd} --genWeights ${genweights} --era era2017 -m ExclusiveZX::RunExclusiveZX --ch 0 --runSysts"
        baseCmd="${baseCmd} -o ${outdir} -q ${queue}"
        
	#python ${baseCmd} --farmappendix ZXData2017BSinglePhoton -i ${datadir} --only 2017B_SinglePhoton;
	#python ${baseCmd} --farmappendix ZXData -i ${datadir} --only ${datajson};
	python ${baseCmd} --farmappendix ZXMC   -i ${mcdir}   --only ${mcjson},${zxjson};
	;;

    CHECKSELINTEG )
        python scripts/checkLocalAnalysisInteg.py ../../../FARM${githash}ZXMC/   analysis/data tree
        #python scripts/checkLocalAnalysisInteg.py ../../../FARM${githash}ZXData/ analysis/data tree
        ;;

    MERGESEL )
	mergeOutputs.py /eos/cms/${outdir} True;
	;;

    PLOTSEL )
        lumiSpecs="--lumiSpecs a:${lptalumi}"
        kFactorList="--procSF #gamma+jets:1.4"        
	commonOpts="-i /eos/cms/${outdir} --puNormSF puwgtctr -l ${lumi} --mcUnc ${lumiUnc} ${kFactorList} ${lumiSpecs}"
	python scripts/plotter.py ${commonOpts} -j ${mcjson},${datajson}    -O /eos/cms/${outdir}/plots/sel -o plotter.root --only mll,pt,eta,met,jets,nvtx,ratevsrun; # --saveLog; 
        python scripts/plotter.py ${commonOpts} -j ${mcjson}    --rawYields --silent --only gen -O /eos/cms/${outdir}/plots/ -o plots/bkg_gen_plotter.root; 
	python scripts/plotter.py ${commonOpts} -j ${zxjson} --rawYields --silent --only gen -O /eos/cms/${outdir}/plots/ -o plots/zx_gen_plotter.root;
        python test/analysis/pps/computeDileptonSelEfficiency.py /eos/cms/${outdir}/plots/       
        cp -v /eos/cms/${outdir}/plots/effsummary_* test/analysis/pps/
	;;

    TESTPREPAREMIX )
        predin=/eos/cms/${outdir}/Chunks
        file=Data13TeV_2017D_DoubleMuon_2.root
        predout=/eos/cms/${anadir}/mixing
        mkdir -p $predout
        
        #run locally
        python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/runExclusiveAnalysis.py --step 0 --jobs 1 \
            --json ${mcjson},${datajson},${signaljson},${signalpostts2json} --RPout ${RPout_json} -o ${predout} -i ${predin} --only ${file} --maxEvents 10000
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
        file=Data13TeV_2017B_DoubleMuon_2.root
        
        predin=/eos/cms/${signaldir}
        file=Z_m_X_960_xangle_120_2017_preTS2.root
        predin=/eos/cms/${signaldir}
        file=Z_m_X_960_xangle_120_2017_postTS2.root
        
        mix_file=/eos/cms/${anadir}/mixing/

        predout=./mixtest${pfix}
        addOpt="--effDir test/analysis/pps"
        python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/runExclusiveAnalysis.py --step 1 --jobs 1 \
            --json ${mcjson},${datajson},${signaljson},${signalpostts2json} --RPout ${RPout_json} -o ${predout} \
            --mix ${mix_file} -i ${predin} --only ${file} ${addOpt} --maxEvents 1000 \
            --allowPix ${ALLOWPIX};

        #predout=./mix1200
        #addOpt="--mixSignal /eos/cms/${anadir}/Z_m_X_1200_xangle_{0}_2017_preTS2_opt_v1_simu_reco.root"
        #python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/runExclusiveAnalysis.py --step 1 --jobs 1 \
        #    --json ${mcjson},${datajson},${signaljson},${signalpostts2json} --RPout ${RPout_json} -o ${predout} --mix ${mix_file} -i ${predin} --only ${file} ${addOpt};
        ;;


    ANA )

        step=1
        predin=/eos/cms/${outdir}/Chunks
        predout=/eos/cms/${anadir}${pfix}
        condor_prep=runana${pfix}_condor.sub
        mix_file=/eos/cms/${anadir}/mixing/
        echo "executable  = ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/pps/wrapAnalysis.sh" > $condor_prep
        echo "output       = ${condor_prep}.out" >> $condor_prep
        echo "error        = ${condor_prep}.err" >> $condor_prep
        echo "log          = ${condor_prep}.log" >> $condor_prep
        #echo "+AccountingGroup = \"group_u_CMST3.all\"" >> $condor_prep
        echo "+JobFlavour = \"tomorrow\"">> $condor_prep
        echo "request_cpus = 4" >> $condor_prep
        echo "arguments    = ${CMSSW_BASE} ${step} ${predout} ${predin} \$(chunk) ${mix_file} ${ALLOWPIX}" >> $condor_prep
        echo "queue chunk matching (${predin}/*.root)" >> $condor_prep
        condor_submit $condor_prep

        #run locally
        #python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/runExclusiveAnalysis.py --step 1 --jobs 8 \
        #    --json ${mcjson},${datajson},${signaljson},${signalpostts2json} --RPout ${RPout_json} -o ${predout} --mix ${mix_file} -i ${predin} --allowPix ${ALLOWPIX};
        
        ;;

    CHECKANA )
        #0-just check
        #1-run locally
        #2-submit to condor
        python test/analysis/pps/checkFinalNtupleInteg.py /eos/cms/${outdir}/Chunks /eos/cms/${anadir}${pfix}/Chunks 2 1 /eos/cms/${anadir}/mixing/ ${ALLOWPIX}
        ;;
    
    ANASIG )

        step=1
        pred_inputs=(/eos/cms/${signaldir} /eos/cms//store/cmst3/group/top/RunIIReReco/2017/vxsimulations_1Jul2021)
        for predin in ${pred_inputs[@]}; do
            predout=/eos/cms/${anadir}${pfix}
            condor_prep=runanasig${pfix}_condor_`basename ${predin}`.sub
            mix_file=/eos/cms/${anadir}/mixing/
            echo "executable  = ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/pps/wrapAnalysis.sh" > $condor_prep
            echo "output      = ${condor_prep}.out" >> $condor_prep
            echo "error       = ${condor_prep}.err" >> $condor_prep
            echo "log         = ${condor_prep}.log" >> $condor_prep
            echo "+AccountingGroup = \"group_u_CMST3.all\"" >> $condor_prep
            echo "+JobFlavour = \"tomorrow\"">> $condor_prep
            echo "request_cpus = 4" >> $condor_prep
            echo "arguments   = ${CMSSW_BASE} ${step} ${predout} ${predin} \$(chunk) ${mix_file} ${ALLOWPIX}" >> $condor_prep
            echo "queue chunk matching (${predin}/*.root)" >> $condor_prep
            condor_submit $condor_prep
            echo "Submitted from $condor_prep"
        done
        
        #run locally
        #python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/runExclusiveAnalysis.py --step 1 --jobs 1 \
        #    --json ${signaljson},${signalpostts2json} --RPout ${RPout_json} --mix ${mix_file} \
        #    --allowPix ${ALLOWPIX} --effDir ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/pps \
        #    -i ${predin} -o anasig/;
        # cp -v anasig/Chunks/*.root /eos/cms/${outdir}/analysis/
        ;;
   
    CHECKANASIG )        
        pred_inputs=(/eos/cms/${signaldir} /eos/cms//store/cmst3/group/top/RunIIReReco/2017/vxsimulations_1Jul2021)
        for predin in ${pred_inputs[@]}; do
            python test/analysis/pps/checkFinalNtupleInteg.py ${predin} /eos/cms/${anadir}${pfix}/Chunks 2 1 /eos/cms/${anadir}/mixing ${ALLOWPIX}
        done
        ;;

    MERGEANA )
        mergeOutputs.py /eos/cms/${anadir}${pfix};
        ;;

    BKGVALIDATION )       

        indirForPlots=/eos/cms/${anadir}${pfix}
        #ptlist=(0 40)        
        ptlist=(40)
        baseOpts="-i ${indirForPlots} --doPerEra --doPerPU --doPerAngle --doPerNch"
        for pt in ${ptlist[@]}; do 
            output=${indirForPlots}/bkg_ptll${pt}
            python test/analysis/pps/doBackgroundValidation.py ${baseOpts} -o ${output} --selCuts "bosonpt>=${pt}&&l1pt>30&&l2pt>20"
        done

        ;;

    PLOTSIGACC )

        indirForPlots=/eos/cms/${anadir}${pfix}

        for m in 600 800 1200 1400; do 
            for i in `seq 2 5`; do 
                python test/analysis/pps/drawSignalPeak.py \
                    ${indirForPlots}/Z_m_X_${m}_xangle_1${i}0_2017_preTS2.root ; 
            done 
        done
        mkdir -p ${indirForPlots}/plots_signal
        mv mmass*sig* ${indirForPlots}/plots_signal
        ;;

    PLOTANA )

        indirForPlots=/eos/cms/${anadir}${pfix};

        lumiSpecs="mmrpin:${ppsLumi},eerpin:${ppsLumi},emrpin:${ppsLumi},a:${lptalumi},arpin:${lptappslumi}"   
        for c in hpur hpur1 hpur2 hpur3 hpur4 hpur0pos hpur0neg hpur1pos hpur1neg hpur2pos hpur2neg; do
            lumiSpecs="${lumiSpecs},mmrpin${c}:${ppsLumi},eerpin${c}:${ppsLumi},emrpin${c}:${ppsLumi},arpin${c}:${lptappslumi}";
        done

	baseOpts="-i ${indirForPlots} --lumiSpecs ${lumiSpecs} --procSF #gamma+jets:1.4 -l ${lumi} --mcUnc ${lumiUnc} ${lumiSpecs} ${kFactorList} --pformats png,pdf,C"
        commonOpts="${baseOpts} -j ${cleanmcjson},${datajson} --signalJson ${plot_signal_json} -O ${indirForPlots}/plots"
        
        plots=xangle_arpinhpur,xangle_eerpinhpur,xangle_mmrpinhpur,xangle_emrpinhpur,xangle_arpinhpur
        python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/plotter.py ${commonOpts} --only ${plots} --strictOnly --saveTeX --rebin 4;

        python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/plotter.py ${commonOpts} --only catcount --saveTeX;

        commonOpts="${baseOpts} -j ${cleanmcjson},${datajson} --signalJson ${plot_signal_json} -O ${indirForPlots}/plots"
        cats=(
            "" 
            "rpinhpur"      
            "rpinhpur0pos"  "rpinhpur1pos"  "rpinhpur2pos"
            "rpinhpur0neg"  "rpinhpur1neg"  "rpinhpur2neg"
            "rpinhpur1"     "rpinhpur2"     "rpinhpur3"     "rpinhpur4"            
        )
        channels=(mm) #(ee mm a) #mm ee a em offz)
        for ch in ${channels[@]}; do
           plots=""

            for c in "${cats[@]}"; do 

                evcat=${ch}${c};
                
                #proton-related plots
                if [[ $evcat == *"rpin"* ]]; then

                    for p in ypp; do # mmass_full mmass ppcount ypp mpp pzpp met njets mpf zjb zj2b nch; do
                        plots="${p}_${evcat},${plots}"
                    done

                    if [[ $evcat == *"pos"* || $evcat == *"neg"* ]]; then
                        for p in ntk csi; do
                            plots="${p}_${evcat},${plots}"
                        done
                    fi            

                #central-kinematics plots
                else
                    continue
                    for p in nextramu extramupt extramueta ptll mll nvtx rho xangle mll mll_full yll etall ptll ptll_high l1eta l1pt l2eta l2pt acopl costhetacs met njets mpf zjb zj2b nch; do
                        plots="${p}_${evcat},${plots}"
                    done
                fi
                
            done

            finalPlotOpts=${commonOpts}
            #if [[ $c == *"a"* ]]; then
            finalPlotOpts="${commonOpts} --normToData"
            #fi

            python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/plotter.py \
                ${finalPlotOpts} --only ${plots} --strictOnly  -o plotter_${ch}.root --saveLog &
        done        
        ;;

    PLOTANAPERERA )

        indirForPlots=/eos/cms/${anadir}${pfix};

        plots=""
        for evcat in a ee mm em; do
            if [[ $evcat == *"a"* ]]; then
                plist=(ptll)
            else
                plist=(ptll mll l1eta l2eta)
            fi
            for p in ${plist[@]}; do 
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
            echo $alumi $eralumi


	    baseOpts="-i ${indirForPlots} --lumiSpecs a:${alumi} --procSF #gamma+jets:1.4 -l ${eralumi} --mcUnc ${lumiUnc} ${lumiSpecs} ${kFactorList}"
            era_json=test/analysis/pps/test_samples_${era}.json;
            commonOpts="${baseOpts} -j ${era_json} --signalJson ${plot_signal_sdz_json} -O ${indirForPlots}/plots${era}"
            python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/plotter.py ${commonOpts} --only ${plots} --strictOnly;
        done
        ;;


    TESTOPTIMSTATANA)

        indirForPlots=/eos/cms/${anadir}${pfix}
        commonOpts="-i ${indirForPlots} --finalStates 22"

        echo "Generating a datacard takes a bit as it'll project the shapes for a given set of cuts"
        echo "You can run locally with python test/analysis/pps/generatedBinnedWorkspace.py and your preferred set of cuts"
        echo "Running with the default values for ${commonOpts} and output @ ppvx_${githash}${pfix}/test"
        
        #python test/analysis/pps/generateBinnedWorkspace.py ${commonOpts} -o ppvx_${githash}${pfix}/test --doBackground    
        python test/analysis/pps/generateBinnedWorkspace.py ${commonOpts} -o ppvx_${githash}${pfix}/test --massList 1200 

        #python test/analysis/pps/generateBinnedWorkspace.py ${commonOpts} -o ppvx_${githash}${pfix}/test --doDataCards 
        
        ;;

    PREPAREOPTIMSTATANA )
        indirForPlots=/eos/cms/${anadir}${pfix}
        python test/analysis/pps/prepareOptimScanCards.py -o ppvx_${githash}${pfix} -i ${indirForPlots}
        ;;

    CHECKOPTIMSTATANA )
        python test/analysis/pps/checkStatAnaInteg.py ppvx_${githash}${pfix}
        ;;

    RUNOPTIMSTATANA )
        
        tagList=(obs exp) #(exp expm1000 obs)
        if [[ $ALLOWPIX == *"1,2"* ]]; then
            tagList=(exp)
        elif [[ $ALLOWPIX == *"-"* ]]; then
            tagList=(exp)
        fi
        for t in ${tagList[@]}; do
            python test/analysis/pps/prepareFinalStatAnalysis.py -i ppvx_${githash}${pfix} -t ${t}
        done
        ;;

    SUMMARIZEOPTIMSTATANA )

        d=ppvx_${githash}${pfix}
        tagList=(obs exp) #(exp expm1000 obs)
        if [[ $ALLOWPIX == *"1,2"* ]]; then
            tagList=(exp)
        elif [[ $ALLOWPIX == *"-"* ]]; then
            tagList=(exp)
        fi
        for t in ${tagList[@]}; do
            python test/analysis/pps/summarizeOptimScanResults.py ${d}/${t}
        done
        echo "You can now use the jupyter notebook to analyze the results of the statistical analysis"
        ;;


    INJECTSIGNAL )

        #afs needs to be used here...
        pfix=_1exc
        indirForPlots=/eos/cms/${anadir}${pfix}
        mu=0.5
        for m in 600 800 1000 1200; do
            baseOpts="-i ${indirForPlots} --injectMass ${m} --injectMu ${mu} --just 0,1,2,3,4,5";
            python test/analysis/pps/prepareOptimScanCards.py -o ppvx_${githash}${pfix}_${m}_${mu}        ${baseOpts};
            python test/analysis/pps/prepareOptimScanCards.py -o ppvx_${githash}${pfix}_inc_${m}_${mu}    ${baseOpts} --xangles 0;
            #python test/analysis/pps/prepareOptimScanCards.py -o ppvx_${githash}${pfix}_signed_${m}_${mu}  ${baseOpts} --signed;
        done

        ;;

    UNBLIND )

        #the real thing
        python test/analysis/pps/prepareOptimScanCards.py -o ppvx_${githash}_obs        --unblind  -i /eos/cms/${anadir} --just 2,10,46;
        python test/analysis/pps/prepareOptimScanCards.py -o ppvx_${githash}_signed_obs --unblind  -i /eos/cms/${anadir}  --just 2,10,46 --signed;
        
        ;;


    WWW )

        indirForPlots=/eos/cms/${anadir}${pfix}
        index=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/index.php
        #cp ${index} ${wwwdir}
       
        #preselection
	#mkdir -p ${wwwdir}/presel
	#cp /eos/cms/${outdir}/plots/*.{png,pdf,dat} ${wwwdir}/presel
	#cp /eos/cms/${outdir}/plots/sel/*.{png,pdf,dat} ${wwwdir}/presel
	#cp $CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/index.php ${wwwdir}/presel
        
        #analysis control plots
        mkdir -p ${wwwdir}/ana
        cp ${indirForPlots}/plots/*.{png,pdf,dat} ${wwwdir}/ana
        cp ${index} ${wwwdir}/ana

        #breakdown per era (control yields etc.)
        for era in B C D E F; do
            odir=${wwwdir}/ana/2017${era};
            mkdir -p ${odir}
            cp ${indirForPlots}/plots${era}/*.{png,pdf,dat} ${odir}
            cp ${index} ${odir};
        done

        #signal-related plots
        #mkdir -p ${wwwdir}/signal
        #cp ${indirForPlots}/plots_signal/*.{png,pdf,dat} ${wwwdir}/signal
        #cp ${index} ${wwwdir}/signal

        #background closure plots
        for pt in 0 40; do
            pdir=${indirForPlots}/bkg_ptll${pt}
            odir=${wwwdir}/bkg/emu_ptll${pt}
            mkdir -p ${odir}
            cp ${pdir}/*.{png,pdf,dat} ${odir}
            cp ${index} ${odir}
        done

        ;;


    #CHECK THIS POINT FWD

    PLOTLOCALSENS )     

        indirForPlots=/eos/cms/${anadir}_1exc  
        odir=${indirForPlots}/localsens
        mkdir -p ${odir};
        for xangle in 120 130 140 150; do
            python test/analysis/pps/estimateLocalSensitivity.py --xangle ${xangle} -i ${indirForPlots} -o ${odir} &            
        done
        ;;



    EXTRA )
        
        scriptDir=test/analysis/pps
        echo ""
        echo "This is a printout of additional scripts used for cross checks or specific plots"
        echo "Scripts can be found in ${scriptDir}"
        echo "[Additional supporting plots]"
        echo "python test/analysis/pps/compareNvtx.py"
        echo ""
        echo "[Additional figures]"
        echo "python test/analysis/pps/drawPPSEfficiency.py"
        echo "for m in 800 1000 1400 1600; do python test/analysis/pps/showMigrationMatrix.py Z_m_X_${m} /eos/cms/store/cmst3/user/psilva/ExclusiveAna/final/2017_unblind_multi/analysis_1exc/; done"
        echo "python test/analysis/pps/checkPPSperEra.py fill"
        echo "python test/analysis/pps/checkPPSperEra.py plot"
        echo "for m in 600 800 1000 1200 1400; do python test/analysis/pps/drawSignalPeak.py /eos/cms//store/cmst3/user/psilva/ExclusiveAna/final/2017_unblind_multi/analysis_1exc/Z_m_X_\${m}_xangle_{0}_2017_{1}.root; done"
        echo "[Display fit shapes/uncertainties]"
        echo "python test/analysis/pps/showUncs.py > statana_zmm_m1000_uncs.dat"
        echo "python test/analysis/pps/showFitShapes.py ppvx_${githash}${pfix}/optim_0  -m 1000 -b z -t mm -u -s"
        echo "python test/analysis/pps/showFitShapes.py ppvx_${githash}${pfix}/optim_25 -m 1000 -b z -t ms -u -s"
        echo "python test/analysis/pps/showFitShapes.py ppvx_${githash}${pfix}/optim_50 -m 1000 -b z -t sm -u -s"
        echo "python test/analysis/pps/showFitShapes.py ppvx_${githash}${pfix}/optim_75 -m 1000 -b z -t ss -u -s"
        echo "[Nuisances post-fit, gof, etc]"
        echo "python test/analysis/pps/doNuisanceReport.py B-only:b:ppvx_2017_unblind_multi_1exc/exp/inc_xangle_nvtx/fitDiagnosticsPPzX.m1000.root S+B:s:ppvx_2017_unblind_multi_1exc/exp/inc_xangle_nvtx/fitDiagnosticsPPzX.m1000.root --tag asimovPPzX"
        echo "python test/analysis/pps/doNuisanceReport.py Obs.:s:ppvx_${githash}${pfix}/obs/inc_xangle_nvtx/fitDiagnosticsPPzX.m1000.root Exp.:s:ppvx_${githash}${pfix}/exp/inc_xangle_nvtx/fitDiagnosticsPPzX.m1000.root --tag PPzX"
        echo "python test/analysis/pps/doNuisanceReport.py Obs.:s:ppvx_${githash}${pfix}/obs/inc_xangle_nvtx/fitDiagnosticsPPgX.m1000.root Exp.:s:ppvx_${githash}${pfix}/exp/inc_xangle_nvtx/fitDiagnosticsPPgX.m1000.root --tag PPgX"
        echo "python test/analysis/pps/showNuisancePerCateg.py ppvx_${githash}${pfix}/obs/inc_xangle_nvtx/fitDiagnosticsPPzX.m1000.root"
        echo "python test/analysis/pps/plotGOF.py ppvx_${githash}${pfix}/obs/inc_xangle_nvtx/higgsCombinePPzX.m1000.gof.GoodnessOfFit.mH1000.root"
        echo "[Postfit plots]"
        echo "python test/analysis/pps/plotPostFitDistributions.py ppvx_2017_unblind_multi_1exc/obs/inc_xangle_nvtx/fitDiagnosticsPPzX.m1000.root"
        echo "python test/analysis/pps/plotPostFitDistributions.py ppvx_2017_unblind_multi_1exc/obs/inc_xangle_nvtx/fitDiagnosticsPPgX.m1000.root"
        echo "[Acceptance plots inputs]"
        echo "python test/analysis/pps/computeFinalAEff.py (takes ~1h to loop over original signal files)"
        echo "python test/analysis/pps/computeFinalAEff.py test/analysis/pps/acc_summary.dat ppvx_2017_unblind_multi_1exc/exp/inc_xangle_nvtx/info.dat"
        ;;
    

    ###                                                                                         ###
    ### THESE OPTIONS WERE USED AT SOME POINT IN THE ANALYSIS BUT HAVE NOT BEEN TESTED RECENTLY ###
    ###                                                                                         ###

    TESTSYNCH )
        
        python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/runLocalAnalysis.py \
            -i /store/cmst3/user/psilva/ExclusiveAna/synch/Data13TeV_2017_DoubleMuon_synch.root \
            -o testsynch.root --genWeights ${genweights} \
            --njobs 1 -q local --debug \
            --era era2017 -m ExclusiveZX::RunExclusiveZX --ch 0 --runSysts;

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








esac
