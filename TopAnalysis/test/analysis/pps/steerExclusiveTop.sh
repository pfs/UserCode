#!/bin/bash

WHAT=$1; 
if [ "$#" -ne 1 ]; then 
    echo "steerExclusiveTop.sh <SEL/MERGE/PLOT/WWW>";
    echo "   SEL        - launches selection jobs to the batch, output will contain summary trees and control plots"; 
    echo "   MERGESEL   - merge output"
    echo "   PLOTSEL    - make plots"
    echo "   WWWSEL     - move plots to web-based are"
    echo "   PREPAREANA - prepare bank of events for event mixing from the summary trees"
    echo "   ANA        - run analysis on the summary trees"
    echo "   PLOTANA    - plot analysis results"
    exit 1; 
fi

#to run locally use local as queue + can add "--njobs 8" to use 8 parallel jobs
queue=tomorrow
#githash=3129835
githash=6816e40
eosdir=/store/cmst3/group/top/RunIIReReco/${githash}
outdir=/store/cmst3/user/psilva/ExclusiveAna/${githash}
signal_json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/signal_samples.json
plot_signal_json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/plot_signal_samples.json
samples_json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/samples.json
jetht_samples_json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/jetht_samples.json
zx_samples_json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/zx_samples.json
RPout_json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/golden_noRP.json
wwwdir=~/www/ExclusiveAna
inputfileTag=MC13TeV_2017_GGH2000toZZ2L2Nu
inputfileTag=MC13TeV_2017_GGToEE_lpair
inputfileTag=Data13TeV_2017F_MuonEG
#inputfileTag=MC13TeV_2017_GJets_HT400to600
inputfileTag=Data13TeV_2017B_DoubleMuon
inputfileTESTSEL=${eosdir}/${inputfileTag}/Chunk_0_ext0.root
lumi=41833
ppsLumi=37500
lumiUnc=0.025


RED='\e[31m'
NC='\e[0m'
case $WHAT in

    TESTSEL )
        
        python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/runLocalAnalysis.py \
            -i ${inputfileTESTSEL} --tag ${inputfileTag} \
            -o testsel.root --genWeights genweights_${githash}.root \
            --njobs 1 -q local --debug \
            --era era2017 -m ExclusiveTop::RunExclusiveTop --ch 0 --runSysts;
        ;;

    SEL )
        baseOpt="-i ${eosdir} --genWeights genweights_${githash}.root"
        baseOpt="${baseOpt} -o ${outdir} -q ${queue} --era era2017 -m ExclusiveTop::RunExclusiveTop --ch 0 --runSysts"
        baseOpt="${baseOpt} --exactonly"        
	python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/runLocalAnalysis.py ${baseOpt} --only ${samples_json};
	#python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/runLocalAnalysis.py ${baseOpt} --only ${zx_samples_json};	
	;;

    MERGESEL )
	mergeOutputs.py /eos/cms/${outdir} True;
	;;

    PLOTSEL )
        kFactors="--procSF MC13TeV_2017_QCDEM_15to20:1.26,MC13TeV_2017_QCDEM_20to30:1.26,MC13TeV_2017_QCDEM_30to50:1.26,MC13TeV_2017_QCDEM_50to80:1.26,MC13TeV_2017_QCDEM_80to120:1.26,MC13TeV_2017_QCDEM_120to170:1.26,MC13TeV_2017_QCDEM_170to300:1.26,MC13TeV_2017_QCDEM_300toInf:1.26,MC13TeV_2017_GJets_HT40to100:1.26,MC13TeV_2017_GJets_HT100to200:1.26,MC13TeV_2017_GJets_HT200to400:1.26,MC13TeV_2017_GJets_HT600toInf:1.26"

	commonOpts="-i /eos/cms/${outdir} --puNormSF puwgtctr -l ${lumi} --mcUnc ${lumiUnc} ${kFactors}"
	python scripts/plotter.py ${commonOpts} -j ${samples_json}    -O plots/sel --only mboson,mtboson,pt,eta,met,jets,nvtx,ratevsrun --saveLog; 
        python scripts/plotter.py ${commonOpts} -j ${samples_json}    --rawYields --silent --only gen -O plots/zx_sel -o bkg_plotter.root ; 
	python scripts/plotter.py ${commonOpts} -j ${zx_samples_json} --rawYields --silent --only gen -O plots/zx_sel;
        python test/analysis/pps/computeDileptonSelEfficiency.py 
        mv *.{png,pdf} plots/zx_sel/
        mv effsummary* plots/
	;;

    WWWSEL )
	mkdir -p ${wwwdir}/presel
	cp plots/sel/*.{png,pdf} ${wwwdir}/presel
	cp plots/zx_sel/*.{png,pdf} ${wwwdir}/presel
	cp $CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/index.php ${wwwdir}/presel
	;;

    PREPAREANA )
        python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/runExclusiveAnalysis.py --step 0 --jobs 8 \
            -i /eos/cms/${outdir}/Chunks \
            --json ${samples_json} --RPout ${RPout_json} -o plots/analysis;               
        python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/collectEventsForMixing.py plots/analysis/Chunks
        ;;

    ANA )
        
        python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/runExclusiveAnalysis.py --step 1 --jobs 8 \
            -i /eos/cms/${outdir}/Chunks \
            --json ${samples_json} --RPout ${RPout_json} -o plots/analysis --mix plots/analysis/Chunks/mixbank.pck;
        
        ;;

    ANASIG )
        
        python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/runExclusiveAnalysis.py --step 1 --jobs 8 \
            --json ${signal_json} --RPout ${RPout_json} -o plots/analysis --mix plots/analysis/Chunks/mixbank.pck \
            -i /eos/cms/store/cmst3/group/top/PPSZX;

        python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/runExclusiveAnalysis.py --step 1 --jobs 8 \
            --json ${zx_samples_json} --RPout ${RPout_json} -o plots/analysis --mix plots/analysis/Chunks/mixbank.pck;

        ;;

    MERGEANA )
        mergeOutputs.py plots/analysis;
        ;;

    PLOTANA )

        
        plots=xangle_ee,xangle_mm,xangle_eeZ,xangle_mmZ,xangle_eehptZ,xangle_mmhptZ
        commonOpts="-i plots/analysis -j ${zx_samples_json}  -l ${ppsLumi} --mcUnc ${lumiUnc} -O plots/analysis/zx_yields" 
        #python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/plotter.py ${commonOpts} --only ${plots} --strictOnly --saveTeX --rebin 4;

        plots=xangle_eeZ,xangle_mmZ,xangle_eehptZ,xangle_mmhptZ,xangle_eehptZelpphighPur,xangle_mmhptZelpphighPur
        commonOpts="-i plots/analysis -j ${samples_json}  -l ${ppsLumi} --mcUnc ${lumiUnc}" # --normToData"
        #python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/plotter.py ${commonOpts} --only ${plots} --strictOnly --saveTeX --rebin 4;
        
        
        plots=mll_ee,mll_mm,mll_em,ptll_eeZ,ptll_mmZ,ptll_em
        plots=${plots},l1pt_eeZ,l1pt_mmZ,l1pt_em,l2pt_eeZ,l2pt_mmZ,l2pt_em
        plots=${plots},acopl_eeZ,acopl_mmZ,acopl_em,l1eta_eeZ,l1eta_mmZ,l1eta_em,l2eta_eeZ,l2eta_mmZ,l2eta_em
        plots=${plots},acopl_eeZelpphighPur140,acopl_mmZelpphighPur140,ptll_eeZelpphighPur140,ptll_mmZelpphighPur140
        plots=${plots},l1pt_eeZelpphighPur140,l1pt_mmZelpphighPur140,l2pt_eeZelpphighPur140,l2pt_mmZelpphighPur140
        plots=${plots},nvtx_eeZ,nvtx_mmZ,nvtx_em,xangle_eeZ,xangle_mmZ,xangle_em
        #python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/plotter.py ${commonOpts} --only ${plots} --strictOnly --signalJson ${plot_signal_json} --saveLog;

        #plots=ntk_eeZpos,ntk_eeZneg,ntk_mmZpos,ntk_mmZneg
        #plots=${plots},ntk_eeZelpp120pos,ntk_eeZelpp130pos,ntk_eeZelpp140pos,ntk_eeZelpp150pos
        #plots=${plots},ntk_eeZelpp120neg,ntk_eeZelpp130neg,ntk_eeZelpp140neg,ntk_eeZelpp150neg
        #plots=${plots},ntk_mmZelpp120pos,ntk_mmZelpp130pos,ntk_mmZelpp140pos,ntk_mmZelpp150pos
        #plots=${plots},ntk_mmZelpp120neg,ntk_mmZelpp130neg,ntk_mmZelpp140neg,ntk_mmZelpp150neg
        ##plots=${plots},csi_eeZelpp120pos,csi_eeZelpp130pos,csi_eeZelpp140pos,csi_eeZelpp150pos
        ##plots=${plots},csi_eeZelpp120neg,csi_eeZelpp130neg,csi_eeZelpp140neg,csi_eeZelpp150neg
        #plots=${plots},csi_eeZelpphighPur120pos,csi_eeZelpphighPur130pos,csi_eeZelpphighPur140pos,csi_eeZelpphighPur150pos
        #plots=${plots},csi_eeZelpphighPur120neg,csi_eeZelpphighPur130neg,csi_eeZelpphighPur140neg,csi_eeZelpphighPur150neg
        ##plots=${plots},csi_mmZelpp120pos,csi_mmZelpp130pos,csi_mmZelpp140pos,csi_mmZelpp150pos
        ##plots=${plots},csi_mmZelpp120neg,csi_mmZelpp130neg,csi_mmZelpp140neg,csi_mmZelpp150neg
        #plots=${plots},csi_mmZelpphighPur120pos,csi_mmZelpphighPur130pos,csi_mmZelpphighPur140pos,csi_mmZelpphighPur150pos
        #plots=${plots},csi_mmZelpphighPur120neg,csi_mmZelpphighPur130neg,csi_mmZelpphighPur140neg,csi_mmZelpphighPur150neg
        #plots=${plots},mpp_eeZelpphighPur120,mpp_eeZelpphighPur130,mpp_eeZelpphighPur140,mpp_eeZelpphighPur150
        #plots=${plots},mmass_eeZelpphighPur120,mmass_eeZelpphighPur130,mmass_eeZelpphighPur140,mmass_eeZelpphighPur150
        #plots=${plots},mpp_mmZelpphighPur120,mpp_mmZelpphighPur130,mpp_mmZelpphighPur140,mpp_mmZelpphighPur150
        #plots=${plots},mmass_mmZelpphighPur120,mmass_mmZelpphighPur130,mmass_mmZelpphighPur140,mmass_mmZelpphighPur150
        #plots=${plots},nch_mmZelpphighPur120,nch_mmZelpphighPur130,nch_mmZelpphighPur140,nch_mmZelpphighPur150
        #plots=${plots},njets_mmZelpphighPur120,njets_mmZelpphighPur130,njets_mmZelpphighPur140,nch_mmZelpphighPur150
        #plots=${plots},nvtx_mmZelpphighPur120,nvtx_mmZelpphighPur130,nvtx_mmZelpphighPur140,nvtx_mmZelpphighPur150
        #plots=${plots},nch_eeZelpphighPur120,nch_eeZelpphighPur130,nch_eeZelpphighPur140,nch_eeZelpphighPur150
        #plots=${plots},njets_eeZelpphighPur120,njets_eeZelpphighPur130,njets_eeZelpphighPur140,njets_eeZelpphighPur150
        #plots=${plots},ntk_eeZelpphighPur120,ntk_eeZelpphighPur130,ntk_eeZelpphighPur140,ntk_eeZelpphighPur150
        
        plots=""
        dists=(nextramu extramupt extramueta nvtx mpp ptll yll acopl)
        cats=(mmZ mmZelpphighPur mmZelpphighPur)
        xangles=(120 130 140 150)
        for d in ${dists[@]}; do
            for c in ${cats[@]}; do
                for a in ${xangles[@]}; do
                    plots=${plots},${d}_${c}${a},${d}_${c}${a}noextramu
                done
            done
        done
       
        #plots=${plots},nextramu_mmZelpphighPur120,nextramu_mmZelpphighPur130,nextramu_mmZelpphighPur140,nextramu_mmZelpphighPur150
        #plots=${plots},extramupt_mmZelpphighPur120,extramupt_mmZelpphighPur130,extramupt_mmZelpphighPur140,extramupt_mmZelpphighPur150
        #plots=${plots},extramueta_mmZelpphighPur120,extramueta_mmZelpphighPur130,extramueta_mmZelpphighPur140,extramueta_mmZelpphighPur150
        #plots=${plots},mpp_mmZelpphighPur120,mpp_mmZelpphighPur130,mpp_mmZelpphighPur140,mpp_mmZelpphighPur150
        #plots=${plots},ptll_mmZelpphighPur120,ptll_mmZelpphighPur130,ptll_mmZelpphighPur140,ptll_mmZelpphighPur150
        #plots=${plots},nvtx_mmZelpphighPur120,nvtx_mmZelpphighPur130,nvtx_mmZelpphighPur140,nvtx_mmZelpphighPur150
        #plots=${plots},ntk_mmZelpp120pos,ntk_mmZelpp130pos,ntk_mmZelpp140pos,ntk_mmZelpp150pos
        #plots=${plots},ntk_mmZelpp120neg,ntk_mmZelpp130neg,ntk_mmZelpp140neg,ntk_mmZelpp150neg
        #plots=${plots},sgny_mmZelpphighPur120,sgny_mmZelpphighPur130,sgny_mmZelpphighPur140,sgny_mmZelpphighPur150
        #plots=${plots},deltaenfwd_mmZelpphighPur120,deltaenfwd_mmZelpphighPur130,deltaenfwd_mmZelpphighPur140,deltaenfwd_mmZelpphighPur150
        #plots=${plots},ptll_mmZelpphighPur120noextramu,ptll_mmZelpphighPur130noextramu,ptll_mmZelpphighPur140noextramu,ptll_mmZelpphighPur150noextramu
        #plots=${plots},nvtx_mmZelpphighPur120noextramu,nvtx_mmZelpphighPur130noextramu,nvtx_mmZelpphighPur140noextramu,nvtx_mmZelpphighPur150noextramu
        #plots=${plots},sgny_mmZelpphighPur120noextramu,sgny_mmZelpphighPur130noextramu,sgny_mmZelpphighPur140noextramu,sgny_mmZelpphighPur150noextramu
        #plots=${plots},deltaenfwd_mmZelpphighPur120noextramu,deltaenfwd_mmZelpphighPur130noextramu,deltaenfwd_mmZelpphighPur140noextramu,deltaenfwd_mmZelpphighPur150noextramu
        #plots=${plots},ntk_mmZelpp120noextramupos,ntk_mmZelpp130noextramupos,ntk_mmZelpp140noextramupos,ntk_mmZelpp150noextramupos
        #plots=${plots},ntk_mmZelpp120noextramuneg,ntk_mmZelpp130noextramuneg,ntk_mmZelpp140noextramuneg,ntk_mmZelpp150noextramuneg
        #plots=${plots},mpp_mmZelpphighPur120noextramu,mpp_mmZelpphighPur130noextramu,mpp_mmZelpphighPur140noextramu,mpp_mmZelpphighPur150noextramu
        #plots=${plots},diff
        python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/plotter.py ${commonOpts} --only ${plots} --strictOnly --signalJson ${signal_json} --saveLog;
        
        #python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/compareBackgroundEstimation.py plots/analysis/plots/plotter.root plots/analysis/plots;
        ;;

    WWWANA )
        mkdir -p ${wwwdir}/ana
        cp plots/analysis/plots/*.{png,pdf,dat} ${wwwdir}/ana
        cp $CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/index.php ${wwwdir}/ana
        ;;

    OPTIMANA )
        for i in 1000 1200 1400; do
            python test/analysis/pps/optimizeSR.py ${i}; 
        done
        python test/analysis/pps/plotLimits.py 1000=optimresults_1000.pck,1200=optimresults_1200.pck,1400=optimresults_1400.pck
        python test//analysis/pps/compareOptimResults.py
        ;;

esac
