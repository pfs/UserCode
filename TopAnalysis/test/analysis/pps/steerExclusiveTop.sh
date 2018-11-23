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
githash=3129835
eosdir=/store/cmst3/group/top/RunIIReReco/${githash}
outdir=/store/cmst3/user/psilva/ExclusiveAna
samples_json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/samples.json
jetht_samples_json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/jetht_samples.json
zx_samples_json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/zx_samples.json
RPout_json=$CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/golden_noRP.json
wwwdir=~/www/ExclusiveAna
inputfileTag=MC13TeV_2017_GGH2000toZZ2L2Nu
#inputfileTag=MC13TeV_2017_GJets_HT400to600
inputfileTESTSEL=${eosdir}/${inputfileTag}/Chunk_0_ext0.root
#inputfileTag=Data13TeV_2017B_DoubleMuon
#inputfileTESTSEL=/store/cmst3/group/top/RunIIReReco/f93b8d8/${inputfileTag}/Chunk_0_ext0.root
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
	python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/runLocalAnalysis.py ${baseOpt} --only ${zx_samples_json};
        #python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/runLocalAnalysis.py ${baseOpt} --only ${jetht_samples_json};
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
        mv *.png plots/zx_sel/
	;;

    WWWSEL )
	mkdir -p ${wwwdir}/presel
	cp plots/sel/*.{png,pdf} ${wwwdir}/presel
	cp plots/zx_sel/*.{png,pdf} ${wwwdir}/presel
	cp $CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/index.php ${wwwdir}/presel
	;;

    PREPAREANA )
        python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/runExclusiveAnalysis.py --step 0 \
            --json ${samples_json} --RPout ${RPout_json} -o plots/analysis;        
        ;;

    ANA )
        python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/runExclusiveAnalysis.py --step 1 \
            --json ${samples_json} --RPout ${RPout_json} -o plots/analysis --mix plots/analysis/evmix.pck; 
	mergeOutputs.py plots/analysis;
        ;;

    PLOTANA )
        commonOpts="-i plots/analysis -j ${samples_json} -l ${ppsLumi} --mcUnc ${lumiUnc} "
        python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/scripts/plotter.py ${commonOpts};
        python $CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/analysis/pps/compareBackgroundEstimation.py plots/analysis/plots/plotter.root plots/analysis/plots;
        ;;

    WWWANA )
        mkdir -p ${wwwdir}/ana
        cp plots/analysis/plots/*.{png,pdf} ${wwwdir}/ana
        cp $CMSSW_BASE/src/TopLJets2015/TopAnalysis/test/index.php ${wwwdir}/ana
        ;;

esac
