#!/bin/bash

WHAT=$1; 
if [ "$#" -lt 1 ]; then 
    echo "steerTOPpPb.sh <SEL>";
    echo "        SEL          - selects data and MC";
    exit 1; 
fi

#suppress batch notifications by mail
export LSB_JOB_REPORT_MAIL=N

queue=8nh
outdir=~/work/LJets-pPb
wwwdir=~/www/LJets-pPb
lumi=27.4


RED='\e[31m'
NC='\e[0m'
case $WHAT in
    MCSEL )
	#queue=local
	for p in PbP pPb; do
	    for ch in 11 13; do
		python scripts/runLocalAnalysis.py \
		    -i /store/cmst3/group/hintt/mverweij/forest_PYQUEN_TTBAR/${p}/v1/merge/ \
		    -o ${outdir}/analysis_${p}_${ch} \
		    --ch ${ch} \
		    --era era5TeV \
		    -m TOP-HIForest::RunToppPb \
		    -q ${queue};
	    done
	done
	;;

    SEL )
	echo -e "[ ${RED} Sending out jobs to batch ${NC} ]"

	commonOpts="--era era5TeV -m TOP-HIForest::RunToppPb"	
	inDir=/store/cmst3/group/top/mverweij/PA8TeV/data/PASingleMuon/crab_FilteredPASingleMuHighPt_PPb_v4/161219_092237/
	a=(`eos ls ${inDir}`)
	for i in ${a[@]}; do
	    python scripts/runLocalAnalysis.py -i ${inDir}/${i}/ -q ${queue} -o ${outdir}/analysis_pPb         ${commonOpts} --ch 13;
	    python scripts/runLocalAnalysis.py -i ${inDir}/${i}/ -q ${queue} -o ${outdir}/analysis_pPb_noniso  ${commonOpts} --ch 1300;
	done

	inDir=/store/cmst3/group/top/mverweij/PA8TeV/data/PASingleMuon/crab_FilteredPASingleMuHighPt_PbP_v6/170105_162741/
	a=(`eos ls ${inDir}`)
        for i in ${a[@]}; do
	    python scripts/runLocalAnalysis.py -i ${inDir}/${i} -q ${queue} -o ${outdir}/analysis_Pbp          ${commonOpts} --ch 13;
	    python scripts/runLocalAnalysis.py -i ${inDir}/${i}/ -q ${queue} -o ${outdir}/analysis_Pbp_noniso  ${commonOpts} --ch 1300;
	done

	;;
    MERGE )
	exit -1
	echo -e "[ ${RED} Merging job output ${NC} ]"
	for i in pPb; do #Pbp
	    ./scripts/mergeOutputs.py ${outdir}/analysis_${i};
	done
	;;
    PLOT )
	exit -1
	echo -e "[ ${RED} Running plotter ${NC} ]"
	a=(pPb Pbp)
	for i in ${a[@]}; do
	    continue
	    python scripts/plotter.py -i ${outdir}/analysis_${i}  -j data/era5TeV/Wsamples.json     --com "5.02 TeV" -l ${lumi} --saveLog --noStack;	
	    mkdir ~/${outdir}/analysis_${i}/wplots;
	    mv ~/${outdir}/analysis_${i}/plots/* ~/${outdir}/analysis_${i}/wplots/;
	    python scripts/plotter.py -i ${outdir}/analysis_${i}  -j data/era5TeV/samples.json      --com "5.02 TeV" -l ${lumi} --saveLog;	
	    python scripts/plotter.py -i ${outdir}/analysis_${i}  -j data/era5TeV/syst_samples.json --com "5.02 TeV" -l ${lumi} -o syst_plotter.root --silent;	
	done
	python test/TopLJAnalysis/combinePlotsForAllCategories.py data/era5TeV/samples.json metpt_0b,metpt_1b,metpt_2b,rankedmjj_0b,rankedmjj_1b,rankedmjj_2b,drjj_0b,drjj_1b,drjj_2b;
	;;

    WWW )
	exit -1
	echo -e "[ ${RED} Moving plots to ${outdir} ${NC} ]"
	a=(mu munoniso e enoniso)
	for i in ${a[@]}; do
	    continue
	    mkdir -p ${wwwdir}/analysis_${i}
	    cp ${outdir}/analysis_${i}/plots/*.{png,pdf} ${wwwdir}/analysis_${i}
	    cp test/index.php ${wwwdir}/analysis_${i}
	done
	mkdir -p ${wwwdir}/combined
	mv plots/* ${wwwdir}/combined
	cp test/index.php ${wwwdir}/combined
	;;

esac