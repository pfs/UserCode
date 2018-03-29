#!/bin/bash

WHAT=$1; 
if [ "$#" -ne 1 ]; then 
    echo "steerVBFVectorBoson.sh <SEL/MERGE/PLOT/WWW>";
    echo "        SEL          - launches selection jobs to the batch, output will contain summary trees and control plots"; 
    echo "        MERGE        - merge output"
    echo "        PLOT         - make plots"
    echo "        WWW          - move plots to web-based are"
    exit 1; 
fi

#to run locally use local as queue + can add "--njobs 8" to use 8 parallel jobs
queue=workday
githash=c29f431 #fbc74ae
eosdir=/store/cmst3/group/top/RunIIFall17/${githash}
fulllumi=41367
vbflumi=7661
lumiUnc=0.025
outdir=${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/VBFVectorBoson
wwwdir=~/www/VBFVectorBoson


RED='\e[31m'
NC='\e[0m'
case $WHAT in

    TESTSEL )
        input=${eosdir}/MC13TeV_DY50toInf/MergedMiniEvents_0_ext0.root
        output=MC13TeV_DY.root
        tag="--tag MC13TeV_DY50toInf"
        #input=${eosdir}/Data13TeV_SingleMuon_2017D/MergedMiniEvents_0_ext0.root
        #output=Data13TeV_SingleMuon_2017D.root
        #input=${eosdir}/MC13TeV_GJets_HT100to200_DR04/MergedMiniEvents_0_ext0.root
        #output=MC13TeV_GJets_HT100to200_DR04.root \
        #input=${eosdir}/Data13TeV_SingleMuon_2017C/MergedMiniEvents_0_ext0.root
        #output=Data13TeV_SingleMuon_2017C.root
	python scripts/runLocalAnalysis.py \
            -i ${input} -o ${output} ${tag} \
            --njobs 1 -q local --debug \
            --era era2017 -m VBFVectorBoson::RunVBFVectorBoson --ch 0 --runSysts;
        ;;
    SEL )
        #--only data/era2017/vbf_samples.json --exactonly \
	python scripts/runLocalAnalysis.py -i ${eosdir} \
            -o ${outdir} \
            -q ${queue} \
            --era era2017 -m VBFVectorBoson::RunVBFVectorBoson --ch 0 --runSysts;
	;;

    MERGE )
	./scripts/mergeOutputs.py ${outdir};
	;;
    PLOT )
	commonOpts="-i ${outdir} --puNormSF puwgtctr -l ${fulllumi}  --saveLog --mcUnc ${lumiUnc} --lumiSpecs VBFA:${vbflumi}"
	python scripts/plotter.py ${commonOpts} -j data/era2017/vbf_samples.json; 
        python scripts/plotter.py ${commonOpts} -j data/era2017/vbf_samples_dr04.json -O ${outdir}/plots_dr04;
        python scripts/plotter.py ${commonOpts} -j data/era2017/gjets_samples.json --noStack -O ${outdir}/plots_gjets;
	;;
    WWW )
        for p in plots plots_dr04 plots_gjets; do
	    mkdir -p ${wwwdir}/${p}
	    cp ${outdir}/${p}/*.{png,pdf} ${wwwdir}/${p};
	    cp test/index.php ${wwwdir}/${p};
            echo "Check plots in ${wwwdir}/${p}"
        done
	;;
esac
