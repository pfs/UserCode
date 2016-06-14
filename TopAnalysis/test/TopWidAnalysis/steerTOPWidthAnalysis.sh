#!/bin/bash

WHAT=$1; 
if [[ "$1" == "" ]]; then 
    echo "steerTOPSensitivityAnalysis.sh <WORKSPACE/SCAN/CLs/TOYS>";
    echo "        WORKSPACE       - create a workspace with the TopHypoTest model";
    echo "        SCAN            - scans the likelihood";
    echo "        CLs             - calculates CLs using higgs_combine";
    echo "        TOYS            - plots toys and saves locally";
    echo "        MERGE_DATACARDS - merge all datacards for one LFS/WID into one WID datacard";
    echo "        <above>_M       - with <above> one of the commands above, run over merged datacards"; 
    echo "        WWW             - moves the analysis plots to an afs-web-based area";

    exit 1; 
fi

queue=8nh
cardsdir=~/CMSSW_7_4_7/src/datacards
outdir=~/CMSSW_7_4_7/src/twid_outputs
wwwdir=~/www/ana2
lumi=2267.84

lfs=(E EE EM MM M)
wid=(0p5w 2p0w 4p0w)

RED='\e[31m'
NC='\e[0m'
            
case $WHAT in
    MERGE_DATACARDS )
        for twid in ${wid[*]}
        do
            mcmd='combineCards.py '
            for tlfs in ${lfs[*]}
            do
                mcmd="${mcmd} ${tlfs}=${cardsdir}/datacard__${twid}_${tlfs}.dat "
            done

            echo $mcmd
            echo " "
            mcmd="${mcmd} > ${cardsdir}/datacard__${twid}.dat"
            
            eval $mcmd
        done
    ;;
    WORKSPACE_M )
        mkdir ${outdir}
            for twid in ${wid[*]}
            do
                echo "Creating workspace for ${twid}" 
	            text2workspace.py ${cardsdir}/datacard__${twid}.dat -P \
                     HiggsAnalysis.CombinedLimit.TopHypoTest:twoHypothesisTest \
                     -m 172.5 --PO verbose --PO altSignal=${twid} --PO muFloating \
                     -o ${outdir}/${twid}.root
            done
    ;;
    SCAN_M )
        cd ${outdir}
            for twid in ${wid[*]}
            do
	             combine ${twid}.root -M MultiDimFit \
                     -m 172.5 -P x --floatOtherPOI=1 --algo=grid --points=200 \
                     -t -1 --expectSignal=1 --setPhysicsModelParameters x=0,r=1 \
                     -n x0_scan_Asimov_${twid}
            done
    ;;
    CLs_M )
        cd ${outdir}
            for twid in ${wid[*]}
            do
                 combine ${twid}.root -M HybridNew --seed 8192 --saveHybridResult \
                      -m 172.5  --testStat=TEV --singlePoint 1 -T 500 -i 2 --fork 6 \
                      --clsAcc 0 --fullBToys  --generateExt=1 --generateNuis=0 \
                      --expectedFromGrid 0.5 -n x_pre-fit_exp_${twid} \
                      &> x_pre-fit_exp__${twid}.log
            done
    ;;
    TOYS_M )
        cd ${outdir}
            for twid in ${wid[*]}
            do
                 cd ${outdir}
                 root -l -q -b higgsCombinex_pre-fit_exp_${twid}.HybridNew.mH172.5.8192.quant0.500.root \
                   "${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/test/plotting/hypoTestResultTreeTopWid.cxx(\"x_pre-fit_exp__${twid}.qvals.root\",172.5,1,\"x\",1000,\"\",\"${twid}\")"
            done
    ;;
    WORKSPACE )
        mkdir ${outdir}
        for tlfs in ${lfs[*]}
        do
            for twid in ${wid[*]}
            do
                echo "Creating workspace for ${twid}" 
	            text2workspace.py ${cardsdir}/datacard__${twid}_${tlfs}.dat -P \
                     HiggsAnalysis.CombinedLimit.TopHypoTest:twoHypothesisTest \
                     -m 172.5 --PO verbose --PO altSignal=${twid} --PO muFloating \
                     -o ${outdir}/${twid}_${tlfs}.root
            done
        done
    ;;
    SCAN )
        cd ${outdir}
        for tlfs in ${lfs[*]}
        do
            for twid in ${wid[*]}
            do
	             combine ${twid}_${tlfs}.root -M MultiDimFit \
                     -m 172.5 -P x --floatOtherPOI=1 --algo=grid --points=200 \
                     -t -1 --expectSignal=1 --setPhysicsModelParameters x=0,r=1 \
                     -n x0_scan_Asimov_${twid}_${tlfs}
            done
        done
    ;;
    CLs )
        cd ${outdir}
        for tlfs in ${lfs[*]}
        do
            for twid in ${wid[*]}
            do
                 combine ${twid}.root -M HybridNew --seed 8192 --saveHybridResult \
                      -m 172.5  --testStat=TEV --singlePoint 1 -T 500 -i 2 --fork 6 \
                      --clsAcc 0 --fullBToys  --generateExt=1 --generateNuis=0 \
                      --expectedFromGrid 0.5 -n x_pre-fit_exp_${twid}_${tlfs} \
                      &> x_pre-fit_exp_${twid}_${tlfs}.log
            done
        done
    ;;
    TOYS )
        cd ${outdir}
        for tlfs in ${lfs[*]}
        do
            for twid in ${wid[*]}
            do
                 cd ${outdir}
                 root -l -q -b higgsCombinex_pre-fit_exp_${twid}_${tlfs}.HybridNew.mH172.5.8192.quant0.500.root \
                   "${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/test/plotting/hypoTestResultTreeTopWid.cxx(\"x_pre-fit_exp_${twid}_${tlfs}.qvals.root\",172.5,1,\"x\",1000,\"${tlfs}\",\"${twid}\")"
            done
        done
    ;;
    WWW )
        mkdir ${www}/
        cp ${outdir}/*.{png,pdf} ${wwwdir}/
	;;
esac
