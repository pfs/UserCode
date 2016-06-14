#!/bin/bash

WHAT=$1; 
if [[ "$1" == "" ]]; then 
    echo "steerTOPWidthAnalysis.sh <SEL/MERGESEL/PLOTSEL/WWWSEL/ANA/MERGE/PLOT/WWW>";
    echo "        SEL        - launches selection jobs to the batch, output will contain summary trees and control plots"; 
    echo "        MERGESEL   - merge the output of the jobs";
    echo "        PLOTSEL    - runs the plotter tool on the selection";
    echo "        WWWSEL     - moves the plots to an afs-web-based area";
    echo "        ANA        - analyze the selected events";
    echo "        MERGE      - merge the output of the analysis jobs";
    echo "        PLOT       - runs the plotter tool on the analysis outputs";
    echo "        WORKSPACE       - create a workspace with the TopHypoTest model";
    echo "        SCAN            - scans the likelihood";
    echo "        CLs             - calculates CLs using higgs_combine";
    echo "        TOYS            - plots toys and saves locally";
    echo "        MERGE_DATACARDS - merge all datacards for one LFS/WID into one WID datacard";
    echo "        QUANTILES       - plots the quantiles chart for all distributions, by lepton final state";
    echo "        <above>_M       - with <above> one of the commands above (past 'PLOT'), run over merged datacards"; 
    echo "        WWW             - moves the analysis plots to an afs-web-based area";

    exit 1; 
fi

queue=8nh
eosdir=/store/cmst3/user/psilva/LJets2015/8c1e7c9
eosdatadir=/store/cmst3/user/psilva/LJets2015/13c9e4c
outdir=~/CMSSW_7_6_3/src/TopLJets2015/TopAnalysis/output
cardsdir=~/CMSSW_7_6_3/src/TopLJets2015/TopAnalysis/datacards
wwwdir=~/www/ana2
lumi=2267.84

lfs=(E EE EM MM M)
wid=(0p5w 2p0w 4p0w)

RED='\e[31m'
NC='\e[0m'

CMSSW_7_4_7dir=~/CMSSW_7_4_7/src/
CMSSW_7_6_3dir=~/CMSSW_7_6_3/src/

case $WHAT in
    SEL )
        cd ${CMSSW_7_6_3dir}
	    #python scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue} -o ${outdir} -m TOPWidth::RunTopWidth --ch 0 --only MC;
	    python scripts/runLocalAnalysis.py -i ${eosdatadir} -q ${queue} -o ${outdir} -m TOPWidth::RunTopWidth --ch 0 --only Data;
	;;
    MERGESEL )
        cd ${CMSSW_7_6_3dir}
	    ./scripts/mergeOutputs.py ${outdir} True;	
	;;
    PLOTSEL )
        cd ${CMSSW_7_6_3dir}
    	python scripts/plotter.py -i ${outdir} --puNormSF puwgtctr  -j data/samples_Run2015.json -l ${lumi};	
	;;
    WWWSEL )
	    mkdir -p ${wwwdir}/sel
    	cp ${outdir}/plots/*.{png,pdf} ${wwwdir}/sel
    	cp test/index.php ${wwwdir}/sel
	;;
    ANA )
        cd ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/
	    python scripts/runTopWidthAnalysis.py -i ${outdir}/Chunks -o ${outdir}/analysis -q 8nh;
	;;
    MERGE )
        cd ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/
	    ./scripts/mergeOutputs.py ${outdir}/analysis;
	;;
    PLOT )
        cd ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/
        python scripts/plotter.py -i ${outdir}/analysis  -j data/samples_Run2015.json -l ${lumi};        
        python scripts/plotter.py -i ${outdir}/analysis  -j data/syst_samples_Run2015.json -l ${lumi} -o syst_plotter.root;        
        ;;
    MERGE_DATACARDS )
        cd ${CMSSW_7_4_7dir}
        for twid in ${wid[*]}
        do
            mcmd="python ${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/scripts/combineCards.py "
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
        cd ${CMSSW_7_4_7dir}
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
	             combine ${outdir}/${twid}.root -M MultiDimFit \
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
                      &> ${outdir}/x_pre-fit_exp__${twid}.log
            done
    ;;
    TOYS_M )
        cd ${CMSSW_7_4_7dir}
        cd ${outdir}
            for twid in ${wid[*]}
            do
                 cd ${outdir}
                 root -l -q -b higgsCombinex_pre-fit_exp_${twid}.HybridNew.mH172.5.8192.quant0.500.root \
                   "${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/test/plotting/hypoTestResultTreeTopWid.cxx(\"x_pre-fit_exp__${twid}.qvals.root\",172.5,1,\"x\",1000,\"\",\"${twid}\")"
            done
    ;;
    WORKSPACE )
        cd ${CMSSW_7_4_7dir}
        mkdir ${outdir}
        for tlfs in ${lfs[*]}
        do
            for twid in ${wid[*]}
            do
                echo "Creating workspace for ${twid}${tlfs}" 
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
	             combine ${outdir}/${twid}_${tlfs}.root -M MultiDimFit \
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
                 combine ${outdir}/${twid}_${tlfs}.root -M HybridNew --seed 8192 --saveHybridResult \
                      -m 172.5  --testStat=TEV --singlePoint 1 -T 500 -i 2 --fork 6 \
                      --clsAcc 0 --fullBToys  --generateExt=1 --generateNuis=0 \
                      --expectedFromGrid 0.5 -n x_pre-fit_exp_${twid}_${tlfs} \
                      &> ${outdir}/x_pre-fit_exp_${twid}_${tlfs}.log
            done
        done
    ;;
    TOYS )
        cd ${CMSSW_7_4_7dir}
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
    QUANTILES )
        cd ${outdir}
        lfsStr=""
        for tlfs in ${lfs[*]}
        do
          if [[ "${lfsStr}" == "" ]];
          then
            lfsStr="${tlfs}"
          else
            lfsStr="${lfsStr},${tlfs}"
          fi
        done
        python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidAnalysis/getQuantilesPlot.py \
               -i ${outdir}/ -o ${outdir}/ --lfs ${lfsStr} 
    ;;
    QUANTILES_M )
        cd ${outdir}
        python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidAnalysis/getQuantilesPlot.py \ 
               -i ${outdir}/ -o ${outdir}/
    ;;
    WWW )
        mkdir -p ${wwwdir}/ana
        cp ${outdir}/analysis/plots/*.{png,pdf} ${wwwdir}/ana        
        cp ${outdir}/*.{png,pdf} ${wwwdir}/ana        
        cp test/index.php ${wwwdir}/ana
	;;
esac
