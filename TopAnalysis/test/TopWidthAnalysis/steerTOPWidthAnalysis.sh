#!/bin/bash

WHAT=$1; 
ERA=$2;
if [ "$#" -ne 2 ]; then 
    echo "steerTOPWidthAnalysis.sh <SEL/MERGESEL/PLOTSEL/WWWSEL/ANA/MERGE/PLOT/WWW>";
    echo "        WORKSPACE       - create a workspace with the TopHypoTest model";
    echo "        SCAN            - scans the likelihood";
    echo "        CLs             - calculates CLs using higgs_combine";
    echo "        TOYS            - plots toys and saves locally";
    echo "        MERGE_DATACARDS - merge all datacards for one LFS/WID into one WID datacard";
    echo "        QUANTILES       - plots the quantiles chart for all distributions, by lepton final state";
    echo "        WWW             - moves the analysis plots to an afs-web-based area";
    echo " "
    echo "        ERA          - era2015/era2016";
    exit 1; 
fi

export LSB_JOB_REPORT_MAIL=N

queue=2nd
outdir=/afs/cern.ch/work/e/ecoleman/public/TopWidth/TopWidth_${ERA}/
cardsdir=${outdir}/datacards
wwwdir=~/www/TopWidth_${ERA}/
CMSSW_7_1_5dir=~/CMSSW_7_1_5/src/
CMSSW_7_4_7dir=~/CMSSW_7_4_7/src/
CMSSW_7_6_3dir=~/CMSSW_8_0_11/src/

unblind=true
nPseudo=1000

lfs=(EE EM MM)
wid=(0p2w 0p4w 0p6w 0p8w 1p0w 1p2w 1p4w 1p6w 1p8w 2p0w 2p2w 2p4w 2p6w 2p8w 3p0w 3p5w 4p0w)
#wid=(0p5w 1p0w 1p5w 2p0w 2p5w 3p0w 3p5w 4p0w)
#dists=(mdrmlb minmlb incmlb sncmlb mt2mlb)
dists=(incmlb)
cat=(1b 2b)
lbCat=(highpt lowpt)

nWids=("0.2" "0.4" "0.6" "0.8" "1.0" "1.2" "1.4" "1.6" "1.8" "2.0" "2.2" "2.4" "2.6" "2.8" "3.0" "3.5" "4.0")
#nMass=("169.5" "172.5" "175.5")
#nWids=("2.2")
nMass=("172.5")

nuisances=(jes jesrate jer pu btag les ltag trigEE trigEM trigMM selEE selEM selMM toppt MEqcdscale PDF Herwig amcnloFxFx Mtop ttPartonShower tWttinterf tWQCDScale DYnorm_thEE DYnorm_thEM DYnorm_thMM) 

RED='\e[31m'
NC='\e[0m'

# Helpers: getting comma-separated lists of variables
distStr=""
lfsStr=""
widStr=""
for dist in ${dists[*]}
do
    if [[ "${distStr}" == "" ]];
    then
        distStr="${dist}"
    else
        distStr="${distStr},${dist}"
    fi
done

for tlfs in ${lfs[*]} ; do
    if [[ "${lfsStr}" == "" ]];
    then
        lfsStr="${tlfs}"
    else
        lfsStr="${lfsStr},${tlfs}"
    fi
done

for twid in ${wid[*]} ; do
    if [[ "${widStr}" == "" ]];
    then
        widStr="${twid}"
    else
        widStr="${widStr},${twid}"
    fi
done
       


case $WHAT in
    SHAPES ) # to get the shapes file / datacards for the full analysis
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh`
        cd ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/
        export PYTHONPATH=$PYTHONPATH:/usr/lib64/python2.6/site-packages/
    
        rm ${cardsdir}/*.root

        for index in `seq 0 17 204` ; do
            min=$index
            max=$(($index+16))

            nohup python test/TopWidthAnalysis/createShapesFromPlotter.py \
                    -s tbart,tW \
                    --dists ${distStr} \
                    -o ${outdir}/datacards/ \
                    -n shapes$index \
                    -i ${outdir}/analysis/plots/plotter.root \
                    --systInput ${outdir}/analysis/plots/syst_plotter.root \
                    --nomorph --min $min --max $max > shapes${index}.txt & 
        done
    ;;
    X4VALIDATION ) # to get the shapes file / datacards for the full analysis
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh`


        export PYTHONPATH=$PYTHONPATH:/usr/lib64/python2.6/site-packages/
        for twid in 4p0w ; do 
            cd ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/
            mkdir ${outdir}/widthx4validation_new
            rm ${outdir}/widthx4validation_new/datacards/*
        for index in `seq 0 17 204` ; do
            min=$index
            max=$(($index+16))

            nohup python test/TopWidthAnalysis/createShapesFromPlotter.py \
                    -s tbart,tW \
                    --dists ${distStr} \
                    -o ${outdir}/widthx4validation_new/datacards/ \
                    -i ${outdir}/analysis/plots/plotter.root \
                    --systInput ${outdir}/analysis/plots/syst_plotter.root \
                    --truth 4p0w \
                    --trx4 widthx4 \
                    --trx4Width 1p0w \
                    --trx4Input ${outdir}/prevAnalysis4/analysis/plots/syst_plotter.root \
                    -n shapes${twid}$index \
                    --min $min --max $max \
                    --nomorph > shapes${twid}${index}.txt &
        done
            # so local jobs don't destroy lxplus (cheap but it works)
            #sleep 900 
        done
    ;;
    X4VALIDATION_SCANS )
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh`

        for swid in 4p0w; do
        
        base=${outdir}/widthx4validation_new/
        # hadd shapes 
        #hadd ${base}/datacards/shapes.root ${base}/datacards/shapes*.root 

        # produce datacards
        # for a given width, merge all
        for dist in ${dists[*]} ; do
        for twid in ${wid[*]} ; do
            allcmd="python ${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/scripts/combineCards.py "
        for tlbCat in ${lbCat[*]} ; do
        for tlfs in ${lfs[*]} ; do
        for tCat in ${cat[*]} ; do
            cardname="${tlbCat}${tlfs}${tCat}=${base}/datacards"
            cardname="${cardname}/datacard__${twid}_${tlbCat}${tlfs}${tCat}_${dist}.dat"
            allcmd="${allcmd} ${cardname} "
        done
        done
        done
            allcmd="${allcmd} > ${base}/datacards/datacard__${twid}_${dist}.dat"
            eval $allcmd
        done
        done

        # workspaces 
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh`
        cd ${base}/
        for dist in ${dists[*]} ; do
        for twid in ${wid[*]} ; do
            
            echo "Creating workspace for ${twid}${dist}" 
            text2workspace.py ${base}/datacards/datacard__${twid}_${dist}.dat -P \
                HiggsAnalysis.CombinedLimit.TopHypoTest:twoHypothesisTest \
                -m 172.5 --PO verbose --PO altSignal=${twid} --PO muFloating \
                -o ${base}/${twid}_${dist}.root 

        done
        done
        
        
        # scans
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh`
        cd ${base}
        for dist in ${dists[*]} ; do
        for twid in ${wid[*]} ; do
            if [[ "${twid}" == "${swid}" ]] ; then

                # x = 0
                cmd="combine ${base}/${twid}_${dist}.root -M MultiDimFit" 
                cmd="${cmd} -m 172.5 --redefineSignalPOIs x --floatOtherPOI=1 --minimizerTolerance 0.0001 --algo=grid --points=50"
                cmd="${cmd} -n x0_fit_exp --setPhysicsModelParameters x=0.0,r=1.0"
                cmd="${cmd} --saveWorkspace --expectSignal=1 -t -1 --robustFit=1"
                
                bsub -q ${queue} \
                    ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh \
                    "${base}/" "${cmd}" 

                # x = 1
                cmd="combine ${base}/${twid}_${dist}.root -M MultiDimFit" 
                cmd="${cmd} -m 172.5 --redefineSignalPOIs x --floatOtherPOI=1 --minimizerTolerance 0.0001 --algo=grid --points=50"
                cmd="${cmd} -n x1_fit_exp --setPhysicsModelParameters x=1.0,r=1.0"
                cmd="${cmd} --saveWorkspace --expectSignal=1 -t -1 --robustFit=1"
                
                bsub -q ${queue} \
                    ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh \
                    "${base}/" "${cmd}" 

                # x from data 
                cmd="combine ${base}/${twid}_${dist}.root -M MultiDimFit" 
                cmd="${cmd} -m 172.5 -P x --minimizerTolerance 0.0001 --algo=grid --points=50"
                cmd="${cmd} -n x_fit${twid}_obs --robustFit=1 --setPhysicsModelParameters x=1.0,r=1.0"
                cmd="${cmd} --saveWorkspace"
                
                bsub -q ${queue} \
                    sh ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh \
                    "${base}/" "${cmd}" 

            fi

            # pre-fit 
            echo "Making CLs for ${twid} ${dist}"
            cmd="combine ${base}/${twid}_${dist}.root -M HybridNew --seed 8192 --saveHybridResult" 
            cmd="${cmd} -m 172.5  --testStat=TEV --singlePoint 1 -T ${nPseudo} -i 2 --fork 6"
            cmd="${cmd} --clsAcc 0 --fullBToys  --generateExt=1 --generateNuis=0"
            cmd="${cmd} --expectedFromGrid 0.5 -n x_pre-fit_exp__${twid}_${dist}"

            bsub -q ${queue} ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh "${base}/" "${cmd}" 
            
            # post-fit expected 
            echo "Making CLs for ${twid} ${dist}"
            cmd="combine ${base}/${twid}_${dist}.root -M HybridNew --seed 8192 --saveHybridResult" 
            cmd="${cmd} -m 172.5  --testStat=TEV --singlePoint 1 -T ${nPseudo} -i 2 --fork 6"
            cmd="${cmd} --clsAcc 0 --fullBToys  --frequentist"
            cmd="${cmd} --expectedFromGrid 0.5 -n x_post-fit_exp__${twid}_${dist}"

            bsub -q ${queue} ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh "${base}/" "${cmd}" 

            # post-fit observed 
            echo "Making CLs for ${twid} ${dist}"
            cmd="combine ${base}/${twid}_${dist}.root -M HybridNew --seed 8192 --saveHybridResult" 
            cmd="${cmd} -m 172.5  --testStat=TEV --singlePoint 1 -T ${nPseudo} -i 2 --fork 6"
            cmd="${cmd} --clsAcc 0 --fullBToys  --frequentist"
            cmd="${cmd} -n x_post-fit_obs__${twid}_${dist}"

            bsub -q ${queue} ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh "${base}/" "${cmd}" 
        done
        done
        done
    ;;
    X4VALIDATION_QUANTILES )
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh` 

        for swid in 4p0w ; do 

        base=${outdir}/widthx4validation_new/

        cd ${base}
        rm statsPlots.root
        for dist in ${dists[*]} ; do
            #for twid in ${wid[*]} ; do
            #    ## pre-fit expected 
            #    rootcmds="${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis"
            #    rootcmds="${rootcmds}/hypoTestResultTreeTopWid.cxx\(\\\"x_pre-fit_exp__${twid}_${dist}.qvals.root"
            #    rootcmds="${rootcmds}\\\",172.5,1,\\\"x\\\",1000,\\\"\\\",\\\"${twid}\\\",\\\"${dist}\\\",${unblind},\\\"pre\\\"\)"

            #    cmd="root -l -q -b"
            #    cmd="${cmd} ${base}"
            #    cmd="${cmd}/higgsCombinex_pre-fit_exp__${twid}_${dist}.HybridNew.mH172.5.8192.quant0.500.root"
            #    cmd="${cmd} ${rootcmds}"

            #    sh ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh \
            #        "${base}" "${cmd}"

            #    ## post-fit expected 
            #    rootcmds="${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis"
            #    rootcmds="${rootcmds}/hypoTestResultTreeTopWid.cxx\(\\\"x_post-fit_exp__${twid}_${dist}.qvals.root"
            #    rootcmds="${rootcmds}\\\",172.5,1,\\\"x\\\",1000,\\\"\\\",\\\"${twid}\\\",\\\"${dist}\\\",${unblind},\\\"post\\\"\)"

            #    cmd="root -l -q -b"
            #    cmd="${cmd} ${base}"
            #    cmd="${cmd}/higgsCombinex_post-fit_exp__${twid}_${dist}.HybridNew.mH172.5.8192.quant0.500.root"
            #    cmd="${cmd} ${rootcmds}"

            #    sh ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh \
            #        "${base}" "${cmd}"

            #    ## post-fit observed 
            #    rootcmds="${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis"
            #    rootcmds="${rootcmds}/hypoTestResultTreeTopWid.cxx\(\\\"x_post-fit_obs__${twid}_${dist}.qvals.root"
            #    rootcmds="${rootcmds}\\\",172.5,1,\\\"x\\\",1000,\\\"\\\",\\\"${twid}\\\",\\\"${dist}\\\",${unblind},\\\"obs\\\"\)"

            #    cmd="root -l -q -b"
            #    cmd="${cmd} ${base}"
            #    cmd="${cmd}/higgsCombinex_post-fit_obs__${twid}_${dist}.HybridNew.mH172.5.8192.root"
            #    cmd="${cmd} ${rootcmds}"

            #    sh ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh \
            #        "${base}" "${cmd}"
            #done
            
            # Quantiles plot with post-fit information 
            python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getQuantilesPlot.py \
                -i ${base}/ -o ${base}/ \
                --wid ${widStr} \
                --dist ${dist}  \
                --prep post \
                --unblind

            # Get CLs plots for pre-fit expectations
            #python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getSeparationTables.py\
            #    -i ${base}/ -o ${base}/ \
            #    --wid ${widStr} \
            #    --prep pre \
            #    --dist ${dist}
            
            # Get CLs plots for post-fit expectations 
            python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getSeparationTables.py\
                -i ${base}/ -o ${base}/ \
                --wid ${widStr} \
                --dist ${dist} \
                --prep post \
                --addPre \
                --splineMin 0.4 --splineMax 7.0 \
                --unblind
            
            # Get CLs plots for post-fit expectations 
            #python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getSeparationTables.py\
            #    -i ${base}/ -o ${base}/ \
            #    --wid ${widStr} \
            #    --dist ${dist} \
            #    --prep obs \
            #    --unblind
        done
        done
    ;;
    X4VALIDATION_FITS )
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh` 

        for swid in ${wid[*]} ; do 

        base=${outdir}/width${swid}validation/

        cd ${base}
        rm statsPlots.root
        for dist in ${dists[*]} ; do

            echo "Limits for ${dist} ||| ${swid}"
            python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getCLsFromFit.py \
                -i ${base}/ \
                --dist ${dist} \
                --prep post \
                --unblind
        done
        done

    ;;
    MORPH ) # to get the morphs for the full analysis
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh`
        cd ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/
        export PYTHONPATH=$PYTHONPATH:/usr/lib64/python2.6/site-packages/
        python test/TopWidthAnalysis/createShapesFromPlotter.py \
                -s tbart,tW \
                --dists ${distStr} \
                -o ${outdir}/datacards/ \
                -i ${outdir}/analysis/plots/plotter.root \
                --systInput ${outdir}/analysis/plots/syst_plotter.root \
                --noshapes 
    ;;
    MORPHALL ) # get shapes/morphs replacing data with a given width template 
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh`
        cd ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/
        export PYTHONPATH=$PYTHONPATH:/usr/lib64/python2.6/site-packages/

        # for each width, generate a truth dataset and produce the proper morph datacards for it
        # --> output in datacards_<wid>/
        for twid in ${wid[*]}; do
           rm -rf ${outdir}/datacards_${twid}
           nohup python test/TopWidthAnalysis/createShapesFromPlotter.py \
                   -s tbart,tW \
                   --dists ${distStr} \
                   -o ${outdir}/datacards_${twid}/ \
                   -i ${outdir}/analysis/plots/plotter.root \
                   --systInput ${outdir}/analysis/plots/syst_plotter.root \
                   --truth ${twid} &
        done
    ;;
    2DBIAS ) # get the bias plot for the 2D scan
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh`

        # Merge datacards
        for twid in ${wid[*]}; do
        for dist in ${dist[*]}; do
            cmd="python ${CMSSW_7_4_7dir}/HiggsAnalysis/CombinedLimit/scripts/combineCards.py "

            for tlbCat in ${lbCat[*]} ; do
            for tlfs in ${lfs[*]} ; do
            for tCat in ${cat[*]} ; do
                cardname="${tlbCat}${tlfs}${tCat}=${outdir}/datacards_${twid}"
                cardname="${cardname}/datacard_widfit__${tlbCat}${tlfs}${tCat}_${dist}.dat"
                cmd="${cmd} ${cardname} "
            done
            done
            done
        done
        done
            cmd="${cmd} > ${outdir}/datacards_${twid}/datacard_widfit__${twid}_${dist}.dat"
            eval `${cmd}`

        # Get workspaces 
        for twid in ${wid[*]}; do 
        for dist in ${dist[*]}; do
            cd ${outdir}/datacards_${twid}
            text2workspace.py ${outdir}/datacards_${twid}/datacard_widfit__${twid}_${dist}.dat \
                -o ${outdir}/datacards_${twid}/workspace_biascheck_${twid}_${dist}.root
        done
        done
       
        # Run combine
        for twid in ${wid[*]}; do 
        for dist in ${dist[*]}; do
            cd ${outdir}/datacards_${twid}
            combine ${outdir}/datacards_${twid}/workspace_biascheck_${twid}_${dist}.root \
                -M MultiDimFit --expectSignal=1 --poi r --poi alpha --poi beta \
                --setPhysicsModelParameterLimits alpha=0,9:beta=0,2 \
                --algo=grid --points=1000
        done
        done
    ;;
    MERGE_SHAPES ) # hadd the shapes files since they are split
        hadd ${outdir}/datacards/shapes.root ${outdir}/datacards/shapes*.root 
    ;;
    MERGE_DATACARDS ) # get all the datacards you could possibly desire for the analysis
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh`

        for dist in ${dists[*]} ; do
        for twid in ${wid[*]} ; do

        # for a given width, merge all
        allcmd="python ${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/scripts/combineCards.py "
        for tlbCat in ${lbCat[*]} ; do

            # for a given width and lbcat, merge all
            lbccmd="python ${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/scripts/combineCards.py "
            for tlfs in ${lfs[*]} ; do

                # for a given width, lfs, and lbcat, merge all
                lfscmd="python ${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/scripts/combineCards.py "
                for tCat in ${cat[*]} ; do
                    cardname="${tlbCat}${tlfs}${tCat}=${cardsdir}"
                    cardname="${cardname}/datacard__${twid}_${tlbCat}${tlfs}${tCat}_${dist}.dat"
                    allcmd="${allcmd} ${cardname} "
                    lbccmd="${lbccmd} ${cardname} "
                    lfscmd="${lfscmd} ${cardname} "
                done

                lfscmd="${lfscmd} > ${cardsdir}/datacard__${twid}_${tlbCat}${tlfs}_${dist}.dat"
            done

            lbccmd="${lbccmd} > ${cardsdir}/datacard__${twid}_${tlbCat}_${dist}.dat"
        done

        allcmd="${allcmd} > ${cardsdir}/datacard__${twid}_${dist}.dat"
            
        eval $allcmd
        eval $lbccmd
        eval $lfscmd

        done
        done
    ;;
    WORKSPACE ) # generate combine workspaces
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh`
        mkdir ${outdir}
        for dist in ${dists[*]} ; do
        for twid in ${wid[*]} ; do
            
            # All datacards
            echo "Creating workspace for ${twid}${dist}" 
            text2workspace.py ${cardsdir}/datacard__${twid}_${dist}.dat -P \
                HiggsAnalysis.CombinedLimit.TopHypoTest:twoHypothesisTest \
                -m 172.5 --PO verbose --PO altSignal=${twid} --PO muFloating \
                -o ${outdir}/${twid}_${dist}.root 
        done
        done
    ;;
    SCAN ) # perform likelihood scans on Asimov datasets 
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh`
        cd ${outdir}
        for dist in ${dists[*]} ; do
        for twid in ${wid[*]} ; do

            # All datacards
            cmd="combine ${outdir}/${twid}_${dist}.root -M MultiDimFit" 
            cmd="${cmd} -m 172.5 -P x --floatOtherPOI=1 --algo=grid --points=200"
            cmd="${cmd} --expectSignal=1 --setPhysicsModelParameters x=0,r=1"
            cmd="${cmd} -n x0_scan_Asimov_${twid}_${dist}"
            if [[ ${unblind} == false ]] ; then 
                echo "Analysis is blinded"
                cmd="${cmd} -t -1"
            fi
            
            bsub -q ${queue} \
                ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh \
                "${outdir}/" "${cmd}" 
        done
        done
    ;;
    CLs ) # get CLs statistics from combine
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh`
        cd ${outdir}
        for dist in ${dists[*]} ; do
        for twid in ${wid[*]} ; do

            # pre-fit 
            echo "Making CLs for ${twid} ${dist}"
            cmd="combine ${outdir}/${twid}_${dist}.root -M HybridNew --seed 8192 --saveHybridResult" 
            cmd="${cmd} -m 172.5  --testStat=TEV --singlePoint 1 -T ${nPseudo} -i 2 --fork 6"
            cmd="${cmd} --clsAcc 0 --fullBToys  --generateExt=1 --generateNuis=0"
            cmd="${cmd} --expectedFromGrid 0.5 -n x_pre-fit_exp__${twid}_${dist}"
            #cmd="${cmd} &> ${outdir}/x_pre-fit_exp__${twid}_${dist}.log"

            bsub -q ${queue} ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh "${outdir}/" "${cmd}" 
            
            # post-fit expected 
            echo "Making CLs for ${twid} ${dist}"
            cmd="combine ${outdir}/${twid}_${dist}.root -M HybridNew --seed 8192 --saveHybridResult" 
            cmd="${cmd} -m 172.5  --testStat=TEV --singlePoint 1 -T ${nPseudo} -i 2 --fork 6"
            cmd="${cmd} --clsAcc 0 --fullBToys  --frequentist"
            cmd="${cmd} --expectedFromGrid 0.5 -n x_post-fit_exp__${twid}_${dist}"
            #cmd="${cmd} &> ${outdir}/x_pre-fit_exp__${twid}_${dist}.log"

            bsub -q ${queue} ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh "${outdir}/" "${cmd}" 

            # post-fit observed 
            echo "Making CLs for ${twid} ${dist}"
            cmd="combine ${outdir}/${twid}_${dist}.root -M HybridNew --seed 8192 --saveHybridResult" 
            cmd="${cmd} -m 172.5  --testStat=TEV --singlePoint 1 -T ${nPseudo} -i 2 --fork 6"
            cmd="${cmd} --clsAcc 0 --fullBToys  --frequentist"
            cmd="${cmd} -n x_post-fit_obs__${twid}_${dist}"
            #cmd="${cmd} &> ${outdir}/x_pre-fit_exp__${twid}_${dist}.log"

            bsub -q ${queue} ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh "${outdir}/" "${cmd}" 
        done
        done
    ;;
    TOYS ) # get toys distributions from the pseudoexperiments
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh`
        cd ${outdir}
        for dist in ${dists[*]} ; do
        for twid in ${wid[*]} ; do

            ## pre-fit expected 
            rootcmds="${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis"
            rootcmds="${rootcmds}/hypoTestResultTreeTopWid.cxx\(\\\"x_pre-fit_exp__${twid}_${dist}.qvals.root"
            rootcmds="${rootcmds}\\\",172.5,1,\\\"x\\\",1000,\\\"\\\",\\\"${twid}\\\",\\\"${dist}\\\",${unblind},\\\"pre\\\"\)"

            cmd=""
            cmd="${cmd}root -l -q -b"
            cmd="${cmd} ${outdir}"
            cmd="${cmd}/higgsCombinex_pre-fit_exp__${twid}_${dist}.HybridNew.mH172.5.8192.quant0.500.root"
            cmd="${cmd} ${rootcmds}"

            #bsub -q ${queue} \
            sh ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh \
                "${outdir}/" "${cmd}"

            ## post-fit expected 
            rootcmds="${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis"
            rootcmds="${rootcmds}/hypoTestResultTreeTopWid.cxx\(\\\"x_post-fit_exp__${twid}_${dist}.qvals.root"
            rootcmds="${rootcmds}\\\",172.5,1,\\\"x\\\",1000,\\\"\\\",\\\"${twid}\\\",\\\"${dist}\\\",${unblind},\\\"post\\\"\)"

            cmd=""
            cmd="${cmd}root -l -q -b"
            cmd="${cmd} ${outdir}"
            cmd="${cmd}/higgsCombinex_post-fit_exp__${twid}_${dist}.HybridNew.mH172.5.8192.quant0.500.root"
            cmd="${cmd} ${rootcmds}"

            #bsub -q ${queue} \
            sh ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh \
                "${outdir}/" "${cmd}"

            ## post-fit observed 
            rootcmds="${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis"
            rootcmds="${rootcmds}/hypoTestResultTreeTopWid.cxx\(\\\"x_post-fit_obs__${twid}_${dist}.qvals.root"
            rootcmds="${rootcmds}\\\",172.5,1,\\\"x\\\",1000,\\\"\\\",\\\"${twid}\\\",\\\"${dist}\\\",${unblind},\\\"obs\\\"\)"

            cmd=""
            cmd="${cmd}root -l -q -b"
            cmd="${cmd} ${outdir}"
            cmd="${cmd}/higgsCombinex_post-fit_obs__${twid}_${dist}.HybridNew.mH172.5.8192.root"
            cmd="${cmd} ${rootcmds}"

            #bsub -q ${queue} \
            sh ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh \
                "${outdir}/" "${cmd}"
        done
        done
    ;;
    QUANTILES ) # plot quantiles distributions of all toys, get CLsPlot
            
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh` 

        cd ${outdir}
        rm statsPlots.root
        for dist in ${dists[*]} ; do

            # Quantiles plot with post-fit information 
            python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getQuantilesPlot.py \
                -i ${outdir}/ -o ${outdir}/ \
                --wid ${widStr} \
                --dist ${dist}  \
                --prep pre

            # Quantiles plot with post-fit information 
            python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getQuantilesPlot.py \
                -i ${outdir}/ -o ${outdir}/ \
                --wid ${widStr} \
                --dist ${dist}  \
                --prep post \
                --unblind

            # Quantiles plot with post-fit information 
            python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getQuantilesPlot.py \
                -i ${outdir}/ -o ${outdir}/ \
                --wid ${widStr} \
                --dist ${dist}  \
                --prep obs \
                --unblind
            
            # Get CLs plots for pre-fit expectations
            python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getSeparationTables.py \
                -i ${outdir}/ -o ${outdir}/ \
                --wid ${widStr} \
                --prep pre \
                --dist ${dist}
            
            # Get CLs plots for post-fit expectations 
            python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getSeparationTables.py \
                -i ${outdir}/ -o ${outdir}/ \
                --wid ${widStr} \
                --dist ${dist} \
                --prep post \
                --addPre \
                --unblind
            
            # Get CLs plots for post-fit expectations 
            python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getSeparationTables.py \
                -i ${outdir}/ -o ${outdir}/ \
                --wid ${widStr} \
                --dist ${dist} \
                --prep obs \
                --unblind

            #python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getCLsFromFit.py \
            #    -i ${outdir}/ \
            #    --dist ${dist} \
            #    --prep pre

            python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getCLsFromFit.py \
                -i ${outdir}/ \
                --dist ${dist} \
                --prep post \
                --unblind
        done
        
       # # Get Separation plots for all dists 
       # python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getSeparationTables.py \
       #     -i ${outdir}/ -o ${outdir}/ \
       #     --dist ${distStr} \
       #     --doAll \
       #     --unblind
    ;;
    NUISANCES ) # run once for each nuisance, doing a full analysis to understand effects of systematics 
        for dist in ${dists[*]} ; do
            systList=""
            i=0
        for syst in ${nuisances[*]} ; do
            if [[ "${syst}" != "${nuisances[0]}" ]]; then
               systList="${systList}," 
            fi
            systList="${systList}${syst}"

            echo "Frozen systematics are now: ${systList[@]}"

            mkdir ${outdir}/nuisanceTurnOn_${dist}_${syst}_${i}
            cd ${outdir}/nuisanceTurnOn_${dist}_${syst}_${i}
            
            # start writing job 
            for twid in ${wid[*]}; do

                # run likelihood scan
                cmd="combine ${outdir}/${twid}_${dist}.root -M MultiDimFit" 
                cmd="${cmd} -m 172.5 -P x --floatOtherPOI=1 --algo=grid --points=200"
                cmd="${cmd} --expectSignal=1 --setPhysicsModelParameters x=0,r=1"
                cmd="${cmd} -n x0_scan_Asimov_${twid}_${dist}"
                cmd="${cmd} --freezeNuisances ${systList}"
                if [[ ${unblind} == false ]] ; then
                  echo "Analysis is blinded"
                  cmd="${cmd} -t -1"   
                fi
                
                # launch job
                bsub -q ${queue} \
                ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh \
                    "${outdir}/nuisanceTurnOn_${dist}_${syst}_${i}/" "${cmd}" 

                # run CLs pre-fit
                cmd="combine ${outdir}/${twid}_${dist}.root -M HybridNew --seed 8192 --saveHybridResult" 
                cmd="${cmd} -m 172.5  --testStat=TEV --singlePoint 1 -T ${nPseudo} -i 2 --fork 6"
                cmd="${cmd} --clsAcc 0 --fullBToys  --generateExt=1 --generateNuis=0"
                cmd="${cmd} --expectedFromGrid 0.5 -n x_pre-fit_exp_${twid}_${dist}"
                cmd="${cmd} --freezeNuisances ${systList}"
                if [[ ${unblind} == false ]] ; then
                  echo "Analysis is blinded"
                  cmd="${cmd} -t -1"   
                fi
                
                # launch job
                bsub -q ${queue} \
                ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh \
                    "${outdir}/nuisanceTurnOn_${dist}_${syst}_${i}/" "${cmd}" 

                # run CLs post-fit exp
                cmd="combine ${outdir}/${twid}_${dist}.root -M HybridNew --seed 8192 --saveHybridResult" 
                cmd="${cmd} -m 172.5  --testStat=TEV --singlePoint 1 -T ${nPseudo} -i 2 --fork 6"
                cmd="${cmd} --clsAcc 0 --fullBToys  --frequentist"
                cmd="${cmd} --expectedFromGrid 0.5 -n x_post-fit_exp_${twid}_${dist}"
                cmd="${cmd} --freezeNuisances ${systList}"
                if [[ ${unblind} == false ]] ; then
                  echo "Analysis is blinded"
                  cmd="${cmd} -t -1"   
                fi
                
                # launch job
                bsub -q ${queue} \
                ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh \
                    "${outdir}/nuisanceTurnOn_${dist}_${syst}_${i}/" "${cmd}" 

                # run CLs post-fit obs
                cmd="combine ${outdir}/${twid}_${dist}.root -M HybridNew --seed 8192 --saveHybridResult" 
                cmd="${cmd} -m 172.5  --testStat=TEV --singlePoint 1 -T ${nPseudo} -i 2 --fork 6"
                cmd="${cmd} --clsAcc 0 --fullBToys  --frequentist"
                cmd="${cmd} -n x_post-fit_obs_${twid}_${dist}"
                cmd="${cmd} --freezeNuisances ${systList}"
                if [[ ${unblind} == false ]] ; then
                  echo "Analysis is blinded"
                  cmd="${cmd} -t -1"   
                fi
                
                # launch job
                bsub -q ${queue} \
                ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh \
                    "${outdir}/nuisanceTurnOn_${dist}_${syst}_${i}/" "${cmd}" 

            done
            
            let "i += 1"
        done    
        done
    ;;
    NUISANCES_TOYS )
        for dist in ${dists[*]} ; do
            i=0
        for syst in ${nuisances[*]} ; do
            cd ${outdir}/nuisanceTurnOn_${dist}_${syst}_${i}
            
            # start writing job 
            for twid in ${wid[*]}; do

                ## pre-fit expected 
                rootcmds="${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis"
                rootcmds="${rootcmds}/hypoTestResultTreeTopWid.cxx\(\\\"x_pre-fit_exp_${twid}_${dist}.qvals.root"
                rootcmds="${rootcmds}\\\",172.5,1,\\\"x\\\",1000,\\\"\\\",\\\"${twid}\\\",\\\"${dist}\\\",${unblind},\\\"pre\\\"\)"

                cmd=""
                cmd="${cmd}root -l -q -b"
                cmd="${cmd} ${outdir}/nuisanceTurnOn_${dist}_${syst}_${i}"
                cmd="${cmd}/higgsCombinex_pre-fit_exp_${twid}_${dist}.HybridNew.mH172.5.8192.quant0.500.root"
                cmd="${cmd} ${rootcmds}"

                #bsub -q ${queue} \
                sh ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh \
                    "${outdir}/nuisanceTurnOn_${dist}_${syst}_${i}/" "${cmd}" 
            

                ## post-fit expected 
                rootcmds="${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis"
                rootcmds="${rootcmds}/hypoTestResultTreeTopWid.cxx\(\\\"x_post-fit_exp_${twid}_${dist}.qvals.root"
                rootcmds="${rootcmds}\\\",172.5,1,\\\"x\\\",1000,\\\"\\\",\\\"${twid}\\\",\\\"${dist}\\\",${unblind},\\\"post\\\"\)"

                cmd=""
                cmd="${cmd}root -l -q -b"
                cmd="${cmd} ${outdir}/nuisanceTurnOn_${dist}_${syst}_${i}"
                cmd="${cmd}/higgsCombinex_post-fit_exp_${twid}_${dist}.HybridNew.mH172.5.8192.quant0.500.root"
                cmd="${cmd} ${rootcmds}"

                #bsub -q ${queue} \
                sh ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh \
                    "${outdir}/nuisanceTurnOn_${dist}_${syst}_${i}/" "${cmd}"
            

                ## post-fit observed 
                rootcmds="${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis"
                rootcmds="${rootcmds}/hypoTestResultTreeTopWid.cxx\(\\\"x_post-fit_obs_${twid}_${dist}.qvals.root"
                rootcmds="${rootcmds}\\\",172.5,1,\\\"x\\\",1000,\\\"\\\",\\\"${twid}\\\",\\\"${dist}\\\",${unblind},\\\"obs\\\"\)"

                cmd=""
                cmd="${cmd}root -l -q -b"
                cmd="${cmd} ${outdir}/nuisanceTurnOn_${dist}_${syst}_${i}"
                cmd="${cmd}/higgsCombinex_post-fit_obs_${twid}_${dist}.HybridNew.mH172.5.8192.root"
                cmd="${cmd} ${rootcmds}"

                #bsub -q ${queue} \
                sh ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh \
                    "${outdir}/nuisanceTurnOn_${dist}_${syst}_${i}/" "${cmd}"
            
            done

            let "i += 1" 
        done
        done
    ;;
    NUISANCE_FITS )
        cd ${CMSSW_7_4_7dir}/
        eval `scramv1 runtime -sh`

        for dist in ${dists[*]} ; do
            i=0
        for syst in ${nuisances[*]} ; do
            dirName=${outdir}/nuisanceTurnOn_${dist}_${syst}_${i}

            python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getSeparationTables.py \
                -i ${dirName}/ \
                -o ${dirName}/ \
                --dist ${dist} \
                --wid ${widStr} \
                --unblind

            let "i += 1" 

        done
        done

        for dist in ${dists[*]} ; do
            i=0
        for syst in ${nuisances[*]} ; do
            dirName=${outdir}/nuisanceTurnOn_${dist}_${syst}_${i}


            echo " "
            echo "$i ||| $dist $syst : "

            python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getCLsFromFit.py \
                -i $dirName \
                --dist ${dist} \
                --unblind
            
            let "i += 1" 

        done
        done
    ;;
    WWW )
        mkdir -p ${wwwdir}/ana_${ERA}
        cp ${outdir}/analysis/plots/*.{png,pdf} ${wwwdir}/ana_${ERA} 
        cp ${outdir}/*.{png,pdf} ${wwwdir}/ana_${ERA} 
        cp test/index.php ${wwwdir}/ana_${ERA}
    ;;
    MAKE_ATGCFILES )
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh`

        python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/createATGCToolInput.py \
            -i ${outdir}/datacards/shapes.root \
            -o ${outdir}/atgccards/ \
            --makeSplit \
            --injWid 1.2 --injMas 175.5
            #--nocards \
    ;;
    MOVE_ATGCFILES )
        cd ${CMSSW_7_1_5dir}
        eval `scramv1 runtime -sh`

        cp $outdir/atgccards/* \
            ${CMSSW_7_1_5dir}/CombinedEWKAnalysis/CommonTools/data/anomalousCoupling
        cp $outdir/atgccards/config_all.dat \
            ${CMSSW_7_1_5dir}/CombinedEWKAnalysis/CommonTools/test/
    ;;
    RUN_ATGCFILES )
        cd ${CMSSW_7_1_5dir}
        eval `scramv1 runtime -sh`

        cd ${CMSSW_7_1_5dir}/CombinedEWKAnalysis/CommonTools/test/

        python buildWorkspace_AC.py \
            --config config_all.dat 

        combineCards.py $(ls aC_[hl]*.txt) > aC_TopMasWid_all.txt

        chString=""
        for tlfs in ${lfs[*]} ; do 
        for tlbCat in ${lbCat[*]} ; do 
        for tCat in ${cat[*]} ; do 
            if [[ "$chString" != "" ]] ; then
                chString="$chString,"
            fi    
        
            chString="$chString$tlbCat${tlfs}${tCat}_incmlb"
        done
        done
        done

        text2workspace.py -m 126 aC_TopMasWid_all.txt -o TopMasWid.root \
            -P CombinedEWKAnalysis.CommonTools.ACModel:par1par2_TH2_Model \
            --PO channels=$chString \
            --PO poi=mt,gam --PO range_mt=169.5,175.5 --PO range_gam=0.2,3.0 

        combine TopMasWid.root -M MultiDimFit -P mt -P gam \
            -t -1 \
            --floatOtherPOIs=1 \
            --algo=grid --points=10000 \
            --minimizerStrategy=2

    ;;
    RUN_ATGCBIAS )
        cd ${CMSSW_7_1_5dir}
        eval `scramv1 runtime -sh`
        
            iWid=0
        for injWid in ${nWids[*]} ; do
            iMas=0
        for injMas in ${nMass[*]} ; do
            chExt="_${iWid}_${iMas}"

        cmd=""

        # make cards and shapes files
        cmd="${cmd} python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/createATGCToolInput.py"
        cmd="${cmd} -i ${outdir}/datacards/shapes.root"
        cmd="${cmd} -o ${outdir}/atgccards/"
        cmd="${cmd} --makeSplit"
        cmd="${cmd} --chNameAdd $chExt"
        cmd="${cmd} --injWid 1.2 --injMas 175.5 ;"

        # move relevant files
        cmd="${cmd} cp $(ls $outdir/atgccards/*${chExt}.{root,dat})"
        cmd="${cmd} ${CMSSW_7_1_5dir}/CombinedEWKAnalysis/CommonTools/data/anomalousCoupling ;"

        cmd="${cmd} cp $outdir/atgccards/config_all${chExt}.dat"
        cmd="${cmd} ${CMSSW_7_1_5dir}/CombinedEWKAnalysis/CommonTools/test/ ;"

        # change to proper directory 
        cmd="${cmd} cd ${CMSSW_7_1_5dir}/CombinedEWKAnalysis/CommonTools/test/ ;"

        # build ATGC workspace 
        cmd="${cmd} python buildWorkspace_AC.py --config config_all${chExt}.dat ;"

        # merge datacards 
        cmd="${cmd} combineCards.py $(ls aC_[hl]*${chExt}.txt) > aC_TopMasWid_all${chExt}.txt ;"

        chString=""
        for tlfs in ${lfs[*]} ; do 
        for tlbCat in ${lbCat[*]} ; do 
        for tCat in ${cat[*]} ; do 
            if [[ "$chString" != "" ]] ; then
                chString="${chString},"
            fi    
        
            chString="${chString}${tlbCat}${tlfs}${tCat}_incmlb${chExt}"
        done
        done
        done

        # build combine workspace
        cmd="${cmd} text2workspace.py -m 126 aC_TopMasWid_all${chExt}.txt -o TopMasWid${chExt}.root" 
        cmd="${cmd} -P CombinedEWKAnalysis.CommonTools.ACModel:par1par2_TH2_Model"
        cmd="${cmd} --PO channels=${chString}"
        cmd="${cmd} --PO poi=mt,gam --PO range_mt=169.5,175.5 --PO range_gam=0.2,3.0 ;"

        # run likelihood scan
        cmd="${cmd} combine TopMasWid${chExt}.root -M MultiDimFit -P mt -P gam"
        cmd="${cmd} --floatOtherPOIs=1"
        cmd="${cmd} --algo=grid --points=10000"
        cmd="${cmd} --minimizerStrategy=2"
        
        # submit to queue
        bsub -q ${queue} \
        sh ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh \
            "${outdir}/" "${cmd}"

            let "iMas += 1"
        done
            let "iWid += 1"
        done

    ;;
    MAKE_SCANCARDS ) 
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh`

        python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/create2DScanDatacards.py \
            -i ${outdir}/datacards/shapes.root \
            -o ${outdir}/scancards/ \
            --nocards \
            --makeSplit

        for dist in ${dists[*]} ; do
            allcmd="python ${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/scripts/combineCards.py "
        for tlbCat in ${lbCat[*]} ; do
        for tlfs in ${lfs[*]} ; do
        for tCat in ${cat[*]} ; do
        for i in $(seq 3 27) ; do
            cardname="bin${i}_${tlbCat}${tlfs}${tCat}_${dist}=${outdir}/scancards/"
            cardname="${cardname}/datacard__1p0w_bin${i}_${tlbCat}${tlfs}${tCat}_${dist}.dat"
            allcmd="${allcmd} ${cardname} "
        done
        done
        done
        done

            allcmd="${allcmd} > ${outdir}/scancards/datacard_fullscan_${dist}.dat"
            eval $allcmd

        done
    ;;
    MAKE_SCANBIASCARDS ) 
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh`

        for injWid in ${nWids[*]} ; do
        for injMas in ${nMass[*]} ; do
            cmd="python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/create2DScanDatacards.py"
            cmd="${cmd} -i ${outdir}/datacards/shapes.root"
            cmd="${cmd} -o ${outdir}/scancards/"
            #cmd="${cmd} --injectMC 4.0,175.5"
            cmd="${cmd} --injectMC ${injWid},${injMas}"
            #cmd="${cmd} --nocards"
            cmd="${cmd} --makeSplit"
            
            bsub -q ${queue} \
            sh ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh \
                "${outdir}/" "${cmd}"
        done
        done

    ;;
    MERGE_SCANBIASCARDS )
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh`

        for injWid in ${nWids[*]} ; do
        for injMas in ${nMass[*]} ; do
        for dist in ${dists[*]} ; do
            allcmd="python ${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/scripts/combineCards.py "
        for tlbCat in ${lbCat[*]} ; do
        for tlfs in ${lfs[*]} ; do
        for tCat in ${cat[*]} ; do
        for i in $(seq 3 27) ; do
            cardname="bin${i}_${tlbCat}${tlfs}${tCat}_${dist}=${outdir}/scancards/"
            cardname="${cardname}/datacard__1p0w_bin${i}_${tlbCat}${tlfs}${tCat}_${dist}_w${injWid}_m${injMas}.dat"
            allcmd="${allcmd} ${cardname} "
        done
        done
        done
        done
            allcmd="${allcmd} > ${outdir}/scancards/datacard_fullscan_${dist}_w${injWid}_m${injMas}.dat"

            newdir=`echo biascards_${dist}_w${injWid}_m${injMas} | tr . p`
            
            cmd="eval $allcmd ;"
            cmd="${cmd} mkdir ${outdir}/scancards/$newdir ;"
            cmd="${cmd} mv ${outdir}/scancards/*__*${dist}_w${injWid}_m${injMas}.dat ${outdir}/scancards/$newdir/"

            bsub -q ${queue} \
            sh ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh \
                "${outdir}/" "${cmd}"
        done
        done
        done
    ;;
    SIMWORKSPACE )
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh`

        for dist in ${dists[*]} ; do
            text2workspace.py ${outdir}/scancards/datacard_fullscan_${dist}.dat -P \
                HiggsAnalysis.CombinedLimit.TopWidMassScan:topWidMassScan \
                -m 172.5 --PO verbose --PO coeffCache=${outdir}/scancards/coeff_cachefile.pck  \
                -o ${outdir}/scancards/fullscan_${dist}.root 
        done
    ;;
    WORKSPACE_SCANBIAS )
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh`

        for injWid in ${nWids[*]} ; do
        for injMas in ${nMass[*]} ; do
        for dist in ${dists[*]} ; do
            cmd="text2workspace.py ${outdir}/scancards/datacard_fullscan_${dist}_w${injWid}_m${injMas}.dat -P"
            cmd="${cmd} HiggsAnalysis.CombinedLimit.TopWidMassScan:topWidMassScan"
            cmd="${cmd} -m 172.5 --PO verbose --PO coeffCache=${outdir}/scancards/coeff_cachefile_w${injWid}_m${injMas}.pck"
            cmd="${cmd} -o ${outdir}/scancards/fullscan_${dist}_w${injWid}_m${injMas}.root"

            bsub -q ${queue} \
            sh ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh \
                "${outdir}/" "${cmd}"
        done
        done
        done
    ;;
    SIMSCAN )
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh`

        for massStep in $(seq 0 9) ; do
        for widStep in $(seq 0 9) ; do
        for dist in ${dists[*]} ; do
            massDn=`python -c "print  $massStep   *(-169.5+175.5)/10+169.5"`
            massUp=`python -c "print ($massStep+1)*(-169.5+175.5)/10+169.5"`
            widDn=`python -c "print  $widStep   *(-0.2648+5.296)/10+0.2648"`
            widUp=`python -c "print ($widStep+1)*(-0.2648+5.296)/10+0.2648"`

            echo "[$massDn,$massUp] : [$widDn,$widUp]"

            cmd="combine ${outdir}/fullscan_${dist}.root -M MultiDimFit --seed 8192" 
            cmd="${cmd} -m 172.5 -t -1 -P x -P m --floatOtherPOIs=1"
            cmd="${cmd} --setPhysicsModelParameterRanges m=$massDn,$massUp:x=$widDn,$widUp"
            cmd="${cmd} --algo=grid --points=20 -n x_m_scan__${dist}_$massStep_$widStep"
        
            bsub -q ${queue} \
            sh ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh \
                "${outdir}/" "${cmd}"
        done
        done
        done
    ;;
    SIMSCAN_BIAS )
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh`

        cd $outdir/scancards/

        for injWid in ${nWids[*]} ; do
        for injMas in ${nMass[*]} ; do
            nBias="w${injWid}_m${injMas}"
        for massStep in $(seq 0 9) ; do
        for widStep in $(seq 0 9) ; do
        for dist in ${dists[*]} ; do
            massDn=`python -c "print  $massStep   *(-169.5+175.5)/10+169.5"`
            massUp=`python -c "print ($massStep+1)*(-169.5+175.5)/10+169.5"`
            widDn=`python -c "print  $widStep   *(-0.2648+5.296)/10+0.2648"`
            widUp=`python -c "print ($widStep+1)*(-0.2648+5.296)/10+0.2648"`

            echo "[$massDn,$massUp] : [$widDn,$widUp]"

            cmd="combine ${outdir}/scancards/fullscan_${dist}_${nBias}.root -M MultiDimFit --seed 8192" 
            cmd="${cmd} -m 172.5 -P x -P m --floatOtherPOIs=1"
            cmd="${cmd} --setPhysicsModelParameterRanges m=$massDn,$massUp:x=$widDn,$widUp"
            cmd="${cmd} --algo=grid --points=20 -n x_m_scan__${dist}_${nBias}_$massStep$widStep"
        
            bsub -q ${queue} \
            sh ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh \
                "${outdir}/" "${cmd}"
        done
        done
        done
        done
        done
    ;;
    SIMSCAN_PARTS_MERGEBIAS )
        for injWid in ${nWids[*]} ; do
        for injMas in ${nMass[*]} ; do
            nBias="w${injWid}_m${injMas}"
            cmd="hadd $outdir/scancards/higgsCombinex_m_scan__${dist}_${nBias}.root"
            cmd="${cmd} $(ls $outdir/higgsCombinex_m_scan__${dist}_${nBias}_*.root)"
            cmd="${cmd} ; rm $(ls $outdir/higgsCombinex_m_scan__${dist}_${nBias}_*.root)"

            bsub -q ${queue} \
            sh ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh \
                "${outdir}/" "${cmd}"
        done
        done
    ;;
esac
