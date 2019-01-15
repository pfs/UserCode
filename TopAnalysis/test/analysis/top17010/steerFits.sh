#!/bin/bash

FITSCRIPT=test/analysis/top17010/createFit.py
OUTDIR=store/TOP17010/fit_results

function run_fit {
    echo "Running fit for $1"
    cd $1
    sh runFit.sh
    cd -
}

function plot_results {
    echo "Plotting results for ${1}"
    python test/analysis/top17010/plotLikelihoodScanResults.py -i ${1} -o ${1}
    python test/analysis/top17010/doNuisanceReport.py -o ${1} \
        Asimov=${1}/nom/fitresults_asimov.root \
        Toys=${1}/nom/fitresults_toys.root
}

#loop over all the available scenarios
a=(`ls em_datacards`)
for i in ${a[@]}; do 

    baseOpts=""
    if [ ${i} == "nom" ]; then
        #baseOpts="${baseOpts} -t 100 -a";
        baseOpts="${baseOpts} -a";
    fi
 
    #emu inclusive
    out=$OUTDIR/em_inc/${i};
    $FITSCRIPT -o ${out} ${baseOpts} em_datacards/${i}/datacard.dat;
    #run_fit ${out};

    #emu+ee+mm inclusive
    out=$OUTDIR/dil_inc/${i};
    $FITSCRIPT -o ${out} ${baseOpts} \
        em=em_datacards/${i}/datacard.dat ee=ee_datacards/${i}/datacard.dat mm=mm_datacards/${i}/datacard.dat;
    #run_fit ${out};

    #highpt+lowpt
    out=$OUTDIR/ptlb_inc/${i};
    $FITSCRIPT -o ${out} ${baseOpts} \
        emhighpt=emhighpt_datacards/${i}/datacard.dat emlowpt=emlowpt_datacards/${i}/datacard.dat \
        mmhighpt=mmhighpt_datacards/${i}/datacard.dat mmlowpt=mmlowpt_datacards/${i}/datacard.dat \
        eehighpt=eehighpt_datacards/${i}/datacard.dat eelowpt=eelowpt_datacards/${i}/datacard.dat
    run_fit ${out};

    #full cat
    out=$OUTDIR/final/${i};
    $FITSCRIPT -o ${out} ${baseOpts} \
        emhighpt2b=emhighpt2b_datacards/${i}/datacard.dat emhighpt1b=emhighpt1b_datacards/${i}/datacard.dat \
        emlowpt2b=emlowpt2b_datacards/${i}/datacard.dat   emlowpt1b=emlowpt1b_datacards/${i}/datacard.dat \
        mmhighpt2b=mmhighpt2b_datacards/${i}/datacard.dat mmhighpt1b=mmhighpt1b_datacards/${i}/datacard.dat \
        mmlowpt2b=mmlowpt2b_datacards/${i}/datacard.dat   mmlowpt1b=mmlowpt1b_datacards/${i}/datacard.dat \
        eehighpt2b=eehighpt2b_datacards/${i}/datacard.dat eehighpt1b=eehighpt1b_datacards/${i}/datacard.dat \
        eelowpt2b=eelowpt2b_datacards/${i}/datacard.dat   eelowpt1b=eelowpt1b_datacards/${i}/datacard.dat
    run_fit ${out};
    
done


#plot likelihood contours, compare postfit-nuisances
for r in em_inc dil_inc ptlb_inc final; do
    plot_results store/TOP17010/fit_results/${r};
done

#comparison between different fit types
python test/analysis/top17010/doNuisanceReport.py \
    -o store/TOP17010/fit_results/ \
    e#mu=store/TOP17010/fit_results/em_inc/nom/fitresults_asimov.root \
    ll=store/TOP17010/fit_results/dil_inc/nom/fitresults_asimov.root

#$FITSCRIPT -o $OUTDIR/em_ptcats -a \
#    emhighpt=emhighpt_datacards/nom/datacard.dat emlowpt=emlowpt_datacards/nom/datacard.dat 
#run_fit $OUTDIR/em_ptcats

#$FITSCRIPT -o $OUTDIR/em_ptbcats -a \
#    emhighpt2b=emhighpt2b_datacards/nom/datacard.dat emhighpt1b=emhighpt1b_datacards/nom/datacard.dat \
#    emlowpt2b=emlowpt2b_datacards/nom/datacard.dat emlowpt1b=emlowpt1b_datacards/nom/datacard.dat 
#run_fit $OUTDIR/em_ptbcats

