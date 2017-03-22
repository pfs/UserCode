#!/bin/bash
echo ""
echo "Running TOP-16-019 analysis."

########################  SETTINGS  #########################

#lfs=$(EE EM MM)
wids=("0p2w") 


########################    SETUP    ########################

cd ${CMSSW_BASE}/src/CombineHarvester/TOP16019/scripts/

case $1 in
MAKE_CARDS)
######################## MAKE CARDS #########################

python makeDataCard.py 

;;
WORKSPACE)
########################  WORKSPACE #########################

getMorphFiles.py
text2workspace.py -b datacards/datacard.txt -o morph.root --default-morphing shape2

;;
COMBINE)
######################## RUN COMBINE ########################

combine -M MaxLikelihoodFit --redefineSignalPOIs gamma morph.root

;;
PLOT)
######################## PLOT OUTPUT ########################
;;

*)
####################### INSTRUCTIONS ########################
if [[ "$1" == "" ]] ; then 
    echo ""
    echo "    Options are:"
    echo "      - MAKE_CARDS: produces datacards"
    echo "      - WORKSPACE: creates the morph workspace" 
    echo "      - COMBINE: runs combine on the full dataset"
    echo "      - PLOT: plots output" 
    echo "    Please enter an option and try again."
    echo ""
    exit 1;
fi
;;
esac
