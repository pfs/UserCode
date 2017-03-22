#!/bin/zsh
echo ""
echo "Running TOP-16-019 analysis."

########################  SETTINGS  #########################

#no emails
export JOB_LSB_REPORT_MAIL=N

queue=1nd
#wids=("0p4w")

wids=("0p2w" "0p4w" "0p6w" "0p8w" "1p0w" 
      "1p2w" "1p4w" "1p6w" "1p8w" "2p0w"
      "2p2w" "2p4w" "2p6w" "2p8w" "3p0w")
ptcs=("highpt" "lowpt")
lfss=("EE" "EM" "MM")
btcs=("1b" "2b")

if [[ "$2" == "TEST" ]] ; then
    wids=("0p4w")
    ptcs=("highpt")
    lfss=("EE")
    btcs=("1b")
fi

#########################  CONFIG  ##########################

function join {local IFS=","; echo "$*"; }

widStr="$(join ${wids[@]})"
lfsStr="$(join ${lfss[@]})"
ptcStr="$(join ${ptcs[@]})"
btcStr="$(join ${btcs[@]})"

echo "${widStr}" 
echo "${lfsStr}" 
echo "${ptcStr}" 
echo "${btcStr}" 

########################    SETUP    ########################

cd ${CMSSW_BASE}/src/CombineHarvester/TOP16019/scripts/

case $1 in
MAKE_FAKES)
######################## MAKE  BIAS #########################

python makeFakeShapes.py    \
    --widStr ${widStr}      \
    --ptcStr ${ptcStr}      \
    --lfsStr ${lfsStr}      \
    --btcStr ${btcStr}

;;
MAKE_CARDS)
######################## MAKE CARDS #########################

for twid in ${wids[*]} ; do
    python makeDataCard.py  --injWid ${twid}        \
                            --injMas 172.5          \
                            --ptcStr ${ptcStr}      \
                            --lfsStr ${lfsStr}      \
                            --btcStr ${btcStr}
done

;;
WORKSPACE)
########################  WORKSPACE #########################

for bwid in ${wids[*]} ; do
    biasDir=datacards/bias_${bwid}_172.5

    ./getMorphFiles.py -i "${biasDir}/*/*.txt" \
                     -o ${biasDir} \
                     --biasWid ${bwid}

    text2workspace.py -b ${biasDir}/morphDatacard.txt \
                      -o ${biasDir}/morph.root \
                      --default-morphing shape2
done

;;
COMBINE)
######################## RUN COMBINE ########################

for bwid in ${wids[*]} ; do
    biasDir=datacards/bias_${bwid}_172.5/
    cd ${biasDir}

    combine -M MaxLikelihoodFit \
            --redefineSignalPOIs gamma \
            morph.root
done

;;
PLOT)
######################## PLOT OUTPUT ########################
;;

JOBS)
########################## LXBATCH ##########################

for bwid in ${wids[*]} ; do
    cmd=""

    cmd="python makeDataCard.py --injWid ${bwid}"
    cmd="${cmd} --injMas 172.5" 
    cmd="${cmd} --ptcStr ${ptcStr}"
    cmd="${cmd} --lfsStr ${lfsStr}"
    cmd="${cmd} --btcStr ${btcStr} ;"

    biasDir=datacards/bias_${bwid}_172.5

    cmd="${cmd} ./getMorphFiles.py -i \\\"${biasDir}/*/*.txt\\\""
    cmd="${cmd}                  -o ${biasDir}"
    cmd="${cmd}                  --biasWid ${bwid} ;"

    cmd="${cmd} text2workspace.py -b ${biasDir}/morphDatacard.txt"
    cmd="${cmd}                   -o ${biasDir}/morph.root"
    cmd="${cmd}                   --default-morphing shape2 ;"

    cmd="${cmd} cd ${biasDir} ;"

    cmd="${cmd} combine -M MaxLikelihoodFit"
    cmd="${cmd} --redefineSignalPOIs gamma"
    cmd="${cmd} -n TopWid1D_w${bwid}_m172.5"
    cmd="${cmd} ${biasDir}/morph.root"

    bsub -q ${queue} \
        "sh ${PWD}/../wrapLXBATCHJob.sh ${PWD} \"${cmd}\""
done

;;
MOVEFILES)

    mkdir biasplots/
    mv higgsCombineTopWid1D_w*_m* biasplots/
    mv mlfitTopWid1D_w*_m*        biasplots/
    mv LSF* biasplots/



;;

*)
####################### INSTRUCTIONS ########################
if [[ "$1" == "" ]] ; then 
    echo ""
    echo "    Options are:"
    echo "      - MAKE_FAKES: produces MC templates to inject as data"
    echo "      - MAKE_CARDS: produces datacards"
    echo "      - WORKSPACE: creates the morph workspace" 
    echo "      - COMBINE: runs combine on the full dataset"
    echo "      - PLOT: plots output" 
    echo "      - JOBS: run as LXBATCH jobs (CARDS-->COMBINE)" 
    echo "    Please enter an option and try again."
    echo ""
    echo "    To enable testmode, enter TEST as a second argument."
    exit 1;
fi
;;
esac
