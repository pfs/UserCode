#!/bin/bash

while getopts "o:e:q:y:d:" opt; do
    case "$opt" in
        o) WHAT=$OPTARG
            ;;
        y) ERA=$OPTARG
            ;;
    esac
done

if [ -z "$WHAT" ]; then 
    echo "steerVBFVectorBoson.sh -o <SEL/MERGE/...>  -y 2016/7/8";
    echo "   {TESTSEL,SEL,CHECK,MERGE,DO,WWW}TRIGEFF - trigger efficiency studies"
    exit 1; 
fi


if [[ ${ERA} == "2016" ]]; then
    year="2016"
    githash=0c522df
    eosdir=/store/cmst3/group/top/RunIIReReco/2016/${githash}  
fi
if [[ ${ERA} == "2017" ]]; then
    year="2017"
    githash=ab05162
    eosdir=/store/cmst3/group/top/RunIIReReco/${githash}
fi
if [[ ${ERA} == "2017ul" ]]; then
    year="2017"
    githash=7120040_ul
    eosdir=/store/cmst3/group/top/RunIIReReco/2017/${githash}
fi
if [[ ${ERA} == "2018" ]]; then
    year="2018"
    githash=219fd1e
    eosdir=/store/cmst3/group/top/RunIIReReco/2018/${githash}
fi

queue=workday
outdir=${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/VBFVectorBoson
wwwdir=/eos/user/p/psilva/www/SMP-19-005/${githash}

echo "Selection adapted to YEAR=${ERA}"
echo "Inputs from ${eosdir}"
echo "Outpits in ${outdir}"

RED='\e[31m'
NC='\e[0m'
case $WHAT in

#######################
### TRIGGER STUDIES ###
#######################

    TESTSELTRIGEFF )
               
        json=test/analysis/vbf/samples_${ERA}.json
        tag=Data13TeV_${year}C_SinglePhoton
        input=${eosdir}/${tag}/Chunk_1_ext0.root
        output=${tag}.root 

	python scripts/runLocalAnalysis.py \
            -i ${input} -o ${output} --tag ${tag} --only ${json} --mvatree\
            --njobs 1 -q local --genWeights genweights_${githash}.root \
            --era era${ERA} -m PhotonTrigEff::RunPhotonTrigEff --ch 0 --runSysts --debug;

        ;;

    SELTRIGEFF )
	python scripts/runLocalAnalysis.py \
	    -i ${eosdir} --only SinglePhoton,EWKAJJ\
            -o ${outdir}/trig/${githash}/${EXTRA} \
            --farmappendix trig${githash} \
            -q ${queue} --genWeights genweights_${githash}.root \
            --era era${ERA} -m PhotonTrigEff::RunPhotonTrigEff --ch 0 --runSysts ${extraOpts};
	;;

    CHECKTRIGEFFINTEG )
        python scripts/checkLocalAnalysisInteg.py ../../../FARMtrig${githash}/
        ;;

    MERGETRIGEFF )
	./scripts/mergeOutputs.py ${outdir}/trig/${githash}/${EXTRA};
	;;

    DOTRIGEFF )
        python test/analysis/computeTriggerEff.py ${outdir}/trig/${githash}/${EXTRA};
        ;;

    WWWTRIGEFF )
        fdir=${wwwdir}/${githash}/${EXTRA}/trigeff
	mkdir -p ${fdir}
	cp ${outdir}/trig/${githash}/${EXTRA}/trigeff/*.{png,pdf} ${fdir}
	cp test/index.php ${fdir};        
	;;

esac
