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
    echo "   {DO,PLOT}GEN - generator level studies"
    echo "   {TESTSEL,SEL,CHECK,MERGE,DO}TRIGEFF - trigger efficiency studies"
    echo "   {TEST,DO,CHECK,MERGE}SEL - selection, histogram filling"
    echo "   WWW - move available plots to a web-based area"
    exit 1; 
fi


if [[ ${ERA} == "2016" ]]; then
    lptlumi=28200
    hptlumi=35900
    year="2016"
    githash=0c522df
    eosdir=/store/cmst3/group/top/RunIIReReco/2016/${githash}  
    githash=1256677
    eosdir=/store/cmst3/group/top/RunIIReReco/2016/${githash}  
fi
if [[ ${ERA} == "2017" ]]; then
    lptlumi=7661
    hptlumi=41367
    year="2017"
    githash=ab05162
    eosdir=/store/cmst3/group/top/RunIIReReco/${githash}
    githash=1256677
    eosdir=/store/cmst3/group/top/RunIIReReco/2017/${githash}
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
#outdir=${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/vbf/${year}_${githash}
outdir=/eos/cms/store/cmst3/group/top/SMP-19-005/${year}_${githash}
wwwdir=/eos/user/p/psilva/www/SMP-19-005/

echo "Selection adapted to YEAR=${ERA}"
echo "Inputs from ${eosdir}"
echo "Outputs in ${outdir}"

RED='\e[31m'
NC='\e[0m'
case $WHAT in

#######################
### GEN STUDIES     ###
#######################
    
    DOGEN )
        python test/analysis/vbf/genLevelStudy.py /eos/cms/${eosdir} ${outdir}/gen 3
        ;;

    PLOTGEN )
       
        opts="--only cat --lumiSpecs cat0:${lptlumi},cat1:${hptlumi}"
	python scripts/plotter.py -i ${outdir}/gen ${opts} -j test/analysis/vbf/signal_${year}.json --noStack
	python scripts/plotter.py -i ${outdir}/gen -l 1    -j test/analysis/vbf/signal_${year}.json --only sel --saveTeX -o sel_plotter.root

        ;;

#######################
### TRIGGER STUDIES ###
#######################

    TESTSELTRIGEFF )
               
        json=test/analysis/vbf/samples_${ERA}.json
        tag=Data13TeV_${year}C_SinglePhoton
        if [[ ${year} == "2018" ]]; then
            tag=Data13TeV_${year}C_EGamma
        fi
        input=${eosdir}/${tag}/Chunk_1_ext0.root
        output=${tag}.root 

	python scripts/runLocalAnalysis.py \
            -i ${input} -o ${output} --tag ${tag} --only ${json} \
            --njobs 1 -q local --genWeights genweights_${githash}.root \
            --era era${ERA} -m PhotonTrigEff::RunPhotonTrigEff --ch 0 --runSysts --debug;

        ;;

    SELTRIGEFF )
        tags=SinglePhoton,EWKAJJ
        if [[ ${year} == "2018" ]]; then
            tags=EGamma,EWKAJJ
        fi
	python scripts/runLocalAnalysis.py \
	    -i ${eosdir} --only ${tags} \
            -o ${outdir}/trig \
            --farmappendix trig${githash} \
            -q ${queue} --genWeights genweights_${githash}.root \
            --era era${ERA} -m PhotonTrigEff::RunPhotonTrigEff --ch 0;
	;;

    CHECKTRIGEFFINTEG )
        python scripts/checkLocalAnalysisInteg.py ../../../FARMtrig${year}${githash}/
        ;;

    MERGETRIGEFF )
	./scripts/mergeOutputs.py ${outdir}/trig;
	;;

    DOTRIGEFF )
        python test/analysis/vbf/computeTriggerEff.py ${outdir}/trig;
        ;;

###################################
## SELECTION AND HISTO FILLING   ##
###################################

    TESTSEL )             
        tag=MC13TeV_${year}_EWKAJJ
        input=${eosdir}/${tag}/Chunk_0_ext0.root
	python scripts/runLocalAnalysis.py \
            -i ${input} -o testsmp19005_sel.root\
            --njobs 1 -q local --genWeights genweights_${githash}.root --tag ${tag}\
            --era era${ERA} -m SMP-19-005::SMP-19-005 --ch 0 --skimtree --runSysts --debug;
        ;;

    DOSEL )             
              
	python scripts/runLocalAnalysis.py \
            -i ${eosdir} -o ${outdir}/analysis \
            --farmappendix analysis${githash} \
            -q ${queue} --genWeights genweights_${githash}.root \
            --era era${ERA} -m SMP-19-005::SMP-19-005 --ch 0 --skimtree --runSysts;
        ;;

    CHECKSEL )
        python scripts/checkLocalAnalysisInteg.py ../../../FARManalysis${githash}/
        ;;
    
    MERGESEL )
	./scripts/mergeOutputs.py ${outdir}/analysis;
	;;

    WWW )
        pdirs=(gen trig)
        for d in ${pdirs[@]}; do           
            fdir=${wwwdir}/${year}_${githash}/${d}
	    mkdir -p ${fdir}
	    cp ${outdir}/${d}/plots/*.{png,pdf,dat} ${fdir}
	    cp test/index.php ${fdir};
        done        
        cp test/index.php ${wwwdir}/${year}_${githash};
        cp test/index.php ${wwwdir};
	;;
 
    *) 
        echo "Operation $WHAT is unknown"
        ;;

esac
