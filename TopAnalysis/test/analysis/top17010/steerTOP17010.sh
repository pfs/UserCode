#!/bin/bash

ERA=2016

while getopts "o:e:q:y:" opt; do
    case "$opt" in
        o) WHAT=$OPTARG
            ;;
        y) ERA=$OPTARG
            ;;
    esac
done

if [ -z "$WHAT" ]; then 
    echo "steerTopWidth.sh -o <SEL/MERGE/...> [-y 2016/7]";
    echo "        RESOL        - preliminary resolution study"
    echo "        TESTSEL      - test selection locally"
    echo "        SEL          - launches selection jobs to the batch, output will contain summary trees and control plots"; 
    echo "        SELSCAN      - launches signal selection jobs to the batch for the mass vs width scan";
    echo "        MERGE        - merge output (if given \"extra\" is appended to the directory)"
    echo "        PLOT         - make plots (if given \"extra\" is appended to the directory)"
    echo "        WWW          - move plots to web-based (if given \"extra\" is appended to the directory)"
    exit 1; 
fi

githash=0c522df
eosdir=/store/cmst3/group/top/RunIIReReco/2016/${githash}
fulllumi=35505
lumiUnc=0.025
if [[ ${ERA} = "2017" ]]; then
    githash=3129835
    eosdir=/store/cmst3/group/top/RunIIReReco/${githash}
    fulllumi=41367
fi
echo "Selection adapted to YEAR=${ERA}, inputs from ${eosdir}"

queue=workday
outdir=${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/top17010
json=test/analysis/top17010/samples_${ERA}.json
syst_json=test/analysis/top17010/syst_samples_${ERA}.json
wwwdir=${HOME}/www/top17010

mkdir -p ${outdir}

RED='\e[31m'
NC='\e[0m'
case $WHAT in

    RESOL )
        tag=MC13TeV_${ERA}_TTJets_m1755
        python test/analysis/top17010/getMlbResolution.py /eos/cms/${eosdir}/${tag}/Chunk_0_ext0.root
        mkdir -p ${outdir}/plots
        mv mlbresol_*.{png,pdf} ${outdir}/plots
        ;;

    TESTSEL )               
        tag=MC13TeV_${ERA}_TTJets
        input=${eosdir}/${tag}/Chunk_0_ext0.root
        output=${tag}.root 
	python scripts/runLocalAnalysis.py \
            -i ${input} -o ${output} --tag ${tag} --only ${json} --flag 0\
            --njobs 1 -q local --genWeights genweights_${githash}.root \
            --era era${ERA} -m TOP17010::TOP17010 --ch 0 --runSysts; --debug;
        ;;

    SEL )
        
	python scripts/runLocalAnalysis.py \
	    -i ${eosdir} --only ${json},${syst_json} --flag 0\
            -o ${outdir}/${githash} \
            --farmappendix ${githash} \
            -q ${queue} --genWeights genweights_${githash}.root \
            --era era${ERA} -m  TOP17010::TOP17010 --ch 0 --runSysts;
	;;


    SELSCAN )
        
        gtList=(0.7 0.75 0.8 0.85 0.9 0.95 1.0 1.05 1.1 1.15 1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.65 1.7 1.75 1.85 1.9 1.95 2.0 2.2 2.4 2.6 2.8 3.0 4.0)
        mtList=(169.5 170 170.5 171 171.5 172 172.5 173 173.5 174 174.5 175 175.5)
        for g in ${gtList[@]}; do
            gidx=`python -c "print int(($g-0.7)/0.01)"`
            for m in ${mtList[@]}; do
                midx=`python -c "print int(($m-169)/0.25)"`
                flag=`python -c "print (($midx<<16)|($gidx))"`

	        python scripts/runLocalAnalysis.py \
	            -i ${eosdir} --only MC13TeV_${ERA}_TTJets --exactonly --flag ${flag} \
                    -o ${outdir}/${githash} \
                    --farmappendix ${githash}SCAN${flag} \
                    -q ${queue} --genWeights genweights_${githash}.root \
                    --era era${ERA} -m  TOP17010::TOP17010 --ch 0 --runSysts;

            done
        done
	;;


    MERGE )
	./scripts/mergeOutputs.py ${outdir}/${githash};
	;;

    PLOT )

	commonOpts="-i ${outdir}/${githash} --puNormSF puwgtctr -l ${fulllumi} --saveLog --mcUnc ${lumiUnc}"
	python scripts/plotter.py ${commonOpts} -j ${json};
        python scripts/plotter.py ${commonOpts} -j ${json}      --only evcount  --saveTeX -o evcount_plotter.root;
        python scripts/plotter.py ${commonOpts} -j ${json}      --only mlb,ptlb --binWid  -o lb_plotter.root;
        python scripts/plotter.py ${commonOpts} -j ${syts_json} --only mlb      --silent  -o syst_plotter.root;

        ;;
    
    WWW )
        pdir=${outdir}/${githash}/plots
        if [ -d ${pdir} ]; then
            fdir=${wwwdir}/${githash}
	    mkdir -p ${fdir}
	    cp ${pdir}/*.{png,pdf,dat} ${fdir};
	    cp test/index.php ${fdir};
            echo "Check plots in ${fdir}"
        fi
        
	;;
esac
