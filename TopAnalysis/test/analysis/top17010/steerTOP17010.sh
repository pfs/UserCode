#!/bin/bash

ERA=2016
STORAGE=""
FITTYPE="em_inc"
COMBINE=`dirname ${CMSSW_BASE}`/CMSSW_10_3_0_pre4
while getopts "o:y:s:f:c:" opt; do
    case "$opt" in
        o) WHAT=$OPTARG
            ;;
        y) ERA=$OPTARG
            ;;
        s) STORAGE=$OPTARG
            ;;
        f) FITTYPE=$OPTARG
            ;;
        c) COMBINE=$OPTARG
            ;;
    esac
done

if [ -z "$WHAT" ]; then 
    echo "steerTopWidth.sh -o <SEL/MERGE/...> [-y 2016/7] [-c higgs_combine] [-f fit_type] [-s storage directory]";
    echo "        PREPARE      - prepare analysis: resolution study + MC2MC corrections"
    echo "        TESTSEL      - test selection locally"
    echo "        SEL          - launches selection jobs to the batch, output will contain summary trees and control plots"; 
    echo "        SELSCAN      - launches signal selection jobs to the batch for the mass vs width scan";
    echo "        MERGE        - merge output"
    echo "        MERGESCAN    - merge output for the mass/width scan and plot local sensitivity"
    echo "        PLOT         - make plots (if given \"extra\" is appended to the directory)"
    echo "        BKG          - performs an estiation of the DY bacgrkound from data"
    echo "        TEMPL        - prepares the ROOT files with the template histograms"
    echo "        DATACARD     - prepares the datacards for combine"
    echo "        FIT          - this submits the fits to condor (additional instructions are printed"
    echo "        WWW          - move plots to web-based (if given \"extra\" is appended to the directory)"
    exit 1; 
fi

githash=0c522df
eosdir=/store/cmst3/group/top/RunIIReReco/2016/${githash}
fulllumi=35882
lumiUnc=0.025
if [[ ${ERA} = "2017" ]]; then
    githash=3129835
    eosdir=/store/cmst3/group/top/RunIIReReco/${githash}
    fulllumi=41367
fi

gtList=(0.7 0.75 0.8 0.85 0.9 0.95 1.0 1.05 1.1 1.15 1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.65 1.7 1.75 1.85 1.9 1.95 2.0 2.2 2.4 2.6 2.8 3.0 4.0)
mtList=(169.5 170 170.5 171 171.5 172 172.5 173 173.5 174 174.5 175 175.5)
dists=em_mlb,ee_mlb,mm_mlb
dists=${dists},emhighpt_mlb,emhighpt1b_mlb,emhighpt2b_mlb
dists=${dists},emlowpt_mlb,emlowpt1b_mlb,emlowpt2b_mlb
dists=${dists},eehighpt_mlb,eehighpt1b_mlb,eehighpt2b_mlb
dists=${dists},eelowpt_mlb,eelowpt1b_mlb,eelowpt2b_mlb
dists=${dists},mmhighpt_mlb,mmhighpt1b_mlb,mmhighpt2b_mlb
dists=${dists},mmlowpt_mlb,mmlowpt1b_mlb,mmlowpt2b_mlb

echo "Selection adapted to YEAR=${ERA}, inputs from ${eosdir}"

queue=workday
outdir=${STORAGE}
if [ -z ${outdir} ]; then
    outdir=${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/top17010
fi

json=test/analysis/top17010/samples_${ERA}.json
syst_json=test/analysis/top17010/syst_samples_${ERA}.json
wwwdir=${HOME}/www/top17010

RED='\e[31m'
NC='\e[0m'
case $WHAT in

    PREPARE )
        
        echo "Deriving baseline expected resolutions"
        tag=MC13TeV_${ERA}_TTJets
        python test/analysis/top17010/getMlbResolution.py /eos/cms/${eosdir}/${tag}/Chunk_0_ext0.root
        mkdir -p ${outdir}/plots
        mv mlbresol_*.{png,pdf} ${outdir}/plots
        
        echo "Computing MC2MC corrections"
        echo "[WARN] currently hardcoded for 2016 samples"
        python test/analysis/top17010/prepareMC2MCCorrections.py;
        ;;

    TESTSEL )               
        tag=MC13TeV_${ERA}_TTJets
        input=${eosdir}/${tag}/Chunk_0_ext0.root
        output=${tag}.root 

        gidx=`python -c "print int((2-0.7)/0.01)"`
        midx=`python -c "print int((172.5-169)/0.25)"`
        flag=`python -c "print (($midx<<16)|($gidx))"`

	python scripts/runLocalAnalysis.py \
            -i ${input} -o ${output} --tag ${tag} --only ${json} --flag ${flag}\
            --njobs 1 -q local --genWeights genweights_${githash}.root \
            --era era${ERA} -m TOP17010::TOP17010 --ch 0 --runSysts --debug;
        ;;

    SEL )
	python scripts/runLocalAnalysis.py \
	    -i ${eosdir} --only ${json},${syst_json} --flag 0 \
            -o ${outdir}/${githash} \
            --farmappendix ${githash} \
            -q ${queue} --genWeights genweights_${githash}.root \
            --era era${ERA} -m  TOP17010::TOP17010 --ch 0 --runSysts;
	;;


    SELSCAN )
        
        for g in ${gtList[@]}; do
            gidx=`python -c "print int(($g-0.7)/0.01)"`
            for m in ${mtList[@]}; do
                midx=`python -c "print int(($m-169)/0.25)"`
                flag=`python -c "print (($midx<<16)|($gidx))"`

	        python scripts/runLocalAnalysis.py \
	            -i ${eosdir} --only MC13TeV_${ERA}_TTJets --exactonly --flag ${flag} \
                    -o ${outdir}/${githash}/scenario${flag} \
                    --farmappendix ${githash}SCAN${flag} \
                    -q ${queue} --genWeights genweights_${githash}.root \
                    --era era${ERA} -m  TOP17010::TOP17010 --ch 0 --runSysts;

            done
        done
	;;


    MERGE )
	./scripts/mergeOutputs.py ${outdir}/${githash};
	;;


    MERGESCAN )
        for g in ${gtList[@]}; do
            gidx=`python -c "print int(($g-0.7)/0.01)"`
            for m in ${mtList[@]}; do
                midx=`python -c "print int(($m-169)/0.25)"`
                flag=`python -c "print (($midx<<16)|($gidx))"`
	        ./scripts/mergeOutputs.py ${outdir}/${githash}/scenario${flag};
            done
        done

        #local sensitivities
        python test/analysis/top17010/estimateLocalSensitivity.py -i ${outdir}/${githash} -o ${outdir}/${githash}/localsens/
	;;

    PLOT )
	commonOpts="-i ${outdir}/${githash} --puNormSF puwgtctr -l ${fulllumi} --saveLog --mcUnc ${lumiUnc}"
	python scripts/plotter.py ${commonOpts} -j ${json};
        python scripts/plotter.py ${commonOpts} -j ${json}      --only evcount  --saveTeX -o evcount_plotter.root;
        python scripts/plotter.py ${commonOpts} -j ${json}      --only mlb,ptlb --binWid  -o lb_plotter.root;
        python scripts/plotter.py ${commonOpts} -j ${syst_json} --only mlb      --silent  -o syst_plotter.root;

        ;;

    BKG )
        python test/analysis/top17010/estimateDY.py -i ${outdir}/${githash}/plots/plotter.root -o ${outdir}/${githash}/plots/;
        ;;

    TEMPL )
        
        inputs=${outdir}/${githash}/plots/plotter.root
        inputs=${inputs},${outdir}/${githash}/plots/syst_plotter.root
        inputs=${inputs},${outdir}/${githash}/plots/plotter_dydata.root
        output=${outdir}/${githash}/templates/
        python test/analysis/top17010/prepareTemplateFiles.py -i ${inputs} -d ${dists} -o ${output} --debug --bbbThr 0.005;

        ;;

    TESTFIT )

        testCat=em
        testDist=${testCat}_mlb
        python test/analysis/top17010/prepareDataCard.py -d ${testDist} \
            -t ${outdir}/${githash}/templates \
            dataDef=sig,${outdir}/${githash}/plots/plotter.root,${testDist}/${testDist}_t#bar{t} \
            -s nom,${outdir}/${githash}/MC13TeV_${ERA}_TTJets.root \
            -o ${outdir}/${githash}/datacards_jerinc \
            --systs test/analysis/top17010/systs_dict_incjer.json
        
        args="${outdir}/${githash}/datacards_jerinc/${testCat}/nom/tbart.datacard.dat"
        python test/analysis/top17010/createFit.py -o ${outdir}/${githash}/fits_jerinc/${testCat}/nom -a -t 50 -c ${COMBINE} --tag tbart ${args}

        sh ${outdir}/${githash}/fits_jerinc/${testCat}/nom/runFit_tbart.sh 

        ;;


    DATACARD )

        python test/analysis/top17010/submitPrepareDataCard.py --dists ${dists} \
            --systs ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/top17010/systs_dict.json \
            --templ ${outdir}/${githash}/templates \
            --nom MC13TeV_${ERA}_TTJets.root \
            -o ${outdir}/${githash}/datacards
        ;;


    FIT )
        
        echo ""
        echo "The fits have to be run with a dedicated script"
        echo "An example can be found test/analysis/top17010/steerFits.sh"
        echo "(Note that the createFit.py script called there needs to be pointed to the location of the installation of the combination tool)"
        echo "If not yet installed please install the combination tool (TensorFlow version) as follows"
        echo ""
        echo "cmsrel CMSSW_10_3_0_pre4"
        echo "cd CMSSW_10_3_0_pre4/src/"
        echo "cmsenv"
        echo "git clone -b tensorflowfit https://github.com/bendavid/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit"
        echo "cd HiggsAnalysis/CombinedLimit"
        echo "scram b -j 8"
        echo "source /afs/cern.ch/user/b/bendavid/work/cmspublic/pythonvenv/tensorflowfit_h5py/bin/activate"
        echo "" 
        echo "Having that said I will now try to run with the following combine location $COMBINE"
        echo "" 

        anchors=(`ls ${outdir}/${githash}/datacards/em`)
        signals=(`ls ${outdir}/${githash}/datacards/em/nom/*datacard.dat`)

        echo "${#anchors[@]} anchors for ${#signals[@]} data/signals... it may take a while to prepare all the scripts to submit"

        echo "Preparing submission of create fit scripts (A)"
        condor_prep=fitprep${FITTYPE}_condor.sub
        echo "executable  = ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/scripts/wrapLocalAnalysisRun.sh" > $condor_prep
        echo "output      = ${condor_prep}.out" >> $condor_prep
        echo "error       = ${condor_prep}.err" >> $condor_prep
        echo "log         = ${condor_prep}.log" >> $condor_prep
        echo "+JobFlavour = \"workday\"" >> $condor_prep        
        for a in ${anchors[@]}; do
            a=`basename ${a}`;

            for s in ${signals[@]}; do                
                out="${outdir}/${githash}/fits/${FITTYPE}/${a}"
                tag=`basename $s | cut -f -1 -d "."`
                
                args=""
                if [ "${FITTYPE}" == "em_inc" ]; then
                    args="${outdir}/${githash}/datacards/em/${a}/${tag}.datacard.dat"
                elif [ "${FITTYPE}" == "dil_inc" ]; then
                    args="em=${outdir}/${githash}/datacards/em/${a}/${tag}.datacard.dat"
                    args="${args} mm=${outdir}/${githash}/datacards/mm/${a}/${tag}.datacard.dat"
                    args="${args} ee=${outdir}/${githash}/datacards/ee/${a}/${tag}.datacard.dat"
                elif [ "${FITTYPE}" == "ptlb_inc" ]; then
                    args="emhighpt=${outdir}/${githash}/datacards/emhighpt/${a}/${tag}.datacard.dat"
                    args="${args} mmhighpt=${outdir}/${githash}/datacards/mmhighpt/${a}/${tag}.datacard.dat"
                    args="${args} eehighpt=${outdir}/${githash}/datacards/eehighpt/${a}/${tag}.datacard.dat"
                    args="${args} emlowpt=${outdir}/${githash}/datacards/emlowpt/${a}/${tag}.datacard.dat"
                    args="${args} mmlowpt=${outdir}/${githash}/datacards/mmlowpt/${a}/${tag}.datacard.dat"
                    args="${args} eelowpt=${outdir}/${githash}/datacards/eelowpt/${a}/${tag}.datacard.dat"
                elif [ "${FITTYPE}" == "final" ]; then
                    args="emhighpt2b=${outdir}/${githash}/datacards/emhighpt2b/${a}/${tag}.datacard.dat"
                    args="${args} mmhighpt2b=${outdir}/${githash}/datacards/mmhighpt2b/${a}/${tag}.datacard.dat"
                    args="${args} eehighpt2b=${outdir}/${githash}/datacards/eehighpt2b/${a}/${tag}.datacard.dat"
                    args="${args} emhighpt1b=${outdir}/${githash}/datacards/emhighpt1b/${a}/${tag}.datacard.dat"
                    args="${args} mmhighpt1b=${outdir}/${githash}/datacards/mmhighpt1b/${a}/${tag}.datacard.dat"
                    args="${args} eehighpt1b=${outdir}/${githash}/datacards/eehighpt1b/${a}/${tag}.datacard.dat"
                    args="${args} emlowpt2b=${outdir}/${githash}/datacards/emlowpt2b/${a}/${tag}.datacard.dat"
                    args="${args} mmlowpt2b=${outdir}/${githash}/datacards/mmlowpt2b/${a}/${tag}.datacard.dat"
                    args="${args} eelowpt2b=${outdir}/${githash}/datacards/eelowpt2b/${a}/${tag}.datacard.dat"
                    args="${args} emlowpt1b=${outdir}/${githash}/datacards/emlowpt1b/${a}/${tag}.datacard.dat"
                    args="${args} mmlowpt1b=${outdir}/${githash}/datacards/mmlowpt1b/${a}/${tag}.datacard.dat"
                    args="${args} eelowpt1b=${outdir}/${githash}/datacards/eelowpt1b/${a}/${tag}.datacard.dat"
                else
                    echo "fit type=${FITTYPE} is not yet implemented... quitting"
                    exit -1
                fi

                extraOpts=""
                if [ "${a}" == "nom" ]; then
                    extraOpts="-a -t 50"
                fi

                echo "arguments   = ${CMSSW_BASE} ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/top17010/createFit.py -o ${out} ${extraOpts} -c ${COMBINE} --tag ${tag} ${args}" >> $condor_prep
                echo "queue 1" >> $condor_prep
            done
        done

        echo "Preparing submission of fits (B)"
        condor_fit=fit${FITTYPE}_condor.sub
        echo "executable  = ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/top17010/runFitWrapper.sh" > $condor_fit
        echo "output      = ${condor_fit}.out" >> $condor_fit
        echo "error       = ${condor_fit}.err" >> $condor_fit
        echo "log         = ${condor_fit}.log" >> $condor_fit
        echo "+JobFlavour = \"workday\"" >> $condor_fit        
        for a in ${anchors[@]}; do
            dir="${outdir}/${githash}/fits/${FITTYPE}/${a}"
            echo "arguments   = ${dir}/runFit_\$(tag).sh" >> $condor_fit
            echo "queue tag from (" >> $condor_fit
            for s in ${signals[@]}; do
                tag=`basename $s | cut -f -1 -d "."`
                echo " ${tag}" >> $condor_fit
            done
            echo ")" >> $condor_fit
        done        

        echo "Submit ${condor_prep}"
        echo "Once done submit ${condor_fit}"

        #echo "Submitting both jobs as a simple DAG (A->B)"
        #condor_dag=fit${FITTYPE}_condor.dag
        #echo "JOB A ${condor_prep}" > $condor_dag
        #echo "JOB B ${condor_fit}" >> $condor_dag
        #echo "PARENT A CHILD B" >> $condor_dag
        #condor_submit_dag $condor_dag

        ;;

    WWW )
        pdir=${outdir}/${githash}
        if [ -d ${pdir} ]; then
            fdir=${wwwdir}/${githash}

            #selection plots
	    mkdir -p ${fdir}
	    mv ${pdir}/plots/*.{png,pdf,dat} ${fdir};
	    cp test/index.php ${fdir};

            #additional plots
            for e in localsens templates; do
                mkdir -p ${fdir}/${e};
	        mv ${pdir}/${e}/*.{png,pdf,dat} ${fdir}/${e};
	        cp test/index.php ${fdir}/${e};
            done

            echo "Check plots in ${fdir}"
        fi
        
	;;
esac
