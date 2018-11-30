#!/bin/bash

WHAT=$1; 
EXTRA=$2
QCD=""

if [ "$#" -lt 1 ]; then 
    echo "steerVBFVectorBoson.sh <SEL/MERGE/PLOT/WWW> [extra]";
    echo "        SEL          - launches selection jobs to the batch, output will contain summary trees and control plots"; 
    echo "        MERGE        - merge output (if given \"extra\" is appended to the directory)"
    echo "        PLOT         - make plots (if given \"extra\" is appended to the directory)"
    echo "        NLOTFACTORS  - compute NLO/LO in DY and test how it looks after applying to gamma+jets"
    echo "        TRIGEFF      - trigger efficiency"
    echo "        BDTTRANSFORM - derive BDT transformation for uniform background yields"
    echo "        WWW          - move plots to web-based (if given \"extra\" is appended to the directory)"
    exit 1; 
fi

if [ "$#" -gt 2 ]; then 
    QCD=$3
fi
#to run locally use local as queue + can add "--njobs 8" to use 8 parallel jobs
queue=workday
githash=3129835
eosdir=/store/cmst3/group/top/RunIIReReco/${githash}
fulllumi=41367
vbflumi=7661
lumiUnc=0.025
outdir=${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/VBFVectorBoson
wwwdir=~/www/VBFVectorBoson
era="2017"

RED='\e[31m'
NC='\e[0m'
case $WHAT in

    TESTSEL )
               
        json=data/era2017/vbf_samples.json
        tag=MC13TeV_2017_EWKAJJ
        #tag=Data13TeV_2017B_DoubleEG
        input=${eosdir}/${tag}/Chunk_0_ext0.root
        output=${tag}.root 

	python scripts/runLocalAnalysis.py \
            -i ${input} -o ${output} --tag ${tag} --only ${json} --mvatree\
            --njobs 1 -q local --genWeights genweights_${githash}.root \
            --era era2017 -m VBFVectorBoson::RunVBFVectorBoson --ch 0 --runSysts --debug;

        #--debug --mvatree \
        ;;

    SEL )
	##### NOTE: There are three options here:
        ### --mvatree: to store trees for BDT training in signal region
        ### --CR     : gives a control region to evaluate fake rates in the photon data samples
        ### --SRfake : gives the distributions of fakes, normalised based on fake rates

        json=data/era2017/vbf_samples.json,data/era2017/vbf_syst_samples.json
	if [[ -z ${EXTRA} ]]; then
	    extraOpts=" --mvatree"
        fi
	python scripts/runLocalAnalysis.py \
	    -i ${eosdir} --only ${json}\
            -o ${outdir}/${githash}/${EXTRA} \
            --farmappendix ${githash} \
            -q ${queue} --genWeights genweights_${githash}.root \
            --era era2017 -m VBFVectorBoson::RunVBFVectorBoson --ch 0 --runSysts ${extraOpts};
	;;


    SELJETHT )
	json=data/era2017/JetHT.json;
	extraOpts=" --CR"
	echo ${QCD} 
	if [[ ${QCD} == "QCDTemp" ]]; then
	    echo 'I do QCD Template photon selection'
	    extraOpts=${extraOpts}" --QCDTemp"
	fi
        if [[ ${QCD} == "SRfake" ]]; then
            echo 'I do SRfake photon selection'
            extraOpts=" --SRfake"
        fi
	python scripts/runLocalAnalysis.py \
	    -i ${eosdir} --only ${json} \
            -o ${outdir}/${githash}/${EXTRA}${QCD} \
            --farmappendix ${githash}${EXTRA}${QCD} \
            -q ${queue} --genWeights genweights_${githash}.root \
            --era era2017 -m VBFVectorBoson::RunVBFVectorBoson --ch 0 --runSysts ${extraOpts} --skipexisting;
	;;


    MERGE )
        gh=${githash}
	./scripts/mergeOutputs.py ${outdir}/${gh}/${EXTRA};
	;;

    PLOT )
	
        json=data/era2017/vbf_samples.json;
	syst_json=data/era2017/vbf_syst_samples.json;
	plotOutDir=${outdir}/${githash}/${EXTRA}/plots/
        kFactors="--procSF MC13TeV_2017_QCDEM_15to20:1.26,MC13TeV_2017_QCDEM_20to30:1.26,MC13TeV_2017_QCDEM_30to50:1.26,MC13TeV_2017_QCDEM_50to80:1.26,MC13TeV_2017_QCDEM_80to120:1.26,MC13TeV_2017_QCDEM_120to170:1.26,MC13TeV_2017_QCDEM_170to300:1.26,MC13TeV_2017_QCDEM_300toInf:1.26,MC13TeV_2017_GJets_HT40to100:1.26,MC13TeV_2017_GJets_HT100to200:1.26,MC13TeV_2017_GJets_HT200to400:1.26,MC13TeV_2017_GJets_HT600toInf:1.26"
	commonOpts="-i ${outdir}/${githash}/${EXTRA} --puNormSF puwgtctr -l ${fulllumi} --saveLog --mcUnc ${lumiUnc} --lumiSpecs LowVPtLowMJJA:${vbflumi},LowVPtHighMJJA:${vbflumi}"
	python scripts/plotter.py ${commonOpts} -j ${json} --only HighMJJ,LowMJJ ${kFactors}
        python scripts/plotter.py ${commonOpts} -j ${json} --only evcount ${kFactors} --saveTeX -o evcout_plotter.root
	python scripts/plotter.py ${commonOpts} -j ${syst_json} ${kFactors} --only HighMJJ,LowMJJ --silent -o syst_plotter.root
        ;;
    
    NLOTFACTORS )
        #transfer factors
        plotOutDir=${outdir}/${githash}/${EXTRA}/plots/
        commonOpts="-p ${plotOutDir}/plotter.root -s ${plotOutDir}/syst_plotter.root -o ${plotOutDir}"
        python test/analysis/computeTransferFactor.py ${commonOpts} --var leadpt     --binList 50,100,150,200,250,300,400,500
        python test/analysis/computeTransferFactor.py ${commonOpts} --var centraleta --binList 0,0.4,0.8,1.2,2,3,5
        python test/analysis/computeTransferFactor.py ${commonOpts} --var forwardeta --binList 0,0.8,1.2,2,3,5
        python test/analysis/computeTransferFactor.py ${commonOpts} --var detajj     --binList 0,0.8,1.6,2.4,3.2,4,5,8
        python test/analysis/computeTransferFactor.py ${commonOpts} --var dphijj     --rebin 2
        python test/analysis/computeTransferFactor.py ${commonOpts} --var mjj        --binList 500,600,700,800,900,1000,1250,1500,2000,3000,4000
        python test/analysis/computeTransferFactor.py ${commonOpts} --var vpt        --binList 75,100,150,200,250,300,350,400,550
        python test/analysis/computeTransferFactor.py ${commonOpts} --var vbfmva     --rebin 5
        ;;

    TRIGEFF )
        #trigger efficiencies
	plotOutDir=${outdir}/${githash}/${EXTRA}/plots/
        python test/analysis/computeVBFTriggerEff.py -p ${plotOutDir}/plotter.root -o ${plotOutDir};
        ;;

    BDTTRANSFORM )
        python test/analysis/VBF_weights/getInverseCDF.py
        ;;

    WWW )

        plotList=(plots)
        gh=${githash}
        for p in ${plotList[@]}; do
            pdir=${outdir}/${gh}/${EXTRA}/${p}
            if [ -d ${pdir} ]; then
                fdir=${wwwdir}/${gh}/${EXTRA}/${p}
	        mkdir -p ${fdir}
	        cp ${pdir}/*.{png,pdf,dat} ${fdir};
	        cp test/index.php ${fdir};
                echo "Check plots in ${fdir} for ${p}"
            fi
        done
	;;
esac
