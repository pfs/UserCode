#!/bin/bash

WHAT=$1; 
EXTRA=$2
QCD=""

if [ "$#" -lt 1 ]; then 
    echo "steerVBFVectorBoson.sh <SEL/MERGE/PLOT/WWW> [extra]";
    echo "        SEL          - launches selection jobs to the batch, output will contain summary trees and control plots"; 
    echo "        MERGE        - merge output (if given \"extra\" is appended to the directory)"
    echo "        PLOT         - make plots (if given \"extra\" is appended to the directory)"
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


RED='\e[31m'
NC='\e[0m'
case $WHAT in

    TESTSEL )
       
        #tag=Data13TeV_2017C_SinglePhoton
        tag=MC13TeV_2017_EWKAJJ
        input=${eosdir}/${tag}/Chunk_0_ext0.root
        output=${tag}.root 

	python scripts/runLocalAnalysis.py \
            -i ${input} -o ${output} --tag ${tag} \
            --njobs 1 -q local --genWeights genweights_${githash}.root \
            --era era2017 -m VBFVectorBoson::RunVBFVectorBoson --ch 0 --runSysts --debug;

        #--debug --mvatree \
        ;;

    SEL )
	if [[ -z ${EXTRA} ]]; then
	    extraOpts=" --mvatree" #" --SRfake"
        fi
	python scripts/runLocalAnalysis.py \
	    -i ${eosdir} \
            -o ${outdir}/${githash}/${EXTRA} \
            --farmappendix ${githash} \
            -q ${queue} --genWeights genweights_${githash}.root \
            --era era2017 -m VBFVectorBoson::RunVBFVectorBoson --ch 0 --runSysts ${extraOpts};
	;;


    SEL2018 )
	python scripts/runLocalAnalysis.py -i /store/cmst3/group/top/RunIISpring18/f3174df\
            -o ${outdir}/2018/raw \
            -q ${queue} \
            --era era2017 -m VBFVectorBoson::RunVBFVectorBoson --ch 0 --runSysts;

	python scripts/runLocalAnalysis.py -i /store/cmst3/group/top/RunIISpring18/hem1516_failure\
            -o ${outdir}/2018/hem1516_failure \
            -q ${queue} \
            --era era2017 -m VBFVectorBoson::RunVBFVectorBoson --ch 0 --runSysts;
	;;

    SELJETHT )
	json=data/era2017/JetHT.json;
	extraOpts=""
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
	    -i ${eosdir} \
            -o ${outdir}/${githash}/${EXTRA}${QCD} \
            --farmappendix ${githash}${EXTRA}${QCD} \
            -q ${queue} --genWeights genweights_${githash}.root \
            --era era2017 -m VBFVectorBoson::RunVBFVectorBoson --ch 0 --only ${json} --runSysts ${extraOpts} --skipexisting;
	;;

    SELWEIGHTED )
	python scripts/runLocalAnalysis.py -i ${eosdir} \
            -o ${outdir}/weighted --flag 1 \
            -q ${queue} \
            --era era2017 -m VBFVectorBoson::RunVBFVectorBoson --ch 0 --runSysts;
        ;;

    MERGE )
        gh=${githash}
        if [[ "${EXTRA}" = *"2018"* ]]; then
            gh=${githash2018}
        fi
	# if [[ "${EXTRA}" = *"JetHT"* ]]; then
	#     gh=${githashJetHT}
	# fi
	./scripts/mergeOutputs.py ${outdir}/${gh}/${EXTRA};
	;;

    PLOT2018)
         commonOpts="-i test/analysis/VBFVectorBoson/2018/hem1516_failure/ --puNormSF puwgtctr -l 1.0"
         python scripts/plotter.py ${commonOpts} -j test/analysis/VBFVectorBoson/2018/hem1516_failure/samples.json --only V1J;
         fdir=${wwwdir}/hem1516_failure;
	 mkdir -p ${fdir}
	 cp test/analysis/VBFVectorBoson/2018/hem1516_failure/plots/*.{png,pdf,dat} ${fdir};
	 cp test/index.php ${fdir};
         echo "Check plots in ${fdir}"
        ;;

    PLOT )

        json=data/era2017/vbf_samples.json;
	syst_json=data/era2017/vbf_syst_samples.json;
        lumi=${fulllumi}        
        gh=${githash}/
	plotOutDir=${outdir}/${githash}/${EXTRA}/plots/
        if [[ "${EXTRA}" = *"2018"* ]]; then
            json=data/era2018/vbf_samples.json;
            lumi=${fulllumi2018}
            vbflumi=${lumi}
            gh=${githash2018}
        fi
        kFactors="--procSF MC13TeV_2017_QCDEM_15to20:1.26,MC13TeV_2017_QCDEM_20to30:1.26,MC13TeV_2017_QCDEM_30to50:1.26,MC13TeV_2017_QCDEM_50to80:1.26,MC13TeV_2017_QCDEM_80to120:1.26,MC13TeV_2017_QCDEM_120to170:1.26,MC13TeV_2017_QCDEM_170to300:1.26,MC13TeV_2017_QCDEM_300toInf:1.26,MC13TeV_2017_GJets_HT40to100:1.26,MC13TeV_2017_GJets_HT100to200:1.26,MC13TeV_2017_GJets_HT200to400:1.26,MC13TeV_2017_GJets_HT600toInf:1.26"
	commonOpts="-i ${outdir}/${gh}/${EXTRA} --puNormSF puwgtctr -l ${lumi} --saveLog --mcUnc ${lumiUnc} --lumiSpecs HighMJJA:${vbflumi},HighMJJLPA:${vbflumi},OfflineVBFA:${fulllumi}"
	commonOpts="-i ${outdir}/${gh}/${EXTRA} --puNormSF puwgtctr --saveLog -l ${lumi} --mcUnc ${lumiUnc} --lumiSpecs HighMJJA:${vbflumi},HighMJJLPA:${vbflumi},LowMJJA:${fulllumi},HighMJJMM:${fulllumi},LowMJJMM:${fulllumi} -O ${plotOutDir}"
	python scripts/plotter.py ${commonOpts} -j ${json} --only HighMJJ,LowMJJ,V1J ${kFactors}
	python scripts/plotter.py ${commonOpts} -j ${syst_json} ${kFactors} --only HighMJJ,LowMJJ,V1J --silent -o syst_plotter.root

        #trigger efficiencies
        python test/analysis/computeVBFTriggerEff.py -p ${plotOutDir}/plotter.root -o ${plotOutDir};

        #transfer factors
	python test/analysis/computeTransferFactor.py -p ${plotOutDir}/plotter.root -s ${plotOutDir}/syst_plotter.root -o ${plotOutDir} --var vbffisher --binList -2,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,2,3;

        # if [[ "${EXTRA}" != *"2018"* ]]; then
        #     python scripts/plotter.py ${commonOpts}  -j data/era2017/vbf_signal_samples.json --only HighPtA_ -O ${outdir}/${githash}/${EXTRA}/plots_signal/ --noStack;
        #     continue
        #     #python scripts/plotter.py ${commonOpts}  -j data/era2017/gjets_samples.json --only HighPtA_mjj -O ${outdir}/${githash}/${EXTRA}/plots_gjets/ --noStack;
        #     trigOpts="-i ${outdir}/${gh}/${EXTRA} --puNormSF puwgtctr -l ${vbflumi}  --saveLog --mcUnc ${lumiUnc} ${kFactors} --only HighPtOfflineVBFA,HighPtVBFA"
        #     #python scripts/plotter.py ${trigOpts} -j data/era2017/vbf_samples_2017F.json -O ${outdir}/${githash}/${EXTRA}/plots_trigger;
        # fi
        ;;
       
    RATIO )
        #python test/analysis/computeVBFRatios.py -t \
        #    -i ${outdir}/${githash}/${EXTRA}/plots_trigger/plotter.root,${outdir}/${githash2018}/raw2018/plots/plotter.root \
        #    --titles "2017","2018" \
        #    -o ${outdir}/${githash}/${EXTRA}/plots_trigger/trigger_ratio_plotter.root

        python test/analysis/computeVBFRatios.py \
            -i ${outdir}/${githash}/${EXTRA}/plots/plotter.root,${outdir}/${githash2018}/raw2018/plots/plotter.root \
            --titles "2017","2018" \
            -o ${outdir}/${githash}/${EXTRA}/plots/ratio_plotter.root
        
	;;

    BDTTRANSFORM )
        python test/analysis/VBF_weights/getInverseCDF.py
        ;;

    WWW )

        plotList=(plots) #plots_signal) # plots plots_trigger plots_signal plots_gjets)
        gh=${githash}
        if [[ "${EXTRA}" = *"2018"* ]]; then
            gh=${githash2018}
        fi

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
