#!/bin/bash

WHAT=$1; 
if [[ "$1" == "" ]]; then 
    echo "steerTOP16006Analysis.sh <SELTEST/SEL/MERGE/PLOT/BKG/FINALPLOT/COMBPLOT/WWW/CinC/SHAPE/GENSTOP> [plotter.root] [signal] [mass]"; 
    echo "        SELTEST - launches test selection jobs to the batch (inclusive, no syst)"; 
    echo "        SEL     - launches selection jobs to the batch"; 
    echo "        MERGE   - merge the output of the jobs";
    echo "        PLOT    - runs the plotter tool"
    echo "        BKG     - runs backgrond estimation script"
    echo "        FINALPLOT - re-runs the plotter tool taking into account the scale factors for some processes derived in the BKG step"    
    echo "        COMBPLOT - combine the plots of different categories"
    echo "        WWW - move the plots to the afs-based web area"
    echo "        CinC    - runs the cut-in-categories analysis from the current plotter files"
    echo "        SHAPE   - similar for the full shape analysis"
    echo "        GENSTOP - interpolates ttbar mass shapes to generate an ad-hoc stealth stop signal with different masses"
    echo "        LIMITS  - runs limit setting on a defined signal (pass as extra parameters the plotter.root file, the signal name, and the mass)"
    exit 1; 
fi

sigplotter=${2}
if [[ "${sigplotter}" == "" ]]; then sigplotter=final_plotter.root; fi

signal=${3}
if [[ "${signal}" == "" ]]; then signal=tbart; fi

mass=${4}
if [[ "${mass}" == "" ]]; then mass=0; fi

queue=8nh
eosdir=/store/cmst3/user/psilva/LJets2015/8c1e7c9
outdir=~/work/LJets2015-arcrev
wwwdir=~/www/LJets2015-arcrev
lumi=2267.84

RED='\e[31m'
NC='\e[0m'

case $WHAT in
    SELTEST)
	echo -e "[ ${RED} Submitting inclusive selection for the signal regions ${NC} without syst ]"
	python scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue} -o ${outdir}/analysis_mu  --ch 13
	python scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue} -o ${outdir}/analysis_e   --ch 11 
	;;
    SEL )
	echo -e "[ ${RED} Submitting the selection for the signal regions ${NC} ]"
	python scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue} --runSysts -o ${outdir}/analysis_muplus   --ch 13   --charge 1 
	python scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue} --runSysts -o ${outdir}/analysis_muminus  --ch 13   --charge -1
	python scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue} --runSysts -o ${outdir}/analysis_eplus   --ch 11   --charge 1
	python scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue} --runSysts -o ${outdir}/analysis_eminus  --ch 11   --charge -1
	
	echo -e "[ ${RED} Submitting the selection for the control regions ${NC} ]"
	python scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue}            -o ${outdir}/analysis_munoniso --ch 1300
	python scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue}            -o ${outdir}/analysis_enoniso --ch 1100
	python scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue} --runSysts -o ${outdir}/analysis_z --ch 21
	;;
    MERGE )
	a=(muplus muminus munoniso eplus eminus enoniso z)
	for i in ${a[@]}; do
	    echo -e "[ ${RED} Merging ${i} ${NC} ]"
	    ./scripts/mergeOutputs.py ${outdir}/analysis_${i};
	done
	;;
    PLOT )
	a=(muplus  muminus munoniso eplus eminus enoniso)	
	for i in ${a[@]}; do
	    echo -e "[ ${RED} Creating plotter for ${i} ${NC} ]";
	    python scripts/plotter.py -i ${outdir}/analysis_${i}/ --puNormSF puwgtctr  -j data/samples_Run2015.json -l ${lumi} --silent;
	done

	a=(muplus muminus eplus eminus)
	for i in ${a[@]}; do
	    echo -e "[ ${RED} Creating syst plotter for ${i} ${NC} ]";
	    python scripts/plotter.py -i ${outdir}/analysis_${i}/ --puNormSF puwgtctr  -j data/syst_samples_Run2015.json -l ${lumi} -o syst_plotter.root --silent;
	done
	;;
    BKG )
	a=(mu e)
	b=(plus minus)
	for i in ${a[@]}; do
	    for j in ${b[@]}; do
		echo -e "[ ${RED} Background estimation for ${i}${j} ${NC} ]"
		python scripts/runQCDEstimation.py     --iso   ${outdir}/analysis_${i}${j}/plots/plotter.root --noniso ${outdir}/analysis_${i}noniso/plots/plotter.root    --out ${outdir}/analysis_${i}${j}/;
		python scripts/getWJetsScaleFactors.py --shape ${outdir}/analysis_${i}${j}/plots/syst_plotter.root --norm   ${outdir}/analysis_${i}${j}/plots/plotter.root --out ${outdir}/analysis_${i}${j}/;
		echo " "
	    done
	done
	;;
    FINALPLOT )
	a=(munoniso muplus muminus)
	a=(enoniso eplus eminus) # z)
	for i in ${a[@]}; do
	    echo -e "[ ${RED} Creating plotter for ${i} ${NC} ]";
	    python scripts/plotter.py -i ${outdir}/analysis_${i}/ \
		--puNormSF puwgtctr \
		-j data/samples_Run2015.json -l ${lumi} \
		--saveLog -o final_plotter.root;
	    #--procSF "W":${outdir}/analysis_${i}/.wjetsscalefactors.pck;
	done
	;;
    COMBPLOT)
	fs=(e mu plus minus all)
	plots=(nvtx lpt leta ht mttbar csv minmlb drlb metpt)
	#plots=(metpt minmlb)
	cats=(1j 2j 3j 4j 1j0t 1j1t 2j0t 2j1t 2j2t 3j0t 3j1t 3j2t 4j0t 4j1t 4j2t)
	for f in ${fs[@]}; do
	    for p in ${plots[@]}; do
		for c in ${cats[@]}; do	
		    python scripts/combinePlotsForAllCategories.py ${p}_${c} ${f};
		done
	    done
	done
	
	plots=(nbtags)
	for f in ${fs[@]}; do
	    for p in ${plots[@]}; do
		python scripts/combinePlotsForAllCategories.py ${p} ${f};
	    done
	done
	;;
    WWW )
	a=(munoniso enoniso muplus muminus eplus eminus z)
	for i in ${a[@]}; do
	    echo -e "[ ${RED} Moving plots for ${i} ${NC} ]"
	    rm -rf ${wwwdir}/${i}
	    mkdir -p ${wwwdir}/${i};
	    cp ${outdir}/analysis_${i}/plots/*.{png,pdf} ${wwwdir}/${i};
	    cp test/index.php ${wwwdir}/${i};
	done
	mkdir -p ${wwwdir}/comb;
	cp plots/*.{png,pdf} ${wwwdir}/comb;
	cp test/index.php ${wwwdir}/comb;
	;;
    CinC )
	echo -e "[ ${RED} Creating datacards ${NC} ]"
	a=(mu e)
	b=(plus minus)
	finalDataCards=""
	minusDataCards=""
	plusDataCards=""
	for i in ${a[@]}; do

	    chDataCards=""
	    for j in ${b[@]}; do
		python scripts/createDataCard.py --signal ${signal} \
		    -i ${outdir}/analysis_${i}${j}/plots/${sigplotter} --systInput ${outdir}/analysis_${i}${j}/plots/syst_plotter.root \
		    -q ${outdir}/analysis_${i}${j}/.qcdscalefactors.pck \
		    -d nbtags -o ${outdir}/analysis_${i}${j}/datacard;

		    #-w ${outdir}/analysis_${i}${j}/.wjetsscalefactors.pck

		cd ${outdir}/analysis_${i}${j}/datacard;
		combineCards.py ${i}${j}1j=datacard_1j.dat ${i}${j}2j=datacard_2j.dat ${i}${j}3j=datacard_3j.dat ${i}${j}4j=datacard_4j.dat > datacard.dat		
		chDataCards="${i}${j}=../../analysis_${i}${j}/datacard/datacard.dat ${chDataCards}"
		if [ "${j}" = "plus" ]; then
		    plusDataCards="${i}${j}=../../analysis_${i}${j}/datacard/datacard.dat ${plusDataCards}"
		else
		    minusDataCards="${i}${j}=../../analysis_${i}${j}/datacard/datacard.dat ${minusDataCards}"
		fi
		cd -
	    done

	    #channel conbination
	    echo "Combining datacards for ch=${i} from ${chDataCards}"
	    mkdir -p ${outdir}/analysis_${i}/datacard
	    cd ${outdir}/analysis_${i}/datacard
	    combineCards.py ${chDataCards} > datacard.dat
	    finalDataCards="${chDataCards} ${finalDataCards}"
	    cd -
	done
	
	echo "Combining datacards for all channels from ${finalDataCards}"
	mkdir -p ${outdir}/analysis/datacard/
	cd ${outdir}/analysis/datacard/
	combineCards.py ${finalDataCards} > datacard.dat
	cd -

	echo "Combining datacards for all + channels from ${plusDataCards}"
	mkdir -p ${outdir}/analysis_plus/datacard/
	cd ${outdir}/analysis_plus/datacard/
	combineCards.py ${plusDataCards} > datacard.dat
	cd -

	echo "Combining datacards for all - channels from ${minusDataCards}"
	mkdir -p ${outdir}/analysis_minus/datacard/
	cd ${outdir}/analysis_minus/datacard/
	combineCards.py ${minusDataCards} > datacard.dat
	cd -

	a=("mu" "e")
	b=("plus" "minus")
	for i in ${a[@]}; do
	    for j in ${b[@]}; do 
		title="#mu"		
		if [ "${i}${j}" = "muplus" ]; then
		    title="#mu^{+}";
		elif [ "${i}${j}" = "muminus" ]; then
		    title="#mu^{-}";
		elif [ "${i}${j}" = "eplus" ]; then
		    title="e^{+}";
		elif [ "${i}${j}" = "eminus" ]; then
		    title="e^{-}";
		fi
		echo -e "[ ${RED} Running the fit for ${title} ${NC} ]"	
		python scripts/fitCrossSection.py "${title}"=${outdir}/analysis_${i}${j}/datacard/datacard.dat -o ${outdir}/analysis_${i}${j}/datacard --unblind; 
	    done
	    
	    #combined per channel
	    title="#mu"
	    if [ "${i}" = "e" ]; then
		    title="e"
	    fi
	    echo -e "[ ${RED} Running the fit for ${title} ${NC} ]"
            python scripts/fitCrossSection.py "${title}"=${outdir}/analysis_${i}/datacard/datacard.dat -o ${outdir}/analysis_${i}/datacard --unblind;

	done

	#combined per charge
	for j in ${b[@]}; do 
	    title="e^{+}/#mu^{+}"
	    if [ "${j}" = "minus" ]; then
                title="e^{-}/#mu^{-}";
            fi
	    echo -e "[ ${RED} Running the fit for ${title} ${NC} ]"
            python scripts/fitCrossSection.py "${title}"=${outdir}/analysis_${j}/datacard/datacard.dat -o ${outdir}/analysis_${j}/datacard --unblind;
	done

	#final combination
	echo -e "[ ${RED} Running the final ${NC} ]"
        python scripts/fitCrossSection.py "e/#mu"=${outdir}/analysis/datacard/datacard.dat -o ${outdir}/analysis/datacard --unblind;

	;;
    SHAPE )
	echo -e "[ ${RED} Creating shape datacards ${NC} ]"
	a=(mu e)
	b=(plus minus)
	finalDataCards=""
	minusDataCards=""
	plusDataCards=""
	for i in ${a[@]}; do

	    chDataCards=""
	    for j in ${b[@]}; do
		python scripts/createDataCard.py --signal ${signal} \
		    -i ${outdir}/analysis_${i}${j}/plots/${sigplotter} --systInput ${outdir}/analysis_${i}${j}/plots/syst_plotter.root \
		    -o  ${outdir}/analysis_${i}${j}/datacard_shape  -q ${outdir}/analysis_${i}${j}/.qcdscalefactors.pck \
		    -d metpt -c 1j0t,2j0t,3j0t,4j0t;
		python scripts/createDataCard.py --signal ${signal} \
		    -i ${outdir}/analysis_${i}${j}/plots/${sigplotter} --systInput ${outdir}/analysis_${i}${j}/plots/syst_plotter.root \
		    -o  ${outdir}/analysis_${i}${j}/datacard_shape  -q ${outdir}/analysis_${i}${j}/.qcdscalefactors.pck \
		    -d minmlb -c 1j1t,2j1t,2j2t,3j1t,3j2t,4j1t,4j2t;

		cd ${outdir}/analysis_${i}${j}/datacard_shape;
		combineCards.py ${i}${j}1j0t=datacard_1j0t.dat \
		    ${i}${j}1j1t=datacard_1j1t.dat \
		    ${i}${j}2j0t=datacard_2j0t.dat \
		    ${i}${j}2j1t=datacard_2j1t.dat \
		    ${i}${j}2j2t=datacard_2j2t.dat \
		    ${i}${j}3j0t=datacard_3j0t.dat \
		    ${i}${j}3j1t=datacard_3j1t.dat \
		    ${i}${j}3j2t=datacard_3j2t.dat \
		    ${i}${j}4j0t=datacard_4j0t.dat \
		    ${i}${j}4j1t=datacard_4j1t.dat \
		    ${i}${j}4j2t=datacard_4j2t.dat > datacard.dat
		chDataCards="${i}${j}=../../analysis_${i}${j}/datacard_shape/datacard.dat ${chDataCards}"
		if [ "${j}" = "plus" ]; then
		    plusDataCards="${i}${j}=../../analysis_${i}${j}/datacard_shape/datacard.dat ${plusDataCards}"
		else
		    minusDataCards="${i}${j}=../../analysis_${i}${j}/datacard_shape/datacard.dat ${minusDataCards}"
		fi
		cd -
	    done

	    #channel conbination
	    echo "Combining datacards for ch=${i} from ${chDataCards}"
	    mkdir -p ${outdir}/analysis_${i}/datacard_shape
	    cd ${outdir}/analysis_${i}/datacard_shape
	    combineCards.py ${chDataCards} > datacard.dat
	    finalDataCards="${chDataCards} ${finalDataCards}"
	    cd -
	done
	
	echo "Combining datacards for all channels from ${finalDataCards}"
	mkdir -p ${outdir}/analysis/datacard_shape/
	cd ${outdir}/analysis/datacard_shape/
	combineCards.py ${finalDataCards} > datacard.dat
	cd -

	echo "Combining datacards for all + channels from ${plusDataCards}"
	mkdir -p ${outdir}/analysis_plus/datacard_shape/
	cd ${outdir}/analysis_plus/datacard_shape/
	combineCards.py ${plusDataCards} > datacard.dat
	cd -

	echo "Combining datacards for all - channels from ${minusDataCards}"
	mkdir -p ${outdir}/analysis_minus/datacard_shape/
	cd ${outdir}/analysis_minus/datacard_shape/
	combineCards.py ${mimusDataCards} > datacard.dat
	cd -

	a=("mu" "e")
	b=("plus" "minus")
	for i in ${a[@]}; do
	    for j in ${b[@]}; do 
		title="#mu"		
		if [ "${i}${j}" = "muplus" ]; then
		    title="#mu^{+}";
		elif [ "${i}${j}" = "muminus" ]; then
		    title="#mu^{-}";
		elif [ "${i}${j}" = "eplus" ]; then
		    title="e^{+}";
		elif [ "${i}${j}" = "eminus" ]; then
		    title="e^{-}";
		fi
		echo -e "[ ${RED} Running the fit for ${title} ${NC} ]"
		python scripts/fitCrossSection.py "${title}"=${outdir}/analysis_${i}${j}/datacard_shape/datacard.dat -o ${outdir}/analysis_${i}${j}/datacard_shape --unblind; 
	    done
	    
	    #combined per channel
	    title="#mu"
	    if [ "${i}" = "e" ]; then
		    title="e"
	    fi
	    echo -e "[ ${RED} Running the fit for ${title} ${NC} ]"
            python scripts/fitCrossSection.py "${title}"=${outdir}/analysis_${i}/datacard_shape/datacard.dat -o ${outdir}/analysis_${i}/datacard_shape --unblind;

	done

	#combined per charge
	for j in ${b[@]}; do 
	    title="e^{+}/#mu^{+}"
	    if [ "${j}" = "minus" ]; then
                title="e^{-}/#mu^{-}";
            fi
	    continue
	    echo -e "[ ${RED} Running the fit for ${title} ${NC} ]"
            python scripts/fitCrossSection.py "${title}"=${outdir}/analysis_${j}/datacard_shape/datacard.dat -o ${outdir}/analysis_${j}/datacard_shape --unblind;
	done

	#final cobmination
	echo -e "[ ${RED} Running the final ${NC} ]"
        python scripts/fitCrossSection.py "e/#mu"=${outdir}/analysis/datacard_shape/datacard.dat -o ${outdir}/analysis/datacard_shape --unblind;

	;;
    
    GENSTOP)
	a=(mu e)
	b=(plus minus)
	for i in ${a[@]}; do
	    for j in ${b[@]}; do
		echo -e "[ ${RED} Emulating stealth stop by interpolation for ${i}${j} ${NC} ]"
		python scripts/emulateStealthStop.py \
		    -i 172.5:${outdir}/analysis_${i}${j}/plots/final_plotter.root,169.5:${outdir}/analysis_${i}${j}/plots/syst_plotter.root,175.5:${outdir}/analysis_${i}${j}/plots/syst_plotter.root \
		    --scan 169.5,175.5,1.0 \
		    -o ${outdir}/analysis_${i}${j}/plots/stealth_stop_plotter.root;

		hadd -f -k ${outdir}/analysis_${i}${j}/plots/sig_plotter.root ${outdir}/analysis_${i}${j}/plots/stealth_stop_plotter.root ${outdir}/analysis_${i}${j}/plots/final_plotter.root;
		echo " ...full signal plotter can be found in ${outdir}/analysis_${i}${j}/plots/sig_plotter.root"
	    done
	done
	;;
  LIMITS )
	echo -e "[ ${RED} Creating datacards ${NC} ]"
	a=(mu e)
	b=(plus minus)
	finalDataCards=""
	for i in ${a[@]}; do	 
	    for j in ${b[@]}; do
		python scripts/createDataCard.py --signal ${signal} --mass ${mass}\
		    -i ${outdir}/analysis_${i}${j}/plots/${sigplotter} --systInput ${outdir}/analysis_${i}${j}/plots/syst_plotter.root \
		    -q ${outdir}/analysis_${i}${j}/.qcdscalefactors.pck -d nbtags -o ${outdir}/analysis_${i}${j}/datacard_limit_${signal}${mass};
		cd ${outdir}/analysis_${i}${j}/datacard_limit_${signal}${mass};
		combineCards.py ${i}${j}1j=datacard_1j.dat ${i}${j}2j=datacard_2j.dat ${i}${j}3j=datacard_3j.dat ${i}${j}4j=datacard_4j.dat > datacard.dat
		finalDataCards="${i}${j}=../../analysis_${i}${j}/datacard_limit_${signal}${mass}/datacard.dat ${finalDataCards}"
		cd -
	    done
	done
	
	echo "Combining datacards for all channels from ${finalDataCards}"
	mkdir -p ${outdir}/analysis/datacard_limit_${signal}${mass}/
	cd ${outdir}/analysis/datacard_limit_${signal}${mass}/
	combineCards.py ${finalDataCards} > datacard.dat
	cd -
	;;



esac