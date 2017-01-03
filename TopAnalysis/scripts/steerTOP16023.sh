#!/bin/bash

WHAT=$1; 
DIST=$2
REBINFACTOR=$3
UNBLIND=$4
if [ "$#" -lt 1 ]; then 
    echo "steerTOP16023.sh <SEL/MERGE/BKG/PLOT/WWW/PREPAREFIT/FIT/DOFITLOOP/SHOWFIT> [DIST [REBINFACTOR [UNBLIND]]]";
    echo "        SEL          - selects data and MC";
    echo "        MERGE        - merge the output of the jobs";
    echo "        BKG          - runs the background estimation from sidebands";
    echo "        PLOT         - runs the plotter tool on the selection";
    echo "        WWW          - moves the plots to an afs-web-based area";    
    echo "        PREPAREFIT   - create datacards for the fit requires extras arguments DIST (distribution name: mjj by default) REBINFACTOR (bins to merge)"
    echo "        FIT          - run the cross section fit (may need a special CMSSW release to use combine)"
    echo "                       requires DIST & REBINFACTOR + if 1 is passed as well it will unblind"
    echo "        SHOWFIT      - show summary plots - requires DIST & REBINFACTOR"
    echo "        FITLOOP      - run all fits"
    exit 1; 
fi

#suppress batch notifications by mail
export LSB_JOB_REPORT_MAIL=N

queue=8nh
sourcedir=/store/cmst3/group/hintt/LJets5TeV/
outdir=~/work/LJets-5TeV
wwwdir=~/www/LJets-5TeV
lumi=27.4
COMBINERELEASE=~/scratch0/CMSSW_7_4_7/src/


RED='\e[31m'
NC='\e[0m'
case $WHAT in
    PREPARE )
	echo -e "[ ${RED} computing total number of eff. events available for analysis ${NC} ]"
	python scripts/produceNormalizationCache.py -i ${sourcedir}                -o data/era5TeV/genweights.root  --HiForest;
	echo -e "[ ${RED} projecting the b-tagging efficiency ${NC} ]"
	python scripts/saveExpectedBtagEff.py       -i ${sourcedir}/MCTTNominal_v2 -o data/era5TeV/expTageff.root   --HiForest; 
	;;
    SEL )
	echo -e "[ ${RED} Sending out jobs to batch ${NC} ]"

	commonOpts="--era era5TeV -m TOP-16-023::RunTop16023"

	#muon channel
	python scripts/runLocalAnalysis.py -i ${sourcedir} -q ${queue} -o ${outdir}/analysis_mu       ${commonOpts} --ch 13   --only MC --runSysts;
	python scripts/runLocalAnalysis.py -i ${sourcedir} -q ${queue} -o ${outdir}/analysis_munoniso ${commonOpts} --ch 1300 --only MC;
	python scripts/runLocalAnalysis.py -i ${sourcedir} -q ${queue} -o ${outdir}/analysis_mu       ${commonOpts} --ch 13   --only FilteredSingleMuHighPt;
	python scripts/runLocalAnalysis.py -i ${sourcedir} -q ${queue} -o ${outdir}/analysis_munoniso ${commonOpts} --ch 1300 --only FilteredSingleMuHighPt;

	#electron channel
        python scripts/runLocalAnalysis.py -i ${sourcedir} -q ${queue} -o ${outdir}/analysis_e           ${commonOpts} --ch 11   --only MC --runSysts;
	python scripts/runLocalAnalysis.py -i ${sourcedir} -q ${queue} -o ${outdir}/analysis_enoniso     ${commonOpts} --ch 1100 --only MC;
	python scripts/runLocalAnalysis.py -i ${sourcedir}    -q ${queue} -o ${outdir}/analysis_e/       ${commonOpts} --ch 11   --only FilteredHighPtPhoton30AndZ;
	python scripts/runLocalAnalysis.py -i ${sourcedir}    -q ${queue} -o ${outdir}/analysis_enoniso/ ${commonOpts} --ch 1100 --only FilteredHighPtPhoton30AndZ;
	;;
    MERGE )
	echo -e "[ ${RED} Merging job output ${NC} ]"
	a=(mu munoniso e enoniso)
	for i in ${a[@]}; do
	    ./scripts/mergeOutputs.py ${outdir}/analysis_${i};
	done
	;;
    BKG )
	echo -e "[ ${RED} Running QCD estimation from non-isolated side-band ${NC} ]"
	a=(mu munoniso e enoniso)
	for i in ${a[@]}; do
	    python scripts/plotter.py -i ${outdir}/analysis_${i}  -j data/era5TeV/samples.json  --com "5.02 TeV"  -l ${lumi} --silent;
	done
	for ch in e mu; do
	    python scripts/runQCDEstimation.py \
		--iso    ${outdir}/analysis_${ch}/plots/plotter.root \
		--noniso ${outdir}/analysis_${ch}noniso/plots/plotter.root \
		--out    ${outdir}/analysis_${ch}/ \
		--sels  ,0b,1b,2b;
	done
	;;
    PLOT )
	echo -e "[ ${RED} Running plotter ${NC} ]"
	a=(mu munoniso e enoniso)
	for i in ${a[@]}; do
	    python scripts/plotter.py -i ${outdir}/analysis_${i}  -j data/era5TeV/Wsamples.json     --com "5.02 TeV" -l ${lumi} --saveLog --noStack;	
	    mkdir ~/${outdir}/analysis_${i}/wplots;
	    mv ~/${outdir}/analysis_${i}/plots/* ~/${outdir}/analysis_${i}/wplots/;
	    python scripts/plotter.py -i ${outdir}/analysis_${i}  -j data/era5TeV/samples.json      --com "5.02 TeV" -l ${lumi} --saveLog;	
	    python scripts/plotter.py -i ${outdir}/analysis_${i}  -j data/era5TeV/syst_samples.json --com "5.02 TeV" -l ${lumi} -o syst_plotter.root --silent;	
	done
	;;

    WWW )
	echo -e "[ ${RED} Moving plots to ${outdir} ${NC} ]"
	a=(mu munoniso e enoniso)
	for i in ${a[@]}; do
	    mkdir -p ${wwwdir}/analysis_${i}
	    cp ${outdir}/analysis_${i}/plots/*.{png,pdf} ${wwwdir}/analysis_${i}
	    cp test/index.php ${wwwdir}/analysis_${i}
	done
	;;
    PREPAREFIT )
	echo -e "[ ${RED} Creating datacards for ${DIST} with rebin=${REBINFACTOR} ${NC} ]"

	FITTAG=${DIST}_${REBINFACTOR}
	for ch in e mu; do	
	    python scripts/createDataCard.py \
		-i ${outdir}/analysis_${ch}/plots/plotter.root \
		--systInput ${outdir}/analysis_${ch}/plots/syst_plotter.root \
		-q ${outdir}/analysis_${ch}/.qcdscalefactors.pck \
		-o ${outdir}/analysis_${ch}/datacard_${FITTAG} \
		--specs TOP-16-015 \
		--signal tbart \
		-d ${DIST} \
		-c 0b,1b,2b \
		--addBinByBin 0.3 \
		--rebin ${REBINFACTOR};
	
	    a=(0b 1b 2b)
	    for i in ${a[@]}; do	
		python scripts/projectShapeUncs.py ${outdir}/analysis_${ch}/datacard_${FITTAG}/shapes_${i}.root btag,othertag,jes,jer;
		python scripts/projectShapeUncs.py ${outdir}/analysis_${ch}/datacard_${FITTAG}/shapes_${i}.root ttPartonShower,Hadronizer,ttFactScale,ttRenScale,ttCombScale;
		python scripts/projectShapeUncs.py ${outdir}/analysis_${ch}/datacard_${FITTAG}/shapes_${i}.root wFactScale,wRenScale,wCombScale W;
	    done
	    mkdir -p ${wwwdir}/shapes_${ch}_${FITTAG}
	    mv ${outdir}/analysis_${ch}/datacard_${FITTAG}/*.{png,pdf} ${wwwdir}/shapes_${ch}_${FITTAG};
	    cp test/index.php ${wwwdir}/shapes_${ch}_${FITTAG};
	done
	;;
    FIT )
	
	FITTAG=${DIST}_${REBINFACTOR}

	#make sure combine is installed
	echo "Sourcing cmsenv for ${COMBINERELEASE}"
	cd ${COMBINERELEASE}/
        eval `scramv1 r -sh`
        cd -
	echo -e "[ ${RED} $CMSSW_BASE will be used - make sure combine is compatible and is installed ${NC} ]"


	for ch in e mu; do
	    echo -e "\t combining datacards/workspace for ${ch} channel";
	    cd ${outdir}/analysis_${ch}/datacard_${FITTAG};
	    combineCards.py ${ch}0b=datacard_0b.dat ${ch}1b=datacard_1b.dat ${ch}2b=datacard_2b.dat > datacard.dat;
	    python ${COMBINERELEASE}/HiggsAnalysis/CombinedLimit/test/systematicsAnalyzer.py datacard.dat --all -m 172.5 -f html > systs.html
	    text2workspace.py datacard.dat -m 0 -o workspace.root
            cd -
	done

	#prepare the combined output
	echo -e "\t combining datacards/workspace for l=combined channel"
	mkdir -p ${outdir}/analysis_l/datacard_${FITTAG};
	cd ${outdir}/analysis_l/datacard_${FITTAG}
	combineCards.py \
	    mu0b=${outdir}/analysis_mu/datacard_${FITTAG}/datacard_0b.dat \
	    mu1b=${outdir}/analysis_mu/datacard_${FITTAG}/datacard_1b.dat \
	    mu2b=${outdir}/analysis_mu/datacard_${FITTAG}/datacard_2b.dat \
	    e0b=${outdir}/analysis_e/datacard_${FITTAG}/datacard_0b.dat \
	    e1b=${outdir}/analysis_e/datacard_${FITTAG}/datacard_1b.dat \
	    e2b=${outdir}/analysis_e/datacard_${FITTAG}/datacard_2b.dat > datacard.dat
	python ${COMBINERELEASE}/HiggsAnalysis/CombinedLimit/test/systematicsAnalyzer.py datacard.dat --all -m 172.5 -f html > systs.html
	text2workspace.py datacard.dat -m 0 -o workspace.root
        cd -

	#now run the fits
	for ch in e mu l; do
	    echo -e "[ ${RED} Running fits for ${ch} channel ${NC} ]"

	    cd ${outdir}/analysis_${ch}/datacard_${FITTAG};

            #expected
	    commonOpts="-t -1 --expectSignal=1 --setPhysicsModelParameterRanges btagRate=-2,2:r=0,2 -m 0 --saveWorkspace";
	    combine workspace.root -M MaxLikelihoodFit ${commonOpts};
	    mv mlfit.root mlfit_exp.root
	    combine workspace.root -M MaxLikelihoodFit ${commonOpts} -S 0;
	    mv mlfit.root mlfit_exp_stat.root
	    combine workspace.root -M MultiDimFit ${commonOpts} --algo=grid --points=100;
	    mv higgsCombineTest.MultiDimFit.mH0.root exp_plr_scan_r.root
	    combine exp_plr_scan_r.root --snapshotName MultiDimFit -M MultiDimFit ${commonOpts} --algo=grid --points=100 -S 0;
	    mv higgsCombineTest.MultiDimFit.mH0.root exp_plr_scan_stat_r.root

	    commonOpts="-t -1 --redefineSignalPOIs btagRate -P btagRate --expectSignal=1 --algo=grid --points=100 --setPhysicsModelParameterRanges btagRate=-2,2:r=0,2 -m 0";	
	    combine workspace.root -M MultiDimFit ${commonOpts};
	    mv higgsCombineTest.MultiDimFit.mH0.root exp_plr_scan_btag.root
	    combine exp_plr_scan_btag.root --snapshotName MultiDimFit -M MultiDimFit ${commonOpts} -S 0;
	    mv higgsCombineTest.MultiDimFit.mH0.root exp_plr_scan_stat_btag.root
	
            combine workspace.root -M MultiDimFit --algo=grid --points=2500 -m 0 -t -1 \
		--redefineSignalPOIs r,btagRate -P r -P btagRate --setPhysicsModelParameterRanges btagRate=-2,2:r=0,2 \
		--expectSignal=1 --saveWorkspace;
	    mv higgsCombineTest.MultiDimFit.mH0.root exp_plr_scan_rvsbtag.root;

    	    #observed...
	    if [ "${UNBLIND}" == "1" ]; then
		echo -e "[ ${RED} will unblind the results now ${NC}]";
		
		commonOpts="--setPhysicsModelParameterRanges btagRate=-2,2:r=0,2 -m 0 --saveWorkspace";
		combine workspace.root -M MaxLikelihoodFit ${commonOpts};
		mv mlfit.root mlfit_obs.root
		combine workspace.root -M MaxLikelihoodFit ${commonOpts} -S 0;
		mv mlfit.root mlfit_obs_stat.root

		combine workspace.root -M MultiDimFit --redefineSignalPOIs btagRate -P btagRate --algo=grid --points=100 ${commonOpts};
		mv higgsCombineTest.MultiDimFit.mH0.root obs_plr_scan_btag.root
		combine obs_plr_scan_btag.root --snapshotName MultiDimFit -M MultiDimFit --redefineSignalPOIs btagRate -P btagRate --algo=grid --points=100 ${commonOpts} -S 0;
		mv higgsCombineTest.MultiDimFit.mH0.root obs_plr_scan_stat_btag.root
		
		combine workspace.root -M MultiDimFit -P r --algo=grid --points=100 ${commonOpts};
		mv higgsCombineTest.MultiDimFit.mH0.root obs_plr_scan_r.root
		combine obs_plr_scan_r.root --snapshotName MultiDimFit -M MultiDimFit -P r --algo=grid --points=100 ${commonOpts} -S 0;
		mv higgsCombineTest.MultiDimFit.mH0.root obs_plr_scan_stat_r.root
		
		combine workspace.root -M MultiDimFit --algo=grid --points=2500 --redefineSignalPOIs r,btagRate -P r -P btagRate  ${commonOpts};
		mv higgsCombineTest.MultiDimFit.mH0.root obs_plr_scan_rvsbtag.root;
		
		#commonOpts="--setPhysicsModelParameterRanges r=0,2 -m 0";
		#combineTool.py -M Impacts -d workspace.root --doInitialFit ${commonOpts} --minimizerTolerance 0.001
		#combineTool.py -M Impacts -d workspace.root --doFits       ${commonOpts} --minimizerTolerance 0.001
		#combineTool.py -M Impacts -d workspace.root ${commonOpts} -o impacts.json
		#plotImpacts.py -i impacts.json -o impacts
	    fi
	    
	    cd -
	done
	;;
    SHOWFIT )

	FITTAG=${DIST}_${REBINFACTOR}
	echo -e "[ ${RED} Fit plots will be made available in ${wwwdir}/fits_${FITTAG} ${NC}]";

	for ch in e mu l; do
	    cardsDir=${outdir}/analysis_${ch}/datacard_${FITTAG}/;
	    python scripts/fitSummaryPlots.py ${ch}=${cardsDir}/datacard.dat --POIs r,btagRate --label "27.4 pb^{-1} (5.02 TeV)" -o ${cardsDir};
	    mkdir -p ${wwwdir}/fits_${ch}_${FITTAG};
	    mv ${cardsDir}/*.{png,pdf,C} ${wwwdir}/fits_${ch}_${FITTAG}/;
	    cp test/index.php ${wwwdir}/fits_${ch}_${FITTAG};
	done
	;;
    FITLOOP )
	
	DISTS=("mjj" "drjj" "rankedmjj" "rankedq70mjj")
	for dist in ${DISTS[@]}; do 
	    REBIN=(4 20)
	    if [ "$dist" == "drjj" ]; then
		REBIN=(2 12)
	    fi
	    STEPS=(PREPAREFIT FIT SHOWFIT)	    
	    for rebin in ${REBIN[@]}; do 
		for step in ${STEPS[@]}; do
		    sh ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/scripts/steerTOP16023.sh ${step} ${dist} ${rebin} 1
		done
	    done
	done

	;;
esac