#!/bin/bash

WHAT=$1;

wwwdir=~/www/HIN-17-002

case $WHAT in

    TESTSEL)

	#for sample in MC8.16TeV_TTbar_pPb_Pohweg MC8.16TeV_WJets_pPb MC8TeV_WJets_pp; do 
	for sample in MC8.16TeV_TTbar_pPb_Pohweg; do
	    for wjjOrder in drjj mjj sumpt; do
		for thadOrder in dr dm2tlep dm2pdg; do
	            python prepareWorkspace.py  -d ${sample} --wjjOrder ${wjjOrder} --thadOrder ${thadOrder};
	            mv plots/${sample}/controlplots.root plots/${sample}/controlplots_${wjjOrder}_${thadOrder}.root;
                    python ./checkMatchingEfficiency.py plots/${sample}/controlplots_${wjjOrder}_${thadOrder}.root #> plots/${sample}/algoeff_${wjjOrder}_${thadOrder}.tex 
		done
	    done

	    rm plots/${sample}/*.{png,pdf}

	    python checkMatchingDistributions.py "min#DeltaR(jj)":plots/${sample}/controlplots_drjj_dm2tlep.root \
	        "max#Sigmap_{T}(j)":plots/${sample}/controlplots_sumpt_dm2tlep.root \
	        "min|M(jj)-M(W)|":plots/${sample}/controlplots_mjj_dm2tlep.root
	    mkdir plots/${sample}/w
	    mv *.{png,pdf} plots/${sample}/w


	    python checkMatchingDistributions.py "min|m_{thad}-m_{tlep}|":plots/${sample}/controlplots_drjj_dm2tlep.root \
		"min|m_{thad}_m_{t}|":plots/${sample}/controlplots_sumpt_dm2pdg.root \
		"min#DeltaR(b,W)":plots/${sample}/controlplots_mjj_dr.root
	    mkdir plots/${sample}/top
	    mv *.{png,pdf} plots/${sample}/top
	done
	;;

    SEL)
	common="--wjjOrder drjj --thadOrder dm2tlep" # --etaRestr
        for sample in MC8.16TeV_TTbar_pPb_tighte MC8.16TeV_DY_pPb MC8.16TeV_WJets_pPb MC8.16TeV_TTbar_pPb MC8.16TeV_TTbar_pPb_Pohweg MC8TeV_WJets_pp MC8.16TeV_TTbar_pPb_hypertighte; do    
	    python prepareWorkspace.py  -d ${sample} ${common};
	    python prepareWorkspace.py  -d ${sample} ${common} --etaRestr;
	done
	;;

    SELDATA )
     	common="--wjjOrder drjj --thadOrder dm2tlep"
        for sample in  Data8.16TeV_pPb_nonsubiso_tighte; do #Data8.16TeV_pPb_nonsubiso_hypertighte Data8.16TeV_pPb_nonsubiso Data8.16TeV_pPb
	    python prepareWorkspace.py  -d ${sample} ${common} --jerProf plots/MC8.16TeV_TTbar_pPb/controlplots.root;
	    python prepareWorkspace.py  -d ${sample} ${common} --jerProf plots/MC8.16TeV_TTbar_pPb/controlplots.root --etaRestr;
        done
	;;

    COMPLOTS )

        python compareControlPlots.py --outDir ${wwwdir}/wttmccomp --ref "t#bar{t}" --shape \
            "W":plots/MC8.16TeV_WJets_pPb/controlplots.root \
            "t#bar{t}":plots/MC8.16TeV_TTbar_pPb_tighte/controlplots.root;
        cp ../../index.php ${wwwdir}/wttmccomp;
        
        python compareControlPlots.py --outDir ${wwwdir}/ttmccomp --ref "PYQUEN (official)" --shape \
            "PYQUEN (official)":plots/MC8.16TeV_TTbar_pPb_tighte/controlplots.root \
            "PYQUEN (pp, private)":plots/controlplots_mcfreeze.root \
            "POWHEG (private)":plots/MC8.16TeV_TTbar_pPb_Pohweg/controlplots.root;
        cp ../../index.php ${wwwdir}/ttmccomp;

        python runQCDestimation.py data:plots/Data8.16TeV_pPb_nonsubiso_tighte/controlplots.root \
            ttbar:plots/MC8.16TeV_TTbar_pPb_tighte/controlplots.root \
            wjets:plots/MC8.16TeV_WJets_pPb/controlplots.root \
            dy:plots/MC8.16TeV_DY_pPb/controlplots.root;
        mkdir ${wwwdir}/qcdest;
        mv *fit_*{pdf,png} ${wwwdir}/qcdest;
        cp ../../index.php ${wwwdir}/qcdest;

        python produceXsecBasedControlPlots.py data:plots/Data8.16TeV_pPb_nonsubiso_tighte/controlplots.root \
            ttbar:plots/MC8.16TeV_TTbar_pPb_tighte/controlplots.root \
            wjets:plots/MC8.16TeV_WJets_pPb/controlplots.root \
            dy:plots/MC8.16TeV_DY_pPb/controlplots.root;
        mkdir ${wwwdir}/controlplots;
        mv *control.{pdf,png} ${wwwdir}/controlplots;
        cp ../../index.php ${wwwdir}/controlplots;

        python generateWshapes.py plots/MC8TeV_WJets_pp/controlplots.root plots/MC8.16TeV_WJets_pPb/controlplots.root;
        mkdir  ${wwwdir}/west;
        mv w_*{pdf,png} ${wwwdir}/west;
        cp ../../index.php ${wwwdir}/west;

        python compareControlPlots.py --outDir ${wwwdir}/wmccomp --ref "W pp (8 TeV)" --shape \
            "W pPb (8.16 TeV)":plots/MC8.16TeV_WJets_pPb/controlplots.root \
            "W pp (8 TeV)":plots/MC8TeV_WJets_pp/controlplots.root
        cp ../../index.php ${wwwdir}/wmccomp;

        python compareControlPlots.py --outDir ${wwwdir}/mccomp-eid --ref "pPb" \
            "pPb":plots/MC8.16TeV_TTbar_pPb_tighte/controlplots.root \
            "pPb (|#eta_{e}|<1.5)":plots/MC8.16TeV_TTbar_pPb_tighte_etarestr/controlplots.root \
            "pPb (loose e iso)":plots/MC8.16TeV_TTbar_pPb/controlplots.root \
            "pPb (tight e)":plots/MC8.16TeV_TTbar_pPb_hypertighte/controlplots.root
        cp ../../index.php ${wwwdir}/mccomp-eid;

        python compareControlPlots.py --outDir ${wwwdir}/datacomp --ref "pPb" --data \
            "pPb":plots/Data8.16TeV_pPb_nonsubiso_tighte/controlplots.root \
            "v2":plots/controlplots_datafreeze.root
        cp ../../index.php ${wwwdir}/datacomp;

        python compareControlPlots.py --outDir ${wwwdir}/datacomp-eid --ref "pPb" --data \
            "pPb":plots/Data8.16TeV_pPb_nonsubiso_tighte/controlplots.root \
            "pPb (|#eta_{e}|<1.5)":plots/Data8.16TeV_pPb_nonsubiso_tighte_etarestr/controlplots.root \
            "pPb (loose e iso)":plots/Data8.16TeV_pPb_nonsubiso/controlplots.root \
            "pPb (tight e)":plots/Data8.16TeV_pPb_nonsubiso_hypertighte/controlplots.root;
        cp ../../index.php ${wwwdir}/datacomp-eid;
        ;;
    
    PARAM)

        for sample in MC8.16TeV_TTbar_pPb_tighte; do #MC8.16TeV_TTbar_pPb_hypertighte MC8.16TeV_TTbar_pPb_tighte_etarestr MC8.16TeV_TTbar_pPb_tighte MC8.16TeV_TTbar_pPb MC8.16TeV_TTbar_pPb_Pohweg; do
	    python parameterizeMCShapes.py workspace_${sample}.root;
	    outDir=plots/$sample/workspace/
	    mkdir ${outDir}
            mv workspace_${sample}.dat ${outDir}/
	    b=(`ls *.{png,pdf}`)
	    for i in ${b[@]}; do
		f=${i/./_ttbar.}
		mv -v ${i} ${outDir}/${f};
	    done
	done
        mkdir -p ${wwwdir}/param;
        cp plots/MC8.16TeV_TTbar_pPb_tighte/workspace/*ttbar.{png,pdf} ${wwwdir}/param;
        cp plots/MC8.16TeV_TTbar_pPb_tighte/workspace/*dat ${wwwdir}/param;
        cp ../../index.php ${wwwdir}/param/;
        
	python parameterizeMCShapes.py workspace_MC8.16TeV_WJets_pPb.root pdf_workspace_MC8.16TeV_TTbar_pPb_tighte.root;
	a=(`ls *.{png,pdf}`)
	for i in ${a[@]}; do
	    f=${i/./_w.}
	    mv -v ${i} plots/MC8.16TeV_WJets_pPb/${f};
	done
        cp plots/MC8.16TeV_WJets_pPb/*_w.{png,pdf}  ${wwwdir}/param/;
	;;

    FITS)
        sigRef=MC8.16TeV_TTbar_pPb_tighte #_Pohweg
        dataRef=Data8.16TeV_pPb_nonsubiso_tighte
        #sigRef=MC8.16TeV_TTbar_pPb_tighte_etarestr
        #dataRef=Data8.16TeV_pPb_nonsubiso_tighte_etarestr
        #sigRef=MC8.16TeV_TTbar_pPb_hypertighte
        #dataRef=Data8.16TeV_pPb_nonsubiso_hypertighte
	for f in `seq 0 0`; do
            for wmodel in `seq 0 0`; do
	        python runDataFit.py --fitType ${f} -s pdf_workspace_${sigRef}.root -i workspace_${dataRef}.root -o plots/${dataRef} --wModel ${wmodel};
                impacts=""
                if [[ "${f}" == "0" && "${wmodel}" == "0" ]]; then
                    impacts="--impacts"
                fi
                impacts=""
	        python runDataFit.py --fitType ${f} -o plots/${dataRef}/ --verbose 9 --finalWorkspace finalworkspace_wmodel${wmodel}.root --wModel ${wmodel} ${impacts}
	    done
        done
	;;

    FINALCOMP )
        python doPRLplots.py -i fit_finalworkspace_wmodel0_0.root
        for var in mjj mthad mtlep; do
            for cat in 1l4j2b 1l4j1b1q 1l4j2q; do
                python mergePLRplotCategs.py $var $cat;
            done
        done

        mkdir ${wwwdir}/final;
        mv *{comb,final}*{pdf,png} ${wwwdir}/final;
        cp ../../index.php ${wwwdir}/final;

	;;

esac
