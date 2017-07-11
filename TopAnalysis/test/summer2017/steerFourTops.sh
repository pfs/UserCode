#!/bin/bash

WHAT=$1; 
if [ "$#" -ne 1 ]; then 
    echo "steerFourTops.sh <SEL/MERGE/PLOT/WWW>";
    echo "        SEL          - launches selection jobs to the batch, output will contain summary trees and control plots"; 
    echo "        MERGE        - merge output"
    echo "        PLOT         - make plots"
    echo "        WWW          - move plots to web-based are"
    exit 1; 
fi

export LSB_JOB_REPORT_MAIL=N
#queue=2nd
queue=condor
githash=b312177
lumi=35922
lumiUnc=0.025
whoami=`whoami`
myletter=${whoami:0:1}
eosdir=/store/cmst3/group/top/ReReco2016/${githash}
outdir=${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/summer2017/FourTops
wwwdir=~/www/FourTops


RED='\e[31m'
NC='\e[0m'
case $WHAT in

    SEL )
        queue=condor
        python scripts/runLocalAnalysis.py -i ${eosdir} \
            --only test/summer2017/4tops_samples.json --exactonly \
            -q ${queue} -o ${outdir} --jobflavour longlunch\
            --era era2016 -m FourTopsAnalyzer::RunFourTopsAnalyzer --ch 0 --runSysts;
    ;;

    SELALL )
        python scripts/runLocalAnalysis.py -i ${eosdir} \
            --only data/era2016/samples-fourtops.json --exactonly \
            -q condor -o ${outdir} --jobflavour workday\
            --era era2016 -m FourTopsAnalyzer::RunFourTopsAnalyzer --ch 0 --runSysts;
    ;;

    MERGE )
    ./scripts/mergeOutputs.py ${outdir};
    ;;

    PLOT )
    commonOpts="-i ${outdir} --puNormSF puwgtctr -j test/summer2017/4tops_samples.json -l ${lumi}  --saveLog --mcUnc ${lumiUnc} --noStack"
    python scripts/plotter.py ${commonOpts}; 
    ;;

    PLOTALL )
    commonOpts="-i ${outdir} --puNormSF puwgtctr -j data/era2016/samples-fourtops.json -l ${lumi}  --saveLog --mcUnc ${lumiUnc}"
    python scripts/plotter.py ${commonOpts}; 
    ;;

    WWW )
    mkdir -p ${wwwdir}/sel
    cp ${outdir}/plots/*.{png,pdf} ${wwwdir}/sel
    cp test/index.php ${wwwdir}/sel
    ;;
esac
