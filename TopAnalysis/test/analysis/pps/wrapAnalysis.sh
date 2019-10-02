#!/bin/bash
HOME=`pwd`
CMSSW=/afs/cern.ch/user/p/psilva/work/CMSSW_9_4_11
PPS=${CMSSW}/src/TopLJets2015/TopAnalysis/test/analysis/pps
samples_json=${PPS}/samples.json,${PPS}/signal_samples.json,${PPS}/signal_samples_postTS2.json
RPout_json=${PPS}/golden_noRP.json
cd ${CMSSW}/src
eval `scram r -sh`
cd ${HOME}
if [ ! -z ${5} ]; then
    extraOpts="--mix ${5}";
fi

echo "python $PPS/runExclusiveAnalysis.py --step $1 --json ${samples_json} --RPout ${RPout_json} -o ./ -i $3 --only $4 ${extraOpts}"
python $PPS/runExclusiveAnalysis.py --step $1 --json ${samples_json} --RPout ${RPout_json} -o ./ -i $3 --only $4 ${extraOpts}

#move to final place
mkdir -p ${2}/Chunks
localout=`basename $4`
localout="${localout%.*}"
a=(`ls Chunks/${localout}.*`)
for i in ${a[@]}; do
    echo $i
    echo "xrdcp -f ${i} root://eoscms/${2/\/eos\/cms/}/Chunks/`basename ${i}`"
    xrdcp -f ${i} root://eoscms/${2/\/eos\/cms/}/Chunks/`basename ${i}`
    rm -v ${i}
done
