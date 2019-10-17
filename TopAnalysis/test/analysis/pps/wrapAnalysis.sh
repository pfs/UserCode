#!/bin/bash
HOME=`pwd`
CMSSW=${1}
PPS=${CMSSW}/src/TopLJets2015/TopAnalysis/test/analysis/pps
samples_json=${PPS}/samples.json,${PPS}/signal_samples.json,${PPS}/signal_samples_postTS2.json
RPout_json=${PPS}/golden_noRP.json
cd ${CMSSW}/src
eval `scram r -sh`
cd ${HOME}
if [ ! -z ${6} ]; then
    extraOpts="--mix ${6}";
fi

echo "python $PPS/runExclusiveAnalysis.py --step $2 --json ${samples_json} --RPout ${RPout_json} -o ./ -i $4 --only $5 ${extraOpts}"
python $PPS/runExclusiveAnalysis.py --step $2 --json ${samples_json} --RPout ${RPout_json} -o ./ -i $4 --only $5 ${extraOpts}

#move to final place
mkdir -p ${3}/Chunks
localout=`basename $5`
localout="${localout%.*}"
a=(`ls Chunks/${localout}.*`)
for i in ${a[@]}; do
    echo $i
    echo "xrdcp -f ${i} root://eoscms/${3/\/eos\/cms/}/Chunks/`basename ${i}`"
    xrdcp -f ${i} root://eoscms/${3/\/eos\/cms/}/Chunks/`basename ${i}`
    rm -v ${i}
done
