#!/bin/bash
HOME=`pwd`
CMSSW=${1}
PPS=${CMSSW}/src/TopLJets2015/TopAnalysis/test/analysis/pps
cd ${CMSSW}/src
eval `scram r -sh`
cd ${HOME}
python $PPS/optimizeSR.py ${2} ${3} ${4}
