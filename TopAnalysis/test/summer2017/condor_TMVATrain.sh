#!/bin/bash
WORKDIR=`pwd`
cd /afs/cern.ch/user/v/vwachira/CMSSW_8_0_28
eval `scram r -sh`
cd ${WORKDIR}
root -l -q /afs/cern.ch/user/v/vwachira/CMSSW_8_0_28/src/TopLJets2015/TopAnalysis/test/summer2017/FourTops/FourTopsTMVA_606.cpp
cp ${WORKDIR}/FourTops/weights/FourTopsML_BDT.weights.xml /afs/cern.ch/user/v/vwachira/CMSSW_8_0_28/src/TopLJets2015/TopAnalysis/test/summer2017/FourTopsML_BDT_606.weights.xml
mv ${WORKDIR}/FourTops/* /afs/cern.ch/user/v/vwachira/CMSSW_8_0_28/src/TopLJets2015/TopAnalysis/test/summer2017/FourTops/