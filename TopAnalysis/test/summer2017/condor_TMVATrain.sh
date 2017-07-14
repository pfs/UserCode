#!/bin/bash
cd /afs/cern.ch/user/v/vwachira/CMSSW_8_0_28
eval `scram r -sh`
cd src/TopLJets2015/TopAnalysis/test/summer2017
root -l -q FourTopsTMVA_606.cpp