#!/bin/bash

#configure environment
WORKDIR=${1}/src
echo "Setting environment from $WORKDIR"
cd $WORKDIR
eval `scram r -sh`
cd -

#run with the arguments passed
shift
echo $@
$*
