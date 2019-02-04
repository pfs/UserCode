#!/bin/bash

for l1pt in 30 40 50; do
    for l2pt in 30 40 50; do
        for bosonpt in 40 50 60 70; do
            python test/analysis/pps/generateWorkspace.py -c "cat==121 && l1pt>${l1pt} && l2pt>${l2pt} && bosonpt>${bosonpt}";
            text2workspace.py shapes-parametric-140murad.datacard.dat
            combine shapes-parametric-140murad.datacard.dat.root -M AsymptoticLimits --expectSignal=1 -t -1;
            mv higgsCombineTest.AsymptoticLimits.mH120.root AsymptoticLimits_${l1pt}_${l2pt}_${bosonpt}.root;
        done
    done
done