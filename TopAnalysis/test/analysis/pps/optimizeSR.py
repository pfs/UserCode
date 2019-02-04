import ROOT
import os
import generateWorkspace
import numpy as np
import sys

sig=sys.argv[1]

task_list=[]
for l1pt in [30,40,50]:
    for l2pt in np.arange(20,l1pt,10):
        for bosonpt in [20,40,50,60,70,100]:
            task_list.append((l1pt,l2pt,bosonpt))

h=ROOT.TH1F('optim',';[p_{T}(l_{1}),p_{T}(l_{2}),p_{T}(Z)];95% CL limit on #sigma/#sigma_{th}',len(task_list),0,len(task_list))
for xbin in range(len(task_list)):
    l1pt,l2pt,bosonpt=task_list[xbin]
    
    #generate new workspace
    generateWorkspace.main(['--sig',sig,'-c',"cat==121 && l1pt>{0} && l2pt>{1} && bosonpt>{2}".format(l1pt,l2pt,bosonpt)])

    #run limits
    os.system('text2workspace.py shapes-parametric-140murad.datacard.dat')
    os.system('combine shapes-parametric-140murad.datacard.dat.root -M AsymptoticLimits -t -1')

    #read 50% quantile 
    fIn=ROOT.TFile.Open('higgsCombineTest.AsymptoticLimits.mH120.root')
    limit=fIn.Get('limit')
    limit.GetEntry(2)
    
    h.GetXaxis().SetBinLabel(xbin+1,'[%d,%d,%d]'%(l1pt,l2pt,bosonpt))
    h.SetBinContent(xbin+1,limit.limit)

fOut=ROOT.TFile.Open('optimresults_%s.root'%sig,'RECREATE')
h.Write()
fOut.Close()
