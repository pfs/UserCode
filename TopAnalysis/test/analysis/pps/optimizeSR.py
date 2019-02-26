import ROOT
import os
import generateWorkspace
import numpy as np
import sys
import pickle

sig=sys.argv[1]
workingpoints={}

task_list=[]
for l1pt in [30,40,50]:
    for l2pt in np.arange(20,l1pt,10):
        for bosonpt in [20,40,50,60,70,80]:
            task_list.append((l1pt,l2pt,bosonpt))

print '%d working points being scanned for %s'%(len(task_list),sig)

h=ROOT.TH1F('optim',';[p_{T}(l_{1}),p_{T}(l_{2}),p_{T}(Z)];95% CL limit on #sigma/#sigma_{th}',len(task_list),0,len(task_list))
hsig={1:h.Clone('optim_sig1'),10:h.Clone('optim_sig10')}
for xbin in range(len(task_list)):
    l1pt,l2pt,bosonpt=task_list[xbin]
    anaCuts='[%d,%d,%d]'%(l1pt,l2pt,bosonpt)
    iwp=[]
    
    #generate new workspace
    generateWorkspace.main(['-d','--sig',sig,'-c',"(cat==121||cat==169) && l1pt>{0} && l2pt>{1} && bosonpt>{2}".format(l1pt,l2pt,bosonpt)])

    allDC=' '.join(['alpha{0}=analysis/shapes-parametric-{0}murad.datacard.dat'.format(xangle) for xangle in [120,130,140,150]])
    combDC='analysis/shapes-parametric.datacard.dat'
    os.system('combineCards.py {0} > {1}'.format(allDC,combDC))
    os.system('text2workspace.py {0}'.format(combDC))

    #run limits and read 50% quantile
    os.system('combine {0}.root -M AsymptoticLimits -t -1'.format(combDC))
    fIn=ROOT.TFile.Open('higgsCombineTest.AsymptoticLimits.mH120.root')
    limit=fIn.Get('limit')
    for i in [0,1,2,3,4]:
        limit.GetEntry(i)
        iwp.append(limit.limit)
    r95=iwp[2]
    h.GetXaxis().SetBinLabel(xbin+1,anaCuts)
    h.SetBinContent(xbin+1,r95)
    fIn.Close()

    #run significance
    for exp in hsig:
        os.system('combine {0}.root -M Significance -t -1 --expectSignal={1}'.format(combDC,exp))
        fIn=ROOT.TFile.Open('higgsCombineTest.Significance.mH120.root')
        limit=fIn.Get('limit')
        limit.GetEntry(0)
        iwp.append(limit.limit)
        fIn.Close()
        hsig[exp].GetXaxis().SetBinLabel(xbin+1,anaCuts)
        hsig[exp].SetBinContent(xbin+1,iwp[-1])


    workingpoints[(l1pt,l2pt,bosonpt)]=iwp

fOut=ROOT.TFile.Open('optimresults_%s.root'%sig,'RECREATE')
h.Write()
for exp in hsig: hsig[exp].Write()
fOut.Close()

with open('optimresults_%s.pck'%sig,'w') as cache:
    pickle.dump(workingpoints,cache,pickle.HIGHEST_PROTOCOL) 
