import ROOT
import os
import generateBinnedWorkspace
import numpy as np
import sys
import pickle

def runWorkspacePacked(args):

    """create the workspace for a set of cuts/categories"""

    itask,m,l1pt,l2pt,bosonpt,categs,baseDir=args
    for xangle in generateBinnedWorkspace.VALIDLHCXANGLES:
        generateBinnedWorkspace.main(['--xangle', '%d'%xangle,
                                      '--sig',    'Z_m_X_{0}_xangle_{1}_2017_preTS2_opt_v1_simu_reco.root'.format(m,xangle),
                                      '--presel', 'cat==169 && l1pt>{0} && l2pt>{1} && bosonpt>{2}'.format(l1pt,l2pt,bosonpt),
                                      '--categs', categs,
                                      '-o',       '%s/optim_task%d'%(baseDir,itask)
                                  ])

#generate a grid search in the parameters/categories of interest
m=int(sys.argv[1])
baseOptim='analysis'
if len(sys.argv)>2:
    baseOptim=sys.argv[2]
baseDir='%s/stat_m%d'%(baseOptim,m)

itask=0
task_list=[]

if baseOptim=='analysis':
    l1pt=30
    l2pt=20
    for bosonpt in [20,30,40,50,60]:
        for categs in ['nvtx>0','nvtx<15,nvtx>=15','nvtx<20,nvtx>=20','nvtx<25,nvtx>=25','nvtx<30,nvtx>=30']:
            itask=len(task_list)+1
            task_list.append((itask,m,l1pt,l2pt,bosonpt,categs,baseDir))

elif baseOptim=='pucatdef':
    l1pt=30
    l2pt=20
    bosonpt=40
    #formVar='(PFPzSumHF-1400)/384+nvtx'
    formVar='TMath::Sqrt(((PFPzSumHF-1400)/384)**2+nvtx**2)'
    for categs in ['nvtx<20,nvtx>=20',
                   '{0}<10,{0}>=20'.format(formVar),
                   '{0}<20,{0}>=30'.format(formVar),
                   '{0}<30,{0}>=40'.format(formVar),
                   '{0}<30,{0}>=50'.format(formVar)
                   ]:

        print categs
        itask=len(task_list)+1
        task_list.append((itask,m,l1pt,l2pt,bosonpt,categs,baseDir))

os.system('mkdir -p %s'%baseDir)
with open('%s/scanpoints.pck'%baseDir,'w') as cache:
    pickle.dump(task_list,cache,pickle.HIGHEST_PROTOCOL) 

import multiprocessing as MP
pool = MP.Pool(8)
print '%d working points being scanned for m=%d'%(len(task_list),m)
pool.map( runWorkspacePacked, task_list)
