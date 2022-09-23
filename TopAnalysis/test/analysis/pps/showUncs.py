import ROOT
import os
import sys
import numpy as np
import argparse

systList=[
    'bkgShapeEMUp',         'bkgShapeEMDown',
    'bkgShapeSingleDiffUp', 'bkgShapeSingleDiffDown',
    'sigShapeEMUp',         'sigShapeEMDown',
    'sigCalibUp',           'sigCalibDown',
    'sigPzModelUp',         'sigPzModelDown',
    'sigPPSEffUp',          'sigPPSEffDown',
]


tableList=[x.replace('Up','') for x in systList[::2]]
tableList += ['mcstats']
           

def printUncertainties(flist,m=1000,bosonTag='zmm'):

    uncVals={}

    for f in flist:
        r=ROOT.TFile.Open(f)
        print(f)
        baseH={'bkg'       : r.Get('bkg_{}'.format(bosonTag)),
               'fidsig'    : r.Get('fidsig_{}_m{}'.format(bosonTag,m)),
               'outfidsig' : r.Get('outfidsig_{}_m{}'.format(bosonTag,m))}
        for k in r.GetListOfKeys():
            kname=k.GetName()
            
            if 'data' in kname : continue
            if 'sig_' in kname:
                sigm=int(kname.split('_')[2].replace('m',''))
                if m!=sigm : continue
                
            isSigFid    = True if 'fidsig' in kname and not 'outfidsig' in kname else False
            isSigOutFid = True if 'outfidsig' in kname else False            
            proc='bkg'
            if isSigFid:
                proc='fidsig'
            elif isSigOutFid:
                proc='outfidsig'
            
            if not proc in uncVals: uncVals[proc]={}
            
            sname       = kname.split('_')[-1]
            if not sname in systList: sname='mcstats'
            sname=sname.replace('Up','')
            sname=sname.replace('Down','')
            if not sname in uncVals[proc]: uncVals[proc][sname]=[]           

            h=k.ReadObj()
            if sname!='mcstats':
                h.Divide(baseH[proc])
                
            vals=[]
            xbin1000=h.GetXaxis().FindBin(1000)
            for xbin in xrange(1,h.GetNbinsX()):
                val=h.GetBinContent(xbin)
                #if proc=='fidsig':
                #    if abs(xbin-xbin1000)>4: continue
                if sname=='mcstats':
                    if val<=0.01 :continue
                    vals.append( h.GetBinError(xbin)/val )
                else:
                    vals.append( h.GetBinContent(xbin)-1 )
            uncVals[proc][sname].append( vals )

        r.Close()
    
    for s in tableList:
        print s,
        for p in ['bkg','fidsig','outfidsig']:
            if not s in uncVals[p]: continue
            
            for i in range(len(flist)):

                q=np.percentile( uncVals[p][s][i], [50,10,90])
                #q=np.percentile( uncVals[p][s][i], [50,16,84])
                
                if i==0: print '\n\t & %10s'%p,
                print ' & %3.3f~~$]%3.3f,%3.3f[$'%(q[0],q[1],q[2]),
                #print ' & %3.3f~~$]%3.3f,%3.3f[$'%(np.mean(np.abs(uncVals[p][s][i])),q[1],q[2]),
                print '\\\\',
        print '\n'

def main(args):

    parser = argparse.ArgumentParser(description='usage: %prog [options]')
    parser.add_argument('-i', '--input',
                        dest='input',   
                        default='ppvx_2017_unblind_multi_1exc/',
                        help='input directory with the optim directories [default: %default]')
    parser.add_argument('--optimPoints',
                        dest='optimPoints',
                        default='0,25,50,75',
                        help='optim points to compare [default %default]')
    parser.add_argument('-b','--boson',
                        dest='boson',
                        default=169,
                        help='boson code [default %default]')
    parser.add_argument('-m',
                        dest='massPoint',
                        default=1000,
                        help='signal point [%default]')
    opt=parser.parse_args(args)


    flist=[]
    for pt in opt.optimPoints.split(','):
        optim_url='{}/optim_{}/shapes_{}.root'.format(opt.input,pt,opt.boson)
        flist.append(optim_url)

    bosonTag='zmm'
    if opt.boson==121 : bosonTag='zee'
    if opt.boson==22  : bosonTag='g'

    printUncertainties(flist,opt.massPoint,bosonTag)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

