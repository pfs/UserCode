import ROOT
import os
import sys
import optparse
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.ndimage import filters
import numpy as np
import re
import pickle

def getTheoryPrediction(x=np.arange(169,176,0.1)):
    
    """top width versus top mass at NLO """

    def nlo_gt(mt):
        GF=1.1663787e-5
        mW=80.385
        aS=0.118
        pi=ROOT.TMath.Pi()
        
        gt  = GF*(mt**3)/(8*pi*ROOT.TMath.Sqrt(2.))
        gt *= (1-(mW/mt)**2)**2
        gt *= (1+2*(mW/mt)**2)
        gt *= (1-(2*aS/(3*pi))*(2*(pi**2)/3-5./2.))

        return gt

    return x,np.array( [nlo_gt(m) for m in x] )

def findLikelihoodMinimum(x,y):

    """scans by brute force the likelihood values or attempts for a second-order polynomial type of fit"""

    #minimum by brute force
    min_idx = np.argmin(y)
    x0      = x[min_idx]
    nll0    = y[min_idx]
    
    #find minimum from the local extremes of a polynomial fit around the min.
    idx2param     = np.where((y < nll0+3))[0]
    y_pol         = np.poly1d( np.polyfit(x[idx2param], y[idx2param], 2) )
    dy_pol        = np.polyder(y_pol)  #first derivative
    pol_extremes  = np.roots(dy_pol)
    x0_pol        = pol_extremes[ np.abs(pol_extremes - x0).argmin() ]
    nll0_pol      = y_pol(x0_pol)

    #use the scan as it is and find the 68%CL scanning around the minimum
    dnll_low,dnll_up=abs(y[0]-nll0-1.),abs(y[-1]-nll0-1.)
    x_low,x_up=x[0],x[-1]
    for i in range(0,min_idx):
        idnll=abs(y[i]-nll0-1.)
        if idnll>dnll_low: continue
        dnll_low=idnll
        x_low=x[i]
    nll_low=nll0+dnll_low
    for i in range(min_idx+1,len(y)):
        idnll=abs(y[i]-nll0-1.)
        if idnll>dnll_up: continue
        dnll_up=idnll
        x_up=x[i]
    nll_up=nll0+dnll_up

    #use the polynomial to find the up/down values of the CI
    baseCoeff     = y_pol.c
    baseCoeff[-1] = baseCoeff[-1]-nll0_pol-1
    yp1_pol       = np.poly1d(baseCoeff)
    x_up_pol,x_low_pol = np.roots( yp1_pol )

    toReturn={
        'brute-force' : [(x_low,nll_low),       (x0,nll0),         (x_up,nll_up)],
        'polyfit'     : [(x_low_pol,nll0_pol+1),(x0_pol,nll0_pol), (x_up_pol,nll0_pol+1)]
    }
    
    return toReturn
                
    


def getScanPoint(inDir,fitTag):

    """read the fit result and return the likelihood value together with the corresponding mtop,width values"""


    #decode mass and width
    tag=os.path.basename(inDir)
    mt,gt=172.5,1.31
    if 'scenario' in tag:
        flag=int(re.findall(r'\d+',tag)[0])
        mask=int('0xffff',16)
        gt = (flag&mask)*0.01+0.7
        mt = (((flag>>16)&mask))*0.25+169        
    else:
        return mt,gt,None

    #read nll from fit
    nll=None
    try:
        url=os.path.join(inDir,'fitresults%s.root'%fitTag)
        if not os.path.isfile(url):
            raise ValueError('%s is missing'%url)
        inF=ROOT.TFile.Open(url)
        if inF.IsZombie() or inF.TestBit(ROOT.TFile.kRecovered):
            raise ValueError('%s probably corrupted'%url)
        tree=inF.Get('fitresults')
        tree.GetEntry(0)
        nll=tree.nllvalfull
        
    except Exception as e:
        shScript=url.replace('/fitresults','/runFit')
        shScript=shScript.replace('.root','.sh')
        return shScript

    return mt,gt,nll
        

def profilePOI(data,outdir,axis=0,sigma=5):

    """ profiles in x and y the POI """

    #raw values
    x=data[:,axis]
    xvar='$m_{t}$' if axis==1 else '$\Gamma_{t}$'
    xtit='%s [GeV]'%xvar
    yvar='$m_{t}$' if axis==0 else '$\Gamma_{t}$'
    ytit='%s [GeV]'%yvar

    ictr=0
    xvals=[]
    llvals=[]
    for xi in np.unique(x):

        rdata=data[data[:,axis]==xi]
        y=rdata[:,0 if axis==1 else 1]
        bounds = [min(y),max(y)]
        z=rdata[:,2]
        
        #check we still have enough points
        if len(y)<2: continue

        #make sure it's sorted correcly
        sortIdx = np.argsort(y)        

        #interpolate to generate equally spaced grid and apply a gaussian filter
        y_unif       = np.arange(bounds[0],bounds[1],0.001*(bounds[1]-bounds[0]))
        z_spline     = interp1d(y[sortIdx],z[sortIdx],kind='cubic',fill_value='extrapolate')
        z_spline_val = z_spline(y_unif)
        z_filt       = filters.gaussian_filter1d(z_spline_val,sigma=sigma)
            
        #minimize likelihood
        minResults = findLikelihoodMinimum(y_unif,2*z_filt)
        bestFitX=minResults['polyfit'][1][0]
        dX_up=minResults['brute-force'][2][0]-minResults['brute-force'][1][0]
        dX_lo=minResults['brute-force'][0][0]-minResults['brute-force'][1][0]
        dX_up=max(dX_up,minResults['polyfit'][2][0]-bestFitX)
        dX_lo=min(dX_lo,minResults['polyfit'][0][0]-bestFitX)
        minLL=minResults['polyfit'][1][1]
        xvals.append(xi)
        llvals.append(minLL)


        #show the likelihood
        plt.clf()
        fig, ax = plt.subplots()
        yp = np.linspace(bounds[0],bounds[1], 100)
        plt.plot(y,      2*z-minLL,      'o',  label='scan points')
        plt.plot(y_unif, 2*z_filt-minLL, '--', label='likelihood')
        plt.xlabel(xtit)
        plt.ylabel(r'$-2\log(\lambda)$')
        plt.ylim(0.,20.0)
        ax.text(0,1.02,'CMS preliminary', transform=ax.transAxes, fontsize=16)
        ax.text(1.0,1.02,r'%s=%3.2f 35.6 fb$^{-1}$ (13 TeV)'%(ytit,xi), transform=ax.transAxes,horizontalalignment='right',fontsize=14)
        ax.text(0.5, 0.94,r'Best fit: %s=$%3.2f^{+%3.2f}_{%3.2f}$ GeV'%(xvar,bestFitX,dX_up,dX_lo), transform=ax.transAxes,horizontalalignment='center',fontsize=14)

        figName='nllprofile_%d_%d'%(axis,ictr)
        for ext in ['.png','.pdf']:
            plt.savefig(os.path.join(outdir,figName+ext))
        print 'Saved likelihood',figName,(xvar,bestFitX,dX_up,dX_lo)
        ictr+=1


    #final profile: notice here we don't need to multiply by a factor of 2
    #the previous step ensures that the raw nll values have already been multiplied by 2 when minimizing
    plt.clf()
    xvals=np.array(xvals)
    llvals=np.array(llvals)
    bounds = [min(xvals),max(xvals)]
    #interpolate to generate equally spaced grid
    #and apply a gaussian filter
    xvals_unif = np.arange(bounds[0],bounds[1],0.0001*(bounds[1]-bounds[0]))
    llvals_spline = interp1d(xvals,llvals,kind='cubic',fill_value='extrapolate')
    llvals_spline_val=llvals_spline(xvals_unif)
    llvals_filt = filters.gaussian_filter1d(llvals_spline_val,sigma=5)

    minResults=findLikelihoodMinimum(xvals_unif,llvals_filt)
    bestFitX=minResults['polyfit'][1][0]
    dX_up=minResults['brute-force'][2][0]-minResults['brute-force'][1][0]
    dX_lo=minResults['brute-force'][0][0]-minResults['brute-force'][1][0]
    dX_up=max(dX_up,minResults['polyfit'][2][0]-bestFitX)
    dX_lo=min(dX_lo,minResults['polyfit'][0][0]-bestFitX)
    minLL=minResults['polyfit'][1][1]

    fig, ax = plt.subplots()
    plt.plot(xvals,      llvals-minLL,      'o',  label='scan points')
    plt.plot(xvals_unif, llvals_filt-minLL, '--', label='interpolated')
    plt.xlabel(ytit)    
    plt.ylabel(r'$-2\Delta\log(\lambda)$')
    plt.ylim(0.,20.0)
    ax.text(0,1.02,'CMS preliminary', transform=ax.transAxes, fontsize=16)
    ax.text(1.0,1.02,r'35.6 fb$^{-1}$ (13 TeV)', transform=ax.transAxes,horizontalalignment='right',fontsize=14)
    ax.text(0.5,0.94,r'Best fit: %s=$%3.2f^{+%3.2f}_{%3.2f}$ GeV'%(ytit,bestFitX,dX_up,dX_lo), transform=ax.transAxes,horizontalalignment='center',fontsize=12)
    
    #ax.legend(framealpha=0.0, fontsize=14, loc='upper left', numpoints=1)    

    figName='finalnllprofile_%d'%axis
    for ext in ['.png','.pdf']:
        plt.savefig(os.path.join(outdir,figName+ext))

    return (bestFitX,dX_up,dX_lo)
    


def doContour(data,
              bestFitX,
              bestFitY,
              outdir,
              myFinalstate,
              mjjcat,
              method='linear',
              levels=[2.30,4.61,9.21],
              levelLabels=['68.3%','90%','99%'],
              linestyles=['-','-','-']):

    """ interpolates the grid to obtain the likelihood contour 
    2 parameter fit levels (see PDG Statistics Table 38.2) """

    #raw values
    x=data[:,0]
    y=data[:,1]
    z=data[:,2]

    #filter for outliers
    #medianz=np.median(z)
    #filtIdx=np.where( abs(z-medianz)>10000)
    #x=x[filtIdx]
    #y=y[filtIdx]
    #z=z[filtIdx]


    #interpolate and find minimum
    xi = np.linspace(169.5, 175.5,100)
    yi = np.linspace(0.7,4.0,100)
    zi = griddata((x, y), z, (xi[None,:], yi[:,None]), method=method)
    minz=zi.min()
    zi=(zi-minz)*2
    z=(z-minz)*2

    fig, ax = plt.subplots()
    cntr=ax.contour(xi, yi, zi, levels=levels, linewidths=1.0, linestyles=linestyles, colors='k')
    cntr.collections[0].set_label('$-2\Delta ln(\lambda)$')

    #add contour level names
    fmt = {}
    for l, s in zip(levels,levelLabels):
        fmt[l] = s
    ax.clabel(cntr, cntr.levels[:], inline=True, fmt=fmt, fontsize=10)

    #add best-fit point
    #plt.plot([bestFitX[0]], [bestFitY[0]], '+', mew=4, markersize=12, color='k',label='Best fit')
    ax.errorbar([bestFitX[0]], [bestFitY[0]], xerr=bestFitX[1:2], yerr=bestFitY[1:2], 
                markersize=12, color='k', fmt='--o',label='Best fit (1D)')

    #0=mumu/1=ee/>2=nothing
    fsTag=""
    if myFinalstate==0 :
        fsTag=" #mu#mu"
    elif myFinalstate==1 :
        fsTag=" ee"
    else :
        print "warning: wrong sfTag input"

    #add information for all  categories
    #0=lowVbf/1=Hightpt/>2=200to500
    #ax.text(0,1.02 ,  "2016 samples with #mu final state", transform=ax.transAxes, fontsize=16)
    if mjjcat == 0 :
        ax.text(0,1.02 ,  "Low" + fsTag, transform=ax.transAxes, fontsize=16)
    elif mjjcat == 1:
        ax.text(0,1.02,"HighVPt  Pt_{z}>200 200<M_{jj}<500GeV" + fsTag,transform=ax.transAxes, fontsize=16)
    else:
           ax.text(0,1.02,"HighVPt  M_{jj}>500GeV" + fsTag,transform=ax.transAxes, fontsize=16) 
        
    #add theory prediction
    theory=getTheoryPrediction()
    plt.plot(theory[0],theory[1],'-',color='lightgray',label='Theory NLO')

    plt.xlabel('$m_{t}$ [GeV]', fontsize=16)
    plt.ylabel('$\Gamma_{t}$ [GeV]', fontsize=16)
    plt.ylim(0.7,4.0)
    plt.xlim(169.5,175.5)
    ax.text(0,1.02,'CMS preliminary', transform=ax.transAxes, fontsize=16)
    ax.text(1.0,1.02,r'35.6 fb$^{-1}$ (13 TeV)', transform=ax.transAxes,horizontalalignment='right',fontsize=14)
    ax.legend(framealpha=0.0, fontsize=14, loc='upper left', numpoints=1)

    for ext in ['png','pdf']:
        plt.savefig(os.path.join(outdir,'nllcontour.%s'%ext))


def main():

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--in',          
                      dest='input',       
                      help='input directory [%default]',  
                      default='/eos/cms/store/cmst3/group/top/TOP17010/final_method1a/0c522df/fits/final',
                      type='string')
    parser.add_option('-o', '--out',          
                      dest='outdir',
                      help='output directory [%default]',  
                      default='fit_results/final',
                      type='string')
    parser.add_option('-t', '--tag',          
                      dest='fitTag',
                      help='fit tag [%default]',  
                      default='_tbart',
                      type='string')
    parser.add_option('--sigma',          
                      dest='filterSigma',
                      help='fiter sigma [%default]',  
                      default=5,
                      type=int)
    parser.add_option('-r', '--recover',
                      dest='recover',
                      help='recover [%default]',  
                      default=False,
                      action='store_true')

    parser.add_option('-m', '--mjjcat',
                      dest='mjjcat',
                      help='enum for mjj category, 0=LVBF, 1=High, 2=Middle [%default]',  
                      default=0,
                      type=int)
    
    parser.add_option( '--fsTag',
                      dest='fsTag',
                      help='enum for final state, 0=mumu, 1=ee, 2=wrong [%default]',  
                      default=0,
                      type=int)
    
    (opt, args) = parser.parse_args()


    os.system('rm -rf {0} && mkdir -p {0}'.format(opt.outdir))

    #build nll scan if the input is a directory
    if os.path.isdir(opt.input):       
        fitres=[]
        toSub=[]
        print 'Scanning available results'
        for f in os.listdir(opt.input):
            scanRes=getScanPoint(inDir=os.path.join(opt.input,f),fitTag=opt.fitTag)
            if isinstance(scanRes,basestring):
                toSub.append(scanRes)
                continue
            elif not scanRes[-1]: 
                continue
            fitres.append( scanRes )

        # treat missing jobs
        if len(toSub)>0:
            condor_file='condor_recover_%s_%s.sub'%( os.path.basename(opt.input),opt.fitTag )
            print 'I have %d missing/corrupted jobs to submit on condor - sub file @ %s'%(len(toSub),condor_file)        
            cmssw=os.environ['CMSSW_BASE']
            with open(condor_file,'w') as condor:            
                condor.write('executable  = %s/src/TopLJets2015/TopAnalysis/test/analysis/top17010/runFitWrapper.sh\n'%cmssw)
                condor.write('output      = datacard_condor.out\n')
                condor.write('error       = datacard_condor.err\n')
                condor.write('log         = datacard_condor.log\n')           
                condor.write('+JobFlavour = "workday"\n')
                for f in toSub:
                    condor.write('arguments  = %s\n'%f)
                    condor.write('queue 1\n')
            if opt.recover:
                print 'Submitting condor file'
                os.system('condor_submit %s'%condor_file)
                os.system('cp -v {0} {1}/{0}'.format(condor_file,opt.outdir))

        #dump current results to a pickle file and update the input parameter to point to file
        opt.input=os.path.join(opt.input,'nllscan_%s.pck'%opt.fitTag)
        with open(opt.input,'w') as cache:
            pickle.dump(fitres,cache,pickle.HIGHEST_PROTOCOL)
        print 'Scan results have been stored in',opt.input
        print 'Next run can parse them directly with -i',opt.input

    #open pickle file with the results
    with open(opt.input,'r') as cache:
        fitres=pickle.load(cache)
    print len(fitres),'fit results are available, using to find best fit points'

    #plot the contour interpolating the available points
    try:
        fitres=np.array(fitres)
        bestFitX=profilePOI(fitres,outdir=opt.outdir,axis=0,sigma=opt.filterSigma)
        bestFitY=profilePOI(fitres,outdir=opt.outdir,axis=1,sigma=opt.filterSigma)
        doContour(fitres,bestFitX,bestFitY,outdir=opt.outdir,opt.fsTag, opt.mjjcat)
    except Exception as e:
        print '<'*50
        print e
        print '<'*50

if __name__ == "__main__":
    sys.exit(main())
