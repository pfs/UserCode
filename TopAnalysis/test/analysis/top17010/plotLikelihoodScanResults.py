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
        print mt,gt,nll,url
        
    except Exception as e:
        shScript=url.replace('/fitresults','/runFit')
        shScript=shScript.replace('.root','.sh')
        return shScript

    return mt,gt,nll
        

def profilePOI(data,outdir,axis=0,sigma=5):

    """ profiles in x and y the POI """

    #raw values
    x=data[:,axis]
    xtit='$m_{t}$ [GeV]' if axis==1 else '$\Gamma_{t}$ [GeV]'
    ytit='$m_{t}$ [GeV]' if axis==0 else '$\Gamma_{t}$ [GeV]'

    #def ll_param(x,x0,aL,aR,bL,cL):
    #    bR=2*(aL-aR)*x0+bL
    #    cR=(aL-aR)*(x0**2)+(bL-bR)*x0+cL
    #    y = np.piecewise(x, 
    #                     [x < x0, x >= x0],
    #                     [lambda x:aL*(x**2)+bL*x+cL, lambda x:aR*(x**2)+bR*x+cR])
    #    return y

    ictr=0
    xvals=[]
    llvals=[]
    for xi in np.unique(x):
        rdata=data[data[:,axis]==xi]
        y=rdata[:,0 if axis==1 else 1]
        bounds = [min(y),max(y)]
        z=rdata[:,2]
        
        #filter for outliers
        #medianz=np.median(z)
        #filtIdx=np.where(abs(z-medianz) >=10000)
        #y=y[filtIdx]
        #z=z[filtIdx]


        #check we still have enough points
        if len(y)<2: continue

        #make sure it's sorted correcly
        sortIdx=np.argsort(y)        

        #interpolate to generate equally spaced grid
        #and apply a gaussian filter
        y_unif = np.arange(bounds[0],bounds[1],0.001*(bounds[1]-bounds[0]))
        z_spline = interp1d(y[sortIdx],z[sortIdx],kind='cubic',fill_value='extrapolate')
        z_spline_val=z_spline(y_unif)
        z_filt = filters.gaussian_filter1d(z_spline_val,sigma=sigma)

        min_idx=np.argmin(z_filt)
        xvals.append(xi)
        llvals.append(z_filt[min_idx])
        
        #custom fit
        #try:
        #    popt, pcov = curve_fit(ll_param, y,z)
        #except:
        #    continue

        #4th pol order fit
        #pcoeff=np.polyfit(y,z,4)
        #p=np.poly1d(pcoeff)        
        #crit_points = [px for px in p.deriv().r if px.imag == 0 and bounds[0] < px.real < bounds[1]]
        #print xi,crit_points

        plt.clf()
        fig, ax = plt.subplots()
        yp = np.linspace(bounds[0],bounds[1], 100)
        plt.plot(y,  2*z,                   'o',  label='scan points')
        #plt.plot(yp, p(yp),               '-',  label='interpolation')
        #plt.plot(yp, ll_param(yp, *popt), '--', label='interpolation')
        plt.plot(y_unif, 2*z_filt,          '--', label='interpolated')
        plt.xlabel(xtit)
        plt.ylabel(r'$-2\log(\lambda)$')
        #plt.ylim(0.,20.0)
        ax.text(0,1.02,'CMS preliminary', transform=ax.transAxes, fontsize=16)
        ax.text(1.0,1.02,r'%s=%3.2f 34.5 fb$^{-1}$ (13 TeV)'%(ytit,xi), transform=ax.transAxes,horizontalalignment='right',fontsize=14)
        ax.legend(framealpha=0.0, fontsize=14, loc='upper left', numpoints=1)        
        #_,_,axymin,axymax = plt.axis()
        #if axymax>axymin+20:
        #    ax.set_ylim([axymin,axymin+20])
        plt.savefig(os.path.join(outdir,'nllprofile_%d_%d.png'%(axis,ictr)))
        plt.savefig(os.path.join(outdir,'nllprofile_%d_%d.pdf'%(axis,ictr)))
        ictr+=1


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

    min_idx=np.argmin(llvals_filt)
    x0=xvals_unif[min_idx]
    ll_0=llvals_filt[min_idx]
    llvals_filt=(llvals_filt-ll_0)*2
    llvals=(llvals-ll_0)*2

    dnll=9999.
    xLow=bounds[0]
    for i in range(0,min_idx):
        idnll=abs(llvals_filt[i]-1)
        if idnll>dnll: continue
        dnll=idnll
        xLow=xvals_unif[i]
    if dnll>0.01: xLow=bounds[0]

    dnll=9999.
    xUp=bounds[1]
    for i in range(min_idx+1,len(llvals_filt)):
        idnll=abs(llvals_filt[i]-1)
        if idnll>dnll: continue
        dnll=idnll
        xUp=xvals_unif[i]
    if dnll>0.01: xUp=bounds[1]


    
    fig, ax = plt.subplots()
    plt.plot(xvals,      llvals,      'o',  label='scan points')
    plt.plot(xvals_unif, llvals_filt, '--', label='interpolated')
    plt.xlabel(ytit)    
    plt.ylabel(r'$-2\Delta\log(\lambda)$')
    ax.text(0,1.02,'CMS preliminary', transform=ax.transAxes, fontsize=16)
    ax.text(1.0,1.02,r'34.5 fb$^{-1}$ (13 TeV)', transform=ax.transAxes,horizontalalignment='right',fontsize=14)
    ax.text(0.95,0.94,r'%s=$%3.2f^{+%3.2f}_{-%3.2f}$ GeV'%(ytit,x0,xUp-x0,x0-xLow), transform=ax.transAxes,horizontalalignment='right',fontsize=12)
    ax.legend(framealpha=0.0, fontsize=14, loc='upper left', numpoints=1)        
    plt.savefig(os.path.join(outdir,'finalnllprofile_%d.png'%(axis)))
    plt.savefig(os.path.join(outdir,'finalnllprofile_%d.pdf'%(axis)))



def doContour(data,outdir,              
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
    bestFitIdx = np.where(zi==minz)
    bestX=xi[ bestFitIdx[1][0] ]
    bestY=yi[ bestFitIdx[0][0] ]
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
    plt.plot([bestX], [bestY], '+', mew=4, markersize=12, color='k',label='Best fit')

    #add theory prediction
    theory=getTheoryPrediction()
    plt.plot(theory[0],theory[1],'-',color='lightgray',label='Theory NLO')

    plt.xlabel('$m_{t}$ [GeV]', fontsize=16)
    plt.ylabel('$\Gamma_{t}$ [GeV]', fontsize=16)
    plt.ylim(0.7,4.0)
    plt.xlim(169.5,175.5)
    ax.text(0,1.02,'CMS preliminary', transform=ax.transAxes, fontsize=16)
    ax.text(1.0,1.02,r'34.5 fb$^{-1}$ (13 TeV)', transform=ax.transAxes,horizontalalignment='right',fontsize=14)
    ax.legend(framealpha=0.0, fontsize=14, loc='upper left', numpoints=1)

    for ext in ['png','pdf']:
        plt.savefig(os.path.join(outdir,'nllcontour.%s'%ext))


def main():

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--in',          
                      dest='input',       
                      help='input directory [%default]',  
                      default='/eos/cms/store/cmst3/group/top/TOP17010/final/0c522df/fits/em_inc',
                      type='string')
    parser.add_option('-o', '--out',          
                      dest='outdir',
                      help='output directory [%default]',  
                      default='fit_results/em_inc',
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
    (opt, args) = parser.parse_args()


    os.system('rm -rf {0} && mkdir -p {0}'.format(opt.outdir))

    #build nll scan
    fitres=[]
    toSub=[]
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

    #plot the contour interpolating the available points
    try:
        fitres=np.array(fitres)
        profilePOI(fitres,outdir=opt.outdir,axis=0,sigma=opt.filterSigma)
        profilePOI(fitres,outdir=opt.outdir,axis=1,sigma=opt.filterSigma)
        doContour(fitres,outdir=opt.outdir)
    except Exception as e:
        print '<'*50
        print e
        print '<'*50

if __name__ == "__main__":
    sys.exit(main())
