import numpy as np
from scipy.ndimage import gaussian_filter1d
import ROOT
import sys

class GFSmoother():

    def __init__(self,h,ntoys=100,sigma=2,truncate=4,fromCDF=True):

        self.x=self.generateToys(h,ntoys)
        self.y=self.smoothToys(self.x,sigma,truncate,fromCDF)
        self.smooth=self.profiledSmoothedToys(h,fromCDF)

    def generateToys(self,h,ntoys):

        """generates n toys taking into account stat unc in bins"""
        
        x=[]
        
        for i in range(ntoys):
            ix=[]

            #loop over the histogram bins
            for xbin in range(h.GetNbinsX()):
                val=h.GetBinContent(xbin+1)
                unc=h.GetBinError(xbin+1)

                #sample the mean taking into account stat fluctuation
                mean=np.random.normal(val,unc)

                #sample this bin
                rand_val=np.random.normal(mean,unc)
                ix.append(rand_val)
            x.append(ix)

        return np.array(x)
    
    def smoothToys(self,x,sigma,truncate,fromCDF):

        """ loops over toys and applies a gaussian filter smooth
            if fromCDF is true, the smoothing is applied on the normalized cumulative
            distribution function """

        nbins=len(x[0])

        #transform to normalized CDF if needed (prepend 0, force 1 at the edge)
        rawx = x
        if fromCDF: 
            rawx = x.cumsum(axis=1) / x.sum(axis=1)[:,None]
            rawx[:,nbins-1] =  1.
            rawx=np.c_[np.zeros(len(rawx)),rawx]

        #apply smoothing
        ntoys = len(rawx)
        y = np.array( [gaussian_filter1d(rawx[i], sigma=sigma, truncate=truncate,mode='nearest') for i in range(ntoys) ] )

        #force the first and last bins of the CDF smoothing
        if fromCDF : y[:,0],y[:,-1]=0.,1.

        return y

    def profiledSmoothedToys(self,h,fromCDF):
        """ gets a copy of the original histogram where the bins filled with the median obtained
            with the results of smoothing the toys and the uncertainties correspond to the the 68% quantiles """

        nbins=h.GetNbinsX()
        norm=h.Integral()
        
        #transform back to PDF if needed (ignore extra aux. bin to set CDF=0)
        smoothedy=self.y
        if fromCDF:
            smoothedy=[]
            ntoys=len(self.y)
            for i in range(ntoys):
                iy = [ self.y[i][j]-self.y[i][j-1] for j in xrange(1,nbins+1) ]
                smoothedy.append(iy)
            smoothedy=np.array(smoothedy)

        #replace contents of the original histogram
        h.Reset('ICE')
        yq=np.percentile(smoothedy,[16,50,84],axis=0)
        for xbin in range(nbins):
            h.SetBinContent(xbin+1, yq[1][xbin])
            h.SetBinError(xbin+1, 0.5*(yq[2][xbin]-yq[0][xbin]))
        h.Scale(norm/h.Integral())

        return h

def testGFSmoother():

    #test case histogram
    y=[(15.9931,4.62851), (57.1474,8.53402), (67.8991,9.29725), (72.7597,9.53358), (55.4123,8.24206), (74.7102,9.70385), (85.5248,10.4569), (62.2573,8.93964), (60.7506,8.72382), (87.9013,10.5906), (94.8423,10.9037), (106.235,11.6979), (76.7571,9.85242), (88.0131,10.6197), (84.174,10.3317), (103.513,11.5849), (87.0253,10.591), (116.627,12.1183), (158.117,14.1105), (603.206,27.8202)]
    nbins=len(y)
    h=ROOT.TH1F('h',';Bin number;Events (a.u.)',nbins,0,nbins)
    for ix in range(nbins):
        h.SetBinContent(ix+1,y[ix][0])
        h.SetBinError(ix+1,y[ix][1])
    rawh=h.Clone('rawh')


    ntoys=100
    gfs=GFSmoother(h,ntoys=100,sigma=1)
    h=gfs.smooth

    import matplotlib.pyplot as plt
    
    #show toys
    fig=plt.figure()
    for ix in range(ntoys):
        plt.plot(gfs.x[ix],color="lightgray")
    plt.errorbar(x=np.linspace(0,nbins-1,nbins), 
                 y=[val for val,_ in y], 
                 yerr=[unc for _,unc in y], 
                 fmt='o')
    plt.xlabel('Bin number')
    plt.ylabel('Counts (a.u.)')
    fig.suptitle('Toys compared to original histogram')
    fig.savefig('gfs_toys.png')

    #show a random interpolation
    plt.clf()
    itoy=np.random.randint(low=0,high=ntoys)
    plt.plot(gfs.y[itoy],color='green',label='Gaussian KDE')
    cdfix = gfs.x[itoy].cumsum() / gfs.x[itoy].sum()
    cdfix[-1]=1
    cdfix=np.insert(cdfix,0,0)
    plt.plot(cdfix,marker='o',color='black',label='Toy CDF')
    plt.xlabel('Bin number')
    plt.ylabel('CDF')
    plt.ylim((0,1))
    leg=plt.legend(loc='upper left',numpoints=1)
    leg.get_frame().set_facecolor('none')
    leg.get_frame().set_linewidth(0.0)
    fig.suptitle('Smoothed toy')
    fig.savefig('gfs_smoothedtoy.png')

    #show final comparison
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetBottomMargin(0.1)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.04)
    rawh.SetMarkerStyle(24)
    rawh.SetLineColor(ROOT.kGray)
    rawh.SetMarkerColor(ROOT.kGray)
    rawh.SetTitle('original')
    rawh.Draw()
    h.SetMarkerStyle(20)
    h.SetLineColor(1)
    h.SetMarkerColor(1)
    h.SetTitle('smoothed from toys')
    h.Draw('e1same')
    leg=c.BuildLegend(0.15,0.9,0.4,0.8)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.04)
    c.SaveAs('gfs_result.png')
    

if __name__ == "__main__":
    sys.exit(testGFSmoother())
