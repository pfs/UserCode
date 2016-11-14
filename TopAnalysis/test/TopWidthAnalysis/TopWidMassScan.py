from HiggsAnalysis.CombinedLimit.PhysicsModel import *
import pickle

### This base class implements signal yields by production and decay mode
### Specific models can be obtained redefining getYieldScale
class TopWidMassScan(PhysicsModel):
    def __init__(self):
        self.poiMap  = []
        self.pois    = {}
        self.verbose = False

        self.minWidth=0.2*1.324
        self.nomWidth=1.324
        self.maxWidth=4.0*1.324

        self.minMass=169.5
        self.nomMass=172.5
        self.maxMass=175.5

        self.coeffCache=""
        self.binCoeffs = {}
    def yankFloat(stringIn):
        return float(''.join([i for i in stringIn if i.isdigit() or i=='.']))
    def setModelBuilder(self, modelBuilder):
        PhysicsModel.setModelBuilder(self, modelBuilder)
        self.modelBuilder.doModelBOnly = False
    def getYieldScale(self,bin,process):
        "Split in production and decay, and call getSignalYieldScale; return 1 for backgrounds "

        if not self.DC.isSignal[process]: return 1

        if self.verbose:
            print '@%s with process=%s'%(bin,process)

        newbin=bin.replace('_incmlb','')
        target = "%(bin)s/%(process)s" % locals()
        scale = '%s_%s' % (newbin,process)

        print "Will scale ", target, " by ", scale
        #self.modelBuilder.var("CMS_th1x").setMin(bin)
        #self.modelBuilder.var("CMS_th1x").setMax(bin)
        return scale;

    def setPhysicsOptions(self,physOptions):
        for po in physOptions:
            if po.startswith("verbose"):
                self.verbose = True
                print 'Verbose is active'
            if po.startswith("widthScale="):
                scale=po.replace('widthScale=','').split(',')
                if(len(scale) < 3) :
                    print "WARNING: Width scale input improperly!"
                    continue
                self.minWidth=scale[0]
                self.nomWidth=scale[1]
                self.maxWidth=scale[2]
            if po.startswith("massScale="):
                scale=po.replace('massScale=','').split(',')
                if(len(scale) < 3) :
                    print "WARNING: Mass scale input improperly!"
                    continue
                self.minMass=scale[0]
                self.nomMass=scale[1]
                self.maxMass=scale[2]
            if po.startswith("coeffCache="):
                self.coeffCache=po.replace('coeffCache=','')

    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""

        with open(self.coeffCache,'r') as cachefile:
          self.binCoeffs=pickle.load(cachefile)

        # width and mass
        self.modelBuilder.doVar("x[%6.6f,%6.6f,%6.6f]"%(self.nomWidth,self.minWidth,self.maxWidth));
        self.modelBuilder.doVar("m[%6.6f,%6.6f,%6.6f]"%(self.nomMass,self.minMass,self.maxMass));
        poi = "x,m"

        # setup a normalization for each linear fit
        for nBin in self.binCoeffs :
            for sigproc in self.binCoeffs[nBin]:
                a=self.binCoeffs[nBin][sigproc][0]
                b=self.binCoeffs[nBin][sigproc][1]
                c=self.binCoeffs[nBin][sigproc][2]
                d=self.binCoeffs[nBin][sigproc][3]
                e=1#self.binCoeffs[nBin][sigproc][4]

                # linear interpolation formula: Taylor series expansion around nominal point
                intrp="(%6.6f*(@0-%6.6f) + %6.6f*(@1-%6.6f) + %6.6f*(@0-%6.6f)*(@1-%6.6f)+%6.6f"
                intrp+=")/%6.6f"
                intrp=intrp%(a,self.nomWidth,b,self.nomMass,c,self.nomWidth,self.nomMass,d,e)

                self.modelBuilder.factory_("expr::bin%i_%s1p0w(\"%s\", x,m)"%(nBin,sigproc,intrp))

        self.modelBuilder.doSet("POI",poi)

topWidMassScan = TopWidMassScan()

