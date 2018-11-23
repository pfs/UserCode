from HiggsAnalysis.CombinedLimit.PhysicsModel import *

### This base class implements signal yields by production and decay mode
### Specific models can be obtained redefining getYieldScale
class TwoHypothesisTest(PhysicsModel):
    def __init__(self):
        self.muAsPOI    = False
        self.muFloating = False
        self.poiMap  = []
        self.pois    = {}
        self.verbose = False
        self.altSignal  = "ALT"
    def setModelBuilder(self, modelBuilder):
        PhysicsModel.setModelBuilder(self, modelBuilder)
        self.modelBuilder.doModelBOnly = False
    def getYieldScale(self,bin,process):
        "Split in production and decay, and call getSignalYieldScale; return 1 for backgrounds "

        if not self.DC.isSignal[process]: return 1

        #check if it's an alternative process
        isAlt = any(tag in process for tag in self.altSignal)
        if self.verbose:
            print '@%s with process=%s (alt signal=%s)'%(bin,process,str(isAlt))

        if self.pois:
            target = "%(bin)s/%(process)s" % locals()
            scale = 1
            for p, l in self.poiMap:
                for el in l:
                    if re.match(el, target):
                        scale = p + self.sigNorms[isAlt]

            print "Will scale ", target, " by ", scale
            return scale;

        else:
            print "Process ", process, " will get norm ", self.sigNorms[isAlt]
            return self.sigNorms[isAlt]

    def setPhysicsOptions(self,physOptions):
        for po in physOptions:
            if po == "muAsPOI":
                self.muAsPOI = True
                self.muFloating = True
                print 'Signal strength will float as  a free parameter and registered as parameter of interest'
            if po == "muFloating":
                self.muFloating = True
                print 'Signal strength will float as  a free parameter (if --PO muAsPOI is given it will become parameter of interest as well)'
            if po.startswith("altSignal="):
                self.altSignal = po.replace('altSignal=','').split(',')
                print 'Alternative signal model defined from the following tag(s)',self.altSignal
            if po.startswith("verbose"):
                self.verbose = True
                print 'Verbose is active'

    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""


        #fraction of nominal and alternative signal models
        self.modelBuilder.doVar("x[0,0,1]");
        poi = "x"


        if self.muFloating:
            if self.pois:
                for pn,pf in self.pois.items():
                    self.modelBuilder.doVar(pf)
                    if self.muAsPOI:
                        print 'Treating %(pn)s as a POI' % locals()
                        poi += ','+pn

                    self.modelBuilder.factory_('expr::%(pn)s_times_not_x("@0*(1-@1)", %(pn)s, x)' % locals())
                    self.modelBuilder.factory_('expr::%(pn)s_times_x("@0*@1", %(pn)s, x)' % locals())

                #signal normalization for alternative and nominal scenarios
                self.sigNorms = { True:'_times_x', False:'_times_not_x' }

            else:
                self.modelBuilder.doVar("r[1,0,4]");
                if self.muAsPOI:
                    print 'Treating r as a POI'
                    poi += ",r"

                self.modelBuilder.factory_("expr::r_times_not_x(\"@0*(1-@1)\", r, x)")
                self.modelBuilder.factory_("expr::r_times_x(\"@0*@1\", r, x)")
                self.sigNorms = { True:'r_times_x', False:'r_times_not_x' }

        else:
            self.modelBuilder.factory_("expr::not_x(\"(1-@0)\", x)")
            self.sigNorms = { True:'x', False:'not_x' }

        self.modelBuilder.doSet("POI",poi)

twoHypothesisTest = TwoHypothesisTest()

