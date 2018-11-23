from HiggsAnalysis.CombinedLimit.PhysicsModel import *
from HiggsAnalysis.CombinedLimit.NuisanceModifier import *

class TopPoleMass(PhysicsModel):
    def __init__(self): 
        
        #analysis-specific parameters
        self.dataCardRefXsec=832
        self.accCorrection='1.0-0.0008*(@0-172.5)'

        #PDF-specific parameters
        self.param='${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/xsecvsmtop.dat'
        self.PDF='NNPDF30_nnlo_as_0118'
        self.xsecParam=None
        self.scaleUnc=None
        self.pdfAsUnc=None
        self.beamEnUnc=None
        self.extrapolUnc='1.015'

    def setModelBuilder(self, modelBuilder):
        PhysicsModel.setModelBuilder(self, modelBuilder)
        self.modelBuilder.doModelBOnly = False

    def getYieldScale(self,bin,process):
        "return 1 for backgrounds, r for signal "
        return 1 if not self.DC.isSignal[process] else 'r'
    
    def setPhysicsOptions(self,physOptions):
        for po in physOptions:
            if 'PDF' in po:
                self.PDF=po.split('=')[1]
                print 'Will use the following PDF=',self.PDF
            if 'param' in po:
                self.param=po.split('=')[1]
                print 'Parameterizations to be used will be read from',self.param
            if 'extrapolUnc' in po:
                self.extrapolUnc=po.split('=')[1]
                print 'Extrapoliation uncertainty',self.extrapolUnc
            if 'refxsec' in po:
                self.dataCardRefXsec=float(po.split('=')[1])
                print 'Switching to ref. xsec in datacards=',self.dataCardRefXsec,'pb'
            if 'accCorrection' in po:
                self.accCorrection=po.split('=')[1]
                print 'Defining acceptance correction as',self.accCorrection

        with open(self.param,'r') as f:
            for line in f:
                line=line.rstrip().split()
                if len(line)!=7 or line[0]=='#' : continue
                if line[0]!=self.PDF : continue
                print 'Parsing parametrizations to be used for',self.PDF
                self.xsecParam=(float(line[1]),float(line[2]),float(line[3]))
                self.scaleUnc=line[4]
                self.pdfAsUnc=line[5]
                self.beamEnUnc=line[6]
                break
                                                                                                            
    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        self.modelBuilder.doVar("mpole[172.5,150,200]");
        self.modelBuilder.factory_('expr::acceptanceCorrection("%s",mpole)'%self.accCorrection)
        self.modelBuilder.factory_('expr::mu("(%3.4f/%3.4f)*TMath::Power(172.5/@0,4)*(1+%3.4f*(@0/172.5-1.0)+%3.4f*TMath::Power(@0/172.5-1.0,2))",mpole)' % (self.xsecParam[0],
                                                                                                                                                             self.dataCardRefXsec,
                                                                                                                                                             self.xsecParam[1],
                                                                                                                                                             self.xsecParam[2] ) )
        self.modelBuilder.factory_('expr::r("@0*@1",acceptanceCorrection,mu)')
        self.modelBuilder.doSet('POI','mpole')

        signals=','.join(self.DC.signals)
        doAddNuisance(self.DC, (signals, '*', 'extrapolUnc',   'lnN', self.extrapolUnc) )
        doAddNuisance(self.DC, (signals, '*', 'mpoleScaleUnc', 'lnU', self.scaleUnc) )
        doAddNuisance(self.DC, (signals, '*', 'pdfAsUnc',      'lnN', self.pdfAsUnc) )
        doAddNuisance(self.DC, (signals, '*', 'beamEnUnc',     'lnN', self.beamEnUnc) )

topPoleMass = TopPoleMass()

