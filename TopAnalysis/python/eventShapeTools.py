#!/usr/bin/env python

import ROOT

class EventShapeTool:

    def __init__(self):
        self.momentumTensor=ROOT.TMatrixDSym(3)
        self.eigenValues=ROOT.TVectorD(3)
        self.sphericity=0
        self.aplanarity=0
        self.C=0
        self.D=0

    def analyseNewEvent(self,inputP4,r=2):
        self.computeMomentumTensor(inputP4,r)
        self.computeEventShapes()

    def computeMomentumTensor(self,inputP4,r=2):
        self.momentumTensor.Zero();
        if len(inputP4)<2: return
        norm=0
        for p4 in inputP4:
            p3vec=p4.Vect()
            pR=ROOT.TMath.Power(p3vec.Mag(),r)
            pRminus2=ROOT.TMath.Power(pR,0.5*r-1)
            norm += pR
            self.momentumTensor[0][0] += pRminus2*p3vec.X()*p3vec.X();
            self.momentumTensor[0][1] += pRminus2*p3vec.X()*p3vec.Y();
            self.momentumTensor[0][2] += pRminus2*p3vec.X()*p3vec.Z();
            self.momentumTensor[1][0] += pRminus2*p3vec.Y()*p3vec.X();
            self.momentumTensor[1][1] += pRminus2*p3vec.Y()*p3vec.Y();
            self.momentumTensor[1][2] += pRminus2*p3vec.Y()*p3vec.Z();
            self.momentumTensor[2][0] += pRminus2*p3vec.Z()*p3vec.X();
            self.momentumTensor[2][1] += pRminus2*p3vec.Z()*p3vec.Y();
            self.momentumTensor[2][2] += pRminus2*p3vec.Z()*p3vec.Z();
        
        for i in xrange(0,3):
            for j in xrange(0,3):
                self.momentumTensor[i][j]=(1./norm)*self.momentumTensor[i][j]

    def computeEventShapes(self):
        if not self.momentumTensor.IsSymmetric() : return
        if self.momentumTensor.NonZeros()==0: return
        self.momentumTensor.EigenVectors(self.eigenValues)
        lambdas=sorted([self.eigenValues(i) for i in xrange(0,3)],reverse=True)
        
        self.sphericity=1.5*(lambdas[1]+lambdas[2])
        self.aplanarity=1.5*lambdas[2]
        self.C=3.*(lambdas[0]*lambdas[1] + lambdas[0]*lambdas[2] + lambdas[1]*lambdas[2])
        self.D=27.*lambdas[0]*lambdas[1]*lambdas[2]
