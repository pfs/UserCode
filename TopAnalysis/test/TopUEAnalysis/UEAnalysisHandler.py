#!/usr/bin/env/python

import ROOT
import pickle
import numpy as np
import array as array


"""
parses the event and counts particles in each region at gen/rec levels
"""
class UEAnalysisHandler:

    def __init__(self,analysisCfg):
        
        with open(analysisCfg,'r') as cachefile:
            self.axes=pickle.load(cachefile)
            self.sliceAxes=pickle.load(cachefile)
            self.histos=pickle.load(cachefile)

    """
    inclusive histogram filling
    """
    def fillInclusive(self,obs,ue,weight,gen_passSel,passSel,sliceVarVals=None):

        if not passSel and not gen_passSel: return

        sliceVar=None
        recSliceShift,genSliceShift=0,0
        try:
            sliceVar,genSliceVal,recSliceVal=sliceVarVals
            genSliceShift=self.getBinForVariable(genSliceVal, self.sliceAxes[ (sliceVar,False) ])-1
            recSliceShift=self.getBinForVariable(recSliceVal, self.sliceAxes[ (sliceVar,True) ])-1
        except:
            pass

        #RECO level counting
        recCts=getattr(ue,'rec_chmult')
        recVal=getattr(ue,'rec_'+obs)
        recBin=self.getBinForVariable( recVal,  self.axes[(obs,'inc',0,True)])-1
        recBin += recSliceShift*(self.axes[(obs,'inc',0,True)].GetNbins())
        if not passSel : recBin=-1
        if recCts>0:
            if sliceVar :
                self.histos[(obs,'inc',True,sliceVar)].Fill(recBin,weight)
            else:
                self.histos[(obs,'inc',True)].Fill(recBin,weight)

        #GEN level counting
        genCts=getattr(ue,'gen_chmult')
        genVal=getattr(ue,'gen_'+obs)
        genBin=self.getBinForVariable(genVal, self.axes[(obs,'inc',0,False)])-1
        genBin += genSliceShift*(self.axes[(obs,'inc',0,False)].GetNbins())
        if not gen_passSel : genBin=-1        
        if genCts>0: 
            if sliceVar:
                self.histos[(obs,'inc',False,sliceVar)].Fill(genBin,weight)
            else:
                self.histos[(obs,'inc',False)].Fill(genBin,weight)
                
        #Migration matrix
        if genCts>0 or recCts>0 : 
            if sliceVar:
                self.histos[(obs,'inc',sliceVar)].Fill(genBin,recBin,weight)
            else:
                self.histos[(obs,'inc')].Fill(genBin,recBin,weight)

    """
    differential histogram filling
    """
    def fillDifferential(self,obs,a,ue,weight,gen_passSel,passSel,sliceVarVals=None):

        if not passSel and not gen_passSel: return

        sliceVar=None
        recSliceShift,genSliceShift=0,0
        try:
            sliceVar,genSliceVal,recSliceVal=sliceVarVals
            genSliceShift=self.getBinForVariable(genSliceVal, self.sliceAxes[ (sliceVar,False) ])-1
            recSliceShift=self.getBinForVariable(recSliceVal, self.sliceAxes[ (sliceVar,True) ])-1
        except:
            pass

        recCts=getattr(ue,'rec_chmult_incWrtTo')[a]
        recCtsMtrx=getattr(ue,'rec_chmult_wrtTo')[a]
        recVal=getattr(ue,'rec_%s_incWrtTo'%obs)[a]
        recValMtrx=getattr(ue,'rec_%s_wrtTo'%obs)[a]

        genCts=getattr(ue,'gen_chmult_wrtTo')[a]        
        genVal=getattr(ue,'gen_%s_wrtTo'%obs)[a]

        #RECO level counting
        recBinOffset=0
        totalBinsRec=0
        for idx_rec in xrange(0,3) : totalBinsRec += self.axes[(obs,a,idx_rec,True)].GetNbins()
        for idx_rec in xrange(0,3):
            
            if recCts[idx_rec]==0: continue
            recBin  = recBinOffset+self.getBinForVariable( recVal[idx_rec],  self.axes[(obs,a,idx_rec,True)])-1
            recBin += recSliceShift*totalBinsRec
            if not passSel : recBin=-1
            if sliceVar:
                self.histos[(obs,a,True,sliceVar)].Fill(recBin,weight)
            else:
                self.histos[(obs,a,True)].Fill(recBin,weight)
            recBinOffset += self.axes[(obs,a,idx_rec,True)].GetNbins()

        
        #MC truth counting
        genBinOffset=0
        totalBinsGen=0
        for idx_gen in xrange(0,3) : totalBinsGen += self.axes[(obs,a,idx_gen,False)].GetNbins()
        for idx_gen in xrange(0,3):
            
            genBin=-1
            if genCts[idx_gen]>0 :
                genBin  = genBinOffset+self.getBinForVariable( genVal[idx_gen], self.axes[(obs,a,idx_gen,False)])-1
                genBin += genSliceShift*totalBinsGen
                if not gen_passSel : genBin=-1
                if sliceVar:
                    self.histos[(obs,a,False,sliceVar)].Fill(genBin,weight)
                else:
                    self.histos[(obs,a,False)].Fill(genBin,weight)
                
            recBinOffset=0
            for idx_rec in xrange(0,3):
            
                recBin=-1
                if recCtsMtrx[idx_gen][idx_rec]>0:
                    recBin  = recBinOffset+self.getBinForVariable( recVal[idx_rec],  self.axes[(obs,a,idx_rec,True)])-1
                    recBin += recSliceShift*totalBinsRec
                    if not passSel : recBin=-1
                if genCts[idx_gen]>0 or recCtsMtrx[idx_gen][idx_rec]>0: 
                    if sliceVar:
                        self.histos[(obs,a,sliceVar)].Fill(genBin,recBin,weight)
                    else:
                        self.histos[(obs,a)].Fill(genBin,recBin,weight)

                recBinOffset += self.axes[(obs,a,idx_rec,True)].GetNbins()

            genBinOffset += self.axes[(obs,a,idx_gen,False)].GetNbins()



    """
    return the most appropriate bin for a given value, taking into account the range available
    """
    def getBinForVariable(self,val,axis):
        xmin,xmax=axis.GetXmin(),axis.GetXmax()
        if val>xmax : return axis.GetNbins()
        if val<xmin : return 0
        return axis.FindBin(val)
