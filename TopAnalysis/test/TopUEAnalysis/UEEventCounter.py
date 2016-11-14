#!/usr/bin/env/python

import ROOT

"""
parses the event and counts particles in each region at gen/rec levels
"""
class UEEventCounter:

    """
    start variables
    """
    def __init__(self,axes=[]):
 
        self.rec_chmult  = 0           #inclusive rec counts
        self.rec_chflux  = 0
        self.rec_chavgpt = 0
        self.gen_chmult  = 0           #inclusive gen counts
        self.gen_chflux  = 0
        self.gen_chavgpt = 0

        self.rec_chmult_incWrtTo={}    #reco counts per reco region
        self.rec_chflux_incWrtTo={}
        self.rec_chavgpt_incWrtTo={}
        self.rec_chmult_wrtTo={}       #migration matrix true-reco regions
        self.rec_chflux_wrtTo={}
        self.rec_chavgpt_wrtTo={}
        self.gen_chmult_wrtTo={}       #gen counts per gen region
        self.gen_chflux_wrtTo={}
        self.gen_chavgpt_wrtTo={}
        for a in axes:
            self.rec_chmult_incWrtTo[a]  = [0]*3
            self.rec_chflux_incWrtTo[a]  = [0]*3
            self.rec_chavgpt_incWrtTo[a] = [0]*3
            self.rec_chmult_wrtTo[a]  = [[0]*3,[0]*3,[0]*3]
            self.rec_chflux_wrtTo[a]  = [[0]*3,[0]*3,[0]*3]
            self.rec_chavgpt_wrtTo[a] = [[0]*3,[0]*3,[0]*3]
            self.gen_chmult_wrtTo[a]  = [0]*3
            self.gen_chflux_wrtTo[a]  = [0]*3
            self.gen_chavgpt_wrtTo[a] = [0]*3
    
    """
    restart counters
    """
    def reset(self):        
        self.rec_chmult  = 0
        self.rec_chflux  = 0
        self.rec_chavgpt = 0
        self.gen_chmult  = 0
        self.gen_chflux  = 0
        self.gen_chavgpt = 0
        for a in self.rec_chmult_wrtTo:
            self.rec_chmult_incWrtTo[a]  = [0]*3
            self.rec_chflux_incWrtTo[a]  = [0]*3
            self.rec_chavgpt_incWrtTo[a] = [0]*3
            self.rec_chmult_wrtTo[a]     = [[0]*3,[0]*3,[0]*3]
            self.rec_chflux_wrtTo[a]     = [[0]*3,[0]*3,[0]*3]
            self.rec_chavgpt_wrtTo[a]    = [[0]*3,[0]*3,[0]*3]
            self.gen_chmult_wrtTo[a]     = [0]*3
            self.gen_chflux_wrtTo[a]     = [0]*3
            self.gen_chavgpt_wrtTo[a]    = [0]*3

    """
    printout the event contents
    """
    def show(self):
        print 'Level # pTsum <pT>'
        print 'RECO %d %3.1f %3.1f'%( self.rec_chmult, self.rec_chflux, self.rec_chavgpt )
        print 'GEN %d %3.1f %3.1f'%( self.gen_chmult, self.gen_chflux, self.gen_chavgpt )
        for a in self.rec_chmult_wrtTo:
            print 'AXIS:',a
            print '\t#'
            print '\t',self.rec_chmult_wrtTo[a]
            print '\t',self.gen_chmult_wrtTo[a]
            print '\tpTsum'
            print '\t',self.rec_chflux_wrtTo[a]
            print '\t',self.gen_chflux_wrtTo[a]
            print '\t<pT>'
            print '\t',self.rec_chavgpt_wrtTo[a]
            print '\t',self.gen_chavgpt_wrtTo[a]
   
    """
    0 - tow(ards), 1 trans(verse), 2 away
    """
    def getRegionFor(self,dphi) :
        if ROOT.TMath.Abs(dphi) < ROOT.TMath.Pi()/3.     : return 0
        elif ROOT.TMath.Abs(dphi) < 2*ROOT.TMath.Pi()/3. : return 1
        return 2
    
    """
    converts index to name
    """
    def getRegionName(self,idx):
        if idx==0   : return 'tow'
        elif idx==1 : return 'tra'
        elif idx==2 : return 'awa'
        return 'inc'

    """
    count the particles in an event
    """
    def count(self,t,varIdx=0,debug=False):
                
        self.reset()

        #reco level
        passSel=((t.passSel>>varIdx)&0x1)
        if passSel:                    
            for n in xrange(0,t.n):
                isInB=(t.isInBFlags[n] & 0x1)
                if isInB : continue
                self.rec_chmult +=1
                self.rec_chflux += t.pt[n]
                for a in self.rec_chmult_wrtTo:
                    
                    phirec=getattr(t,a)[varIdx]
                    idxrec=self.getRegionFor( ROOT.TVector2.Phi_mpi_pi(t.phi[n]-phirec) )
                    
                    phigen=getattr(t,'gen_'+a)
                    idxgen=self.getRegionFor( ROOT.TVector2.Phi_mpi_pi(t.phi[n]-phigen) )

                    self.rec_chmult_incWrtTo[a][idxrec] +=1
                    self.rec_chflux_incWrtTo[a][idxrec] += t.pt[n]                    
                    self.rec_chmult_wrtTo[a][idxgen][idxrec] +=1
                    self.rec_chflux_wrtTo[a][idxgen][idxrec] += t.pt[n]


            #average pt
            self.rec_chavgpt = self.rec_chflux/self.rec_chmult if self.rec_chmult>0 else 0.
            for a in self.rec_chmult_wrtTo:
                for k in xrange(0,len(self.rec_chmult_wrtTo[a])):
                    for l in xrange(0,len(self.rec_chmult_wrtTo[a][k])):
                        ncounted=self.rec_chmult_wrtTo[a][k][l]
                        if ncounted==0 : continue
                        self.rec_chavgpt_wrtTo[a][k][l]=self.rec_chflux_wrtTo[a][k][l]/ncounted
                for k in xrange(0,len(self.rec_chmult_incWrtTo[a])):
                    ncounted=self.rec_chmult_incWrtTo[a][k]
                    if ncounted==0 : continue
                    self.rec_chavgpt_incWrtTo[a][k]=self.rec_chflux_incWrtTo[a][k]/ncounted

        #gen level
        passSel=(t.gen_passSel&0x1)
        if passSel:
            for n in xrange(0,t.gen_n):
                self.gen_chmult +=1
                self.gen_chflux += t.gen_pt[n]
                for a in self.gen_chmult_wrtTo:
                                        
                    phigen=getattr(t,'gen_'+a)
                    idxgen=self.getRegionFor( ROOT.TVector2.Phi_mpi_pi(t.gen_phi[n]-phigen) )
                    self.gen_chmult_wrtTo[a][idxgen] +=1
                    self.gen_chflux_wrtTo[a][idxgen] += t.gen_pt[n]

            #average pt
            self.gen_chavgpt = self.gen_chflux/self.gen_chmult if self.gen_chmult>0 else 0.
            for a in self.gen_chmult_wrtTo:
                for k in xrange(0,len(self.gen_chmult_wrtTo[a])):
                    ncounted=self.gen_chmult_wrtTo[a][k]
                    if ncounted==0 : continue
                    self.gen_chavgpt_wrtTo[a][k]=self.gen_chflux_wrtTo[a][k]/ncounted
        
        if debug : self.show()
