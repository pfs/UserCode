import FWCore.ParameterSet.Config as cms
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
# EGM corrections :  https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaMiniAODV2

def customizeEGM(process,era):

    egmEra='2017-Nov17ReReco'
    if '2016' in era: egmEra='2016-Legacy'
    setupEgammaPostRecoSeq(process,runVID=True,era=egmEra)

    process.egammaPostReco=cms.Path(process.egammaPostRecoSeq)
