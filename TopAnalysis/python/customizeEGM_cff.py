import FWCore.ParameterSet.Config as cms
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
# EGM corrections :  https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaMiniAODV2

def customizeEGM(process,era,runWithAOD=False):

    if '2016' in era: 
        egmEra='2016-Legacy'
        runEnergyCorrections=False
    if '2017' in era:
        egmEra='2017-Nov17ReReco'   
        runEnergyCorrections=False
    if '2018' in era: 
        egmEra="2018-Prompt"
        runEnergyCorrections=True

    setupEgammaPostRecoSeq(process,
                           isMiniAOD=True,
                           era=egmEra,
                           runEnergyCorrections=runEnergyCorrections)

    if runWithAOD:
        print 'Adapting e/g sources to AOD'
        process.electronMVAValueMapProducer.src = cms.InputTag("")
        process.photonMVAValueMapProducer.src = cms.InputTag("")
        process.photonIDValueMapProducer.src = cms.InputTag("")
    
    process.egammaPostReco=cms.Path(process.egammaPostRecoSeq)
