import FWCore.ParameterSet.Config as cms
from EgammaUser.EgammaPostRecoTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
# EGM corrections : https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPostRecoRecipes#106X

def customizeEGM(process,era,runWithAOD=False):

    if '2016' in era: 
        egmEra='2016-Legacy'
        runEnergyCorrections=True
    if '2017' in era:
        egmEra='2017-UL'   
        runEnergyCorrections=True
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
