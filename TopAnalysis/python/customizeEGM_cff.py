import FWCore.ParameterSet.Config as cms
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq

# EGM corrections : https://hypernews.cern.ch/HyperNews/CMS/get/egamma/2405/1/1/1.html

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

    if runWithAOD: #2 tier solution
        print 'Adapting e/g sources to AOD'
        process.electronMVAValueMapProducer.src = cms.InputTag("")
        process.photonMVAValueMapProducer.src = cms.InputTag("")
        process.photonIDValueMapProducer.src = cms.InputTag("")

    process.egammaPostReco=cms.Path(process.egammaPostRecoSeq)
