import FWCore.ParameterSet.Config as cms

corrjetproducer = cms.EDProducer(
    'CorrJetProducer',
    src          = cms.InputTag(''),
    jecLevel     = cms.string(''),
    instanceName = cms.string('')
    )
