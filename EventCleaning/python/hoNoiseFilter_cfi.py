import FWCore.ParameterSet.Config as cms

hoNoiseFilter = cms.EDFilter(
    "HONoiseFilter",
    patJetsInputTag   = cms.InputTag(''),
    jetPtMin          = cms.double(30),
    jetEtaMax         = cms.double(5),
    maxHOEfrac        = cms.double(0.4),
    taggingMode       = cms.bool(False)
)
