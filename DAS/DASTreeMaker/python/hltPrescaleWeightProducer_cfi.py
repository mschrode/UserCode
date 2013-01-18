import FWCore.ParameterSet.Config as cms

hltPrescaleWeightProducer = cms.EDProducer(
    'HLTPrescaleWeightProducer',
    HLTName = cms.string('')
    )
