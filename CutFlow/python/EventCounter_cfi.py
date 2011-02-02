import FWCore.ParameterSet.Config as cms
  
EventCounter = cms.EDAnalyzer('EventCounter',
                              FilterName = cms.string('GenericEventCounter')
                              )
