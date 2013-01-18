# $Id: dasTreeMaker_cfi.py,v 1.2 2013/01/10 17:18:13 mschrode Exp $

import FWCore.ParameterSet.Config as cms

dasTreeMaker = cms.EDAnalyzer(
    'DASTreeMaker',
    MCdata        = cms.bool(False),
    sampleID      = cms.int32(0),
    isSUSY        = cms.bool(False),
    evtWgtTag     = cms.InputTag("weight"), # If non-empty, use this weight variable in the event
    evtWgt        = cms.double(1.),         # If evtWgtTag is empty and evtWgt > 0 use evtWgt
    genjets       = cms.InputTag("ak5GenJets"),  
    genmet        = cms.InputTag("genMetCalo"),
    vertex        = cms.InputTag("offlinePrimaryVertices"),  
    jets          = cms.InputTag("allLayer1JetsIC5"),
    PATmet        = cms.InputTag("allLayer1METsIC5"),
    muons         = cms.InputTag("allLayer1Muons"),
    electrons     = cms.InputTag("allLayer1Electrons"),
    PATphotons    = cms.InputTag("allLayer1Photons"),
    PFRhoTag      = cms.InputTag("kt6PFJetsForGammaIso","rho"),       
    OutFile       = cms.string('ra2DAStree.root')
    )
