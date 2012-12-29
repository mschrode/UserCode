# $Id: $

import FWCore.ParameterSet.Config as cms

dasTreeMaker = cms.EDAnalyzer(
    'DASTreeMaker',
    MCdata        = cms.bool(False),
    sampleID      = cms.int32(0),
    isSUSY        = cms.bool(False),
    evtWgt        = cms.double(1.),                          
    evtWgtTag     = cms.InputTag("weight"),  
    genjets       = cms.InputTag("ak5GenJets"),  
    genmet        = cms.InputTag("genMetCalo"),
    vertex        = cms.InputTag("offlinePrimaryVertices"),  
    PATjets       = cms.InputTag("allLayer1JetsIC5"),
    PATmet        = cms.InputTag("allLayer1METsIC5"),
    PATmuons      = cms.InputTag("allLayer1Muons"),
    muID          = cms.string('GlobalMuonPromptTight'),
    PATelectrons  = cms.InputTag("allLayer1Electrons"),
    eleID         = cms.string('eidTight'),
    PATphotons    = cms.InputTag("allLayer1Photons"),
    PFRhoTag      = cms.InputTag("kt6PFJetsForGammaIso","rho"),       
    OutFile       = cms.string('ra2DAStree.root')
    )
