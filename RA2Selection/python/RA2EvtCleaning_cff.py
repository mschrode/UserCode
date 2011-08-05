# $Id: RA2EvtCleaning_cff.py,v 1.2 2011/06/11 13:43:09 mschrode Exp $
# Definition of the filters at
# https://twiki.cern.ch/twiki/bin/view/CMS/SusyRA2NJetsInData2011

import FWCore.ParameterSet.Config as cms


# Vertex filter
from SandBox.Skims.vertex_cfi import goodVerticesRA2, oneGoodVertex

# PKAM filter
from SandBox.Skims.noscraping_cfi import noscraping

# Tracking failure filter
from SandBox.Skims.trackingFailureFilter_cfi import trackingFailureFilter
trackingFailureFilterOnAOD = trackingFailureFilter.clone(
    JetSource = cms.InputTag('ak5PFJets')
    )

# Bad PF muons
from RecoParticleFlow.PostProcessing.greedyMuonPFCandidateFilter_cfi import greedyMuonPFCandidateFilter
from RecoParticleFlow.PostProcessing.inconsistentMuonPFCandidateFilter_cfi import inconsistentMuonPFCandidateFilter
inconsistentMuons = inconsistentMuonPFCandidateFilter.clone()
greedyMuons = greedyMuonPFCandidateFilter.clone()
goodPFEventsSequence = cms.Sequence(
    ~inconsistentMuons +
    ~greedyMuons  
    )

# HBHENoise
from SandBox.Skims.HBHENoiseFilter_cff import HBHENoiseFilter

# Beam halo filter
from SandBox.Skims.beamHaloFilter_cfi import beamHaloFilter

# Dead ECAL cell TP filter
from JetMETAnalysis.ecalDeadCellTools.RA2TPfilter_cff import ecalDeadCellTPfilter




# Cut flow counter
from UserCode.CutFlow.EventCounter_cfi import EventCounter

EvtCntTotal              =  EventCounter.clone(FilterName = cms.string("Total"))
EvtCntHLT                =  EventCounter.clone(FilterName = cms.string("HLT"))
EvtCntOneGoodVertex      =  EventCounter.clone(FilterName = cms.string("OneGoodVertex"))
EvtCntNoScraping         =  EventCounter.clone(FilterName = cms.string("NoScraping"))
EvtCntBeamHalo           =  EventCounter.clone(FilterName = cms.string("BeamHalo"))
EvtCntBadPFMuon          =  EventCounter.clone(FilterName = cms.string("BadPFMuon"))
EvtCntTrackingFailure    =  EventCounter.clone(FilterName = cms.string("TrackingFailure"))
EvtCntHBHENoise          =  EventCounter.clone(FilterName = cms.string("HBHENoise"))
EvtCntTP                 =  EventCounter.clone(FilterName = cms.string("TP"))
EvtCntBE                 =  EventCounter.clone(FilterName = cms.string("BE"))


# RA2 cleaning sequence
ra2Cleaning = cms.Sequence(
    EvtCntTotal *
    EvtCntHLT *
    goodVerticesRA2 * oneGoodVertex * EvtCntOneGoodVertex *
    noscraping * EvtCntNoScraping *
    beamHaloFilter * EvtCntBeamHalo *
    goodPFEventsSequence * EvtCntBadPFMuon *
    trackingFailureFilterOnAOD * EvtCntTrackingFailure *
    HBHENoiseFilter * EvtCntHBHENoise *
    ecalDeadCellTPfilter * EvtCntTP *
    EvtCntBE
    )
