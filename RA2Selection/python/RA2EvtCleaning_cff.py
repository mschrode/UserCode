# $Id: RA2EvtCleaning_cff.py,v 1.3 2011/08/05 09:05:14 mschrode Exp $
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

# ECAL endcap noise filter
from SandBox.Skims.eeNoiseFilter_cfi import eeNoiseFilter

# Lepton veto
from SandBox.Skims.RA2Leptons_cff import *




# Cut flow counter
from UserCode.CutFlow.EventCounter_cfi import EventCounter

EvtCntTotal              =  EventCounter.clone(FilterName = cms.string("Total"))
EvtCntHLT                =  EventCounter.clone(FilterName = cms.string("HLT"))
EvtCntOneGoodVertex      =  EventCounter.clone(FilterName = cms.string("OneGoodVertex"))
EvtCntHBHENoise          =  EventCounter.clone(FilterName = cms.string("HBHENoise"))
EvtCntBeamHalo           =  EventCounter.clone(FilterName = cms.string("BeamHalo"))
EvtCntEENoise            =  EventCounter.clone(FilterName = cms.string("EENoise"))
EvtCntNoScraping         =  EventCounter.clone(FilterName = cms.string("NoScraping"))
EvtCntBadPFMuon          =  EventCounter.clone(FilterName = cms.string("BadPFMuon"))
EvtCntTrackingFailure    =  EventCounter.clone(FilterName = cms.string("TrackingFailure"))
EvtCntTP                 =  EventCounter.clone(FilterName = cms.string("TP"))
EvtCntBE                 =  EventCounter.clone(FilterName = cms.string("BE"))
EvtCntLeptonVeto         =  EventCounter.clone(FilterName = cms.string("LeptonVeto"))



# RA2 cleaning sequence
ra2Cleaning = cms.Sequence(
    EvtCntTotal *
    EvtCntHLT *
    goodVerticesRA2 * oneGoodVertex * EvtCntOneGoodVertex *
    HBHENoiseFilter * EvtCntHBHENoise *
    beamHaloFilter * EvtCntBeamHalo *
    eeNoiseFilter * EvtCntEENoise *
    noscraping * EvtCntNoScraping *
    goodPFEventsSequence * EvtCntBadPFMuon *
    trackingFailureFilterOnAOD * EvtCntTrackingFailure *
    ecalDeadCellTPfilter * EvtCntTP *
    EvtCntBE *
    ra2PFElectrons * ra2PFElectronVeto *
    ra2PFMuons * ra2PFMuonVeto *
    EvtCntLeptonVeto
    )
