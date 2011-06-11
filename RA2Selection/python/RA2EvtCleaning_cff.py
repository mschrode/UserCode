# $Id: $
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


# RA2 cleaning sequence
ra2Cleaning = cms.Sequence(
    goodVerticesRA2 *
    oneGoodVertex *
    noscraping *
    beamHaloFilter *
    goodPFEventsSequence *
    trackingFailureFilterOnAOD *
    HBHENoiseFilter
    )
