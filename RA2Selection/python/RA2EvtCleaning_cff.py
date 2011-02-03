# RA2 cleaning sequences

import FWCore.ParameterSet.Config as cms

from SandBox.Skims.vertex_cfi import goodVertices
from SandBox.Skims.vertex_cfi import oneGoodVertex
from SandBox.Skims.noscraping_cfi import noscraping
from SandBox.Skims.HBHENoiseFilter_cff import HBHENoiseFilter
from SandBox.Skims.singleEventFilter_cfi import singleEventFilter
vetoHCALNoiseEv = singleEventFilter.clone(Run = 149063, Lumi = 19, Event = 21076423)
from SandBox.Skims.eeNoiseFilter_cfi import eeNoiseFilter
from SandBox.Skims.beamHaloFilter_cfi import beamHaloFilter
from SandBox.Skims.trackingFailureFilter_cfi import trackingFailureFilter
trackingFailureFilter.JetSource = cms.InputTag('patJetsPF')


# One EventCounter per filter
from UserCode.CutFlow.EventCounter_cfi import EventCounter

EvtCntTotal              =  EventCounter.clone(FilterName = cms.string("Total"))
EvtCntHLT                =  EventCounter.clone(FilterName = cms.string("HLT"))
EvtCntOneGoodVertex      =  EventCounter.clone(FilterName = cms.string("OneGoodVertex"))
EvtCntNoScraping         =  EventCounter.clone(FilterName = cms.string("NoScraping"))
EvtCntHBHENoise          =  EventCounter.clone(FilterName = cms.string("HBHENoise"))
EvtCntEENoise            =  EventCounter.clone(FilterName = cms.string("EENoise"))
EvtCntBeamHalo           =  EventCounter.clone(FilterName = cms.string("BeamHalo"))
EvtCntTrackingFailure    =  EventCounter.clone(FilterName = cms.string("TrackingFailure"))
EvtCntTPBE               =  EventCounter.clone(FilterName = cms.string("TPBE"))
EvtCntInconsistentPFMuon =  EventCounter.clone(FilterName = cms.string("InconsistentPFMuon"))
EvtCntBadPFMuon          =  EventCounter.clone(FilterName = cms.string("BadPFMuon"))
EvtCntGreedyPFMuon       =  EventCounter.clone(FilterName = cms.string("GreedyPFMuon"))
EvtCntPFElectron         =  EventCounter.clone(FilterName = cms.string("PFElectron"))
EvtCntPFMuon             =  EventCounter.clone(FilterName = cms.string("PFMuon"))



cleaningStd = cms.Sequence(
    goodVertices *
    oneGoodVertex * EvtCntOneGoodVertex
    )

cleaningStdData = cms.Sequence(
    noscraping * EvtCntNoScraping *
    HBHENoiseFilter *
    ~vetoHCALNoiseEv * EvtCntHBHENoise *
    eeNoiseFilter * EvtCntEENoise *
    beamHaloFilter * EvtCntBeamHalo *
    trackingFailureFilter * EvtCntTrackingFailure
    )

from UserCode.EcalDeadCellEventFilter.RA2EcalDeadCell_cff import *
ecalDeadCellTPonlyFilter = ecalDeadCellFilter.clone( filterSelector = 1 )
ecalDeadCellBEonlyFilter = ecalDeadCellFilter.clone( filterSelector = 2 )

deadEcalCellEvtCleaningTPBE = cms.Sequence(
    ecalAnomalousFilter *
    ecalDeadCellTPonlyFilter *
    ecalDeadCellBEonlyFilter *
    EvtCntTPBE
    )

# For MC Fall10
deadEcalCellEvtCleaningBE = ecalAnomalousFilter.clone()
deadEcalCellEvtCleaningBE.FilterAlgo = cms.untracked.string("FilterMode")
deadEcalCellEvtCleaningBE.cutBoundEnergyDeadCellsEB = cms.untracked.double(10)
deadEcalCellEvtCleaningBE.cutBoundEnergyDeadCellsEE = cms.untracked.double(10)
deadEcalCellEvtCleaningBE.limitDeadCellToChannelStatusEB = cms.vint32(12,14)
deadEcalCellEvtCleaningBE.limitDeadCellToChannelStatusEE = cms.vint32(12,14)



from SandBox.Skims.badPFMET_cfi import badPFMET
from SandBox.Skims.muonPFCandidateProducer_cfi import muonPFCandidateProducer

cleaningPFProducts = cms.Sequence(
  badPFMET *
  muonPFCandidateProducer
)


from SandBox.Skims.badPFMuonFilter_cfi import badPFMuonFilter
badPFMuonFilter.PFCandSource = cms.InputTag('muonPFCandidateProducer')
from RecoParticleFlow.PostProcessing.greedyMuonPFCandidateFilter_cfi import greedyMuonPFCandidateFilter
greedyMuonPFCandidateFilter.PFCandidates = cms.InputTag('muonPFCandidateProducer')
from RecoParticleFlow.PostProcessing.inconsistentMuonPFCandidateFilter_cfi import inconsistentMuonPFCandidateFilter
inconsistentMuonPFCandidateFilter.PFCandidates = cms.InputTag('muonPFCandidateProducer')

cleaningPF = cms.Sequence(
  ~inconsistentMuonPFCandidateFilter * EvtCntInconsistentPFMuon *
  badPFMuonFilter * EvtCntBadPFMuon *
  ~greedyMuonPFCandidateFilter * EvtCntGreedyPFMuon
)


from SandBox.Skims.RA2Leptons_cff import *

leptonVetoPF = cms.Sequence(
    ra2PFMuons *
    ra2PFMuonVeto * EvtCntPFMuon *
    ra2PFElectrons *
    ra2PFElectronVeto *EvtCntPFElectron
    )


