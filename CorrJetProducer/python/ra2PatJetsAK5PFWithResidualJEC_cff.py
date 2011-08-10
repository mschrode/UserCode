import FWCore.ParameterSet.Config as cms

from JetMETCorrections.Configuration.DefaultJEC_cff import *
from UserCode.CorrJetProducer.corrjetproducer_cfi import corrjetproducer

ra2PatJetsAK5PFWithResidualJEC = corrjetproducer.clone(
    src          = cms.InputTag('patJetsAK5PF'),
    jecLevel     = cms.string('ak5PFResidual')
    )
