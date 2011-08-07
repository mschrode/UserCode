# $Id: $
#
# Process setup for Kalibri Ntuple including RA2 selection
# See https://twiki.cern.ch/twiki/bin/view/CMS/SusyRA2NJetsInData2011
# for a definition of the RA2 selection

import FWCore.ParameterSet.Config as cms

def runKalibriNtupleWithRA2Selection(process,
                                     globalTag,
                                     isData=True,
                                     vetoBEFilterKey='vetoBEFilter',
                                     reportEveryEvt=5000,
                                     testFileName='/store/data/Run2011A/Jet/AOD/May10ReReco-v1/0000/8CE5E293-227C-E011-B58F-0024E86E8D4C.root',
                                     numProcessedEvt=100):


    # ---- Configuration ----------------------------------------------------------
    process.load('Configuration/StandardSequences/Services_cff')
    process.load('Configuration.StandardSequences.MagneticField_38T_cff')
    process.load('Configuration.StandardSequences.Reconstruction_cff')
    process.load('Configuration.StandardSequences.EndOfProcess_cff')
    process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
    process.GlobalTag.globaltag = globalTag


    # ---- Message logger ---------------------------------------------------------
    process.load('FWCore.MessageService.MessageLogger_cfi')
    process.MessageLogger.cerr.threshold             = 'INFO'
    process.MessageLogger.cerr.FwkReport.reportEvery = reportEveryEvt


    # ---- Input ------------------------------------------------------------------
    process.source = cms.Source(
        "PoolSource",
        fileNames = cms.untracked.vstring(testFileName)
    )
    process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(numProcessedEvt)
    )
    process.options = cms.untracked.PSet(
        Rethrow = cms.untracked.vstring('ProductNotFound'),
        wantSummary = cms.untracked.bool(True)
    )


    # ---- HLT --------------------------------------------------------------------
    process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')
    process.hltHighLevel.HLTPaths = cms.vstring('HLT_DiJetAve*')
    process.hltHighLevel.andOr = cms.bool(True)
    process.hltHighLevel.throw = cms.bool(False)


    # ---- SUSYPAT specifics ------------------------------------------------------
    from PhysicsTools.Configuration.SUSY_pattuple_cff import addDefaultSUSYPAT, getSUSY_pattuple_outputCommands

    hltMenu = 'HLT'
    if not isData:
        hltMenu = 'REDIGI'

    jetMetCorr = ['L1FastJet', 'L2Relative', 'L3Absolute']
    if isData:
        jetMetCorr.append('L2L3Residual')
        
    addDefaultSUSYPAT(process,mcInfo=not isData,HLTMenu=hltMenu,jetMetCorrections=jetMetCorr,
                      mcVersion='',theJetNames=['AK5PF'],doSusyTopProjection=True)


    # ---- RA2 cleaning ----------------------------------------------------------
    process.load('UserCode.RA2Selection.RA2EvtCleaning_cff')
    if isData:
        process.ra2Cleaning.replace( process.EvtCntHLT,
                                     process.hltHighLevel * process.EvtCntHLT )

    if vetoBEFilterKey in process.sequences_():
        process.ra2Cleaning.replace( process.EvtCntBE,
                                     process.sequences_()[vetoBEFilterKey] * process.EvtCntBE )

    process.TFileService = cms.Service(
        "TFileService",
        fileName = cms.string("CutFlow.root")
        )


    # ---- JEC --------------------------------------------------------------------
    process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
    process.load('RecoJets.Configuration.RecoJPTJets_cff')
    process.kt6CaloJets.doRhoFastjet = True
    process.kt6CaloJets.doAreaFastjet = True
    process.ak5CaloJets.doAreaFastjet = True
    process.kt6PFJets.doRhoFastjet = True
    process.kt6PFJets.doAreaFastjet = True
    process.ak5PFJets.doAreaFastjet = True


    # ---- Kalibri ntuple maker ---------------------------------------------------
    process.load("Calibration.CalibTreeMaker.CalibTreeMaker_cfi")
    process.calibTreeMakerCalo.WriteStableGenParticles = False
    process.calibTreeMakerCalo.WriteStableGenParticles = False
    process.calibTreeMakerCalo.OutputFile              = 'ak5Calo.root'
    process.calibTreeMakerCalo.NJet_Jets               = 'ak5CaloJets'
    process.calibTreeMakerCalo.NJet_JetIDs             = 'ak5JetID'
    process.calibTreeMakerCalo.NJet_GenJets            = 'ak5GenJets'
    process.calibTreeMakerCalo.NJetZSPJets             = 'ZSPJetCorJetAntiKt5'
    process.calibTreeMakerCalo.NJet_L1JetCorrector     = cms.string('ak5CaloL1Offset')
    process.calibTreeMakerCalo.NJet_L2JetCorrector     = cms.string('ak5CaloL2Relative')
    process.calibTreeMakerCalo.NJet_L3JetCorrector     = cms.string('ak5CaloL3Absolute')
    process.calibTreeMakerCalo.NJet_JPTZSPCorrector    = cms.string('JetPlusTrackZSPCorrectorAntiKt5')
    if isData:
        process.calibTreeMakerCalo.NJet_L1L2L3JetCorrector = cms.string('ak5CaloL1L2L3Residual')
    else:
        process.calibTreeMakerCalo.NJet_L1L2L3JetCorrector = cms.string('ak5CaloL1L2L3')
    process.calibTreeMakerCalo.NJet_L1L2L3L4JWJetCorrector = process.calibTreeMakerCalo.NJet_L1L2L3JetCorrector
    process.calibTreeMakerCalo.NJet_L2L3JetCorrectorJPT = cms.string('ak5JPTL2L3')
    process.calibTreeMakerCalo.NJetConeSize            = 0.5
 

    process.calibTreeMakerPF.WriteStableGenParticles     = False
    process.calibTreeMakerPF.OutputFile                  = 'ak5PF.root'
    process.calibTreeMakerPF.NJet_Jets                   = 'ak5PFJets'
    process.calibTreeMakerPF.NJet_JetIDs                 = ''
    process.calibTreeMakerPF.NJet_GenJets                = 'ak5GenJets'
    process.calibTreeMakerPF.NJetZSPJets                 = 'ZSPJetCorJetAntiKt5'
    process.calibTreeMakerPF.NJet_L1JetCorrector         = cms.string('ak5PFL1Offset')
    process.calibTreeMakerPF.NJet_L2JetCorrector         = cms.string('ak5PFL2Relative')
    process.calibTreeMakerPF.NJet_L3JetCorrector         = cms.string('ak5PFL3Absolute')
    process.calibTreeMakerPF.NJet_JPTZSPCorrector        = cms.string('JetPlusTrackZSPCorrectorAntiKt5')
    if isData:
        process.calibTreeMakerPF.NJet_L1L2L3JetCorrector     = cms.string('ak5PFL1L2L3Residual')
    else:
        process.calibTreeMakerPF.NJet_L1L2L3JetCorrector     = cms.string('ak5PFL1L2L3')
    process.calibTreeMakerPF.NJet_L1L2L3L4JWJetCorrector = process.calibTreeMakerPF.NJet_L1L2L3JetCorrector
    process.calibTreeMakerPF.NJet_L2L3JetCorrectorJPT    = cms.string('ak5JPTL2L3')
    process.calibTreeMakerPF.NJetConeSize                = 0.5


    process.calibTreeMakerAK5FastCalo = process.calibTreeMakerCalo.clone(
        OutputFile                  = cms.string('ak5FastCalo.root'),
        NJet_L1JetCorrector         = cms.string('ak5CaloL1Fastjet')
        )
    if isData:
        process.calibTreeMakerAK5FastCalo.NJet_L1L2L3JetCorrector     = cms.string('ak5CaloL1FastL2L3Residual')
        process.calibTreeMakerAK5FastCalo.NJet_L1L2L3L4JWJetCorrector = cms.string('ak5CaloL1FastL2L3Residual')
    else:
        process.calibTreeMakerAK5FastCalo.NJet_L1L2L3JetCorrector     = cms.string('ak5CaloL1FastL2L3')
        process.calibTreeMakerAK5FastCalo.NJet_L1L2L3L4JWJetCorrector = cms.string('ak5CaloL1FastL2L3')


    process.calibTreeMakerAK5FastPF = process.calibTreeMakerPF.clone(
        OutputFile                  = cms.string('ak5FastPF.root'),
        NJet_L1JetCorrector         = cms.string('ak5PFL1Fastjet')
        )
    if isData:
        process.calibTreeMakerAK5FastPF.NJet_L1L2L3JetCorrector     = cms.string('ak5PFL1FastL2L3Residual')
        process.calibTreeMakerAK5FastPF.NJet_L1L2L3L4JWJetCorrector = cms.string('ak5PFL1FastL2L3Residual')
    else:
        process.calibTreeMakerAK5FastPF.NJet_L1L2L3JetCorrector     = cms.string('ak5PFL1FastL2L3')
        process.calibTreeMakerAK5FastPF.NJet_L1L2L3L4JWJetCorrector = cms.string('ak5PFL1FastL2L3')



    # ---- The paths  -------------------------------------------------------------
    process.pKalibriNtupleMaker = cms.Path(
        process.susyPatDefaultSequence *
        process.ra2Cleaning *
        process.recoJets *
        process.recoPFJets *
        process.calibTreeMakerAK5FastCalo *
        process.calibTreeMakerAK5FastPF
        )

    process.schedule = cms.Schedule(process.pKalibriNtupleMaker)


