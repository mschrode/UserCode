# $Id: makeDASTree_cfg.py,v 1.5 2013/01/11 14:25:09 mschrode Exp $


## --- GLOBAL PARAMETERS -----------------------------------------------------

# Read parameters
from RA2Classic.Utils.CommandLineParams import CommandLineParams
parameters = CommandLineParams()

dataSet_   = parameters.value("data_set","")
globalTag_ = parameters.value("global_tag","")+"::All"
isMC_      = parameters.value("is_mc",True)
isSUSY_    = parameters.value("is_susy",False)
lumi_      = parameters.value("lumi",5295)
hltPath_   = parameters.value("hlt_path","none")


print "***** SETUP ************************************"
print "    dataSet_ : "+dataSet_
print "    hltPath_ : "+hltPath_
print "  globalTag_ : "+globalTag_
print "       isMC_ : "+str(isMC_)
print "     isSUSY_ : "+str(isSUSY_)
print "       lumi_ : "+str(lumi_)
print "************************************************"



## --- GENERAL SETUP ---------------------------------------------------------
from PhysicsTools.PatAlgos.patTemplate_cfg import *
process.outpath.remove(process.out)     # We only want the ROOT tree

process.GlobalTag.globaltag = globalTag_

#-- Message Logger -----------------------------------------------------------
process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.PATSummaryTables = cms.untracked.PSet(
    limit = cms.untracked.int32(10),
    reportEvery = cms.untracked.int32(1)
    )
process.MessageLogger.cerr.FwkReport.reportEvery = 5000


#-- Input Source -------------------------------------------------------------
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(dataSet_)
    )
process.maxEvents.input = -100



## --- SUSY PAT --------------------------------------------------------------
from PhysicsTools.Configuration.SUSY_pattuple_cff import addDefaultSUSYPAT, getSUSY_pattuple_outputCommands

hltMenu = 'HLT'

theJetColls = ['AK5PF']

jetMetCorr = ['L1FastJet', 'L2Relative', 'L3Absolute']
if not isMC_:
    jetMetCorr.append('L2L3Residual')  

# Due to problem in production of LM samples: same event number appears multiple times
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

addDefaultSUSYPAT(process,
                  mcInfo=isMC_,
                  HLTMenu=hltMenu,
                  jetMetCorrections=jetMetCorr,
                  mcVersion='',
                  theJetNames=theJetColls,
                  doSusyTopProjection=False)

# Remove the PAT cleaning and filtering sequences
process.patDefaultSequence.remove(process.selectedPatCandidates)
process.patDefaultSequence.remove(process.cleanPatCandidates)
process.patDefaultSequence.remove(process.countPatCandidates)

# Disable embedment so that lepton cleaning method works
process.patJetsAK5PF.embedCaloTowers = False
process.patJetsAK5PF.embedPFCandidates = False
process.patJetsPF.embedCaloTowers = False
process.patJetsPF.embedPFCandidates = False

#-- Adjust collections to use PFNoPU jets ------------------------------------
    
# do not use Z-mass window for PU subtraction
# such that JEC works properly
process.pfPileUpPF.checkClosestZVertex = cms.bool(False)

# do not remove muons and electrons from the jet clustering input
process.pfIsolatedElectronsPF.isolationCut = -1.0
process.pfIsolatedMuonsPF.isolationCut = -1.0

# do not remove taus from the jet collection
process.pfTausPF.discriminators = cms.VPSet()
process.pfUnclusteredTausPF = process.pfTausPF.clone(
    cut = cms.string("pt < 0")
    )
process.pfTauSequencePF.replace(process.pfTausPF, process.pfTausPF+ process.pfUnclusteredTausPF)
process.pfNoTauPF.topCollection = "pfUnclusteredTausPF"

# make loose clones of the original electron collection
process.pfRelaxedElectronsPF = process.pfIsolatedElectronsPF.clone()
process.pfRelaxedElectronsPF.isolationCut = 9999.0
process.pfElectronsFromVertexPF.dzCut = 9999.0
process.pfElectronsFromVertexPF.d0Cut = 9999.0
process.pfSelectedElectronsPF.cut = ""
process.patElectronsPF.pfElectronSource  = "pfRelaxedElectronsPF"
process.pfElectronSequencePF.replace(process.pfIsolatedElectronsPF,
                                     process.pfIsolatedElectronsPF +
                                     process.pfRelaxedElectronsPF)

# make loose clones of the original muon collection
process.pfRelaxedMuonsPF = process.pfIsolatedMuonsPF.clone()
process.pfRelaxedMuonsPF.isolationCut = 9999.0
process.pfMuonsFromVertexPF.dzCut = 9999.0
process.pfMuonsFromVertexPF.d0Cut = 9999.0
process.pfSelectedMuonsPF.cut = ""
process.patMuonsPF.pfMuonSource  = "pfRelaxedMuonsPF"
process.pfMuonSequencePF.replace(process.pfIsolatedMuonsPF,
                                 process.pfIsolatedMuonsPF +
                                 process.pfRelaxedMuonsPF)


# overwrite default output content
from SandBox.Skims.RA2Content_cff import getRA2PATOutput
process.out.outputCommands = getRA2PATOutput(process)
process.out.dropMetaData = cms.untracked.string('DROPPED')



## --- HLT -------------------------------------------------------------------

process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')
process.hltHighLevel.HLTPaths = cms.vstring(hltPath_+"*")
process.hltHighLevel.andOr = cms.bool(True)
process.hltHighLevel.throw = cms.bool(False)

process.hltSelection = cms.Sequence(
    process.hltHighLevel
    )
if isMC_:
    process.hltSelection.remove(process.hltHighLevel)



## --- RA2 OBJECTS AND CLEANING ----------------------------------------------

process.load('SandBox.Skims.RA2Objects_cff')
process.load('SandBox.Skims.RA2Selection_cff')
process.load('SandBox.Skims.RA2Cleaning_cff')
process.load("SandBox.Skims.provInfoMuons_cfi")
process.load("SandBox.Skims.provInfoElectrons_cfi")

# Set filters in 'filter' mode (rejects events)
process.trackingFailureFilter.taggingMode = False
process.beamHaloFilter.taggingMode        = False
process.inconsistentMuons.taggingMode     = False
process.greedyMuons.taggingMode           = False
process.ra2EcalTPFilter.taggingMode       = False
process.ra2EcalBEFilter.taggingMode       = False
process.hcalLaserEventFilter.taggingMode  = False
process.eeBadScFilter.taggingMode         = False
process.ecalLaserCorrFilter.taggingMode   = False

# Replace eeNoiseFilter by ra2PBNR filter
process.ra2NoiseCleaning.remove(process.eeNoiseFilter)
process.ra2PostCleaning += process.ra2PBNR

# This is the 'tagging' mode of the HBHENoiseFilter
# but we want filtering
process.ra2NoiseCleaning.remove(process.HBHENoiseFilterRA2)

# Tracking-POG Cleaning (decisions have inverted meaning
# compared to the other filters!)
process.manystripclus53X.taggedMode         = cms.untracked.bool(False)
process.toomanystripclus53X.taggedMode      = cms.untracked.bool(False)
process.logErrorTooManyClusters.taggedMode  = cms.untracked.bool(False)
process.trackingPOGCleaningFilters = cms.Sequence(
    ~process.manystripclus53X *
    ~process.toomanystripclus53X *
    ~process.logErrorTooManyClusters
    )
process.ra2PostCleaning += process.trackingPOGCleaningFilters

# The final RA2 object + cleaning sequence
process.cleanpatseq = cms.Sequence(
        process.susyPatDefaultSequence *
        process.ra2StdCleaning *
        process.ra2Objects *
        process.provInfoMuons *
        process.provInfoElectrons *
        process.ra2PostCleaning
        )



## --- HT & MHT PRESELECTION -------------------------------------------------
process.htPFchsFilter.MinHT   = cms.double(300)
process.mhtPFchsFilter.MinMHT = cms.double(200)



## --- RHO FOR PHOTON ISOLATION ----------------------------------------------
from RecoJets.Configuration.RecoPFJets_cff import *
from RecoJets.Configuration.RecoJets_cff import *
process.kt6PFJetsForGammaIso = kt4PFJets.clone(
    rParam       = 0.6,
    doRhoFastjet = True,
    voronoiRfact = 0.9,
    Rho_EtaMax   = 2.5,
    )
process.kt6CaloJetsForGammaIso = kt4CaloJets.clone(
    rParam       = 0.6,
    doRhoFastjet = True,
    voronoiRfact = 0.9,
    Rho_EtaMax   = 2.5,
    )
process.calculateRhoForGamma = cms.Sequence(
    process.kt6PFJetsForGammaIso +
    process.kt6CaloJetsForGammaIso
)



## --- WEIGHT PRODUCER -------------------------------------------------------
from RA2Classic.WeightProducer.getWeightProducer_cff import getWeightProducer
process.WeightProducer = getWeightProducer(dataSet_)
if isMC_:
    process.WeightProducer.Lumi = cms.double(lumi_)
    process.WeightProducer.FileNamePUDataDistribution = cms.string("RA2Classic/WeightProducer/data/DataPileupHistogram_RA2Summer12_190456-196531_AB.root")

from DAS.DASTreeMaker.hltPrescaleWeightProducer_cfi import hltPrescaleWeightProducer
process.HLTPrescaleWeight = hltPrescaleWeightProducer.clone(
    HLTName = cms.string(hltPath_)
    )



## --- DAS TREE MAKER --------------------------------------------------------
from DAS.DASTreeMaker.dasTreeMaker_cfi import dasTreeMaker
from DAS.DASTreeMaker.getSampleID_cff import getSampleID
process.dasTree = dasTreeMaker.clone(
    MCdata        = cms.bool(isMC_),
    isSUSY        = cms.bool(isSUSY_),
    sampleID      = cms.int32(getSampleID(dataSet_)),
    evtWgt        = cms.double(-1.),
    evtWgtTag     = cms.InputTag('WeightProducer:weight'),
#    evtWgtTag     = cms.InputTag('HLTPrescaleWeight'),  
    genjets       = cms.InputTag("ak5GenJets"),  
    genmet        = cms.InputTag("genMetCalo"),
    vertex        = cms.InputTag("offlinePrimaryVertices"),  
    jets          = cms.InputTag("patJetsPF"),
    PATmet        = cms.InputTag("patMETsPF"),
    muons         = cms.InputTag("patMuonsPFIDIso"),
    muID          = cms.string('GlobalMuonPromptTight'),
    electrons     = cms.InputTag("patElectronsIDIso"),
    eleID         = cms.string('eidTight'),
    PATphotons    = cms.InputTag("patPhotons"),
    PFRhoTag      = cms.InputTag("kt6PFJetsForGammaIso","rho"),       
    OutFile       = cms.string('RA2DASTree.root')
    )

#process.dump = cms.EDAnalyzer("EventContentAnalyzer")

process.ppfchs = cms.Path(
    process.hltSelection *
    process.cleanpatseq *
    #process.dump *
    process.calculateRhoForGamma *
    #process.htPFchsFilter *
    #process.mhtPFchsFilter *
    #process.HLTPrescaleWeight *
    process.WeightProducer *
    process.dasTree
    )


process.schedule = cms.Schedule(process.ppfchs)
