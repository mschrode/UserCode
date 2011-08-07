## $Id: $
##
## --- Configuration file for Kalibri Ntuples including RA2 selection -----------------
##
## For MC
##


import FWCore.ParameterSet.Config as cms


# ---- Load the PAT defaults -------------------------------------------------------
from PhysicsTools.PatAlgos.patTemplate_cfg import *


# ---- Load the BE filter lists ----------------------------------------------------
process.load("PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_15to30_cfi")
process.load("PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_30to50_cfi")
process.load("PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_50to80_cfi")
process.load("PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_80to120_cfi")
process.load("PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_120to170_cfi")
process.load("PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_170to300_cfi")
process.load("PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_300to470_cfi")
process.load("PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_470to600_cfi")
process.load("PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_600to800_cfi")
process.load("PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_800to1000_cfi")
process.load("PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_1000to1400_cfi")
process.load("PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_1400to1800_cfi")
process.load("PhysicsTools.EcalAnomalousEventFilter.BEFilter_QCD_Pythia_Summer11_1800_cfi")
process.vetoBEFilter = cms.Sequence(
    process.vetoBEFilterQCDPythia15To30 *
    process.vetoBEFilterQCDPythia30To50 *
    process.vetoBEFilterQCDPythia50To80 *
    process.vetoBEFilterQCDPythia80To120 *
    process.vetoBEFilterQCDPythia120To170 *
    process.vetoBEFilterQCDPythia170To300 *
    process.vetoBEFilterQCDPythia300To470 *
    process.vetoBEFilterQCDPythia470To600 *
    process.vetoBEFilterQCDPythia600To800 *
    process.vetoBEFilterQCDPythia800To1000 *
    process.vetoBEFilterQCDPythia1000To1400 *
    process.vetoBEFilterQCDPythia1400To1800 *
    process.vetoBEFilterQCDPythia1800
    )


# ---- Run the process -------------------------------------------------------------
from UserCode.RA2Selection.KalibriNtupleRA2Sequence_cff import runKalibriNtupleWithRA2Selection
runKalibriNtupleWithRA2Selection(process,
                                 globalTag="START42_V13::All",
                                 isData=False,
                                 vetoBEFilterKey='vetoBEFilter',
                                 reportEveryEvt=5000,
                                 testFileName='/store/mc/Summer11/QCD_Pt-15to3000_TuneD6T_Flat_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v1/0000/6AD764B7-D878-E011-9AB5-E41F131817C4.root')
