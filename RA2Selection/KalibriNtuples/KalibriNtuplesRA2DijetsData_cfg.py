## $Id: $
##
## --- Configuration file for Kalibri Ntuples including RA2 selection -----------------
##
## For data
##


import FWCore.ParameterSet.Config as cms


# ---- Load the PAT defaults -------------------------------------------------------
from PhysicsTools.PatAlgos.patTemplate_cfg import *


# ---- Load the BE filter lists ----------------------------------------------------
process.load("PhysicsTools.EcalAnomalousEventFilter.BEFilter_Jet_Run2011A_May10ReReco_v1_cfi")
process.load("PhysicsTools.EcalAnomalousEventFilter.BEFilter_Jet_Run2011A_PromptReco_v4_cfi")
process.vetoBEFilter = cms.Sequence(
    process.vetoBEFilterJetRun2011AMay10ReRecov1 *
    process.vetoBEFilterJetRun2011Av4
    )


# ---- Run the process -------------------------------------------------------------
from UserCode.RA2Selection.KalibriNtupleRA2Sequence_cff import runKalibriNtupleWithRA2Selection
runKalibriNtupleWithRA2Selection(process,
                                 globalTag="GR_R_42_V19::All",
                                 isData=True,
                                 vetoBEFilterKey='vetoBEFilter',
                                 reportEveryEvt=5000,
                                 testFileName='/store/data/Run2011A/Jet/AOD/May10ReReco-v1/0000/8CE5E293-227C-E011-B58F-0024E86E8D4C.root')
