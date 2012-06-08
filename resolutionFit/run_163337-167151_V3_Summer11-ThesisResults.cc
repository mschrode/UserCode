// $Id: run_163337-167151_V3_Summer11-ThesisResults.cc,v 1.7 2012/06/07 21:10:55 mschrode Exp $

#include <iostream>

#include "TError.h"
#include "TString.h"

#include "CommanderCool.h"
#include "FitResult.h"
#include "OutputManager.h"
#include "Parameters.h"
#include "ResolutionFunction.h"


using namespace resolutionFit;

int main(int argc, char *argv[]) {

  util::StyleSettings::setStyleNoteNoTitle();
  
  gErrorIgnoreLevel = 1001;
  
  //  Parameters* par = new Parameters("~/scratch/ThesisPlots/Res_163337-167151","config/Analysis2011/Binning/BinningAdmin2011_v2.cfg",0);
  Parameters* par = new Parameters("Res_163337-167151","config/Analysis2011/Binning/BinningAdmin2011_v2.cfg",0);
  par->setJetProperties(JetProperties::AK5,JetProperties::PF);
  //par->setOutMode(OutputManager::EPSSingleFilesPlusROOT);
  par->setOutMode(OutputManager::EPSSingleFiles);
  par->setLumi(855.);
  par->setPtSoftAbsMin(10.);
  //par->setNEtaBinsUser(1);
  
  TString pathToSrc = "/afs/naf.desy.de/user/m/mschrode/UserCode/mschrode/";
  TString pathToConfig = pathToSrc+"resolutionFit/config/Analysis2011/";
  TString pathToFitResultsData = pathToSrc+"results/Analysis2011/Run2011A_163337-167151";
  TString pathToFitResultsMC = pathToSrc+"results/Analysis2011/QCD_Pt-15to3000_TuneZ2_Flat_7TeV-pythia6_Summer11-PU_S3_START42_V11-v2";

  FitResult::Type nominalFitResType = FitResult::FullMaxLikeAbs;

  CommanderCool* cmd = new CommanderCool(par);
  
  // Samples and FitResult types
  TString jetTypeStr;
  if( par->jetType() == JetProperties::Calo ) jetTypeStr = "Calo";
  else if( par->jetType() == JetProperties::JPT ) jetTypeStr = "JPT";
  else if( par->jetType() == JetProperties::PF )  jetTypeStr = "PF";
  
  // MC truth resolution
  cmd->setMCTruthResolution(pathToConfig+"Parameters_MCTruthResolution_Summer11_PythiaZ2_L1FastJet_ThesisResults.txt",ResolutionFunction::ModifiedNSC);
  
  
  // Particle level imbalance
  //cmd->fitPLI("CMS Simulation ",pathToFitResultsMC+"/ResFit_PtGenAveBins_MCSummer11_V3_"+jetTypeStr+"_L1FastJet.root",ResolutionFunction::ModifiedNSC);
  cmd->setPLI(pathToConfig+"Parameters_PLI_Summer11_PythiaZ2.txt",ResolutionFunction::ModifiedNSC);
  
  // Samples
  TString idData = "Data";
  TString idMC = "CMS Simulation";

  // Summer11 primary result
  cmd->addDataSample(idData,pathToFitResultsData+"/ResFit_PtAveBins_Data_163337-167151_V3_"+jetTypeStr+"_L1FastJet.root");
  cmd->addMCSample(idMC,pathToFitResultsMC+"/ResFit_PtAveBins_MCSummer11_V3_"+jetTypeStr+"_L1FastJet_Nominal.root");
  cmd->addMCSample("JES Down",pathToFitResultsMC+"/ResFit_PtAveBins_MCSummer11_V3_"+jetTypeStr+"_L1FastJet_JESDown.root");
  cmd->addMCSample("JES Up",pathToFitResultsMC+"/ResFit_PtAveBins_MCSummer11_V3_"+jetTypeStr+"_L1FastJet_JESUp.root");
  cmd->addMCSample("Spec Herwig",pathToFitResultsMC+"/ResFit_PtAveBins_MCSummer11_V3_"+jetTypeStr+"_L1FastJet_SpectrumHerwigpp.root");
  cmd->addMCSample("PU Up",pathToFitResultsMC+"/ResFit_PtAveBins_MCSummer11_V3_"+jetTypeStr+"_L1FastJet_PUUp.root");
       
  cmd->addFitResult(nominalFitResType);

  // Systematic uncertainties
  cmd->addUncertaintyFromVariedSample("PU",1.,idMC,nominalFitResType,"PU Up",kGreen-2);
  cmd->addUncertaintyFromVariedSample("Spectrum",1.,idMC,nominalFitResType,"Spec Herwig",kCyan-7);
  cmd->addUncertaintyFromVariedSample("JES",1.,idMC,nominalFitResType,"JES Down","JES Up",kBlue);
  cmd->addPLIUncertainty(idMC,nominalFitResType,kRed+4);
  cmd->addExtrapolationUncertainty(idMC,nominalFitResType,kRed+2);

  
  // Samples to be compared
  cmd->compareSamples(idData,idMC);
  cmd->fitKValues(nominalFitResType);

  // Output
  cmd->printSetup();
  cmd->makeAllPlots();
  cmd->printResult();
  cmd->printPLI();

  delete cmd;
  delete par;
  
  return 0;
}
