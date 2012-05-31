// $Id: $

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
  
  Parameters* par = new Parameters("Res_163337-167151","config/Analysis2011/Binning/BinningAdmin2011_v2.cfg",0);
  par->setJetProperties(JetProperties::AK5,JetProperties::PF);
  //par->setOutMode(OutputManager::PSAllInOne);
  par->setOutMode(OutputManager::EPSSingleFiles);
  //par->setOutMode(OutputManager::EPSSingleFilesPlusROOT);
  par->setLumi(855.);
  par->setPtSoftAbsMin(10.);
  par->setNEtaBinsUser(1);
  
  TString pathToHome = "/afs/naf.desy.de/user/m/mschrode/";
  TString pathToConfig = pathToHome+"UserCode/mschrode/resolutionFit/config/Analysis2011/";

  // Updated Summer11 control plots (inclusive pt spectra)
  TString pathToFitResultsData = pathToHome+"results/ResolutionFit/NewControlPlots";
  TString pathToFitResultsMC = pathToHome+"results/ResolutionFit/NewControlPlots";

  FitResult::Type nominalFitResType = FitResult::PtAsym;

  CommanderCool* cmd = new CommanderCool(par);
  
  // Samples and FitResult types
  TString jetTypeStr;
  if( par->jetType() == JetProperties::Calo ) jetTypeStr = "Calo";
  else if( par->jetType() == JetProperties::JPT ) jetTypeStr = "JPT";
  else if( par->jetType() == JetProperties::PF )  jetTypeStr = "PF";
  
  // MC truth resolution
  cmd->setMCTruthResolution(pathToConfig+"Parameters_MCTruthResolution_Summer11_PythiaZ2_L1FastJet_NumPUMay10ReReco_v2.txt",ResolutionFunction::ModifiedNSC);
  
  
  // Particle level imbalance
  cmd->setPLI(pathToConfig+"Parameters_PLI_Summer11_PythiaZ2.txt",ResolutionFunction::ModifiedNSC);
  
  // Samples
  TString idData = "Data";
  TString idMC = "CMS Simulation";

  // Summer11: updated pt control plots
  cmd->addDataSample(idData,pathToFitResultsData+"/ResFit_PtAveBins_Data_163337-167151_V3_"+jetTypeStr+"_L1FastJet.root");
  cmd->addMCSample(idMC,pathToFitResultsMC+"/ResFit_PtAveBins_MCSummer11_V3_"+jetTypeStr+"_L1FastJet_Nominal.root");
  
  cmd->addFitResult(nominalFitResType);

  // Samples to be compared
  cmd->compareSamples(idData,idMC);
  cmd->fitKValues(nominalFitResType);

  // Output
  cmd->printSetup();
  cmd->makeAllPlots();

  delete cmd;
  delete par;

  return 0;
}
