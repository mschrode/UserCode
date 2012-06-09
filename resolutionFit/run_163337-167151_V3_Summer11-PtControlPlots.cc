// $Id: run_163337-167151_V3_Summer11-PtControlPlots.cc,v 1.3 2012/06/05 22:44:46 mschrode Exp $

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
  
  Parameters* par = new Parameters("Res_163337-167151_PtControlPlots","config/Analysis2011/Binning/BinningAdmin2011_v2.cfg",0);
  par->setJetProperties(JetProperties::AK5,JetProperties::PF);
  par->setOutMode(OutputManager::EPSSingleFiles);
  par->setLumi(855.);
  par->setPtSoftAbsMin(10.);
  par->setNEtaBinsUser(1);
  
  TString pathToSrc = "/afs/naf.desy.de/user/m/mschrode/UserCode/mschrode/";
  TString pathToConfig = pathToSrc+"resolutionFit/config/Analysis2011/";

  // Updated Summer11 control plots (inclusive pt spectra)
  TString pathToFitResultsData = pathToSrc+"results/Analysis2011/ControlPlots";
  TString pathToFitResultsMC = pathToSrc+"results/Analysis2011/ControlPlots";

  FitResult::Type nominalFitResType = FitResult::PtAsym;

  CommanderCool* cmd = new CommanderCool(par);
  
  // Samples and FitResult types
  TString jetTypeStr;
  if( par->jetType() == JetProperties::Calo ) jetTypeStr = "Calo";
  else if( par->jetType() == JetProperties::JPT ) jetTypeStr = "JPT";
  else if( par->jetType() == JetProperties::PF )  jetTypeStr = "PF";
  
  // MC truth resolution
  cmd->setMCTruthResolution(pathToConfig+"Parameters_MCTruthResolution_Summer11_PythiaZ2_L1FastJet_ThesisResults.txt",ResolutionFunction::ModifiedNSC);
  
  
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
