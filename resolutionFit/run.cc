// $Id: $

#include <iostream>

#include "TString.h"

#include "../util/StyleSettings.h"
#include "CommanderCool.h"
#include "FitResult.h"
#include "Parameters.h"
#include "ResolutionFunction.h"


using namespace resolutionFit;

int main(int argc, char *argv[]) {

  util::StyleSettings::presentationNoTitle();

  Parameters* par = new Parameters("Test","testFiles/BinningAdmin4.cfg",0);
  par->setJetProperties(JetProperties::AK5,JetProperties::PF);
  
  CommanderCool* cmd = new CommanderCool(par);
  cmd->setMCTruthResolution(ResolutionFunction::ModifiedNSC);
  cmd->setPLI(ResolutionFunction::NSC);
  cmd->addMCSample("PYTHIA MC","testFiles/ResFitThres_Calo_MCFall10");
  cmd->addFitResult(FitResult::FullMaxLike);
  cmd->printSetup();
  cmd->makeAllPlots();
  cmd->printResult();

  delete cmd;
  delete par;

  return 0;
}
