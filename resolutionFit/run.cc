// $Id: run.cc,v 1.42 2011/02/15 18:16:55 mschrode Exp $

#include <iostream>

#include "TString.h"

#include "CommanderCool.h"
#include "FitResult.h"
#include "Parameters.h"
#include "ResolutionFunction.h"


using namespace resolutionFit;

int main(int argc, char *argv[]) {

   util::StyleSettings::paperNoTitle();

   Parameters* par = new Parameters("Test","testFiles/BinningAdmin4.cfg",0);
   par->setJetProperties(JetProperties::AK5,JetProperties::PF);

   CommanderCool* cmd = new CommanderCool(par);
   cmd->setMCTruthResolution(ResolutionFunction::ModifiedNSC);
   cmd->setPLI(ResolutionFunction::NSC);
   cmd->addMCSample("PYTHIA MC","testFiles/ResFitThres_Calo_MCFall10");
   cmd->addFitResult(FitResult::FullMaxLikeRel);
   cmd->printSetup();
   cmd->makeAllPlots();
   cmd->printResult();

   delete cmd;
   delete par;

  return 0;
}
