// $Id: run.cc,v 1.43 2011/02/17 13:42:32 mschrode Exp $

#include <iostream>

#include "TString.h"

#include "CommanderCool.h"
#include "FitResult.h"
#include "Parameters.h"
#include "ResolutionFunction.h"


using namespace resolutionFit;

int main(int argc, char *argv[]) {

   util::StyleSettings::paperNoTitle();

   Parameters* par = new Parameters("Test","config/BinningAdmin4.cfg",1);
   par->setJetProperties(JetProperties::AK5,JetProperties::Calo);

   CommanderCool* cmd = new CommanderCool(par);
   cmd->setMCTruthResolution("config/Parameters_MCTruthResolution.txt",ResolutionFunction::ModifiedNSC);
   cmd->setPLI("config/Parameters_PLI.txt",ResolutionFunction::NSC);
   cmd->addMCSample("PYTHIA MC","input/ResFitThres_Calo_MCFall10");
   cmd->addDataSample("Data","input/ResFitThres_Calo_Data");
   cmd->addFitResult(FitResult::FullMaxLikeRel);
   cmd->compareSamples("Data","PYTHIA MC");
   cmd->printSetup();
   cmd->makeAllPlots();
   cmd->printResult();

   delete cmd;
   delete par;

  return 0;
}
