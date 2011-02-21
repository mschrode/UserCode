// $Id: run.cc,v 1.44 2011/02/18 18:42:22 mschrode Exp $

#include <iostream>

#include "TString.h"

#include "CommanderCool.h"
#include "FitResult.h"
#include "Parameters.h"
#include "ResolutionFunction.h"


using namespace resolutionFit;

int main(int argc, char *argv[]) {

   util::StyleSettings::paperNoTitle();

   Parameters* par = new Parameters("Test","config/BinningAdminPt3.cfg",0);
   par->setJetProperties(JetProperties::AK5,JetProperties::Calo);

   CommanderCool* cmd = new CommanderCool(par);
   cmd->setMCTruthResolution("config/Parameters_MCTruthResolution.txt",ResolutionFunction::ModifiedNSC);
   cmd->setPLI("config/Parameters_PLI.txt",ResolutionFunction::NSC);
   cmd->addMCSample("PYTHIA MC","~/results/ResolutionFit/Pt3Cuts/ResFitThres_Calo_MCFall10");
   //cmd->addMCSample("PYTHIA MC","input/ResFitThres_Calo_MCFall10");
   //cmd->addDataSample("Data","input/ResFitThres_Calo_Data");
   //   cmd->compareSamples("Data","PYTHIA MC");
   cmd->addFitResult(FitResult::FullMaxLikeRel);
   cmd->addFitResult(FitResult::PtAsym);
   cmd->printSetup();
   cmd->makeAllPlots();
   cmd->printResult();

   delete cmd;
   delete par;

  return 0;
}
