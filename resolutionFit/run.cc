// $Id: run.cc,v 1.45 2011/02/21 18:25:46 mschrode Exp $

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

   // MC truth resolution
   cmd->setMCTruthResolution("config/Parameters_MCTruthResolution.txt",ResolutionFunction::ModifiedNSC);

   // Particle level imbalance
   cmd->setPLI("config/Parameters_PLI.txt",ResolutionFunction::NSC);

   // Samples and FitResult types
   cmd->addMCSample("PYTHIA MC","~/results/ResolutionFit/Pt3Cuts/ResFitThres_Calo_MCFall10");
   //cmd->addMCSample("PYTHIA MC","input/ResFitThres_Calo_MCFall10");
   //cmd->addDataSample("Data","input/ResFitThres_Calo_Data");
   cmd->addFitResult(FitResult::FullMaxLikeRel);
   //cmd->addFitResult(FitResult::PtAsym);

   // Samples to be compared
   //   cmd->compareSamples("Data","PYTHIA MC");

   // Systematic uncertainties
   cmd->addExtrapolationUncertainty("PYTHIA MC",FitResult::FullMaxLikeRel,7);
   cmd->addPLIUncertainty("PYTHIA MC",FitResult::FullMaxLikeRel,8);

    // Output
   cmd->printSetup();
   cmd->makeAllPlots();
   cmd->printResult();

   delete cmd;
   delete par;

  return 0;
}
