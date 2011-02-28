// $Id: run.cc,v 1.48 2011/02/26 19:01:13 mschrode Exp $

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
   cmd->fitPLI("MC truth","~/results/ResolutionFit/Pt3Cuts/ResFit_PLI",ResolutionFunction::NSC);
   //cmd->setPLI("config/Parameters_PLI.txt",ResolutionFunction::NSC);

   // Samples and FitResult types
   cmd->addMCSample("PYTHIA MC","~/results/ResolutionFit/Pt3Cuts/ResFitThres_Calo_MCFall10");

//    cmd->addMCSample("PYTHIA MC","~/results/ResolutionFit/Note/ResFitThres_Calo_MCFall10");
//    cmd->addDataSample("Data","~/results/ResolutionFit/Note/ResFitThres_Calo_Data");

//    cmd->addMCSample("PYTHIA MC JES-","~/results/ResolutionFit/Note/ResFitThres_Calo_MCFall10_JESDown10");
//    cmd->addMCSample("PYTHIA MC JES+","~/results/ResolutionFit/Note/ResFitThres_Calo_MCFall10_JESUp10");
//    cmd->addMCSample("PYTHIA MC Spec-","~/results/ResolutionFit/Note/ResFitThres_Calo_MCFall10_SpectrumDown");
//    cmd->addMCSample("PYTHIA MC Spec+","~/results/ResolutionFit/Note/ResFitThres_Calo_MCFall10_SpectrumUp");

   cmd->addFitResult(FitResult::FullMaxLikeRel);
   //cmd->addFitResult(FitResult::PtAsym);

   // Samples to be compared
   //cmd->compareSamples("Data","PYTHIA MC");

   // Systematic uncertainties
//    cmd->addUncertaintyFromVariedSample("JEC",1.,"PYTHIA MC",FitResult::FullMaxLikeRel,"PYTHIA MC JES-","PYTHIA MC JES+",14);
//    cmd->addMCClosureUncertainty("PYTHIA MC",FitResult::FullMaxLikeRel,46);
//    cmd->addExtrapolationUncertainty("PYTHIA MC",FitResult::FullMaxLikeRel,7);
//    cmd->addUncertaintyFromVariedSample("Spectrum",1.,"PYTHIA MC",FitResult::FullMaxLikeRel,"PYTHIA MC Spec-","PYTHIA MC Spec+",38);
//    cmd->addPLIUncertainty("PYTHIA MC",FitResult::FullMaxLikeRel,8);

    // Output
   cmd->printSetup();
   cmd->makeAllPlots();
   cmd->printResult();

   delete cmd;
   delete par;

  return 0;
}
