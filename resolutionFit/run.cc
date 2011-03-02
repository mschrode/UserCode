// $Id: run.cc,v 1.50 2011/03/01 16:52:42 mschrode Exp $

#include <iostream>

#include "TString.h"

#include "CommanderCool.h"
#include "FitResult.h"
#include "Parameters.h"
#include "ResolutionFunction.h"


using namespace resolutionFit;

int main(int argc, char *argv[]) {

   util::StyleSettings::paperNoTitle();

   Parameters* par = new Parameters("Test","config/BinningAdmin4.cfg",0);
   par->setJetProperties(JetProperties::AK5,JetProperties::Calo);

   TString pathToFitResults = "~/results/ResolutionFit/Pt3Cuts/";
   FitResult::Type nominalFitResType = FitResult::MaxLikeKSoftRel;

   CommanderCool* cmd = new CommanderCool(par);

   // MC truth resolution
   cmd->setMCTruthResolution("config/Parameters_MCTruthResolution.txt",ResolutionFunction::ModifiedNSC);

   // Particle level imbalance
   //cmd->fitPLI("MC truth",pathToFitResults+"ResFit_PLI",ResolutionFunction::NSC);
   cmd->setPLI("config/Parameters_PLI.txt",ResolutionFunction::NSC);

   // Samples and FitResult types
//    cmd->addMCSample("PYTHIA MC",pathToFitResults+"ResFitThres_Calo_MCFall10");
//    cmd->addDataSample("Data",pathToFitResults+"ResFitThres_Calo_Data");


   cmd->addMCSample("PYTHIA MC","~/results/ResolutionFit/Note/ResFitThres_Calo_MCFall10");
   cmd->addDataSample("Data","~/results/ResolutionFit/Note/ResFitThres_Calo_Data");

   cmd->addMCSample("PYTHIA MC JES-","~/results/ResolutionFit/Note/ResFitThres_Calo_MCFall10_JESDown10");
   cmd->addMCSample("PYTHIA MC JES+","~/results/ResolutionFit/Note/ResFitThres_Calo_MCFall10_JESUp10");
   cmd->addMCSample("PYTHIA MC Spec-","~/results/ResolutionFit/Note/ResFitThres_Calo_MCFall10_SpectrumDown");
   cmd->addMCSample("PYTHIA MC Spec+","~/results/ResolutionFit/Note/ResFitThres_Calo_MCFall10_SpectrumUp");

   cmd->addFitResult(nominalFitResType);

   // Systematic uncertainties
   cmd->addUncertaintyFromVariedSample("JEC",1.,"PYTHIA MC",nominalFitResType,"PYTHIA MC JES-","PYTHIA MC JES+",14);
   //cmd->addMCClosureUncertainty("PYTHIA MC",nominalFitResType,46);
   cmd->addExtrapolationUncertainty("PYTHIA MC",nominalFitResType,7);
   cmd->addUncertaintyFromVariedSample("Spectrum",1.,"PYTHIA MC",nominalFitResType,"PYTHIA MC Spec-","PYTHIA MC Spec+",38);
   cmd->addPLIUncertainty("PYTHIA MC",nominalFitResType,8);

   // Samples to be compared
   cmd->compareSamples("Data","PYTHIA MC");
   cmd->fitKValues(nominalFitResType);


    // Output
   cmd->printSetup();
   cmd->makeAllPlots();
   cmd->printResult();

   delete cmd;
   delete par;

   std::cout << "\n\n\n****************************************************" << std::endl;
   std::cout << "Plots to be implemented:" << std::endl;
   std::cout << "  extrapolation slope vs pt --> kSoft fits required, less wiggles?" << std::endl;
   std::cout << "  wiggles: LVMINI problem? --> ptasym fits (trigger thresholds?)" << std::endl; 
   std::cout << "  Pt3Rel spectra" << std::endl;
   std::cout << "****************************************************\n\n\n" << std::endl;

  return 0;
}
