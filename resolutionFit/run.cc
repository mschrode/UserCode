// $Id: run.cc,v 1.51 2011/03/02 11:55:51 mschrode Exp $

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
   par->setJetProperties(JetProperties::AK5,JetProperties::PF);

   TString pathToFitResults = "~/results/ResolutionFit/Pt3Cuts/";
   //TString pathToFitResults = "~/results/ResolutionFit/BinsHauke/";

   FitResult::Type nominalFitResType = FitResult::FullMaxLikeRel;


   CommanderCool* cmd = new CommanderCool(par);

   // MC truth resolution
   cmd->setMCTruthResolution("config/Parameters_MCTruthResolution.txt",ResolutionFunction::ModifiedNSC);

   // Particle level imbalance
   cmd->fitPLI("MC truth",pathToFitResults+"ResFit_PLI",ResolutionFunction::NSC);
   //cmd->setPLI("config/Parameters_PLI.txt",ResolutionFunction::NSC);

   // Samples and FitResult types
   TString jetTypeStr;
   if( par->jetType() == JetProperties::Calo ) jetTypeStr = "Calo";
   else if( par->jetType() == JetProperties::JPT ) jetTypeStr = "JPT";
   else if( par->jetType() == JetProperties::PF )  jetTypeStr = "PF";

   cmd->addMCSample("PYTHIA MC",pathToFitResults+"ResFitThres_"+jetTypeStr+"_MCFall10");
   cmd->addDataSample("Data",pathToFitResults+"ResFitThres_"+jetTypeStr+"_Data");
//    cmd->addMCSample("JES Down",pathToFitResults+"ResFitThres_"+jetTypeStr+"_MCFall10_JESDown");
//    cmd->addMCSample("JES Up",pathToFitResults+"ResFitThres_"+jetTypeStr+"_MCFall10_JESUp");
//    cmd->addMCSample("Spec Down",pathToFitResults+"ResFitThres_"+jetTypeStr+"_MCFall10_SpectrumDown");
//    cmd->addMCSample("Spec Up",pathToFitResults+"ResFitThres_"+jetTypeStr+"_MCFall10_SpectrumUp");

   // Dec22ReReco
   //cmd->addDataSample("Dec 22","~/results/ResolutionFit/Dec22ReReco/ResFitThres_"+jetTypeStr+"_DataDec22");


   // Hauke's binning: Eta0 < 0.5
//    cmd->addMCSample("PYTHIA MC",pathToFitResults+"ResFitThres_BarrelHauke_"+jetTypeStr+"_MCFall10");
//    cmd->addDataSample("Data",pathToFitResults+"ResFitThres_BarrelHauke_"+jetTypeStr+"_Data");


//    cmd->addMCSample("PYTHIA MC","~/results/ResolutionFit/Note/ResFitThres_"+jetTypeStr+"_MCFall10");
//    cmd->addDataSample("Data","~/results/ResolutionFit/Note/ResFitThres_"+jetTypeStr+"_Data");

//    cmd->addMCSample("PYTHIA MC JES-","~/results/ResolutionFit/Note/ResFitThres_"+jetTypeStr+"_MCFall10_JESDown5");
//    cmd->addMCSample("PYTHIA MC JES+","~/results/ResolutionFit/Note/ResFitThres_"+jetTypeStr+"_MCFall10_JESUp5");
//    cmd->addMCSample("PYTHIA MC Spec-","~/results/ResolutionFit/Note/ResFitThres_"+jetTypeStr+"_MCFall10_SpectrumDown");
//    cmd->addMCSample("PYTHIA MC Spec+","~/results/ResolutionFit/Note/ResFitThres_"+jetTypeStr+"_MCFall10_SpectrumUp");

   cmd->addFitResult(nominalFitResType);

   // Systematic uncertainties
   //cmd->addUncertaintyFromVariedSample("JEC",1.,"PYTHIA MC",nominalFitResType,"JES Down","JES Up",14);
   //cmd->addMCClosureUncertainty("PYTHIA MC",nominalFitResType,46);
   cmd->addExtrapolationUncertainty("PYTHIA MC",nominalFitResType,7);
   //cmd->addUncertaintyFromVariedSample("Spectrum",1.,"PYTHIA MC",nominalFitResType,"Spec Down","Spec Up",38);
   cmd->addPLIUncertainty("PYTHIA MC",nominalFitResType,8);

   // Samples to be compared
   cmd->compareSamples("Data","PYTHIA MC");
   //cmd->compareSamples("Dec 22","PYTHIA MC");
   cmd->fitKValues(nominalFitResType);


    // Output
   cmd->printSetup();
   cmd->makeAllPlots();
   cmd->printResult();

   delete cmd;
   delete par;

   std::cout << "\n\n\n****************************************************" << std::endl;
   std::cout << "Plots to be implemented:" << std::endl;
   std::cout << "\n  SET TITLE_\n" << std::endl;
   std::cout << "\n  Print PLI\n" << std::endl;
   std::cout << "  wiggles: LVMINI problem? --> ptasym fits (trigger thresholds?)" << std::endl; 
   std::cout << "  Pt3Rel spectra" << std::endl;
   std::cout << "****************************************************\n\n\n" << std::endl;

  return 0;
}
