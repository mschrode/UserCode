// $Id: run.cc,v 1.55 2011/06/07 18:23:31 mschrode Exp $

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

  //util::StyleSettings::paper();
  util::StyleSettings::presentationNoTitle();
//   gStyle->SetPadTopMargin(0.06);
//   gStyle->SetPadRightMargin(0.05);
  gErrorIgnoreLevel = 1001;

  Parameters* par = new Parameters("Res_May10ReReco","config/BinningAdmin.cfg",0);
  par->setNEtaBinsUser(4);
  par->setJetProperties(JetProperties::AK5,JetProperties::PF);
  par->setOutMode(OutputManager::PSAllInOne);
  //par->setOutMode(OutputManager::EPSSingleFiles);

  TString pathToHome = "/afs/naf.desy.de/user/m/mschrode/";
  TString pathToFitResults = pathToHome+"results/ResolutionFit/";
  TString pathToConfig = pathToHome+"UserCode/mschrode/resolutionFit/config/";
  TString spectrumBaseName = "~/Kalibri/input/Kalibri_DijetSpectrum_PythiaD6T_Summer11_PtSoft015";

  FitResult::Type nominalFitResType = FitResult::FullMaxLikeAbs;


  CommanderCool* cmd = new CommanderCool(par);

  // Samples and FitResult types
  TString jetTypeStr;
  if( par->jetType() == JetProperties::Calo ) jetTypeStr = "Calo";
  else if( par->jetType() == JetProperties::JPT ) jetTypeStr = "JPT";
  else if( par->jetType() == JetProperties::PF )  jetTypeStr = "PF";
  
  // MC truth resolution
  cmd->setMCTruthResolution(pathToConfig+"Analysis2011/Parameters_MCTruthResolution_Summer11_PythiaD6T_L1FastJet_NumPUMay10ReReco_v2.txt",ResolutionFunction::ModifiedNSC);
  //cmd->setMCTruthResolution(pathToConfig+"Analysis2011/Parameters_MCTruthResolution_Summer11_PythiaD6T_L1FastJet_NumPUGreater9_v2.txt",ResolutionFunction::ModifiedNSC);
  //cmd->setMCTruthResolution(pathToConfig+"Analysis2011/Parameters_MCTruthResolution_Summer11_PythiaD6T_L1FastJet_NumPULess5_v2.txt",ResolutionFunction::ModifiedNSC);
  

  // Particle level imbalance
  //cmd->fitPLI("PLI",pathToFitResults+"QCD_Pt-15to3000_TuneD6T_Flat_7TeV-pythia6_Summer11-PU_S3_START42_V11-v2/mvec10/ResFit_MCSummer11_"+jetTypeStr+"_L1FastJet",ResolutionFunction::ModifiedNSC);
  cmd->setPLI(pathToConfig+"Parameters_PLI.txt",ResolutionFunction::ModifiedNSC);
  
  // Samples
  cmd->addMCSample("MC Summer11",pathToFitResults+"QCD_Pt-15to3000_TuneD6T_Flat_7TeV-pythia6_Summer11-PU_S3_START42_V11-v2/ResFit_MCSummer11_PUMay10ReReco_"+jetTypeStr+"_L1FastJet",spectrumBaseName);
  //cmd->addMCSample("N(PU) #leq 4",pathToFitResults+"QCD_Pt-15to3000_TuneD6T_Flat_7TeV-pythia6_Summer11-PU_S3_START42_V11-v2/ResFit_MCSummer11_NumPULess5_"+jetTypeStr+"_L1FastJet",spectrumBaseName);
  //cmd->addMCSample("N(PU) #geq 10",pathToFitResults+"QCD_Pt-15to3000_TuneD6T_Flat_7TeV-pythia6_Summer11-PU_S3_START42_V11-v2/ResFit_MCSummer11_NumPUGreater9_"+jetTypeStr+"_L1FastJet",spectrumBaseName);
  //cmd->addDataSample("Data May10",pathToFitResults+"Run2011A-May10ReReco_v2/ResFit_Data_"+jetTypeStr+"_L1FastJet",spectrumBaseName);

  cmd->addFitResult(nominalFitResType);

  // Systematic uncertainties
  cmd->addExtrapolationUncertainty("MC Summer11",nominalFitResType,7);
  cmd->addPLIUncertainty("MC Summer11",nominalFitResType,8);

  
  
  // Samples to be compared
//   cmd->compareSamples("Data May10","MC Summer11");
//   cmd->fitKValues(nominalFitResType);

  //cmd->compareSamples("N(PU) #leq 4","Inclusive");
  //cmd->compareSamples("N(PU) #geq 10","Incl");
  //cmd->compareSamples("N(PU) #geq 10","N(PU) #leq 4");

  // Output
  cmd->printSetup();
  cmd->makeAllPlots();
  cmd->printResult();
  //cmd->printRatioToCombinationFormat();
  
  delete cmd;
  delete par;
  
  std::cout << "\n\n\n****************************************************" << std::endl;
  std::cout << "Plots to be implemented:" << std::endl;
  std::cout << "\n  Print PLI\n" << std::endl;
  std::cout << "  Pt3Rel spectra" << std::endl;
  std::cout << "****************************************************\n\n\n" << std::endl;
  
  return 0;
}
