// $Id: run.cc,v 1.58 2011/07/18 09:36:47 mschrode Exp $

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

//   util::StyleSettings::paper();
//   gStyle->SetPadTopMargin(0.06);
//   gStyle->SetPadRightMargin(0.05);
  util::StyleSettings::presentationNoTitle();

  gErrorIgnoreLevel = 1001;

  Parameters* par = new Parameters("Res_163337-167151","config/BinningAdmin.cfg",0);
  //par->setNEtaBinsUser(1);
  par->setJetProperties(JetProperties::AK5,JetProperties::PF);
  //par->setOutMode(OutputManager::PSAllInOne);
  par->setOutMode(OutputManager::EPSSingleFiles);

  TString pathToHome = "/afs/naf.desy.de/user/m/mschrode/";
  TString pathToConfig = pathToHome+"UserCode/mschrode/resolutionFit/config/";
  TString pathToFitResultsData = pathToHome+"results/ResolutionFit/Run2011A_163337-167151";
  TString pathToFitResultsMC = pathToHome+"results/ResolutionFit/QCD_Pt-15to3000_TuneZ2_Flat_7TeV-pythia6_Summer11-PU_S3_START42_V11-v2";


  FitResult::Type nominalFitResType = FitResult::FullMaxLikeAbs;


  CommanderCool* cmd = new CommanderCool(par);

  // Samples and FitResult types
  TString jetTypeStr;
  if( par->jetType() == JetProperties::Calo ) jetTypeStr = "Calo";
  else if( par->jetType() == JetProperties::JPT ) jetTypeStr = "JPT";
  else if( par->jetType() == JetProperties::PF )  jetTypeStr = "PF";
  
  // MC truth resolution
  cmd->setMCTruthResolution(pathToConfig+"Analysis2011/Parameters_MCTruthResolution_Summer11_PythiaZ2_L1FastJet_NumPUMay10ReReco_v2.txt",ResolutionFunction::ModifiedNSC);
  

  // Particle level imbalance
  cmd->fitPLI("PLI",pathToFitResultsMC+"/ResFit_PtGenAveBins_MCSummer11_"+jetTypeStr+"_L1FastJet.root",ResolutionFunction::ModifiedNSC);
  //cmd->setPLI(pathToConfig+"Parameters_PLI.txt",ResolutionFunction::ModifiedNSC);
  
  // Samples
  TString idData = "Data";
  TString idMC = "PYTHIA Z2";
  cmd->addDataSample(idData,pathToFitResultsData+"/ResFit_PtAveBins_Data_163337-167151_"+jetTypeStr+"_L1FastJet.root");
  cmd->addMCSample(idMC,pathToFitResultsMC+"/ResFit_PtAveBins_MCSummer11_"+jetTypeStr+"_L1FastJet_Nominal.root");
  cmd->addMCSample("JES Down",pathToFitResultsMC+"/ResFit_PtAveBins_MCSummer11_"+jetTypeStr+"_L1FastJet_JESDown.root");
  cmd->addMCSample("JES Up",pathToFitResultsMC+"/ResFit_PtAveBins_MCSummer11_"+jetTypeStr+"_L1FastJet_JESUp.root");
  cmd->addMCSample("Spec Herwig",pathToFitResultsMC+"/ResFit_PtAveBins_MCSummer11_"+jetTypeStr+"_L1FastJet_SpectrumHerwigpp.root");
  cmd->addMCSample("PU Up",pathToFitResultsMC+"/ResFit_PtAveBins_MCSummer11_"+jetTypeStr+"_L1FastJet_PUUp.root");

  cmd->addFitResult(nominalFitResType);

  // Systematic uncertainties
  cmd->addExtrapolationUncertainty(idMC,nominalFitResType,kCyan+2);
  //cmd->addMCClosureUncertainty(idMC,nominalFitResType,46);
  cmd->addUncertaintyFromVariedSample("JES",1.,idMC,nominalFitResType,"JES Down","JES Up",kBlue-9);
  cmd->addPLIUncertainty(idMC,nominalFitResType,kGreen-1);
  cmd->addUncertaintyFromVariedSample("Spectrum",1.,idMC,nominalFitResType,"Spec Herwig",kBlue-2);
  cmd->addUncertaintyFromVariedSample("PU",1.,idMC,nominalFitResType,"PU Up",kBlue+1);
  
  // Samples to be compared
  cmd->compareSamples(idData,idMC);
  cmd->fitKValues(nominalFitResType);

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
