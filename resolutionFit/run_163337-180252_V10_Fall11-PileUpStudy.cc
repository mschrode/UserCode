// $Id:  $

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

  util::StyleSettings::setStyleNoteNoTitle();
  
  gErrorIgnoreLevel = 1001;
  
  Parameters* par = new Parameters("PileUpStudy_163337-180252","config/Analysis2011/Binning/BinningAdmin2011_v5.cfg",0);
  //par->setNEtaBinsUser(1);
  par->setJetProperties(JetProperties::AK5,JetProperties::PF);
  //par->setOutMode(OutputManager::PSAllInOne);
  par->setOutMode(OutputManager::EPSSingleFiles);
  //par->setOutMode(OutputManager::EPSSingleFilesPlusROOT);
  par->setLumi(4598.);
  par->setPtSoftAbsMin(10.);
  
  TString pathToHome = "/afs/naf.desy.de/user/m/mschrode/";
  TString pathToConfig = pathToHome+"UserCode/mschrode/resolutionFit/config/Analysis2011/";

  // 2011 pile-up study (asymmetry)
  TString pathToFitResultsData = pathToHome+"/results/ResolutionFit/PileUpStudy";
  TString pathToFitResultsMC = pathToHome+"/results/ResolutionFit/PileUpStudy";

  CommanderCool* cmd = new CommanderCool(par);
  
  // MC truth resolution
  cmd->setMCTruthResolution(pathToConfig+"Parameters_MCTruthResolution_Summer11_PythiaZ2_L1FastJet_NumPUMay10ReReco_v2.txt",ResolutionFunction::ModifiedNSC);
  
  // Particle level imbalance
  cmd->setPLI(pathToConfig+"Parameters_PLI_Summer11_PythiaZ2.txt",ResolutionFunction::ModifiedNSC);
  
  // Samples
  cmd->addDataSample("N_{vtx} #leq 5",pathToFitResultsData+"/PileUpStudy_Data2011_PF_L1FastJet_V10_NVtx00-05.root");
  cmd->addDataSample("N_{vtx} #geq 9",pathToFitResultsData+"/PileUpStudy_Data2011_PF_L1FastJet_V10_NVtx09-99.root");
//   cmd->addMCSample("N_{vtx} #leq 5",pathToFitResultsData+"/PileUpStudy_MCFall11_PF_L1FastJet_V10_NVtx00-05.root");
//   cmd->addMCSample("N_{vtx} #geq 9",pathToFitResultsData+"/PileUpStudy_MCFall11_PF_L1FastJet_V10_NVtx09-99.root");

  // Samples to be compared
  cmd->compareSamples("N_{vtx} #leq 5","N_{vtx} #geq 9");

  // Output
  cmd->printSetup();
  cmd->makeAllPlots();

  delete cmd;
  delete par;

  return 0;
}
