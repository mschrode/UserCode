// $Id: run_163337-167151_V3_Summer11-PileUpStudy.cc,v 1.3 2012/06/01 18:40:32 mschrode Exp $

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
  
  Parameters* par = new Parameters("Res_163337-167151_PileUpStudy","config/Analysis2011/Binning/BinningAdmin2011_v2.cfg",0);
  par->setJetProperties(JetProperties::AK5,JetProperties::PF);
  par->setOutMode(OutputManager::EPSSingleFilesPlusROOT);
  par->setNEtaBinsUser(1);
  par->setLumi(855.);
  par->setPtSoftAbsMin(10.);
  
  TString pathToSrc = "/afs/naf.desy.de/user/m/mschrode/UserCode/mschrode/";
  TString pathToConfig = pathToSrc+"resolutionFit/config/Analysis2011/";
  TString pathToFitResultsMC = pathToSrc+"results/Analysis2011/QCD_Pt-15to3000_TuneZ2_Flat_7TeV-pythia6_Summer11-PU_S3_START42_V11-v2";

  FitResult::Type nominalFitResType = FitResult::FullMaxLikeAbs;

  CommanderCool* cmd = new CommanderCool(par);
  
  // Samples and FitResult types
  
  // MC truth resolution
  cmd->setMCTruthResolution(pathToConfig+"Parameters_MCTruthResolution_Summer11_PythiaZ2_L1FastJet_ThesisResults.txt",ResolutionFunction::ModifiedNSC);
  
  
  // Particle level imbalance
  cmd->setPLI(pathToConfig+"Parameters_PLI_Summer11_PythiaZ2.txt",ResolutionFunction::ModifiedNSC);
  
  // Summer11 pile-up study
  cmd->addMCSample("N_{PU} #leq 5",pathToFitResultsMC+"/ResFit_PtAveBins_MCSummer11Unweighted_V3_PF_L1FastJet_NPU00-05.root");
  cmd->addMCSample("N_{PU} #geq 6",pathToFitResultsMC+"/ResFit_PtAveBins_MCSummer11Unweighted_V3_PF_L1FastJet_NPU06-99.root");
       
  cmd->addFitResult(nominalFitResType);

  // Samples to be compared
  cmd->compareSamples("N_{PU} #leq 5","N_{PU} #geq 6");
  cmd->fitKValues(nominalFitResType);

  // Output
  cmd->printSetup();
  cmd->makeAllPlots();
  cmd->printResult();
  cmd->printPLI();
  //cmd->printRatioToCombinationFormat();

  delete cmd;
  delete par;
  
  return 0;
}
