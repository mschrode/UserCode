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
  
  Parameters* par = new Parameters("Res_163337-167151","config/Analysis2011/Binning/BinningAdmin2011_v2.cfg",0);
  par->setJetProperties(JetProperties::AK5,JetProperties::PF);
  //par->setOutMode(OutputManager::PSAllInOne);
  par->setOutMode(OutputManager::EPSSingleFiles);
  //par->setOutMode(OutputManager::EPSSingleFilesPlusROOT);
  par->setLumi(855.);
  par->setPtSoftAbsMin(10.);
  //par->setNEtaBinsUser(1);
  
  TString pathToHome = "/afs/naf.desy.de/user/m/mschrode/";
  TString pathToConfig = pathToHome+"UserCode/mschrode/resolutionFit/config/Analysis2011/";
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
  cmd->setMCTruthResolution(pathToConfig+"Parameters_MCTruthResolution_Summer11_PythiaZ2_L1FastJet_NumPUMay10ReReco_v2.txt",ResolutionFunction::ModifiedNSC);
  
  
  // Particle level imbalance
  //cmd->fitPLI("CMS Simulation ",pathToFitResultsMC+"/ResFit_PtGenAveBins_MCSummer11_V3_"+jetTypeStr+"_L1FastJet.root",ResolutionFunction::ModifiedNSC);
  cmd->setPLI(pathToConfig+"Parameters_PLI_Summer11_PythiaZ2.txt",ResolutionFunction::ModifiedNSC);
  
  // Samples
  TString idData = "Data";
  TString idMC = "CMS Simulation";

  // Summer11 primary result
  cmd->addDataSample(idData,pathToFitResultsData+"/ResFit_PtAveBins_Data_163337-167151_V3_"+jetTypeStr+"_L1FastJet.root");
  cmd->addMCSample(idMC,pathToFitResultsMC+"/ResFit_PtAveBins_MCSummer11_V3_"+jetTypeStr+"_L1FastJet_Nominal.root");
  cmd->addMCSample("JES Down",pathToFitResultsMC+"/ResFit_PtAveBins_MCSummer11_V3_"+jetTypeStr+"_L1FastJet_JESDown.root");
  cmd->addMCSample("JES Up",pathToFitResultsMC+"/ResFit_PtAveBins_MCSummer11_V3_"+jetTypeStr+"_L1FastJet_JESUp.root");
  cmd->addMCSample("PU Up",pathToFitResultsMC+"/ResFit_PtAveBins_MCSummer11_V3_"+jetTypeStr+"_L1FastJet_PUUp.root");
  cmd->addMCSample("Spec Herwig",pathToFitResultsMC+"/ResFit_PtAveBins_MCSummer11_V3_"+jetTypeStr+"_L1FastJet_SpectrumHerwigpp.root");
       
  cmd->addFitResult(nominalFitResType);

  // Systematic uncertainties
//   cmd->addExtrapolationUncertainty(idMC,nominalFitResType,kBlue-4);
//   cmd->addPLIUncertainty(idMC,nominalFitResType,kBlue+2);
//   cmd->addUncertaintyFromVariedSample("PU",1.,idMC,nominalFitResType,"PU Up",kRed+1);
//   cmd->addUncertaintyFromVariedSample("Spectrum",1.,idMC,nominalFitResType,"Spec Herwig",kOrange-1);
//   cmd->addUncertaintyFromVariedSample("JES",1.,idMC,nominalFitResType,"JES Down","JES Up",kRed+3);

  cmd->addExtrapolationUncertainty(idMC,nominalFitResType,kRed+2);
  cmd->addPLIUncertainty(idMC,nominalFitResType,kRed+4);
  cmd->addUncertaintyFromVariedSample("PU",1.,idMC,nominalFitResType,"PU Up",kBlue);
  cmd->addUncertaintyFromVariedSample("Spectrum",1.,idMC,nominalFitResType,"Spec Herwig",kBlue-10);
  cmd->addUncertaintyFromVariedSample("JES",1.,idMC,nominalFitResType,"JES Down","JES Up",kBlue+3);


  // Samples to be compared
  cmd->compareSamples(idData,idMC);
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
