// $Id: run.cc,v 1.63 2012/02/04 21:51:49 mschrode Exp $

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

  //util::StyleSettings::setStylePAS();
  util::StyleSettings::setStyleNoteNoTitle();
  //util::StyleSettings::setStylePresentationNoTitle();
  
  gErrorIgnoreLevel = 1001;
  
  Parameters* par = new Parameters("Res_163337-167151","config/Analysis2011/Binning/BinningAdmin2011_v2.cfg",0);
  //Parameters* par = new Parameters("/afs/naf.desy.de/user/m/mschrode/lustre/RESULT_PLOTS/Resolution_SUS-11-004/Res_163337-167151","config/Analysis2011/Binning/BinningAdmin2011_v2.cfg",0);
  //Parameters* par = new Parameters("ResSUS-11-004","config/Analysis2011/Binning/BinningAdmin2011_v2.cfg",0);
  //Parameters* par = new Parameters("ResNewBinning","config/Analysis2011/Binning/BinningAdmin2011_v3.cfg",0);
  //Parameters* par = new Parameters("PileUpStudy","config/Analysis2011/Binning/BinningAdmin2011_v4.cfg",0);
  par->setNEtaBinsUser(1);
  par->setJetProperties(JetProperties::AK5,JetProperties::PF);
  //par->setOutMode(OutputManager::PSAllInOne);
  par->setOutMode(OutputManager::EPSSingleFiles);
  //par->setOutMode(OutputManager::EPSSingleFilesPlusROOT);
  par->setLumi(838.);
  par->setPtSoftAbsMin(10.);
  
     TString pathToHome = "/afs/naf.desy.de/user/m/mschrode/";
     TString pathToConfig = pathToHome+"UserCode/mschrode/resolutionFit/config/Analysis2011/";

   //   TString pathToFitResultsData = pathToHome+"results/ResolutionFit/Run2011A_V5_163337-167151";
   //   TString pathToFitResultsMC = pathToHome+"results/ResolutionFit/QCD_Pt-15to3000_TuneZ2_Flat_7TeV-pythia6_Summer11-PU_S3_START42_V11-v2_V5";

     //TString pathToFitResultsData = pathToHome+"results/ResolutionFit/PileUpStudy";
     //TString pathToFitResultsMC = pathToHome+"results/ResolutionFit/PileUpStudy/QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Summer11-PU_S3_START42_V11-v2_V5_PU-2011_to_172255_LP_LumiScale";

     // Summ11 results
//      TString pathToFitResultsData = pathToHome+"results/ResolutionFit/Run2011A_163337-167151";
//      TString pathToFitResultsMC = pathToHome+"results/ResolutionFit/QCD_Pt-15to3000_TuneZ2_Flat_7TeV-pythia6_Summer11-PU_S3_START42_V11-v2";
//      // Summer11 control plots (inclusive pt spectra)
     TString pathToFitResultsData = pathToHome+"results/ResolutionFit/NewControlPlots";
     TString pathToFitResultsMC = pathToHome+"results/ResolutionFit/NewControlPlots";


     //FitResult::Type nominalFitResType = FitResult::FullMaxLikeAbs;
     FitResult::Type nominalFitResType = FitResult::PtAsym;


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

   //   cmd->addDataSample("2011A-PromptReco-v4",pathToFitResultsData+"/ResFit_PtAveBinsV4_JetRun2011A-PromptReco-v4_V5_NoBEFilter_"+jetTypeStr+"_L1FastJet.root");
   //   cmd->addDataSample("2011B-PromptReco-v1",pathToFitResultsData+"/ResFit_PtAveBinsV4_JetRun2011B-PromptReco-v1_V5_NoBEFilter_"+jetTypeStr+"_L1FastJet.root");
     //cmd->addMCSample("Summer11 QCDFlat",pathToFitResultsMC+"/ResFit_PtAveBinsV3_MCSummer11_"+jetTypeStr+"_L1FastJet.root");

   //   cmd->addDataSample(idData,pathToFitResultsData+"/ResFit_PtAveBins_Data_163337-167151_"+jetTypeStr+"_L1FastJet_V5.root");
   //   cmd->addMCSample(idMC,pathToFitResultsMC+"/ResFit_PtAveBins_MCSummer11Flat_"+jetTypeStr+"_L1FastJet_V5.root");

     cmd->addDataSample(idData,pathToFitResultsData+"/ResFit_PtAveBins_Data_163337-167151_V3_"+jetTypeStr+"_L1FastJet.root");
     cmd->addMCSample(idMC,pathToFitResultsMC+"/ResFit_PtAveBins_MCSummer11_V3_"+jetTypeStr+"_L1FastJet_Nominal.root");
 //     cmd->addMCSample("JES Down",pathToFitResultsMC+"/ResFit_PtAveBins_MCSummer11_V3_"+jetTypeStr+"_L1FastJet_JESDown.root");
//      cmd->addMCSample("JES Up",pathToFitResultsMC+"/ResFit_PtAveBins_MCSummer11_V3_"+jetTypeStr+"_L1FastJet_JESUp.root");
//      cmd->addMCSample("Spec Herwig",pathToFitResultsMC+"/ResFit_PtAveBins_MCSummer11_V3_"+jetTypeStr+"_L1FastJet_SpectrumHerwigpp.root");
//      cmd->addMCSample("PU Up",pathToFitResultsMC+"/ResFit_PtAveBins_MCSummer11_V3_"+jetTypeStr+"_L1FastJet_PUUp.root");

     cmd->addFitResult(nominalFitResType);

     // Systematic uncertainties
//      cmd->addExtrapolationUncertainty(idMC,nominalFitResType,kCyan+2);
//      //cmd->addMCClosureUncertainty(idMC,nominalFitResType,46);
//      cmd->addUncertaintyFromVariedSample("JES",1.,idMC,nominalFitResType,"JES Down","JES Up",kBlue-9);
//      cmd->addPLIUncertainty(idMC,nominalFitResType,kGreen-1);
//      cmd->addUncertaintyFromVariedSample("Spectrum",1.,idMC,nominalFitResType,"Spec Herwig",kBlue-2);
//      cmd->addUncertaintyFromVariedSample("PU",1.,idMC,nominalFitResType,"PU Up",kBlue+1);

     // Samples to be compared
     cmd->compareSamples(idData,idMC);
     cmd->fitKValues(nominalFitResType);

   //  cmd->compareSamples("2011B-PromptReco-v1","2011A-PromptReco-v4");

     // Output
     cmd->printSetup();
     cmd->makeAllPlots();
     cmd->printResult();
     cmd->printPLI();
     //cmd->printRatioToCombinationFormat();

     delete cmd;
     delete par;

     return 0;


// // ===== 2010 Results ================================================================
//     util::StyleSettings::setStyleNoteNoTitle();

//     gErrorIgnoreLevel = 1001;

//     Parameters* par = new Parameters("Res2010_Pt3paraCuts","config/Analysis2010/Binning2010/BinningAdmin4.cfg",0);
//     //par->setNEtaBinsUser(1);
//     par->setJetProperties(JetProperties::AK5,JetProperties::PF);
//     //par->setOutMode(OutputManager::PSAllInOne);
//     par->setOutMode(OutputManager::EPSSingleFiles);
//     //par->setOutMode(OutputManager::EPSSingleFilesPlusROOT);
//     par->setLumi(36.);
//     //    par->setPtSoftAbsMin(6.);
//     par->setPtSoftAbsMin(1.);

//     TString pathToHome = "/afs/naf.desy.de/user/m/mschrode/";
//     TString pathToConfig = pathToHome+"UserCode/mschrode/resolutionFit/config/Analysis2010/";

// //     TString pathToFitResultsData = pathToHome+"results/ResolutionFit/Analysis2010/Pt3Cuts";
// //     TString pathToFitResultsMC = pathToHome+"results/ResolutionFit/Analysis2010/Pt3Cuts";
//     TString pathToFitResultsData = pathToHome+"results/ResolutionFit/Analysis2010/Pt3paraCuts";
//     TString pathToFitResultsMC = pathToHome+"results/ResolutionFit/Analysis2010/Pt3paraCuts";

//     FitResult::Type nominalFitResType = FitResult::FullMaxLikeRel;

//     CommanderCool* cmd = new CommanderCool(par);

//     // Samples and FitResult types
//     TString jetTypeStr;
//     if( par->jetType() == JetProperties::Calo ) jetTypeStr = "Calo";
//     else if( par->jetType() == JetProperties::JPT ) jetTypeStr = "JPT";
//     else if( par->jetType() == JetProperties::PF )  jetTypeStr = "PF";

//     // MC truth resolution
//     cmd->setMCTruthResolution(pathToConfig+"Parameters_MCTruthResolution.txt",ResolutionFunction::ModifiedNSC);


//     // Particle level imbalance
//     //cmd->fitPLI("PLI",pathToFitResultsMC+"/ResFit2010_PLI_MCFall10.root",ResolutionFunction::ModifiedNSC);
//     cmd->setPLI(pathToConfig+"Parameters_PLI_Fall10.txt",ResolutionFunction::ModifiedNSC);

//     // Samples
//     TString idData = "Data";
//     TString idMC = "CMS Simulation";

// //     cmd->addDataSample(idData,pathToFitResultsData+"/ResFit2010_Pt3Cuts_"+jetTypeStr+"_Data.root");
// //     cmd->addMCSample(idMC,pathToFitResultsMC+"/ResFit2010_Pt3Cuts_"+jetTypeStr+"_MCFall10.root");
// //     cmd->addMCSample("JES Down",pathToFitResultsMC+"/ResFit2010_Pt3Cuts_"+jetTypeStr+"_MCFall10_JESDown.root");
// //     cmd->addMCSample("JES Up",pathToFitResultsMC+"/ResFit2010_Pt3Cuts_"+jetTypeStr+"_MCFall10_JESUp.root");
// //     cmd->addMCSample("Spectrum Down",pathToFitResultsMC+"/ResFit2010_Pt3Cuts_"+jetTypeStr+"_MCFall10_SpectrumDown.root");
// //     cmd->addMCSample("Spectrum Up",pathToFitResultsMC+"/ResFit2010_Pt3Cuts_"+jetTypeStr+"_MCFall10_SpectrumUp.root");

//     cmd->addDataSample(idData,pathToFitResultsData+"/ResFit2010_Pt3paraCuts_"+jetTypeStr+"_Data.root");
//     cmd->addMCSample(idMC,pathToFitResultsMC+"/ResFit2010_Pt3paraCuts_"+jetTypeStr+"_MCFall10.root");
//     cmd->addMCSample("JES Down",pathToFitResultsMC+"/ResFit2010_Pt3paraCuts_"+jetTypeStr+"_MCFall10_JESDown5.root");
//     cmd->addMCSample("JES Up",pathToFitResultsMC+"/ResFit2010_Pt3paraCuts_"+jetTypeStr+"_MCFall10_JESUp5.root");
//     cmd->addMCSample("Spectrum Down",pathToFitResultsMC+"/ResFit2010_Pt3paraCuts_"+jetTypeStr+"_MCFall10_SpectrumDown.root");
//     cmd->addMCSample("Spectrum Up",pathToFitResultsMC+"/ResFit2010_Pt3paraCuts_"+jetTypeStr+"_MCFall10_SpectrumUp.root");

//     cmd->addFitResult(nominalFitResType);

//     // Systematic uncertainties
//     cmd->addExtrapolationUncertainty(idMC,nominalFitResType,kCyan+2);
//     //cmd->addMCClosureUncertainty(idMC,nominalFitResType,46);
//     cmd->addUncertaintyFromVariedSample("JES",1.,idMC,nominalFitResType,"JES Down","JES Up",kBlue-9);
//     cmd->addUncertaintyFromVariedSample("Spectrum",1.,idMC,nominalFitResType,"Spectrum Down","Spectrum Up",kBlue-2);
//     cmd->addPLIUncertainty(idMC,nominalFitResType,kGreen-1);

//     // Samples to be compared
//     cmd->compareSamples(idData,idMC);
//     cmd->fitKValues(nominalFitResType);

//     // Output
//     cmd->printSetup();
//     cmd->makeAllPlots();
//     cmd->printResult();
//     //cmd->printRatioToCombinationFormat();

//     delete cmd;
//     delete par;

    return 0;
}
