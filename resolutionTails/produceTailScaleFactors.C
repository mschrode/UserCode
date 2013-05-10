// $Id: produceTailScaleFactors.C,v 1.2 2013/05/10 13:22:19 mschrode Exp $

#include <iostream>
#include <vector>

#include "TString.h"

#include "CoreScaleFactors.h"
#include "Output.h"
#include "FinalResultProducer.cc"
#include "FitParameters.h"
#include "ScaleFactorProducer.cc"
#include "Style.h"
#include "Uncertainty.h"
#include "TError.h"
#include "TROOT.h"

#include "../sampleTools/BinningAdmin.h"
#define UTILS_AS_HEADER_FILE
#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"


using namespace resolutionTails;

// Main steering method
void produceTailScaleFactors() {
  //// ----- SPECIFY PARAMETER VALUES ------------------------------------------------

  // Fit parameters
  FitParameters* fitPars = new FitParameters();
  fitPars->setNSigTailStart(2.);
  fitPars->setNSigTailEnd(10000.);
  fitPars->setNSigCore(2.);
  fitPars->setFixDataShape(false);
  fitPars->setMinPt3Data(0.);


  // Binning
  sampleTools::BinningAdmin* adm = new sampleTools::BinningAdmin("BinningAdmin2011.cfg");


  // Core scale factors
  CoreScaleFactors csf(CoreScaleFactors::Run2011AB_42X_S6);
  

  // Input files
  const double lumi = 3000;
//   const TString srcData = "../results/Analysis2012/DiJetRun2012ABC-ReReco13JulAB-PromptRecoC_190456-202305/";
//   const TString srcMC   = "../results/Analysis2012/QCD_Pt-15to3000_TuneZ2_Flat_8TeV_pythia6-Summer12_DR53X-PU_S10_START53_V7A-v1/";
  
//   const TString fileNameData    = srcData+"ResTails_PtAveBins_Data2012_PF_L1FastJet_REBINNED.root";
//   const TString fileNameMC      = srcMC  +"ResTails_PtAveBins_MCPythiaSummer12_S10_ReReco13JulAB-PromptRecoC_PF_L1FastJet_REBINNED.root";
//   const TString fileNameMCPUDn  = srcMC  +"";
//   const TString fileNameMCPUUp  = srcMC  +"";
//   const TString fileNameMCTruth = srcMC  +"ResTails_PtAveBins_MCPythiaSummer12_S10_ReReco13JulAB-PromptRecoC_PF_L1FastJet_REBINNED.root";
//   const TString fileNameTailWindow = "";

  // Thesis
  const TString srcData = "../results/Analysis2011/Run2011_163337-180252_V10/";
  const TString srcMC   = "../results/Analysis2011/QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Fall11-PU_S6_START42_V14B-v1_V10/";
  
  const TString fileNameData    = srcData+"ResTails_PtAveBins_Data2011_PF_L1FastJet_V10_REBINNED.root";
  const TString fileNameMC      = srcMC  +"ResTails_PtAveBins_MCFall11_PF_L1FastJet_V10_REBINNED.root";
  const TString fileNameMCPUDn  = srcMC  +"ResTails_PtAveBins_MCFall11_PF_L1FastJet_V10_PUDn_REBINNED.root";
  const TString fileNameMCPUUp  = srcMC  +"ResTails_PtAveBins_MCFall11_PF_L1FastJet_V10_PUUp_REBINNED.root";
  const TString fileNameMCTruth = srcMC  +"ResTails_PtGenAveBins_MCFall11_PF_L1FastJet_V10_REBINNED.root";
  const TString fileNameTailWindow = "";


  // Output labels
  const TString uid          = "TEST";
  const TString jetAlgo      = "PF";


  // Systematic variations
  std::vector<Uncertainty::SystematicVariation> variations;
  variations.push_back(Uncertainty::Nominal);
  variations.push_back(Uncertainty::CoreUp);
  variations.push_back(Uncertainty::CoreDn);
  variations.push_back(Uncertainty::Extrapolation);
  variations.push_back(Uncertainty::Closure);
  variations.push_back(Uncertainty::PUUp);
  variations.push_back(Uncertainty::PUDn);



  // Style parameters
  Style* style = new Style();
  style->showTitle(false);
  style->setTitle("CMS preliminary, L = "+util::StyleSettings::luminosity(lumi)+",  #sqrt{s} = 8 TeV");
  style->setLabelLumi("L = "+util::StyleSettings::luminosity(lumi));
  style->setLabelData(util::LabelFactory::data(util::StyleSettings::luminosity(lumi)));
  style->setLabelMC(util::LabelFactory::mc());
  style->setLabelMCSmear(util::LabelFactory::mc()+" (Corrected #sigma_{A})^{#color[10]{A}}");

  Output::DEBUG = false;


  //// -------------------------------------------------------------------------------




  //// ----- RUN ANALYSIS ------------------------------------------------------------

  // Set style
  if( style->showTitle() ) {
    util::StyleSettings::setStylePAS();
  } else {
    util::StyleSettings::setStylePresentationNoTitle();
  }
  gErrorIgnoreLevel = 1001;

  // Determine scale factors for nominal setup and systematic variations
  TString fileNameNomResult;
  std::vector<Uncertainty::SystematicVariation>::const_iterator vit = variations.begin();
  for(;vit != variations.end(); ++vit) {
    
    std::cout << "\n***** Deriving scale factors (" << Uncertainty::name(*vit) << " variation)" << std::endl;

    // Create output object
    TString outLabel = "Tail_"+uid+"_Sig"+style->nameWindow(fitPars->nSigTailStart(),fitPars->nSigTailEnd())+"_"+jetAlgo;
    if( *vit != Uncertainty::Nominal )
      outLabel += "_"+Uncertainty::id(*vit);
    Output* out = new Output(outLabel,true,true,*vit != Uncertainty::Nominal);
    if( *vit == Uncertainty::Nominal ) fileNameNomResult = outLabel; // store for final results

    // Setup analysis
    ScaleFactorProducer* sfp = new ScaleFactorProducer(fileNameData,fileNameMC,fileNameMCTruth,fileNameTailWindow,fitPars,csf(*vit),adm,out,style);
    sfp->makeScaleFactors((*vit != Uncertainty::Extrapolation),(*vit == Uncertainty::Closure));
    sfp->makeControlPlots();
    sfp->writeWindowBorders();
    sfp->writeExtrapolation();
    if( *vit == Uncertainty::Closure ) sfp->writeMCClosure();

    // Clean up
    delete sfp;
    delete out;
  } // End of loop over variations


  // Combine variations to final result
  FinalResultProducer* frp = new FinalResultProducer(fileNameNomResult,fitPars,adm,variations,style);
  frp->makeScaleFactorPlots();
  frp->makeUncertaintyPlots();
  frp->print();


  // Clean up
  delete frp;
  delete fitPars;
  delete style;
  delete adm;
}
