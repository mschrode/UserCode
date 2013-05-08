// $Id: $

#include <iostream>
#include <vector>

#include "TString.h"

#include "CoreScaleFactors.h"
#include "Output.h"
#include "FinalResultProducer.cc"
#include "FitParameters.h"
#include "ScaleFactorProducer.cc"
#include "Style.h"
#include "SystematicVariation.h"
#include "TError.h"
#include "TROOT.h"

#include "../sampleTools/BinningAdmin.h"
#define UTILS_AS_HEADER_FILE
#include "../util/LabelFactory.h"
#include "../util/StyleSettings.h"




// Main steering method
void produceTailScaleFactors() {
  //// ----- SPECIFY PARAMETER VALUES ------------------------------------------------

  // Fit parameters
  resolutionTails::FitParameters* fitPars = new resolutionTails::FitParameters();
  fitPars->setNSigTailStart(2.);
  fitPars->setNSigTailEnd(10000.);
  fitPars->setNSigCore(2.);
  fitPars->setFixDataShape(false);
  fitPars->setMinPt3Data(0.);


  // Binning
  sampleTools::BinningAdmin* adm = new sampleTools::BinningAdmin("../resolutionFit/config/Analysis2012/Binning/BinningAdmin2012_v1_mergedPtBins.cfg");


  // Core scale factors
  resolutionTails::CoreScaleFactors csf(resolutionTails::CoreScaleFactors::Run2012ABCReReco53X);
  

  // Input files
  const double lumi = 3000;
  const TString srcData = "../results/Analysis2012/DiJetRun2012ABC-ReReco13JulAB-PromptRecoC_190456-202305/";
  const TString srcMC   = "../results/Analysis2012/QCD_Pt-15to3000_TuneZ2_Flat_8TeV_pythia6-Summer12_DR53X-PU_S10_START53_V7A-v1/";
  
  const TString fileNameData    = srcData+"ResTails_PtAveBins_Data2012_PF_L1FastJet_REBINNED.root";
  const TString fileNameMC      = srcMC  +"ResTails_PtAveBins_MCPythiaSummer12_S10_ReReco13JulAB-PromptRecoC_PF_L1FastJet_REBINNED.root";
  const TString fileNameMCPUDn  = srcMC  +"";
  const TString fileNameMCPUUp  = srcMC  +"";
  const TString fileNameMCTruth = srcMC  +"ResTails_PtAveBins_MCPythiaSummer12_S10_ReReco13JulAB-PromptRecoC_PF_L1FastJet_REBINNED.root";
  const TString fileNameTailWindow = "";


  // Output labels
  const TString uid          = "TEST";
  const TString jetAlgo      = "PF";


  // Style parameters
  resolutionTails::Style* style = new resolutionTails::Style();
  style->showTitle(false);
  style->setTitle("CMS preliminary, L = "+util::StyleSettings::luminosity(lumi)+",  #sqrt{s} = 8 TeV");
  style->setLabelLumi("L = "+util::StyleSettings::luminosity(lumi));
  style->setLabelData(util::LabelFactory::data(util::StyleSettings::luminosity(lumi)));
  style->setLabelMC(util::LabelFactory::mc());
  style->setLabelMCSmear(util::LabelFactory::mc()+" (Corrected #sigma_{A})^{#color[10]{A}}");

  resolutionTails::Output::DEBUG = false;


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
  std::vector<resolutionTails::SystematicVariation*> variations = resolutionTails::SystematicVariation::create();
  // Loop over systematic variations (incl. nominal analysis)
  std::vector<resolutionTails::SystematicVariation*>::const_iterator vit = variations.begin();
  for(;vit != variations.end(); ++vit) {
    
    std::cout << "\n***** Deriving scale factors (" << (*vit)->name() << " variation)" << std::endl;

    // Create output object
    TString outLabel = "Tail_"+uid+"_Sig"+style->nameWindow(fitPars->nSigTailStart(),fitPars->nSigTailEnd())+"_"+jetAlgo;
    if( !((*vit)->isNominal()) ) outLabel += "_"+(*vit)->label();
    resolutionTails::Output* out = new resolutionTails::Output(outLabel,true,true,!(*vit)->isNominal());
    if( (*vit)->isNominal() ) fileNameNomResult = outLabel; // store for final results

    // Setup analysis
    resolutionTails::ScaleFactorProducer* sfp = new resolutionTails::ScaleFactorProducer(fileNameData,fileNameMC,fileNameMCTruth,fileNameTailWindow,fitPars,csf(*vit),adm,out,style);
    sfp->makeScaleFactors(!((*vit)->isExtrapolation()),(*vit)->isClosure());
    sfp->makeControlPlots();
    sfp->writeWindowBorders();
    sfp->writeExtrapolation();
    if( (*vit)->isClosure() ) sfp->writeMCClosure();

    // Clean up
    delete sfp;
    delete out;
  } // End of loop over variations


  // Combine variations to final result
  resolutionTails::FinalResultProducer* frp = new resolutionTails::FinalResultProducer(fileNameNomResult,fitPars,adm,variations,style);
  frp->makeScaleFactorPlots();
  frp->makeUncertaintyPlots();
  frp->print();


  // Clean up
  vit = variations.begin();
  for(;vit != variations.end(); ++vit) {
    delete *vit;
  }
  delete frp;
  delete fitPars;
  delete style;
  delete adm;
}
