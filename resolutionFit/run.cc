// $Id: run.cc,v 1.6 2010/05/04 19:20:08 mschrode Exp $

#include <cassert>
#include <iostream>
#include <vector>

#include "TString.h"

#include "PtBin.h"
#include "FittedResolution.h"


int main(int argc, char *argv[]) {

  if( argc > 1 ) {
    TString outNamePrefix;
    TString baseNameStdSel;
    TString baseMCStat;
    std::vector<TString> baseCutVars;
    std::vector<double> ptCuts;
    std::vector<TString> baseSystUp;
    std::vector<TString> baseSystDown;
    std::vector<TString> labelSyst;
    std::vector<double> ptBinEdges;
    int start; 
    int end;
    std::vector<double> trueResPar;


    if( atoi(argv[1]) == 0 ) {

      baseCutVars.clear();
      ptCuts.clear();
      baseSystUp.clear();
      baseSystDown.clear();
      labelSyst.clear();
      ptBinEdges.clear();
      trueResPar.clear();

      outNamePrefix = "ResFit_QCD_Gauss_Eta0_";
      baseNameStdSel = "~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Eta0_";
      baseMCStat = "~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Eta0_Flat_";
 
      baseCutVars.push_back("~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Eta0_Rel3rdJet006_");
      baseCutVars.push_back("~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Eta0_Rel3rdJet008_");
      baseCutVars.push_back("~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Eta0_");
      baseCutVars.push_back("~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Eta0_Rel3rdJet012_");
      baseCutVars.push_back("~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Eta0_Rel3rdJet015_");
      baseCutVars.push_back("~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Eta0_Rel3rdJet020_");

      ptCuts.push_back(0.06);
      ptCuts.push_back(0.08);
      ptCuts.push_back(0.1);
      ptCuts.push_back(0.12);
      ptCuts.push_back(0.15);
      ptCuts.push_back(0.20);

      ptBinEdges.push_back(80.);
      ptBinEdges.push_back(100.);
      ptBinEdges.push_back(120.);
      ptBinEdges.push_back(140.);
      ptBinEdges.push_back(170.);
      ptBinEdges.push_back(200.);
      ptBinEdges.push_back(250.);
      ptBinEdges.push_back(300.);
      ptBinEdges.push_back(400.);
      ptBinEdges.push_back(600.);
      ptBinEdges.push_back(1000.);

      start = 0; 
      end  = 9;

      trueResPar.push_back(1.88);
      trueResPar.push_back(1.205);
      trueResPar.push_back(0.0342);

    } else if( atoi(argv[1]) == 1 ) {

      baseCutVars.clear();
      ptCuts.clear();
      baseSystUp.clear();
      baseSystDown.clear();
      labelSyst.clear();
      ptBinEdges.clear();
      trueResPar.clear();

      outNamePrefix = "ResFit_QCD_Gauss_Eta1_";
      baseNameStdSel = "~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Eta1_";
      baseMCStat = "~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Eta1_Flat_";
 
      baseCutVars.push_back("~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Eta1_Rel3rdJet006_");
      baseCutVars.push_back("~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Eta1_Rel3rdJet008_");
      baseCutVars.push_back("~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Eta1_");
      baseCutVars.push_back("~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Eta1_Rel3rdJet012_");
      baseCutVars.push_back("~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Eta1_Rel3rdJet015_");
      baseCutVars.push_back("~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Eta1_Rel3rdJet020_");

      ptCuts.push_back(0.06);
      ptCuts.push_back(0.08);
      ptCuts.push_back(0.1);
      ptCuts.push_back(0.12);
      ptCuts.push_back(0.15);
      ptCuts.push_back(0.20);

      ptBinEdges.push_back(60.);
      ptBinEdges.push_back(80.);
      ptBinEdges.push_back(100.);
      ptBinEdges.push_back(120.);
      ptBinEdges.push_back(150.);
      ptBinEdges.push_back(200.);
      ptBinEdges.push_back(300.);
      ptBinEdges.push_back(500.);
      ptBinEdges.push_back(700.);

      start = 0; 
      end  = 7;

      trueResPar.push_back(2.51);
      trueResPar.push_back(0.968);
      trueResPar.push_back(0.0483);

    } else {
      std::cerr << "ERROR: '" << argv[1] << "' is not a valid argument.\n";
      exit(1);
    }

    // Check input
    assert( ptBinEdges.size() >= (2+end-start) );
    assert( baseSystUp.size() == baseSystDown.size() );

    // Create pt bins
    std::vector<resolutionFit::PtBin*> ptBins;
    for(int bin = start; bin <= end; bin++) {
      // File names for cut variations
      std::vector<TString> nameCutVars = baseCutVars;
      for(size_t i = 0; i < nameCutVars.size(); i++) {
	nameCutVars[i] += bin;
	nameCutVars[i] += "/jsResponse.root";
      }
      // File names for MC stat uncertainties
      TString nameMCStat = baseMCStat;
      if( nameMCStat != "" ) {
	nameMCStat += bin;
	nameMCStat += "/jsResponse.root";
      }
      // File names for systematic uncertainties
      std::vector<TString> nameSystUp = baseSystUp;
      std::vector<TString> nameSystDown = baseSystDown;
      for(size_t i = 0; i < nameSystUp.size(); i++) {
	nameSystUp[i] += bin;
	nameSystUp[i] += "/jsResponse.root";
	nameSystDown[i] += bin;
	nameSystDown[i] += "/jsResponse.root";
      }
      TString nameStdSel = baseNameStdSel;
      nameStdSel += bin;
      nameStdSel += "/jsResponse.root";
      
      int idx = bin - start;
      ptBins.push_back(new resolutionFit::PtBin(nameStdSel,nameCutVars,ptCuts,nameMCStat,
						nameSystUp,nameSystDown,labelSyst,
						ptBinEdges[idx],ptBinEdges[idx+1]));
    }  
    resolutionFit::FittedResolution *fit = new resolutionFit::FittedResolution(ptBins,trueResPar,outNamePrefix);
    // Plot extrapolation per bin
    fit->plotExtrapolation();
    
    // Plot extrapolated resolution
    fit->plotResolution();
    fit->plotResolutionBins();
    fit->plotSpectra();
    fit->plotSystematicUncertainties();
    
    // Print
    fit->print();

  } else {
    std::cerr << "ERROR: Wrong number of arguments.\n";
    std::cerr << "Usage: 'run <Eta bin number>'\n";
    exit(1);
  }
    

  return 0;
}
