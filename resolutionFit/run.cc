// $Id: run.cc,v 1.5 2010/04/01 11:00:02 mschrode Exp $

#include <cassert>
#include <vector>

#include "TString.h"

#include "PtBin.h"
#include "FittedResolution.h"

int main() {

  TString baseNameStdSel = "~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_";

  std::vector<TString> baseCutVars;
  baseCutVars.push_back("~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Rel3rdPt006_");
  baseCutVars.push_back("~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Rel3rdPt008_");
  baseCutVars.push_back("~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_");
  baseCutVars.push_back("~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Rel3rdPt012_");
  baseCutVars.push_back("~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Rel3rdPt020_");

  std::vector<double> ptCuts;
  ptCuts.push_back(0.06);
  ptCuts.push_back(0.08);
  ptCuts.push_back(0.1);
  ptCuts.push_back(0.12);
//   ptCuts.push_back(0.15);
  ptCuts.push_back(0.20);

  TString baseMCStat = "~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Flat_";

  std::vector<TString> baseSystUp;
//   baseSystUp.push_back("~/results/ResolutionFit/Gauss/PtBins/SigmaUp30/gauss_SigmaUp30_StdSel_");
//   baseSystUp.push_back("~/results/ResolutionFit/Gauss/PtBins/SigmaDown30/gauss_SigmaDown30_StdSel_");
//   baseSystUp.push_back("~/results/ResolutionFit/Gauss/PtBins/SlopeUp50/gauss_SlopeUp50_StdSel_");
//   baseSystUp.push_back("~/results/ResolutionFit/Gauss/PtBins/SlopeDown50/gauss_SlopeDown50_StdSel_");

  std::vector<TString> baseSystDown;
//   baseSystDown.push_back("~/results/ResolutionFit/Gauss/PtBins/SigmaUp30/gauss_SigmaUp30_StdSel_");
//   baseSystDown.push_back("~/results/ResolutionFit/Gauss/PtBins/SigmaDown30/gauss_SigmaDown30_StdSel_");
//   baseSystDown.push_back("~/results/ResolutionFit/Gauss/PtBins/SlopeUp50/gauss_SlopeUp50_StdSel_");
//   baseSystDown.push_back("~/results/ResolutionFit/Gauss/PtBins/SlopeDown50/gauss_SlopeDown50_StdSel_");

  std::vector<TString> labelSyst;
//   labelSyst.push_back("#sigma_{MC} +30%");
//   labelSyst.push_back("#sigma_{MC} -30%");
//   labelSyst.push_back("Spectrum +50%");
//   labelSyst.push_back("Spectrum -50%");

  std::vector<double> ptBinEdges;
  //  ptBinEdges.push_back(80.);
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

  std::vector<resolutionFit::PtBin*> ptBins;
  int start = 0; 
  int end  = 8;

  assert( ptBinEdges.size() >= (2+end-start) );
  assert( baseSystUp.size() == baseSystDown.size() );
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
  std::vector<double> trueResPar;
  trueResPar.push_back(0.);
  trueResPar.push_back(1.15);
  trueResPar.push_back(0.04);
  resolutionFit::FittedResolution *fit = new resolutionFit::FittedResolution(ptBins,trueResPar,"resFit_QCD_Gauss_");
  // Plot extrapolation per bin
  fit->plotExtrapolation();

  // Plot extrapolated resolution
  fit->plotResolution();
  fit->plotResolutionBins();
  fit->plotSpectra();
  fit->plotSystematicUncertainties();

  // Print
  fit->print();
  

  return 0;
}
