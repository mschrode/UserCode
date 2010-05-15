// $Id: run.cc,v 1.8 2010/05/14 09:05:39 mschrode Exp $

#include <cassert>
#include <iostream>
#include <vector>

#include "TString.h"

#include "Parameters.h"
#include "PtBin.h"
#include "FittedResolution.h"


int main(int argc, char *argv[]) {

  if( argc > 1 ) {
    resolutionFit::Parameters *par = 0;
    TString outNamePrefix = "ResFit_Spring10QCDFlat_Gauss_";
    std::vector<double> ptBinEdges;
    int verbosity = 1;
    if( argc > 2 ) verbosity = atoi(argv[2]);

    if( atoi(argv[1]) == 0 ) {
      std::cout << "Setting up parameters for eta bin 0" << std::endl;
      ptBinEdges.clear();
      ptBinEdges.push_back(80.);
      ptBinEdges.push_back(100.);
      ptBinEdges.push_back(120.);
      ptBinEdges.push_back(140.);
      ptBinEdges.push_back(170.);
      ptBinEdges.push_back(200.);
      ptBinEdges.push_back(250.);
      ptBinEdges.push_back(300.);
      ptBinEdges.push_back(350.);
      ptBinEdges.push_back(400.);
      ptBinEdges.push_back(500.);
      ptBinEdges.push_back(600.);
      ptBinEdges.push_back(800.);
      ptBinEdges.push_back(1000.);

      par = new resolutionFit::Parameters(0.,1.2,
					  "~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Eta0_",
					  ptBinEdges,0,12,outNamePrefix+"Eta0_",verbosity);
      par->setTrueGaussResPar(1.88,1.205,0.0342);
	
      par->addPt3Cut(0.06,"~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Eta0_Rel3rdJet006_");
      par->addPt3Cut(0.08,"~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Eta0_Rel3rdJet008_");
      par->addPt3Cut(0.10,"~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Eta0_");
      par->addPt3Cut(0.12,"~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Eta0_Rel3rdJet012_");
      par->addPt3Cut(0.15,"~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Eta0_Rel3rdJet015_");
      par->addPt3Cut(0.20,"~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Eta0_Rel3rdJet020_");

      par->addFileBaseNameMCStat("~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Eta0_Flat_");
      
    } else if( atoi(argv[1]) == 1 ) {
      std::cout << "Setting up parameters for eta bin 1" << std::endl;
      ptBinEdges.clear();
      ptBinEdges.push_back(60.);
      ptBinEdges.push_back(80.);
      ptBinEdges.push_back(100.);
      ptBinEdges.push_back(120.);
      ptBinEdges.push_back(150.);
      ptBinEdges.push_back(200.);
      ptBinEdges.push_back(300.);
      ptBinEdges.push_back(500.);
      ptBinEdges.push_back(700.);

      par = new resolutionFit::Parameters(1.2,2.6,
					  "~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Eta1_",
					  ptBinEdges,0,7,outNamePrefix+"Eta1_",verbosity);
      par->setTrueGaussResPar(2.51,0.968,0.0483);
	
      par->addPt3Cut(0.06,"~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Eta1_Rel3rdJet006_");
      par->addPt3Cut(0.08,"~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Eta1_Rel3rdJet008_");
      par->addPt3Cut(0.10,"~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Eta1_");
      par->addPt3Cut(0.12,"~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Eta1_Rel3rdJet012_");
      par->addPt3Cut(0.15,"~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Eta1_Rel3rdJet015_");
      par->addPt3Cut(0.20,"~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Eta1_Rel3rdJet020_");

      par->addFileBaseNameMCStat("~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Eta1_Flat_");

    } else if( atoi(argv[1]) == 2 ) {
      std::cout << "Setting up parameters for eta bin 2" << std::endl;
      ptBinEdges.clear();
      ptBinEdges.push_back(80.);
      ptBinEdges.push_back(100.);
      ptBinEdges.push_back(150.);

      par = new resolutionFit::Parameters(2.6,3.2,
					  "~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Eta2_",
					  ptBinEdges,1,3,outNamePrefix+"Eta2_",verbosity);
      par->setTrueGaussResPar(1.68,0.757,0.0348);
	
      par->addPt3Cut(0.06,"~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Eta2_Rel3rdJet006_");
      par->addPt3Cut(0.08,"~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Eta2_Rel3rdJet008_");
      par->addPt3Cut(0.10,"~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Eta2_");
      par->addPt3Cut(0.12,"~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Eta2_Rel3rdJet012_");
      par->addPt3Cut(0.15,"~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Eta2_Rel3rdJet015_");
      par->addPt3Cut(0.20,"~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Eta2_Rel3rdJet020_");

      par->addFileBaseNameMCStat("~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_Eta2_Flat_");

    } else {
      std::cerr << "ERROR: '" << argv[1] << "' is not a valid argument.\n";
      exit(1);
    }

    // Create pt bins
    std::cout << "Creating pt bins" << std::endl;
    std::vector<resolutionFit::PtBin*> ptBins;
    for(int i = 0;  i < par->nPtBins(); ++i) {
      ptBins.push_back(new resolutionFit::PtBin(par->createPtBinParameters(i)));
    }  
    std::cout << "Creating fit" << std::endl;
    resolutionFit::FittedResolution *fit = new resolutionFit::FittedResolution(ptBins,par);

    // Plot extrapolation per bin
    fit->plotExtrapolation();
    
    // Plot extrapolated resolution
    fit->plotResolution();
    fit->plotResolutionBins();
    fit->plotPtAsymmetryBins();
    fit->plotSpectra();
    fit->plotSystematicUncertainties();
    
    // Print
    fit->print();

    // Clean up
    std::cout << "Cleaning up" << std::endl;
    delete fit;
    for(std::vector<resolutionFit::PtBin*>::iterator it = ptBins.begin();
	it != ptBins.end(); it++) {
      delete *it;
    }
    delete par;

  } else {
    std::cerr << "ERROR: Wrong number of arguments.\n";
    std::cerr << "Usage: 'run <Eta bin number> [<verbosity>]'\n";
    exit(1);
  }
    

  return 0;
}
