// $Id: run.cc,v 1.11 2010/05/26 16:43:44 mschrode Exp $

#include <cassert>
#include <iostream>
#include <vector>

#include "TString.h"

#include "Parameters.h"
#include "PtBin.h"
#include "FittedResolution.h"


int main(int argc, char *argv[]) {

  if( argc > 2 ) {
    resolutionFit::Parameters *par = 0;
    TString outNamePrefix;
    TString inNamePrefix;
    std::vector<double> ptBinEdges;
    int verbosity = 1;
    if( argc > 3 ) verbosity = atoi(argv[3]);
    
    TString respType = argv[1];
    int etaBin = atoi(argv[2]);
    if( respType == "Gauss" ) {
      outNamePrefix = "ResFit_Spring10QCDFlat_Gauss_";
      inNamePrefix = "~/results/ResolutionFit/Gauss/resolutionSpring10_Gauss_";

      // outNamePrefix = "ResFit_Spring10QCDFlat_GaussDown30It0_";
      // inNamePrefix = "~/results/ResolutionFit/GaussDown30/Iteration0/resolutionSpring10_Gauss_";

      // outNamePrefix = "ResFit_Spring10QCDFlat_GaussUp30It0_";
      // inNamePrefix = "~/results/ResolutionFit/GaussUp30/Iteration0/resolutionSpring10_Gauss_";
      // outNamePrefix = "ResFit_Spring10QCDFlat_GaussUp30It1_";
      // inNamePrefix = "~/results/ResolutionFit/GaussUp30/Iteration1/resolutionSpring10_Gauss_";

      // outNamePrefix = "ResFit_Spring10QCDFlat_GaussPar1Up30It0_";
      // inNamePrefix = "~/results/ResolutionFit/GaussPar1Up30/Iteration0/resolutionSpring10_Gauss_";


      if( etaBin == 0 ) {
	std::cout << "Setting up parameters for eta bin " << etaBin << std::endl;
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
      
	inNamePrefix += "Eta0_";

	par = new resolutionFit::Parameters(0.,1.2,inNamePrefix,ptBinEdges,0,12,outNamePrefix+"Eta0_",
					    resolutionFit::Parameters::Gauss,verbosity);
	par->setTrueGaussResPar(1.88,1.205,0.0342);
	
	par->addPt3Cut(0.06,inNamePrefix+"Rel3rdJet006_");
	par->addPt3Cut(0.08,inNamePrefix+"Rel3rdJet008_");
	par->addPt3Cut(0.10,inNamePrefix);
	par->addPt3Cut(0.12,inNamePrefix+"Rel3rdJet012_");
	par->addPt3Cut(0.15,inNamePrefix+"Rel3rdJet015_");
	par->addPt3Cut(0.20,inNamePrefix+"Rel3rdJet020_");

	par->addFileBaseNameMCStat(inNamePrefix+"Flat_");

	//par->addMCTruthBins(4,30,100,0.05);
	//par->fitExtrapolatedSigma(true);

	//       par->addStartOffset(0.7);
	par->fitRatio(true);
     
      } else if( etaBin == 1 ) {
	std::cout << "Setting up parameters for eta bin " << etaBin << std::endl;
	ptBinEdges.clear();
	//       ptBinEdges.push_back(60.);
	//       ptBinEdges.push_back(80.);
	ptBinEdges.push_back(100.);
	ptBinEdges.push_back(120.);
	ptBinEdges.push_back(150.);
	ptBinEdges.push_back(200.);
	ptBinEdges.push_back(300.);
	ptBinEdges.push_back(500.);
	ptBinEdges.push_back(700.);

	inNamePrefix += "Eta1_";

	par = new resolutionFit::Parameters(1.2,2.6,inNamePrefix,ptBinEdges,2,7,outNamePrefix+"Eta1_",
					    resolutionFit::Parameters::Gauss,verbosity);
	par->setTrueGaussResPar(2.51,0.968,0.0483);
	
	par->addPt3Cut(0.06,inNamePrefix+"Rel3rdJet006_");
	par->addPt3Cut(0.08,inNamePrefix+"Rel3rdJet008_");
	par->addPt3Cut(0.10,inNamePrefix);
	par->addPt3Cut(0.12,inNamePrefix+"Rel3rdJet012_");
	par->addPt3Cut(0.15,inNamePrefix+"Rel3rdJet015_");
	par->addPt3Cut(0.20,inNamePrefix+"Rel3rdJet020_");

	par->addFileBaseNameMCStat(inNamePrefix+"Flat_");

	//par->addMCTruthBins(4,30,100,0.05);
	par->fitExtrapolatedSigma(true);

      } else {
	std::cerr << "ERROR: '" << etaBin << "' is not a valid eta bin for '" << respType << "' response.\n";
	exit(1);
      }
    } else if( respType == "CrystalBall" ) {

      outNamePrefix = "ResFit_Spring10QCDFlat_CB_";
      inNamePrefix = "~/results/ResolutionFit/CrystalBall/Iteration1/resolutionSpring10_CB_";

      if( etaBin == 0 ) {
	std::cout << "Setting up parameters for eta bin " << etaBin << std::endl;
	ptBinEdges.clear();
	ptBinEdges.push_back(120.);
	ptBinEdges.push_back(140.);
      
	inNamePrefix += "Eta0_";

	par = new resolutionFit::Parameters(0.,1.2,inNamePrefix,ptBinEdges,2,2,outNamePrefix+"Eta0_",
					    resolutionFit::Parameters::CrystalBall,verbosity);
	par->setTrueGaussResPar(1.88,1.205,0.0342);
	
	par->addPt3Cut(0.06,inNamePrefix+"Rel3rdJet006_");
	par->addPt3Cut(0.08,inNamePrefix+"Rel3rdJet008_");
	par->addPt3Cut(0.10,inNamePrefix);
	par->addPt3Cut(0.15,inNamePrefix+"Rel3rdJet015_");
	par->addPt3Cut(0.20,inNamePrefix+"Rel3rdJet020_");
      } else {
	std::cerr << "ERROR: '" << etaBin << "' is not a valid eta bin for '" << respType << "' response.\n";
	exit(1);
      }
    } else {
      std::cerr << "ERROR: '" << respType << "' is not a valid response function.\n";
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

    // Plots
    fit->plotExtrapolation();
    if( respType == "Gauss" ) fit->plotResolution();
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
    std::cerr << "Usage: 'run <Response function> <Eta bin number> [<verbosity>]'\n";
    exit(1);
  }
    

  return 0;
}
