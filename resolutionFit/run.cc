// $Id: run.cc,v 1.37 2010/09/22 08:35:12 mschrode Exp $

#ifndef RUN_RESOLUTION_FIT
#define RUN_RESOLUTION_FIT

#include <cassert>
#include <iostream>
#include <vector>

#include "TError.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"

#include "FittedResolution.h"
#include "Parameters.h"
#include "PtBin.h"
#include "ResponseFunction.h"

#include "../util/StyleSettings.h"

int main(int argc, char *argv[]) {

  if( argc > 2 ) {
    // Set style
    util::StyleSettings::presentationNoTitle();
    //util::StyleSettings::paperNoTitle();
    //util::StyleSettings::cms();

    gErrorIgnoreLevel = 1001;        // Do not print ROOT message if eps file has been created
    gSystem->Load("libHistPainter"); // Remove 'Warning in <TClass::TClass>: no dictionary for class TPaletteAxis is available'
    gROOT->ForceStyle();


    resolutionFit::Parameters *par = 0;
    TString outNamePrefix;
    TString inNamePrefix;
    TString specNamePrefix;
    std::vector<double> ptBinEdges;
    int verbosity = 1;
    if( argc > 3 ) verbosity = atoi(argv[3]);
    
    TString respType = argv[1];
    int etaBin = atoi(argv[2]);
    if( respType == "Gauss" ) {
       inNamePrefix = "~/results/ResolutionFit/Gauss/Res_Calo_Eta00-11_Pt200_274E-1pb_PSoft";
       outNamePrefix = "Test_";

      if( etaBin == 0 ) {
	std::cout << "Setting up parameters for eta bin " << etaBin << std::endl;
	ptBinEdges.clear();
	ptBinEdges.push_back(200.);
  	ptBinEdges.push_back(250.);
     	ptBinEdges.push_back(300.);
      	ptBinEdges.push_back(350.);
      	ptBinEdges.push_back(400.);
      	ptBinEdges.push_back(500.);
      	ptBinEdges.push_back(1000.);
	
	par = new resolutionFit::Parameters(0.,1.1,2.7,inNamePrefix+"10/jsResponse.root",ptBinEdges,outNamePrefix+"Eta00-13_",resolutionFit::ResponseFunction::Gauss,resolutionFit::FitModeMaxLikeFull,resolutionFit::BinPtAve,verbosity);
	//par->fitPtGenAsym(true);
	par->setParPtGenAsym(2.54877,0.149045,0.0109168); // 2.5sigma

      	par->isData(true);
      	par->setLumi(2.2);
	//      	par->isData(false);
	//      	par->setLumi(-1.);

	// Combined input files:
	// - fitted hists
	// - pt bin edges
	// - mc truth resolution
	par->setTrueGaussResPar(3.04407,1.16007,0.0348195); // MC truth Eta 0. - 1.1

   	par->addPt3Threshold(resolutionFit::Pp3Rel,0.06,inNamePrefix+"06/jsResponse.root");
   	par->addPt3Threshold(resolutionFit::Pp3Rel,0.08,inNamePrefix+"08/jsResponse.root");
   	par->addPt3Threshold(resolutionFit::Pp3Rel,0.10,inNamePrefix+"10/jsResponse.root");
   	par->addPt3Threshold(resolutionFit::Pp3Rel,0.12,inNamePrefix+"12/jsResponse.root");
    	par->addPt3Threshold(resolutionFit::Pp3Rel,0.15,inNamePrefix+"15/jsResponse.root");

	//par->addFileBaseNameMCStat(inNamePrefix+"Flat_");

	//par->addMCTruthBins(4,30,100,0.05);
	//par->fitExtrapolatedSigma(true);

	//par->fitRatio(true);

	//par->addFileBaseNameMCClosure(inNamePrefix+"10_");

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
    std::cout << "Creating fit (" << std::flush;
    if( par->pt3Bins() ) std::cout << "exclusive pt3 bins)" << std::endl;
    else std::cout << "accumulative pt3 bins)" << std::endl;
    resolutionFit::FittedResolution *fit = new resolutionFit::FittedResolution(ptBins,par);

    // Plots
     fit->plotExtrapolation();
     fit->plotResolution();
     fit->plotPtAsymmetry();
     fit->plotSpectra();
     //fit->plotAdditionalJetActivity();
     //fit->plotControlDistributions();
     fit->plotMCClosure();
    
    // Print
    fit->print();
    //fit->printPoints();
    //fit->createSlides();

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
#endif
