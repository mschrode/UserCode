// $Id: run.cc,v 1.38 2010/11/11 12:57:04 mschrode Exp $

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
      inNamePrefix = "~/results/ResolutionFit/Gauss/Pp3Cuts/ResFit_PF_Data_Eta0_PtSoft";
      outNamePrefix = "Data_PF_";

//        inNamePrefix = "~/results/ResolutionFit/Gauss/Extrapolation7/Res_Fall10_PF/Res_Fall10_LowStat_PF_Eta00-11_PSoft";
//        outNamePrefix = "Res_Fall10_LowStat_PF_";


      if( etaBin == 0 ) {
	std::cout << "Setting up parameters for eta bin " << etaBin << std::endl;
	ptBinEdges.clear();
  	ptBinEdges.push_back(90.);
  	ptBinEdges.push_back(100.);
  	ptBinEdges.push_back(120.);
  	ptBinEdges.push_back(150.);
  	ptBinEdges.push_back(170.);
  	ptBinEdges.push_back(200.);
    	ptBinEdges.push_back(250.);
      	ptBinEdges.push_back(300.);
       	ptBinEdges.push_back(350.);
       	ptBinEdges.push_back(400.);
       	ptBinEdges.push_back(500.);
       	ptBinEdges.push_back(1000.);


// 	ptBinEdges.push_back(350.);
// 	ptBinEdges.push_back(400.);
// 	ptBinEdges.push_back(500.);
// 	ptBinEdges.push_back(600.);
// 	ptBinEdges.push_back(800.);
	
	par = new resolutionFit::Parameters(0.,1.1,2.7,inNamePrefix+"3.root",ptBinEdges,outNamePrefix,resolutionFit::ResponseFunction::Gauss,resolutionFit::FitModeMaxLikeFull,resolutionFit::BinPtAve,verbosity);
	//par->fitPtGenAsym(true);
	par->setParPtGenAsym(2.54877,0.149045,0.0109168); // 2.5sigma

 	par->isData(true);
 	par->setLumi(32.7);
//   	par->isData(false);
//   	par->setLumi(-1.);

	// Combined input files:
	// - fitted hists
	// - pt bin edges
	// - mc truth resolution
	//par->setTrueGaussResPar(3.04407,1.16007,0.0348195); // Spring10 MC truth Eta 0. - 1.3
	//par->setTrueGaussResPar(0.,1.23999,0.0362056); // Fall10 Calo MC truth Eta 0. - 1.1
	par->setTrueGaussResPar(0.,0.830036,0.0389373); // Fall10 PF MC truth Eta 0. - 1.1

  	par->addPt3Threshold(resolutionFit::Pt3Rel,0.04,inNamePrefix+"0.root");
     	par->addPt3Threshold(resolutionFit::Pt3Rel,0.06,inNamePrefix+"1.root");
     	par->addPt3Threshold(resolutionFit::Pt3Rel,0.08,inNamePrefix+"2.root");
     	par->addPt3Threshold(resolutionFit::Pt3Rel,0.10,inNamePrefix+"3.root");
     	par->addPt3Threshold(resolutionFit::Pt3Rel,0.12,inNamePrefix+"4.root");
      	par->addPt3Threshold(resolutionFit::Pt3Rel,0.15,inNamePrefix+"5.root");
  	//par->addPt3Threshold(resolutionFit::Pt3Rel,0.20,inNamePrefix+"20/jsResponse.root");

//  	par->addPt3Threshold(resolutionFit::Pt3Rel,0.08,inNamePrefix+"04.root");
//     	par->addPt3Threshold(resolutionFit::Pt3Rel,0.12,inNamePrefix+"06.root");
//     	par->addPt3Threshold(resolutionFit::Pt3Rel,0.16,inNamePrefix+"08.root");
//     	par->addPt3Threshold(resolutionFit::Pt3Rel,0.20,inNamePrefix+"10.root");
//     	par->addPt3Threshold(resolutionFit::Pt3Rel,0.24,inNamePrefix+"12.root");
//      	par->addPt3Threshold(resolutionFit::Pt3Rel,0.30,inNamePrefix+"15.root");
// 	par->addPt3Threshold(resolutionFit::Pt3Rel,0.60,inNamePrefix+"20.root");


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
