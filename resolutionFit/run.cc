// $Id: run.cc,v 1.23 2010/08/09 12:43:35 mschrode Exp $

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
      outNamePrefix = "ResFit_Spring10QCDDiJet_SimpleGauss_PtGenCuts_PtGenJ3_SJCut_PtGenOrdered_";
      //inNamePrefix = "~/results/ResolutionFit/Spring10QCDDiJet_SimpleGauss_PtGenCuts_PtGenOrdered/Res_Spring10QCDDiJet_SimpleGauss_PtGenCuts_";
      inNamePrefix = "~/results/ResolutionFit/Spring10QCDDiJet_SimpleGauss_PtGenCuts_PtGenJ3_SJCut_PtGenOrdered_Cauchy/Res_Spring10QCDDiJet_SimpleGauss_PtGenCuts_";

      //      specNamePrefix = "~/results/ResolutionFit/Spring10QCDDiJet_PtGenSpectrum0030-3500_Eta00-10_Pt3";


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

	par = new resolutionFit::Parameters(0.,1.,inNamePrefix,ptBinEdges,0,12,outNamePrefix+"Eta0_",resolutionFit::ResponseFunction::Gauss,resolutionFit::FitModeMaxLikeSimple,resolutionFit::RefPtGen,verbosity);
	par->setTrueGaussResPar(3.249,1.0954,0.0457); // Spring10 QCDDiJet
	//par->setTrueGaussResPar(3.82289,1.05984,0.0379744); // Spring10 QCDDiJet Pt3 Core
	//	par->setTrueGaussResPar(3.18079,1.06712,0.0571787); // Spring10 QCDDiJet Pt3 All
	par->setLumi(10);
	
	//par->addPt3Threshold(0.04,inNamePrefix+"Rel3rdJet004_");//,specNamePrefix+"-004.root");
	par->addPt3Threshold(0.06,inNamePrefix+"Rel3rdJet006_");//,specNamePrefix+"-006.root");
	par->addPt3Threshold(0.08,inNamePrefix+"Rel3rdJet008_");//,specNamePrefix+"-008.root");
	par->addPt3Threshold(0.10,inNamePrefix);//,specNamePrefix+"-010.root");
	par->addPt3Threshold(0.12,inNamePrefix+"Rel3rdJet012_");//,specNamePrefix+"-012.root");
	par->addPt3Threshold(0.15,inNamePrefix+"Rel3rdJet015_");//,specNamePrefix+"-015.root");
	par->addPt3Threshold(0.20,inNamePrefix+"Rel3rdJet020_");//,specNamePrefix+"-020.root");

// 	par->addPt3Bin(0.00,0.06,0.03622,inNamePrefix+"Rel3rdJet00-06_",specNamePrefix+"00-06.root");
// 	par->addPt3Bin(0.06,0.08,0.06992,inNamePrefix+"Rel3rdJet06-08_",specNamePrefix+"06-08.root");
// 	par->addPt3Bin(0.08,0.10,0.08988,inNamePrefix+"Rel3rdJet08-10_",specNamePrefix+"08-10.root");
// 	par->addPt3Bin(0.10,0.12,0.10983,inNamePrefix+"Rel3rdJet10-12_",specNamePrefix+"10-12.root");
// 	par->addPt3Bin(0.12,0.15,0.13464,inNamePrefix+"Rel3rdJet12-15_",specNamePrefix+"12-15.root");
// 	par->addPt3Bin(0.15,0.20,0.17396,inNamePrefix+"Rel3rdJet15-20_",specNamePrefix+"15-20.root");

	//par->addFileBaseNameMCStat(inNamePrefix+"Flat_");

	//par->addMCTruthBins(4,30,100,0.05);
	//par->fitExtrapolatedSigma(true);

	//	par->addStartOffset(1.005);
	//par->fitRatio(true);

	//	par->addFileBaseNameMCClosure("~/results/ResolutionFit/MCClosure_CaloOrdered_Rel3rdJet010/Res_Spring10QCDDiJet_Gauss_Eta0_");
	par->addFileBaseNameMCClosure(inNamePrefix);
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
    //fit->plotResolutionBins();
    fit->plotPtAsymmetryBins();
    fit->plotPtAsymmetryAndResponseWidth();
    //    fit->plotPtGenAsymmetry();
    //fit->plotAsymmetrySlopes();
    fit->plotSpectra();
    //fit->plotSystematicUncertainties();
    fit->plotMCClosure();
    //fit->plotCrystalBallTest();
    
    // Print
    fit->print();
    fit->createSlides();

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
