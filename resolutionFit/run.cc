// $Id: run.cc,v 1.30 2010/08/28 19:34:19 mschrode Exp $

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
      outNamePrefix = "TestExtrapolation2_";
      inNamePrefix = "~/results/ResolutionFit/TestExtrapolation_2/Res_Spring10QCDDiJet_SimpleGauss_PtAveCuts_Eta0_PpRel";
      //      specNamePrefix = "~/results/ResolutionFit/Spring10QCDDiJet_PtGenSpectrum0030-3500_Eta00-10_Pt3";

      //      outNamePrefix = "ResFit_132440-143328_DiJetAve15U_SimpleGauss_PpCuts_";
      //      inNamePrefix = "~/results/ResolutionFit/JetMET_Run2010A-PromptReco-v4_DCSONLY_132440-143328_DiJetAve15U_SimpleGauss_Pt3Rel_Eta00-13/Res_Data_SimpleGauss_PtAveCuts_Eta0_Pt3Rel";

      if( etaBin == 0 ) {
	std::cout << "Setting up parameters for eta bin " << etaBin << std::endl;
	ptBinEdges.clear();
//       	ptBinEdges.push_back(40.);
//   	ptBinEdges.push_back(60.);
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
     
	par = new resolutionFit::Parameters(0.,1.3,2.7,inNamePrefix+"10_",ptBinEdges,2,14,outNamePrefix+"Eta00-13_",resolutionFit::ResponseFunction::Gauss,resolutionFit::FitModeMaxLikeSimple,resolutionFit::BinPtAve,verbosity);
	par->isData(false);
	par->setLumi(-1.);

	
	//par->setTrueGaussResPar(3.249,1.0954,0.0457); // Spring10 QCDDiJet
	//par->setTrueGaussResPar(3.82289,1.05984,0.0379744); // Spring10 QCDDiJet Pt3 Core
	//par->setTrueGaussResPar(3.18079,1.06712,0.0571787); // Spring10 QCDDiJet Pt3 All

	//par->setTrueGaussResPar(3.60458,1.12229,0.0363994); // DeltaPhi
	//par->setTrueGaussResPar(3.14604,1.1096,0.0377418); // DeltaPhi, pt3Gen10, ptSGen10
	//par->setTrueGaussResPar(3.78212,1.07669,0.0380769); // DeltaPhi, pt3RelGen10, ptSGen10
	//par->setTrueGaussResPar(3.74661,1.08124,0.0377255); // DeltaPhi, pt3RelGen10

	par->setTrueGaussResPar(3.04407,1.16007,0.0348195); // MC truth Eta 0. - 1.3
	//par->setTrueGaussResPar(2.94481,0.815236,0.0578183); // MC truth Eta 1.3 - 3

	par->addPt3Threshold(resolutionFit::Pt3Rel,0.04,inNamePrefix+"04_");
	par->addPt3Threshold(resolutionFit::Pt3Rel,0.06,inNamePrefix+"06_");
	par->addPt3Threshold(resolutionFit::Pt3Rel,0.08,inNamePrefix+"08_");
	par->addPt3Threshold(resolutionFit::Pt3Rel,0.10,inNamePrefix+"10_");
	par->addPt3Threshold(resolutionFit::Pt3Rel,0.12,inNamePrefix+"12_");
 	par->addPt3Threshold(resolutionFit::Pt3Rel,0.15,inNamePrefix+"15_");
	//par->addPt3Threshold(resolutionFit::Pt3Rel,0.20,inNamePrefix+"Pt3Rel20_");
	
// 	par->addPt3Threshold(resolutionFit::Pt3Abs,6.,inNamePrefix+"Pt3Gen06_");
// 	par->addPt3Threshold(resolutionFit::Pt3Abs,8.,inNamePrefix+"Pt3Gen08_");
// 	par->addPt3Threshold(resolutionFit::Pt3Abs,10.,inNamePrefix+"Pt3Gen10_");
// 	par->addPt3Threshold(resolutionFit::Pt3Abs,12.,inNamePrefix+"Pt3Gen12_");
//	par->addPt3Threshold(resolutionFit::Pt3Abs,15.,inNamePrefix+"Pt3Gen15_");
// 	par->addPt3Threshold(resolutionFit::Pt3Abs,20.,inNamePrefix+"Pt3Gen20_");
//	par->addPt3Threshold(resolutionFit::Pt3Abs,40.,inNamePrefix+"Pt3Gen40_");

	//par->addFileBaseNameMCStat(inNamePrefix+"Flat_");

	//par->addMCTruthBins(4,30,100,0.05);
	//par->fitExtrapolatedSigma(true);

	//	par->addStartOffset(1.005);
	//par->fitRatio(true);

	//	par->addFileBaseNameMCClosure("~/results/ResolutionFit/MCClosure_CaloOrdered_Rel3rdJet010/Res_Spring10QCDDiJet_Gauss_Eta0_");
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
     fit->plotControlDistributions();
     fit->plotMCClosure();
    
    // Print
    fit->print();
    //    fit->printPoints();
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
