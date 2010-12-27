// $Id: run.cc,v 1.39 2010/12/02 14:32:16 mschrode Exp $

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
    //util::StyleSettings::presentationNoTitle();
    util::StyleSettings::paperNoTitle();
    //util::StyleSettings::cms();

    gErrorIgnoreLevel = 1001;        // Do not print ROOT message if eps file has been created
    gSystem->Load("libHistPainter"); // Remove 'Warning in <TClass::TClass>: no dictionary for class TPaletteAxis is available'
    gROOT->ForceStyle();


    resolutionFit::Parameters *par = 0;
    TString outNamePrefix;
    TString inNamePrefix;
    TString specNamePrefix;
    double etaMin;
    double etaMax;
    std::vector<double> ptBinEdges;
    int verbosity = 1;
    if( argc > 3 ) verbosity = atoi(argv[3]);
    
    TString respType = argv[1];
    int etaBin = atoi(argv[2]);
    if( respType == "Gauss" ) {
      bool isData = true;
      inNamePrefix = "~/results/ResolutionFit/Note/ResFitThres_PF_Data";
      outNamePrefix = "Data_PF";

//       inNamePrefix = "~/results/ResolutionFit/Note/ResFitThres_PF_MCFall10";
//       outNamePrefix = "Closure_PF";
      //outNamePrefix = "VariationExtrapolation_PF";
      //outNamePrefix = "VariationPLIUp_PF";
      //outNamePrefix = "VariationJESUp_PF";
      //outNamePrefix = "VariationSpectrumDown_PF";


      std::cout << "Setting up parameters for eta bin " << etaBin << std::endl;
      if( etaBin == 0 ) {
	etaMin = 0;
	etaMax = 1.1;
	inNamePrefix += "_Eta0_PtSoft";
	outNamePrefix += "_Eta0_";

	ptBinEdges.clear();
	ptBinEdges.push_back(43.);
	ptBinEdges.push_back(55.);
	ptBinEdges.push_back(70.);
	ptBinEdges.push_back(85.);
   	ptBinEdges.push_back(100.);
   	ptBinEdges.push_back(115.);
   	ptBinEdges.push_back(130.);
   	ptBinEdges.push_back(150.);
   	ptBinEdges.push_back(170.);
     	ptBinEdges.push_back(190.);
       	ptBinEdges.push_back(220.);
	ptBinEdges.push_back(250.);
	ptBinEdges.push_back(300.);
	ptBinEdges.push_back(350.);
 	ptBinEdges.push_back(400.);
    	ptBinEdges.push_back(500.);
	ptBinEdges.push_back(1000.);

      } else if( etaBin == 1 ) {
	etaMin = 1.1;
	etaMax = 1.7;
	inNamePrefix += "_Eta1_PtSoft";
	outNamePrefix += "_Eta1_";
	ptBinEdges.clear();
	ptBinEdges.push_back(43.);
	ptBinEdges.push_back(55.);
	ptBinEdges.push_back(70.);
	ptBinEdges.push_back(85.);
   	ptBinEdges.push_back(100.);
   	ptBinEdges.push_back(115.);
   	ptBinEdges.push_back(130.);
   	ptBinEdges.push_back(170.);
       	ptBinEdges.push_back(220.);
	ptBinEdges.push_back(300.);
    	ptBinEdges.push_back(500.);
 	ptBinEdges.push_back(1000.);

      } else if( etaBin == 2 ) {
	etaMin = 1.7;
	etaMax = 2.3;
	inNamePrefix += "_Eta2_PtSoft";
	outNamePrefix += "_Eta2_";
	ptBinEdges.clear();
	ptBinEdges.push_back(43.);
	ptBinEdges.push_back(55.);
	ptBinEdges.push_back(70.);
	ptBinEdges.push_back(85.);
   	ptBinEdges.push_back(100.);
   	ptBinEdges.push_back(115.);
   	ptBinEdges.push_back(130.);
   	ptBinEdges.push_back(170.);
       	ptBinEdges.push_back(220.);
 	ptBinEdges.push_back(1000.);

      } else if( etaBin == 3 ) {
	etaMin = 2.3;
	etaMax = 5.0;
	inNamePrefix += "_Eta3_PtSoft";
	outNamePrefix += "_Eta3_";
	ptBinEdges.clear();
	ptBinEdges.push_back(43.);
	ptBinEdges.push_back(55.);
	ptBinEdges.push_back(70.);
	ptBinEdges.push_back(85.);
   	ptBinEdges.push_back(100.);
   	ptBinEdges.push_back(115.);
   	ptBinEdges.push_back(130.);
   	ptBinEdges.push_back(170.);
       	ptBinEdges.push_back(220.);
 	//ptBinEdges.push_back(1000.);

      } else {
	std::cerr << "ERROR: '" << etaBin << "' is not a valid eta bin for '" << respType << "' response.\n";
	exit(1);
      }


      par = new resolutionFit::Parameters(etaMin,etaMax,2.7,inNamePrefix+"1.root",ptBinEdges,outNamePrefix,resolutionFit::ResponseFunction::Gauss,resolutionFit::FitModeMaxLikeFull,resolutionFit::BinPtAve,verbosity);
      //par->fitPtGenAsym(true);
      par->setParPtGenAsym(2.54877,0.149045,0.0109168); // 2.5sigma
      
      if( isData ) {
	par->isData(true);
	par->setLumi(32.7);
      } else {
	par->isData(false);
	par->setLumi(-1.);
      }

      // Combined input files:
      // - fitted hists
      // - pt bin edges
      // - mc truth resolution
      //par->setTrueGaussResPar(3.04407,1.16007,0.0348195); // Spring10 MC truth Eta 0. - 1.3
      //par->setTrueGaussResPar(0.,1.23999,0.0362056); // Fall10 Calo MC truth Eta 0. - 1.1
      //par->setTrueGaussResPar(0.,0.830036,0.0389373); // Fall10 PF MC truth Eta 0. - 1.1

      if( etaBin == 0 )      par->setTrueGaussResPar(-1.18591,0.40573,0.,0.352069);    // PF Eta 0 - 1.1
      else if( etaBin == 1 ) par->setTrueGaussResPar(-1.49592,0.627368,0,0.218365);    // PF Eta 1.1 - 1.7
      else if( etaBin == 2 ) par->setTrueGaussResPar(-1.51996,0.766658,0,0.0228755);   // PF eta 1.7 - 2.3
      else if( etaBin == 3 ) par->setTrueGaussResPar(-0.336561,0.572859,0,0.144683);   // PF Eta 2.3 - 5.0
    
//       if( etaBin == 0 )      par->setTrueGaussResPar(3.8663,0.728714,0.,0.224013);    // Calo Eta 0 - 1.1
//       else if( etaBin == 1 ) par->setTrueGaussResPar(3.96546,0.836469,0.,0.200146);    // Calo Eta 1.1 - 1.7
//       else if( etaBin == 2 ) par->setTrueGaussResPar(3.22132,0.731672,0.,0.159447);    // Calo Eta 1.7 - 2.3
//       else if( etaBin == 3 ) par->setTrueGaussResPar(-3.09025,1.08265,0,-0.0753305);   // Calo Eta 2.3 - 5.0

//       par->addPt3Threshold(resolutionFit::Pt3Rel,0.04,inNamePrefix+"0.root");
//       par->addPt3Threshold(resolutionFit::Pt3Rel,0.06,inNamePrefix+"1.root");
//       par->addPt3Threshold(resolutionFit::Pt3Rel,0.08,inNamePrefix+"2.root");
//       par->addPt3Threshold(resolutionFit::Pt3Rel,0.10,inNamePrefix+"3.root");
//       par->addPt3Threshold(resolutionFit::Pt3Rel,0.12,inNamePrefix+"4.root");
//       par->addPt3Threshold(resolutionFit::Pt3Rel,0.15,inNamePrefix+"5.root");

      //       par->addPt3Threshold(resolutionFit::Pt3Rel,0.02,inNamePrefix+"0.root");
      par->addPt3Threshold(resolutionFit::Pt3Rel,0.04,inNamePrefix+"1.root");
      par->addPt3Threshold(resolutionFit::Pt3Rel,0.06,inNamePrefix+"2.root");
      par->addPt3Threshold(resolutionFit::Pt3Rel,0.08,inNamePrefix+"3.root");
      par->addPt3Threshold(resolutionFit::Pt3Rel,0.10,inNamePrefix+"4.root");
      par->addPt3Threshold(resolutionFit::Pt3Rel,0.12,inNamePrefix+"5.root");
      par->addPt3Threshold(resolutionFit::Pt3Rel,0.15,inNamePrefix+"6.root");


// For systematic variation
      //par->setSystExtrapolation("high");
      //par->setSystScalePli(1.25);


      //par->addFileBaseNameMCStat(inNamePrefix+"Flat_");

      //par->addMCTruthBins(4,30,100,0.05);
      //par->fitExtrapolatedSigma(true);

      //par->fitRatio(true);

      //par->addFileBaseNameMCClosure(inNamePrefix+"10_");

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
     //fit->plotPtAsymmetry();
     fit->plotSpectra();
     //fit->plotAdditionalJetActivity();
     //fit->plotControlDistributions();
     fit->plotMCClosure();
     fit->writeRootOutput();
    
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
