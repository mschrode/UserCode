// $Id: $

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "TROOT.h"
#include "TString.h"

#include "BinningAdmin.h"
#include "../util/utils.h"
#include "../util/ConfigParser.h"


//! Write Kalibri config files from skeleton config
//!
//! Per eta, pt, and ptSoft bin, one config file is
//! written. The relevant steering parameters (e.g.
//! cuts on pt) are adjusted to the current bin but
//! other parameters are taken from the skeleton.
void writeConfigFiles(const TString &binningConfig, const TString &skeleton, const TString ntuplePrefix, const TString &outFilePrefix) {

  sampleTools::BinningAdmin binAdmin(binningConfig);
  binAdmin.printBinning();

  // Loop over eta, pt, and ptSoft bins
  for(unsigned int etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {
    gROOT->ProcessLine(".! mkdir "+outFilePrefix+"_Eta"+util::toTString(etaBin));
    for(unsigned int ptBin = 0; ptBin < binAdmin.nPtBins(etaBin); ++ptBin) {
      for(unsigned int ptSoftBin = 0; ptSoftBin < binAdmin.nPtSoftBins(); ++ptSoftBin) {
	// Open input file
	// (need to do this here since unable to return to
	// beginning of stream each iteration...)
	std::ifstream inFile;
	inFile.open(skeleton.Data());
	if( !inFile.is_open() ) {
	  std::cerr << "ERROR opening file '" << skeleton << "'" << std::endl;
	  exit(1);
	}

	// Prepare output file
	TString outFileName = outFilePrefix+"_Eta"+util::toTString(etaBin)+"_Pt"+util::toTString(ptBin)+"_PtSoft"+util::toTString(ptSoftBin)+".cfg";
	std::cout << "Writing file " << outFileName << "... " << std::flush;
	std::ofstream outFile(outFileName.Data());
	
	// Write output
	if( inFile.is_open() && outFile.is_open() ) {
	  std::string line = "";
	  while( !inFile.eof() ) {
	    getline(inFile,line);
	    if( util::ConfigParser::noComment(line) ) {
	      size_t pos;
	      std::string tag = util::ConfigParser::getTag(line,"=",pos);
 	      if( tag == "Et min cut on dijet" )
 		outFile << tag << " = " << binAdmin.ptMin(etaBin,ptBin) << std::endl;
 	      else if( tag == "Et max cut on dijet" )
 		outFile << tag << " = " << binAdmin.ptMax(etaBin,ptBin) << std::endl;
 	      else if( tag == "Eta min cut on jet" )
 		outFile << tag << " = " << binAdmin.etaMin(etaBin) << std::endl;
 	      else if( tag == "Eta max cut on jet" )
 		outFile << tag << " = " << binAdmin.etaMax(etaBin) << std::endl;
	      else if( tag == "Min cut on relative Soft Jet Et" )
		outFile << tag << " = " << binAdmin.ptSoftMin(ptSoftBin) << std::endl;
	      else if( tag == "Max cut on relative Soft Jet Et" )
		outFile << tag << " = " << binAdmin.ptSoftMax(ptSoftBin) << std::endl;
	      else if( tag == "Di-Jet input file" ) 
		outFile << tag << " = " << (ntuplePrefix+"_Eta"+util::toTString(etaBin)+"_Pt"+util::toTString(ptBin)+".root") << std::endl;
	      else if( tag == "create plots" )
		outFile << tag << " = true" << std::endl;
	      else if( tag == "plots save as eps" )
		outFile << tag << " = false" << std::endl;
	      else if( tag == "create parameter error plots" )
		outFile << tag << " = true" << std::endl;
	      else if( tag == "create response plots" )
		outFile << tag << " = true" << std::endl;
	      else
		outFile << line << std::endl;
	    } else {
	      outFile << line << std::endl;
	    }
	  }
	  outFile.close();
	  inFile.close();  
	  std::cout << "done" << std::endl;
	} else {
	  std::cerr << "ERROR opening file '" << outFileName << "' for writing\n";
	  exit(1);
	}
      } // End of loop over PtSoft bins
    } // End of loop over Pt bins

    gROOT->ProcessLine(".! mv "+outFilePrefix+"_Eta"+util::toTString(etaBin)+"*.cfg "+outFilePrefix+"_Eta"+util::toTString(etaBin));

  } // End of loop over Eta bins
}
