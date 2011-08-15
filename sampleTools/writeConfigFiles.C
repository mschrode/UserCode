// $Id: writeConfigFiles.C,v 1.4 2011/06/23 17:56:15 mschrode Exp $

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "TROOT.h"
#include "TString.h"

#include "BinningAdmin.h"
#include "../util/utils.h"
#include "../util/ConfigParser.h"


TString OUT_FILE_PREFIX_;


TString outputDir(bool isMC, int i) {
  TString out = OUT_FILE_PREFIX_;
  if( isMC) out += "_"+util::toTString(i);

  return out;
}



TString outputDir(bool isMC, unsigned int etaBin, unsigned int ptBin, unsigned int ptSoftBin) {
  int idx = 0;
  if( isMC ) {
    if( etaBin == 0 || etaBin == 1 ) {
      if( ptBin <= 2 ) idx = 0;
      else if( ptBin <= 5 ) idx = 1;
      else if( ptBin <= 7 ) idx = 2;
      else if( ptBin <= 9 ) idx = 3;
      else if( ptBin == 10 ) idx = 4;
      else if( ptBin == 11 ) idx = 5;
      else if( ptBin == 12 ) idx = 6;
      else if( ptBin == 13 ) idx = 7;
      else if( ptBin == 14 ) idx = 8;
      else if( ptBin == 15 && ptSoftBin%2 == 0 ) idx = 9;
      else idx = 10;
    } else if( etaBin == 2 ) {
      if( ptBin <= 4 ) idx = 2;
      else if( ptBin <= 7 ) idx = 3;
      else idx = 4;
    } else if( etaBin == 3 ) {
      if( ptBin <= 4 ) idx = 5;
      else if( ptBin <= 7 ) idx = 6;
      else idx = 7;
    } else if( etaBin == 4 ) {
      idx = 8;
    }
  }

  return outputDir(isMC,idx);
}


//! Write Kalibri config files from skeleton config
//!
//! Per eta, pt, and ptSoft bin, one config file is
//! written. The relevant steering parameters (e.g.
//! cuts on pt) are adjusted to the current bin but
//! other parameters are taken from the skeleton.
void writeConfigFiles(const TString &binningConfig, const TString &skeleton, const TString ntuplePrefix, const TString &spectrumFileName, const TString &outFilePrefix, bool isMC) {
  

  // Create output directories
  OUT_FILE_PREFIX_ = outFilePrefix;
  for(int i = 0; i < (isMC?11:1); ++i) {
    gROOT->ProcessLine(".! mkdir "+outputDir(isMC,i));
    gROOT->ProcessLine(".! mkdir "+outputDir(isMC,i)+"/done");
  }

  sampleTools::BinningAdmin binAdmin(binningConfig);
  binAdmin.printBinning();

  // Loop over eta, pt, and ptSoft bins
  for(unsigned int etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {
    for(unsigned int ptBin = 0; ptBin < binAdmin.nPtBins(etaBin); ++ptBin) {
      for(unsigned int ptSoftBin = 0; ptSoftBin < binAdmin.nPtSoftBins(); ++ptSoftBin) {
	// Open input file
	// (need to do this here since unable to return to
	// beginning of stream each iteration...)
	std::ifstream inFile;
	inFile.open((util::absolutePath(skeleton)).Data());
	if( !inFile.is_open() ) {
	  std::cerr << "ERROR opening file '" << skeleton << "'" << std::endl;
	  exit(1);
	}

	// Prepare i/o files
	TString outFileName = OUT_FILE_PREFIX_+"_Eta"+util::toTString(etaBin)+"_Pt"+util::toTString(ptBin)+"_PtSoft"+util::toTString(ptSoftBin)+".cfg";
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
	      if( tag == "jet spectrum" ) 
		outFile << tag << " = " << spectrumFileName << "; hPtGen_Eta" << util::toTString(etaBin) << "_PtSoft" << util::toTString(ptSoftBin) << std::endl;
 	      else if( tag == "Et min cut on dijet" )
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
	      else if( tag == "Di-Jet weight" ) {
		if( isMC ) outFile << tag << " = -1 # Use weight from ntuple" << std::endl;
		else outFile << tag << " = 1 # Data has event weight 1" << std::endl;
	      } 
	      else if( tag == "Di-Jet weight relative to ntuple weight" ) {
		if( isMC ) outFile << tag << " = " << binAdmin.hltLumi(etaBin,ptBin) << std::endl;
		else outFile << "# " << tag << " = 1." << std::endl;
	      }
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
	gROOT->ProcessLine(".! mv "+outFileName+" "+outputDir(isMC,etaBin,ptBin,ptSoftBin));
      } // End of loop over PtSoft bins
    } // End of loop over Pt bins
  } // End of loop over Eta bins
}
