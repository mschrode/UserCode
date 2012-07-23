// $Id: combineBins.C,v 1.8 2012/03/01 19:34:50 mschrode Exp $

//! Combine Kalibri output histograms from different
//! eta, pt, and ptSoft bins


#include <cassert>
#include <iostream>
#include <vector>

#include "TDirectory.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TROOT.h"
#include "TString.h"

#include "BinningAdmin.h"
#define UTILS_AS_HEADER_FILE
#include "../util/utils.h"
#include "../util/HistOps.h"


std::vector<TString> namesBinnedHists_;



// -------------------------------------------------------------------------------------
bool getBinnedHists(const TString &fileName, const std::vector<TString> &histNames, const TString &newNameSuffix, util::HistVec &hists) {
  bool result = true;

  TFile file(fileName,"READ");
  if( file.IsZombie() ) {
    std::cerr << "  ERROR opening file '" << fileName << "'" << std::endl;
    result = false;
  } else {
    for(std::vector<TString>::const_iterator hNameIt = histNames.begin();
	hNameIt != histNames.end(); ++hNameIt) {
      TH1 *h = 0;
      file.GetObject(*hNameIt,h);
      if( h ) {
	h->SetDirectory(0);
	TString name = *hNameIt;
	if( name.Length() > 1 && name(name.Length()-2,name.Length()-1) == "_0" )
	  name = name(0,name.Length()-2);
	h->SetName(name+newNameSuffix);
	hists.push_back(h);
      } else {
	std::cerr << "  ERROR getting histogram '" << *hNameIt << "' from file '" << fileName << "'" << std::endl;
	result = false;
	break;
      }
    }
    file.Close();
  }
    
  return result;
}



// -------------------------------------------------------------------------------------
bool setHistNames(const TString &type) {
  bool result = true;

  namesBinnedHists_.clear();

  if( type == "ResolutionFit" || type == "default") {
    namesBinnedHists_.push_back("hPtGen_0");
    namesBinnedHists_.push_back("hPtGenJet1_0");
    namesBinnedHists_.push_back("hPtAveAbs_0");
    namesBinnedHists_.push_back("hTruthPDF_0");
    namesBinnedHists_.push_back("hUnderlyingTruthPDF_0");
    namesBinnedHists_.push_back("hRespMeas_0");
    namesBinnedHists_.push_back("hRespMeasAbs_0");
    namesBinnedHists_.push_back("hRespFit_0");
    namesBinnedHists_.push_back("hPtAsym_0");
    namesBinnedHists_.push_back("hPtAbsAsym_0");
    namesBinnedHists_.push_back("hFitPtAsym_0");
    namesBinnedHists_.push_back("hPtGenAsym_0");
    namesBinnedHists_.push_back("hPtJet1_0");
    namesBinnedHists_.push_back("hPtJet2_0");
    namesBinnedHists_.push_back("hPtJet3_0");
    namesBinnedHists_.push_back("hPtJet4_0");
    namesBinnedHists_.push_back("hPJet3_0");
    namesBinnedHists_.push_back("hPJet3Rel_0");
    namesBinnedHists_.push_back("hPJet3GenRel_0");
    namesBinnedHists_.push_back("hPSJ_0");
    namesBinnedHists_.push_back("hPSJRel_0");
    namesBinnedHists_.push_back("hPSJGenRel_0");
    namesBinnedHists_.push_back("hEta_0");
    namesBinnedHists_.push_back("hDeltaPhi12_0");
    namesBinnedHists_.push_back("hDeltaPtJet12_0"); 
    namesBinnedHists_.push_back("hNumVtx_0"); 
    namesBinnedHists_.push_back("hNumPU_0"); 
    namesBinnedHists_.push_back("hWeights_0"); 
    namesBinnedHists_.push_back("hPredPtAsym_0");
    namesBinnedHists_.push_back("hPredPtAve_0");
    namesBinnedHists_.push_back("hPredPtJet1_0");
    namesBinnedHists_.push_back("hPredPtJet2_0");
    namesBinnedHists_.push_back("hPredPtAsymStart_0");
    namesBinnedHists_.push_back("hPredPtAveStart_0");
    namesBinnedHists_.push_back("hPredPtJet1Start_0");
    namesBinnedHists_.push_back("hPredPtJet2Start_0");
    namesBinnedHists_.push_back("hPtAveCombined");
    namesBinnedHists_.push_back("hAbsoluteParameters");
  } else {
    std::cerr << "Unknown type '" << type << "'" << std::endl;
    result = false;
  }

  return result;
}




// -------------------------------------------------------------------------------------
void combineBins(const TString &binningConfig, const TString &dirPrefix, const TString &type = "default") {
  

  // Prepare parameters
  std::cout << "Preparing parameters" << std::endl;

  const TString absFileName = util::absolutePath(dirPrefix);
  const TString relFileName = util::fileName(dirPrefix);
  sampleTools::BinningAdmin admin(binningConfig);
  std::vector<TString> missingFiles;

  // Prepare histograms
  std::cout << "Preparing histograms" << std::endl;
  setHistNames(type);

  // The output file containing all histograms of interest
  TFile* outFile  = new TFile(relFileName+".root","RECREATE");

  // Loop over bins and combine root files
  std::cout << "Looping over bins" << std::endl;
  for(unsigned int etaBin = 0; etaBin < admin.nEtaBins(); ++etaBin) {
    //for(unsigned int etaBin = 0; etaBin < 1; ++etaBin) {
    TDirectory* dirEtaBin = outFile->mkdir("Eta"+util::toTString(etaBin));
      
    for(unsigned int ptBin = 0; ptBin < admin.nPtBins(etaBin); ++ptBin) {
      TDirectory* dirPtBin = dirEtaBin->mkdir("Pt"+util::toTString(ptBin));

      std::cout << "  Combining files for eta bin " << etaBin << " and pt bin " << ptBin << "... " << std::flush;
      for(unsigned int ptSoftBin = 0; ptSoftBin < admin.nPtSoftBins(); ++ptSoftBin) {
	TDirectory* outDir = dirPtBin->mkdir("PtSoft"+util::toTString(ptSoftBin));

       	TString newNameSuffix = "_Eta"+util::toTString(etaBin)+"_Pt"+util::toTString(ptBin)+"_PtSoft"+util::toTString(ptSoftBin);
 	TString inFileName = absFileName+newNameSuffix+"/jsResponse.root";
 	util::HistVec binnedHists;
  	if( getBinnedHists(inFileName,namesBinnedHists_,newNameSuffix,binnedHists) ) {
  	  for(util::HistItConst hIt = binnedHists.begin(); hIt != binnedHists.end(); ++hIt) {
   	    outDir->WriteTObject(*hIt);
  	  }
  	} else {
	  missingFiles.push_back("Eta "+util::toTString(etaBin)+"   Pt "+util::toTString(ptBin)+"\tPtSoft "+util::toTString(ptSoftBin));
	}
      } // End of loop over ptSoft bins
      
      std::cout << "ok" << std::endl;
    } // End of loop over pt bins
  } // End of loop over eta bins
  outFile->Close();
  delete outFile;

  if( missingFiles.size() > 0 ) {
    std::cerr << "ERROR: There were missing input files" << std::endl;
    std::cerr << "  Deleting combined file... " << std::flush;
    gROOT->ProcessLine(".! rm "+relFileName+".root");
    std::cerr << "ok" << std::endl;
    std::cerr << "\nThe input files in the following bins were missing:" << std::endl;
    for(std::vector<TString>::const_iterator it = missingFiles.begin();
	it != missingFiles.end(); ++it) {
      std::cerr << "  " << *it << std::endl;
    }
  } else {
    std::cout << "Successfully combined files to " << relFileName+".root" << std::endl;
  }
}
