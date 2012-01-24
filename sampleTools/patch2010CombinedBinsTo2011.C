// $Id: $

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
	Ssiz_t pos = name.First('_');
	name = name(0,pos);
	h->SetName(name+newNameSuffix);
	hists.push_back(h);

//  	// Fix bug in Kalibri's likelihood definition
//  	if( name.Contains("hAbsoluteParameters") ) {
//  	  for(int bin = 1; bin <= h->GetNbinsX(); ++bin) {
//  	    h->SetBinError(bin,sqrt(2.)*h->GetBinError(bin));
//  	  }
//  	}
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
bool setHistNames() {
  namesBinnedHists_.clear();
  namesBinnedHists_.push_back("hPtGen");
  namesBinnedHists_.push_back("hPtGenJet1");
  namesBinnedHists_.push_back("hPtAveAbs");
  namesBinnedHists_.push_back("hTruthPDF");
  //  namesBinnedHists_.push_back("hUnderlyingTruthPDF");
  namesBinnedHists_.push_back("hRespMeas");
  //  namesBinnedHists_.push_back("hRespMeasAbs");
  namesBinnedHists_.push_back("hRespFit");
  namesBinnedHists_.push_back("hPtAsym");
  namesBinnedHists_.push_back("hPtAbsAsym");
  namesBinnedHists_.push_back("hFitPtAsym");
  namesBinnedHists_.push_back("hPtGenAsym");
  namesBinnedHists_.push_back("hPtJet1");
  namesBinnedHists_.push_back("hPtJet2");
  namesBinnedHists_.push_back("hPtJet3");
  namesBinnedHists_.push_back("hPtJet4");
  namesBinnedHists_.push_back("hPJet3");
  namesBinnedHists_.push_back("hPJet3Rel");
  namesBinnedHists_.push_back("hPJet3GenRel");
  namesBinnedHists_.push_back("hPSJ");
  namesBinnedHists_.push_back("hPSJRel");
  namesBinnedHists_.push_back("hPSJGenRel");
  namesBinnedHists_.push_back("hEta");
  namesBinnedHists_.push_back("hDeltaPhi12");
  namesBinnedHists_.push_back("hDeltaPtJet12"); 
  //  namesBinnedHists_.push_back("hNumVtx"); 
  //  namesBinnedHists_.push_back("hNumPU"); 
  //  namesBinnedHists_.push_back("hWeights"); 
  //  namesBinnedHists_.push_back("hPredPtAsym");
  //  namesBinnedHists_.push_back("hPredPtAve");
//   namesBinnedHists_.push_back("hPredPtJet1");
//   namesBinnedHists_.push_back("hPredPtJet2");
//   namesBinnedHists_.push_back("hPredPtAsymStart");
//   namesBinnedHists_.push_back("hPredPtAveStart");
//   namesBinnedHists_.push_back("hPredPtJet1Start");
//   namesBinnedHists_.push_back("hPredPtJet2Start");
//  namesBinnedHists_.push_back("hPtAveCombined");
  namesBinnedHists_.push_back("hAbsoluteParameters");

  return true;
}




// -------------------------------------------------------------------------------------
void rebin(const TString &binningConfig, const TString &oldFilePrefix, const TString &newFile) {
  

  // Prepare parameters
  std::cout << "Preparing parameters" << std::endl;
  const TString absFileName = util::absolutePath(oldFilePrefix);
  const TString relFileName = util::extractFileName(oldFilePrefix);
  sampleTools::BinningAdmin admin(binningConfig);

  // Prepare histograms
  std::cout << "Preparing histograms" << std::endl;
  setHistNames();

  // The output file containing all histograms of interest
  TFile* outFile  = new TFile(newFile+".root","RECREATE");

  // Loop over bins and combine root files
  std::cout << "Looping over bins" << std::endl;
  for(unsigned int etaBin = 0; etaBin < admin.nEtaBins(); ++etaBin) {
    TDirectory* dirEtaBin = outFile->mkdir("Eta"+util::toTString(etaBin));
      
    for(unsigned int ptBin = 0; ptBin < admin.nPtBins(etaBin); ++ptBin) {
      TDirectory* dirPtBin = dirEtaBin->mkdir("Pt"+util::toTString(ptBin));

      std::cout << "  Combining files for eta bin " << etaBin << " and pt bin " << ptBin << "... " << std::flush;
      for(unsigned int ptSoftBin = 1; ptSoftBin < admin.nPtSoftBins(); ++ptSoftBin) {
	TDirectory* outDir = dirPtBin->mkdir("PtSoft"+util::toTString(ptSoftBin-1));

      	TString newNameSuffix = "_Eta"+util::toTString(etaBin)+"_Pt"+util::toTString(ptBin)+"_PtSoft"+util::toTString(ptSoftBin-1);
 	TString inFileName = absFileName+"_Eta"+util::toTString(etaBin)+"_PtSoft"+util::toTString(ptSoftBin)+".root";
 	util::HistVec binnedHists;
	std::vector<TString> namesBinnedHists;
	for(std::vector<TString>::const_iterator it = namesBinnedHists_.begin();
	    it != namesBinnedHists_.end(); ++it) {
	  namesBinnedHists.push_back((*it)+"_Eta"+util::toTString(etaBin)+"_Pt"+util::toTString(ptBin));
	}
  	if( getBinnedHists(inFileName,namesBinnedHists,newNameSuffix,binnedHists) ) {
  	  for(util::HistItConst hIt = binnedHists.begin(); hIt != binnedHists.end(); ++hIt) {
   	    outDir->WriteTObject(*hIt);
  	  }
  	} else {
	  std::cerr << "Missing input file: Eta "+util::toTString(etaBin)+"   Pt "+util::toTString(ptBin)+"\tPtSoft "+util::toTString(ptSoftBin) << std::endl;
	  exit(0);
	}
      } // End of loop over ptSoft bins
      
      std::cout << "ok" << std::endl;
    } // End of loop over pt bins
  } // End of loop over eta bins
  outFile->Close();
  delete outFile;

  std::cout << "Successfully combined files to " << relFileName+".root" << std::endl;
}



void patch2010CombinedBinsTo2011() {
  TString path = "/afs/naf.desy.de/user/m/mschrode/lustre/Analysis2010/ResolutionFit_Fall10/Note/";
  TString adm = "../resolutionFit/config/Analysis2010/Binning2010/BinningAdmin4.cfg";

  std::vector<TString> oldFileNameKey;
  oldFileNameKey.push_back("Data");
  oldFileNameKey.push_back("MCFall10");
  oldFileNameKey.push_back("MCFall10_JESDown5");
  oldFileNameKey.push_back("MCFall10_JESUp5");
  oldFileNameKey.push_back("MCFall10_SpectrumDown");
  oldFileNameKey.push_back("MCFall10_SpectrumUp");
  
  for(std::vector<TString>::const_iterator it = oldFileNameKey.begin();
      it != oldFileNameKey.end(); ++it) {
    rebin(adm,(path+"ResFitThres_PF_"+*it),("ResFit2010_Pt3paraCuts_PF_"+*it));
  }
}
