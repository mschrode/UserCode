// $Id: mergePtBins.C,v 1.1 2011/06/23 17:49:52 mschrode Exp $

#include <cassert>
#include <iostream>
#include <vector>

#include "TFile.h"
#include "TH1.h"
#include "TString.h"

#include "BinningAdmin.h"
#define UTILS_AS_HEADER_FILE
#include "../util/utils.h"
#include "../util/FileOps.h"
#include "../util/HistOps.h"


const TString TMP_ID_DATA = "TMPIDDATA";
const TString TMP_ID_MC = "TMPIDMC";
std::vector<TString> HIST_NAMES;



// -------------------------------------------------------------------------------------
unsigned int readHistsFromFile(unsigned int etaBin, unsigned int ptBin, unsigned int ptSoftBin, const TString &fileName, const TString &suffix, std::vector<TH1*> &hists) {
  std::vector<TString> histNames;
  TString absFileName = util::absolutePath(fileName);
  for(std::vector<TString>::const_iterator it = HIST_NAMES.begin();
      it != HIST_NAMES.end(); ++it) {
    histNames.push_back(absFileName+":///Eta"+util::toTString(etaBin)+"/Pt"+util::toTString(ptBin)+"/PtSoft"+util::toTString(ptSoftBin)+"/"+(*it)+"_Eta"+util::toTString(etaBin)+"_Pt"+util::toTString(ptBin)+"_PtSoft"+util::toTString(ptSoftBin));
  }
  hists = util::FileOps::readHistVec(absFileName,histNames,suffix);

  return histNames.size();
}


// -------------------------------------------------------------------------------------
unsigned int setNewNames(unsigned int etaBin, unsigned int ptBin, unsigned int ptSoftBin, const TString &suffix, std::vector<TH1*> &hists) {
  unsigned int n = 0;
  std::vector<TString>::const_iterator nIt = HIST_NAMES.begin();
  std::vector<TH1*>::iterator hIt = hists.begin();
  for(;nIt != HIST_NAMES.end(); ++nIt, ++hIt, ++n) {
    TString name = (*nIt)+"_Eta"+util::toTString(etaBin)+"_Pt"+util::toTString(ptBin)+"_PtSoft"+util::toTString(ptSoftBin)+suffix;
    (*hIt)->SetName(name);
  }

  return n;
}


// -------------------------------------------------------------------------------------
void writeHistsToFile(const TString &outFileName, const TString &strippedSuffix, std::vector<TH1*> hists) {
  // Strip suffix from hist names
  for(std::vector<TH1*>::iterator hIt = hists.begin(); hIt != hists.end(); ++hIt) {
    TString name = (*hIt)->GetName();
    name.ReplaceAll(strippedSuffix,"");
    (*hIt)->SetName(name);
  }  
  // Write to file
  TFile f(outFileName,"RECREATE");
  for(std::vector<TH1*>::iterator hIt = hists.begin(); hIt != hists.end(); ++hIt) {
    f.WriteTObject(*hIt);
  }
  f.Close();
}



// Scale histograms in hMC. The scale factor is obtained by requiring
// the first histogram in hMC to have the same integral (area) as the
// first histogram in hData.
// -------------------------------------------------------------------------------------
double scaleMCtoData(const std::vector<TH1*> &hData, std::vector<TH1*> &hMC) {
  double scale = 1.;
  if( hData.size() ) {
    double normData = hData.front()->Integral("width");
    double normMC = hMC.front()->Integral("width");
    if( normMC > 0. ) {
      scale = normData/normMC;
      for(std::vector<TH1*>::iterator it = hMC.begin();
	  it != hMC.end(); ++it) {
	(*it)->Scale(scale);
      }
    }
  }
    
  return scale;
}



// Add histograms in h2 to the corresponding histograms in h1.
// Note: this will increase the number of histogram entries
// (TH1::GetEntries()) of the histograms in h1 by the number of
// entries of the added histograms.
// -------------------------------------------------------------------------------------
void addHists(std::vector<TH1*> &h1, const std::vector<TH1*> &h2) {
  std::vector<TH1*>::iterator h1It = h1.begin();
  std::vector<TH1*>::const_iterator h2It = h2.begin();
  for(; h1It != h1.end(); ++h1It, ++h2It) {
    (*h1It)->Add(*h2It);
  }
}



// -------------------------------------------------------------------------------------
void mergeBins(const TString &fileNameData, const TString &fileNameMC, const TString &cfgBinAdmOld, const TString &cfgBinAdmNew) {

  std::cout << "Preparing script" << std::endl;

  // Names of output files
  TString outFileNameData = util::fileName(fileNameData);
  outFileNameData.ReplaceAll(".root","_REBINNED.root");
  TString outFileNameMC = util::fileName(fileNameMC);
  outFileNameMC.ReplaceAll(".root","_REBINNED.root");


  // These histograms will be merged
  HIST_NAMES.push_back("hPtAveCombined"); // Scale factor MC/Data is determined from first histogram
  HIST_NAMES.push_back("hPtAbsAsym");
  HIST_NAMES.push_back("hRespMeasAbs");


  // Load the old and the new binning
  sampleTools::BinningAdmin admOld(cfgBinAdmOld);
  sampleTools::BinningAdmin admNew(cfgBinAdmNew);
  assert( admOld.nPtSoftBins() == admNew.nPtSoftBins() );
  assert( admOld.nEtaBins() == admNew.nEtaBins() );

  // These will contain the new histograms in this ptSoft bin
  std::vector<TH1*> hNewData;
  std::vector<TH1*> hNewMC;  

  std::cout << "Combining histograms" << std::endl;

  // Loop over ptSoft bins
  for(unsigned int ptSoftBin = 0; ptSoftBin < admNew.nPtSoftBins(); ++ptSoftBin) {
    std::cout << "  PtSoft bin " << ptSoftBin << std::endl;

    // Loop over eta bins
    for(unsigned int etaBin = 0; etaBin < admNew.nEtaBins(); ++etaBin) {
      std::cout << "    Eta bin " << etaBin << std::endl;

      // pt bin index for original binning
      unsigned int ptBinOrig = 0;

      // Loop over new pt bins
      for(unsigned int ptBin = 0; ptBin < admNew.nPtBins(etaBin); ++ptBin) {
	std::cout << "      Processing new pt bin " << ptBin << std::endl;
	std::cout << "        Adding original pt bin " << ptBinOrig << std::endl;

	// Check bin borders
	assert( admNew.ptMin(etaBin,ptBin) == admOld.ptMin(etaBin,ptBinOrig) );
	assert( admNew.ptMax(etaBin,ptBin) >= admOld.ptMax(etaBin,ptBinOrig) );

	// Get original histograms from file
	std::vector<TH1*> hData;
	std::vector<TH1*> hMC;
	readHistsFromFile(etaBin,ptBinOrig,ptSoftBin,fileNameData,"DATA",hData);
	readHistsFromFile(etaBin,ptBinOrig,ptSoftBin,fileNameMC,"MC",hMC);

	// Set new names corresponding to this ptBin
	setNewNames(etaBin,ptBin,ptSoftBin,TMP_ID_DATA,hData);
	setNewNames(etaBin,ptBin,ptSoftBin,TMP_ID_MC,hMC);
	
	// Scale MC to data hists
	scaleMCtoData(hData,hMC);


	// Depending on new binning, add histograms from old binning
	while( admNew.ptMax(etaBin,ptBin) > admOld.ptMax(etaBin,ptBinOrig) ) {
	  // Raise count for original pt bin
	  ++ptBinOrig;
	  std::cout << "        Adding original pt bin " << ptBinOrig << std::endl;
	  assert( admNew.ptMax(etaBin,ptBin) >= admOld.ptMax(etaBin,ptBinOrig) );
	  
	  // Get original histograms from file
	  std::vector<TH1*> hDataAdded;
	  std::vector<TH1*> hMCAdded;
	  readHistsFromFile(etaBin,ptBinOrig,ptSoftBin,fileNameData,"DATA_ADDED",hDataAdded);
	  readHistsFromFile(etaBin,ptBinOrig,ptSoftBin,fileNameMC,"MC_ADDED",hMCAdded);

	  // Scale MC to data hists
	  scaleMCtoData(hDataAdded,hMCAdded);

	  // Add to other histograms in this new pt bin
	  addHists(hData,hDataAdded);
	  addHists(hMC,hMCAdded);	  

	  // Clean up
	  for(std::vector<TH1*>::iterator hIt = hDataAdded.begin(); hIt != hDataAdded.end(); ++hIt) {
	    delete *hIt;
	  }
	  for(std::vector<TH1*>::iterator hIt = hMCAdded.begin(); hIt != hMCAdded.end(); ++hIt) {
	    delete *hIt;
	  }
	}


	// Raise count for original pt bin
	++ptBinOrig;

	// Add to list of new hisograms
	hNewData.insert(hNewData.end(),hData.begin(),hData.end());
	hNewMC.insert(hNewMC.end(),hMC.begin(),hMC.end());

      } // End of loop over pt bins
    } // End of loop over eta bins
  } // End of loop over ptSoft bins

  // Write new histograms to file
  writeHistsToFile(outFileNameData,TMP_ID_DATA,hNewData);
  writeHistsToFile(outFileNameMC,TMP_ID_MC,hNewMC);
  
  // Clean up
  for(std::vector<TH1*>::iterator hIt = hNewData.begin(); hIt != hNewData.end(); ++hIt) {
    delete *hIt;
  }
  hNewData.clear();
  for(std::vector<TH1*>::iterator hIt = hNewMC.begin(); hIt != hNewMC.end(); ++hIt) {
    delete *hIt;
  }
  hNewMC.clear();
}
