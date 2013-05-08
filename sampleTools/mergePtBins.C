// $Id: $

#include <cassert>
#include <iostream>
#include <vector>

#include "TFile.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TString.h"

#include "BinningAdmin.h"
#include "Parameters.h"
#define UTILS_AS_HEADER_FILE
#include "../util/utils.h"
#include "../util/FileOps.h"
#include "../util/HistOps.h"



namespace sampleTools {
  const TString TMP_ID_DATA = "TMPIDDATA";
  const TString TMP_ID_MC = "TMPIDMC";


  unsigned int readHistsFromFile(unsigned int etaBin, unsigned int ptBin, unsigned int ptSoftBin, const TString &fileName, const std::vector<TString> &histNames, const TString &suffix, std::vector<TH1*> &hists);
  unsigned int setNewNames(unsigned int etaBin, unsigned int ptBin, unsigned int ptSoftBin, const std::vector<TString> &histNames, const TString &suffix, std::vector<TH1*> &hists);
  void writeHistsToFile(const sampleTools::BinningAdmin &binAdmin, const TString &outFileName, const TString &strippedSuffix, std::vector< std::vector< std::vector< std::vector<TH1*> > > > &hists);
  double scaleMCtoData(const std::vector<TH1*> &hData, std::vector<TH1*> &hMC);
  void addHists(std::vector<TH1*> &h1, const std::vector<TH1*> &h2);




  // -------------------------------------------------------------------------------------
  void mergePtBins() {

    std::cout << "Preparing script" << std::endl;

    // These histograms will be merged
    std::vector<TString> histNames;
    histNames.push_back(Parameters::hNamePtAve);  // Scale factor MC/Data is determined from this histogram
    histNames.push_back(Parameters::hNamePtAsymAbs);
    histNames.push_back(Parameters::hNameResp);


    // Load the old and the new binning
    sampleTools::BinningAdmin admOld(Parameters::configAdm,Parameters::configBin);
    sampleTools::BinningAdmin admNew(Parameters::configAdm,Parameters::configBinMerged);
    assert( admOld.nPtSoftBins() == admNew.nPtSoftBins() );
    assert( admOld.nEtaBins() == admNew.nEtaBins() );

    // These will contain the new histograms in (ptsoft,eta,pt) bins
    std::vector< std::vector< std::vector< std::vector<TH1*> > > > hNewData;
    std::vector< std::vector< std::vector< std::vector<TH1*> > > > hNewMC;  

    std::cout << "Combining histograms" << std::endl;

    // Loop over ptSoft bins
    for(unsigned int ptSoftBin = 0; ptSoftBin < admNew.nPtSoftBins(); ++ptSoftBin) {
      std::cout << "  PtSoft bin " << ptSoftBin << std::endl;

      std::vector< std::vector< std::vector<TH1*> > > hNewDataInPtSoftBin;
      std::vector< std::vector< std::vector<TH1*> > > hNewMCInPtSoftBin;  

      // Loop over eta bins
      for(unsigned int etaBin = 0; etaBin < admNew.nEtaBins(); ++etaBin) {
	std::cout << "    Eta bin " << etaBin << std::endl;

	std::vector< std::vector<TH1*> > hNewDataInEtaBin;
	std::vector< std::vector<TH1*> > hNewMCInEtaBin;  

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
	  readHistsFromFile(etaBin,ptBinOrig,ptSoftBin,Parameters::histFileData,histNames,"DATA",hData);
	  readHistsFromFile(etaBin,ptBinOrig,ptSoftBin,Parameters::histFileMC,histNames,"MC",hMC);

	  // Set new names corresponding to this ptBin
	  setNewNames(etaBin,ptBin,ptSoftBin,histNames,TMP_ID_DATA,hData);
	  setNewNames(etaBin,ptBin,ptSoftBin,histNames,TMP_ID_MC,hMC);
	
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
	    readHistsFromFile(etaBin,ptBinOrig,ptSoftBin,Parameters::histFileData,histNames,"DATA_ADDED",hDataAdded);
	    readHistsFromFile(etaBin,ptBinOrig,ptSoftBin,Parameters::histFileMC,histNames,"MC_ADDED",hMCAdded);

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
	  hNewDataInEtaBin.push_back(hData);
	  hNewMCInEtaBin.push_back(hMC);

	} // End of loop over pt bins

	hNewDataInPtSoftBin.push_back(hNewDataInEtaBin);
	hNewMCInPtSoftBin.push_back(hNewMCInEtaBin);

      } // End of loop over eta bins

      hNewData.push_back(hNewDataInPtSoftBin);
      hNewMC.push_back(hNewMCInPtSoftBin);

    } // End of loop over ptSoft bins

    // Write new histograms to file
    TString outFileNameData = util::fileName(Parameters::histFileData);
    outFileNameData.ReplaceAll(".root","_REBINNED.root");
    writeHistsToFile(admNew,outFileNameData,TMP_ID_DATA,hNewData);
    TString outFileNameMC = util::fileName(Parameters::histFileMC);
    outFileNameMC.ReplaceAll(".root","_REBINNED.root");
    writeHistsToFile(admNew,outFileNameMC,TMP_ID_MC,hNewMC);

  
    // Clean up
    for(std::vector< std::vector< std::vector< std::vector<TH1*> > > >::iterator hIt1 = hNewData.begin(); hIt1 != hNewData.end(); ++hIt1) {
      for(std::vector< std::vector< std::vector<TH1*> > >::iterator hIt2 = hIt1->begin(); hIt2 != hIt1->end(); ++hIt2) {
	for(std::vector< std::vector<TH1*> >::iterator hIt3 = hIt2->begin(); hIt3 != hIt2->end(); ++hIt3) {
	  for(std::vector<TH1*>::iterator hIt4 = hIt3->begin(); hIt4 != hIt3->end(); ++hIt4) {
	    delete *hIt4;
	  }
	}
      }
    }
    for(std::vector< std::vector< std::vector< std::vector<TH1*> > > >::iterator hIt1 = hNewMC.begin(); hIt1 != hNewMC.end(); ++hIt1) {
      for(std::vector< std::vector< std::vector<TH1*> > >::iterator hIt2 = hIt1->begin(); hIt2 != hIt1->end(); ++hIt2) {
	for(std::vector< std::vector<TH1*> >::iterator hIt3 = hIt2->begin(); hIt3 != hIt2->end(); ++hIt3) {
	  for(std::vector<TH1*>::iterator hIt4 = hIt3->begin(); hIt4 != hIt3->end(); ++hIt4) {
	    delete *hIt4;
	  }
	}
      }
    }
  }



  // -------------------------------------------------------------------------------------
  unsigned int readHistsFromFile(unsigned int etaBin, unsigned int ptBin, unsigned int ptSoftBin, const TString &fileName, const std::vector<TString> &histNames, const TString &suffix, std::vector<TH1*> &hists) {
    std::vector<TString> fullHistNames;
    TString absFileName = util::absolutePath(fileName);
    for(std::vector<TString>::const_iterator it = histNames.begin();
	it != histNames.end(); ++it) {
      fullHistNames.push_back(absFileName+":///Eta"+util::toTString(etaBin)+"/Pt"+util::toTString(ptBin)+"/PtSoft"+util::toTString(ptSoftBin)+"/"+(*it)+"_Eta"+util::toTString(etaBin)+"_Pt"+util::toTString(ptBin)+"_PtSoft"+util::toTString(ptSoftBin));
    }
    hists = util::FileOps::readHistVec(absFileName,fullHistNames,suffix);

    return fullHistNames.size();
  }


  // -------------------------------------------------------------------------------------
  unsigned int setNewNames(unsigned int etaBin, unsigned int ptBin, unsigned int ptSoftBin, const std::vector<TString> &histNames, const TString &suffix, std::vector<TH1*> &hists) {
    unsigned int n = 0;
    std::vector<TString>::const_iterator nIt = histNames.begin();
    std::vector<TH1*>::iterator hIt = hists.begin();
    for(;nIt != histNames.end(); ++nIt, ++hIt, ++n) {
      TString name = (*nIt)+"_Eta"+util::toTString(etaBin)+"_Pt"+util::toTString(ptBin)+"_PtSoft"+util::toTString(ptSoftBin)+suffix;
      (*hIt)->SetName(name);
    }

    return n;
  }


  // -------------------------------------------------------------------------------------
  void writeHistsToFile(const sampleTools::BinningAdmin &binAdmin, const TString &outFileName, const TString &strippedSuffix, std::vector< std::vector< std::vector< std::vector<TH1*> > > > &hists) {
    TFile* resultFile = new TFile(outFileName,"RECREATE");

    for(unsigned int etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {
      TDirectory* dirEtaBin = resultFile->mkdir("Eta"+util::toTString(etaBin));
      
      for(unsigned int ptBin = 0; ptBin < binAdmin.nPtBins(etaBin); ++ptBin) {
	TDirectory* dirPtBin = dirEtaBin->mkdir("Pt"+util::toTString(ptBin));

	for(unsigned int pt3RelBin = 0; pt3RelBin < binAdmin.nPtSoftBins(); ++pt3RelBin) {
	  TDirectory* outDir = dirPtBin->mkdir("PtSoft"+util::toTString(pt3RelBin));

	  for(std::vector<TH1*>::iterator hIt = hists.at(pt3RelBin).at(etaBin).at(ptBin).begin(); hIt != hists.at(pt3RelBin).at(etaBin).at(ptBin).end(); ++hIt) {
	    TString name = (*hIt)->GetName();
	    name.ReplaceAll(strippedSuffix,"");
	    (*hIt)->SetName(name);
	    outDir->WriteTObject(*hIt);
	  }
	} // End of loop over pt3Rel bins
      } // End of loop over pt bins
    } // End of loop over eta bins

    resultFile->Close();
    delete resultFile;
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
}

