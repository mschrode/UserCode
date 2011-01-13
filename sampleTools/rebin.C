// $Id: $

#include <cassert>
#include <iostream>
#include <vector>

#include "TFile.h"
#include "TH1.h"
#include "TString.h"

#include "BinningAdmin.h"
#include "../util/utils.h"
#include "../util/HistOps.h"
#include "../util/FileOps.h"


void rebin(const TString &inFile, const TString &outFile, const TString &binAdmCfg, const TString &binRebinCfg = "fromentries", bool rescaleWithLumi = false) {
  
  std::cout << "Entering 'rebin'" << std::endl;

  bool verbose = false;
  bool rebinFromEntries = false;

  sampleTools::BinningAdmin* admin = new sampleTools::BinningAdmin(binAdmCfg);
  std::vector<int> nBins;

  sampleTools::Binning* newBinning = 0;
  std::vector<int> nCombBins;

  if( binRebinCfg == "fromentries" ) {
    std::cout << "Rebinning from entries" << std::endl;
    rebinFromEntries = true;
    nBins.push_back(3);
    nBins.push_back(2);
    nBins.push_back(1);
    nBins.push_back(1);
  } else {
    std::cout << "Rebinning from config" << std::endl;
    newBinning = new sampleTools::Binning(binRebinCfg);
  }

  unsigned int nEtaBins = admin->nEtaBins();
  for(unsigned int etaBin = 0; etaBin < nEtaBins; ++etaBin) {

    if( !rebinFromEntries ) {
      nBins.clear();
      assert( newBinning->ptMin(etaBin) == admin->ptMin(etaBin) );
      assert( newBinning->ptMax(etaBin) == admin->ptMax(etaBin) );

      std::cout << "Computing number of bins to be combined  " << std::flush;
      
      unsigned int combBinIdx = 0;
      int numCombBins = 0;
      for(unsigned int i = 0; i < admin->nPtBins(etaBin); ++i) {
	numCombBins++;
	if( newBinning->ptMax(etaBin,combBinIdx) == admin->ptMax(etaBin,i) ) {
	  nCombBins.push_back(numCombBins);
	  numCombBins = 0;
	  combBinIdx++;
	}
      }

      std::cout << "ok" << std::endl;
    }
    unsigned int nPtSoftBins = admin->nPtSoftBins();
    for(unsigned int ptSoftBin = 0; ptSoftBin < nPtSoftBins; ++ptSoftBin) {
      TString fileNameIn = inFile+"_Eta"+util::toTString(etaBin)+"_PtSoft"+util::toTString(ptSoftBin)+".root";
      TString fileNameOut = outFile+"_Eta"+util::toTString(etaBin)+"_PtSoft"+util::toTString(ptSoftBin)+".root";

      // Read histograms from file
      TString histName = "hPtAbsAsym_Eta"+util::toTString(etaBin)+"_Pt";
      //      TString histName = "hRespSymAbs_";
      //TString histName = "hRespMeas_";
      util::HistVec histOrig = util::FileOps::readHistVec(fileNameIn,histName,histName+"tmp");
      if( rescaleWithLumi ) {
	for(unsigned int i = 0; i < histOrig.size(); ++i) {
	  histOrig.at(i)->Scale(admin->hltLumi(etaBin,i));
	  std::cout << "Rescaling with lumi factor " << i << ": " << admin->hltLumi(etaBin,i) << std::endl;
	}
      }

      // Combine hists
      if( rebinFromEntries && !verbose ) std::cout << admin->etaMin(etaBin) << "  " << admin->etaMax(etaBin) << "  " << nBins.at(etaBin)+1 << "    " << admin->ptMin(etaBin,0) << std::flush;

      util::HistVec newHists;
      double nTotalEntries = 0;
      if( rebinFromEntries ) {
 	// Determine what bins to combine
	if( nBins.size() > 1 ) {
	  for(util::HistIt it = histOrig.begin(); it != histOrig.end(); ++it) {
	    nTotalEntries += (*it)->GetEntries();
	  }
	  nTotalEntries /= nBins.at(etaBin);
	  
	  nCombBins.clear();
	  double nEntries = 0;
	  int nComb = 0;
	  if( verbose ) std::cout << "\n\nCombining bins " << std::flush;
	  for(unsigned int i = 0; i < histOrig.size(); ++i) {
	    nEntries += histOrig.at(i)->GetEntries();	  
	    nComb++;
	    if( nEntries >= nTotalEntries ) {
	      if( verbose ) std::cout << "  (" << i << ") " << nComb << std::flush;
	      nCombBins.push_back(nComb);
	      nEntries = 0;
	      nComb = 0;
	    }
	  }
	  if( verbose ) std::cout << "  " << nComb << std::flush;
	  nCombBins.push_back(nComb);
	  if( verbose ) std::cout << std::endl;
	} else {
	  nCombBins.push_back(histOrig.size());
	}
	if( nCombBins.back() == 0 ) {
	  nCombBins.at(nCombBins.size()-2) = nCombBins.at(nCombBins.size()-2)-1;
	  nCombBins.at(nCombBins.size()-1) = 1;
	}
	for(size_t i = 0; i < nCombBins.size(); ++i) {
	  if( verbose ) std::cout << i << ": " << nCombBins.at(i) << std::endl;
	}
	if( verbose ) std::cout << std::endl;
	  
      }

      // Combine hists according to nCombBins
      int nCombHists = 0;
      int nNewHists = 0;
      TH1* hTmp = static_cast<TH1D*>(histOrig.at(0)->Clone("hTmp"));
      hTmp->Reset();
      if( verbose ) std::cout << "Combining pt bins " << std::flush;
      //      std::cout << "Adding entries" << std::endl;
      for(size_t i = 0; i < histOrig.size(); ++i) {
 	if( verbose ) std::cout << "  " << i << std::flush;
 	hTmp->Add(histOrig.at(i));
	//	std::cout << "  " << i << ": " << histOrig.at(i)->GetEntries() << " (" << histOrig.at(i)->Integral() << ")" << std::endl;
 	nCombHists++;
	if( nCombHists == nCombBins.at(nNewHists) ) {
	  if( rebinFromEntries ) std::cout << "  " << admin->ptMax(etaBin,i) << std::flush;
	  newHists.push_back(static_cast<TH1D*>(hTmp->Clone(histName+util::toTString(nNewHists))));
	  //	  std::cout << "  >> " << newHists.back()->GetEntries() << " (" << newHists.back()->Integral() << ")" << std::endl;
	  nNewHists++;
	  nCombHists = 0;
	  hTmp->Reset();
	  if( verbose ) std::cout << "  ok" << std::endl;
	  if( verbose ) std::cout << " and combining pt bins " << std::flush;
	}
      }
      if( rebinFromEntries || verbose ) std::cout << std::endl;
      delete hTmp;
    
      
      // Write combined hists to output file
      TFile oFile(fileNameOut,"RECREATE");
      //TFile oFile(fileNameOut,"UPDATE");
      for(util::HistIt it = newHists.begin(); it != newHists.end(); ++it) {
 	oFile.WriteTObject(*it);
      }   
      oFile.Close();
      for(util::HistIt it = newHists.begin(); it != newHists.end(); ++it) {
 	delete *it;
      }
      for(util::HistIt it = histOrig.begin(); it != histOrig.end(); ++it) {
	delete *it;
      }
    }
  }

  delete admin;
  if( newBinning ) delete newBinning;
}
