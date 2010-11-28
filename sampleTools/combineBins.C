// $Id: combineBins.C,v 1.1 2010/11/28 13:50:03 mschrode Exp $

//! Combine Kalibri output histograms from different
//! eta, pt, and ptSoft bins


#include <cassert>
#include <iostream>
#include <vector>

#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TString.h"

#include "BinningAdmin.h"
#include "../util/utils.h"
#include "../util/HistOps.h"


std::vector<TString> namesBinnedHists_;
std::vector<TString> namesCombinedHists_;



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
bool getCombinedHists(const TString &fileName, util::HistVec &hists) {
  bool result = true;

  TFile file(fileName,"READ");
  if( file.IsZombie() ) {
    std::cerr << "  ERROR opening file '" << fileName << "'" << std::endl;
    result = false;
  } else {
    for(unsigned int i = 0; i < hists.size(); ++i) {
      TH1 *hNew = 0;
      file.GetObject(hists[i]->GetName(),hNew);
      if( hNew ) {
	hNew->SetDirectory(0);

	double hNewMin = 0;
	double hNewMax = 0;	
	util::HistOps::getBinBorders(hNew,hNewMin,hNewMax);
	double hOldMin = 0;
	double hOldMax = 0;	
	util::HistOps::getBinBorders(hists[i],hOldMin,hOldMax);
	if( hNewMin == hOldMax ) {
	  TH1 *hCombi = new TH1D("tmp",hists[i]->GetTitle(),
				 hists[i]->GetNbinsX()+hNew->GetNbinsX(),hOldMin,hNewMax);
	  hCombi->SetXTitle(hists[i]->GetXaxis()->GetTitle());
	  hCombi->SetYTitle(hists[i]->GetYaxis()->GetTitle());
	  hCombi->SetMarkerStyle(hists[i]->GetMarkerStyle());
	  int combiBin = 1;
	  for(int bin = 1; bin <= hists[i]->GetNbinsX(); ++bin, ++combiBin) {
	    hCombi->SetBinContent(combiBin,hists[i]->GetBinContent(bin));
	    hCombi->SetBinError(combiBin,hists[i]->GetBinError(bin));
	  }
	  for(int bin = 1; bin <= hNew->GetNbinsX(); ++bin, ++combiBin) {
	    hCombi->SetBinContent(combiBin,hNew->GetBinContent(bin));
	    hCombi->SetBinError(combiBin,hNew->GetBinError(bin));
	  }
	  hCombi->SetDirectory(0); // This NEEDS to be done... don't understand ROOT!

	  TString name = hists[i]->GetName();
	  delete hists[i];
	  hists[i] = hCombi;
	  hists[i]->SetName(name);
	} else {
	  std::cerr << "  WARNING: bin borders of '" << hNew->GetName() << "' do not match '" << hists[i]->GetName() << "'" << std::endl;
	  result = false;
	}
      } else {
	std::cerr << "  ERROR getting histogram '" << hists[i]->GetName() << "' from file '" << fileName << "'" << std::endl;
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
  namesCombinedHists_.clear();

  if( type == "ResolutionFit" || type == "default") {
    namesBinnedHists_.push_back("hPtGen_0");
    namesBinnedHists_.push_back("hPtGenJet1_0");
    namesBinnedHists_.push_back("hPtAveAbs_0");
    namesBinnedHists_.push_back("hTruthPDF_0");
    namesBinnedHists_.push_back("hRespMeas_0");
    namesBinnedHists_.push_back("hRespFit_0");
    namesBinnedHists_.push_back("hPtAsym_0");
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
    namesBinnedHists_.push_back("hAbsoluteParameters");

    namesCombinedHists_.push_back("hPtAveAbs");  
  } else {
    std::cerr << "Unknown type '" << type << "'" << std::endl;
    result = false;
  }

  return result;
}




// -------------------------------------------------------------------------------------
void combineBins(const TString &binningConfig, const TString &namePrefix, const TString &type = "default") {
  

  // Prepare parameters
  std::cout << "Preparing parameters" << std::endl;

  Ssiz_t fileNameBegin = namePrefix.Last('/');
  if( fileNameBegin < 0 ) fileNameBegin = 0;
  else fileNameBegin++;
  const TString path = namePrefix(0,fileNameBegin);
  const TString file = namePrefix(fileNameBegin,namePrefix.Length()-fileNameBegin);

  sampleTools::BinningAdmin admin(binningConfig);



  // Prepare histograms
  std::cout << "Preparing histograms" << std::endl;
  setHistNames(type);


  // Loop over bins and combine root files
  std::cout << "Looping over bins" << std::endl;
  for(unsigned int ptSoftBin = 0; ptSoftBin < admin.nPtSoftBins(); ++ptSoftBin) {
    for(unsigned int etaBin = 0; etaBin < 1; ++etaBin) {//admin.nEtaBins(); ++etaBin) {
      std::cout << "  Combining files for ptSoft bin " << ptSoftBin << " and eta bin " << etaBin << "... " << std::flush;
      TFile outFile(file+"_Eta"+util::toTString(etaBin)+"_PtSoft"+util::toTString(ptSoftBin)+".root","RECREATE");

      util::HistVec combinedHists;
      bool goodCombination = true;
      for(unsigned int ptBin = 0; ptBin < 3; ++ptBin) {//admin.nPtBins(etaBin); ++ptBin) {
	util::HistVec binnedHists;
	TString newNameSuffix = "_Eta"+util::toTString(etaBin)+"_Pt"+util::toTString(ptBin);
	TString inFileName = path+file+newNameSuffix+"_PtSoft3/jsResponse.root";
 	if( getBinnedHists(inFileName,namesBinnedHists_,newNameSuffix,binnedHists) ) {
 	  for(util::HistItConst hIt = binnedHists.begin(); hIt != binnedHists.end(); ++hIt) {
 	    outFile.WriteTObject(*hIt);
 	  }
 	}
	if( ptBin == 0 ) {
	  goodCombination = getBinnedHists(inFileName,namesCombinedHists_,"",combinedHists);
	} else if( goodCombination ) {
	  goodCombination = getCombinedHists(inFileName,combinedHists);
	}
      } // End of loop over pt bins

      if( goodCombination ) {
	for(util::HistItConst hIt = combinedHists.begin(); hIt != combinedHists.end(); ++hIt) {
	  outFile.WriteTObject(*hIt);
	}
      }

      outFile.Close();
      std::cout << "ok" << std::endl;
    } // End of loop over eta bins
  } // End of loop over ptSoft bins
}
