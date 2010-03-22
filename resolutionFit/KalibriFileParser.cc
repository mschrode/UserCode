// $Id: $

#include "KalibriFileParser.h"

#include <iostream>

#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"

namespace resolutionFit {

  //! Parse the ROOT file of name \p fileName.
  //!
  //! The parameter \p verbose defines the level
  //! of verbosity
  //!  - 0: no output
  //!  - 1: some useful information (default)
  //!  - 2: a lot of information for debugging
  // --------------------------------------------
  KalibriFileParser::KalibriFileParser(const TString &fileName, int verbose)
    : verbose_(verbose) {
    // Histograms to be read from file
    hists_["hPtGen"] = 0;
    hists_["hPtDijet"] = 0;
    hists_["hTruthPDF"] = 0;
    hists_["hRespMeas_0"] = 0;
    hists_["hRespFit_0"] = 0;

    // Parse file
    if( parse(fileName) ) exit(-1);

    // Calculate mean pt values
    meanPtGen_ = 0.;
    meanPtDijet_ = 0.;
    meanPdfPtTrue_ = 0.;
    setMeanPt();
  }


  // --------------------------------------------
  KalibriFileParser::~KalibriFileParser() {
    for(std::map<TString,TH1F*>::iterator it = hists_.begin();
	it != hists_.end(); it++) {
      if( it->second ) delete it->second;
    }
  }



  //! The following histograms are accessible via
  //! specification of \p name
  //!  - 'hPtGen': ptGen spectrum
  //!  - 'hPtDijet': ptDijet spectrum
  //!  - 'hTruthPDF': ptTrue pdf assumed by the fit
  //!  - 'hRespMeas_0': response distribution pt / ptGen
  //!  - 'hRespFit_0': fitted resolution
  //! The returned histogram is created newly with
  //! the name \p newName; deletion has to be taken
  //! care of by the calling instance!
  // --------------------------------------------
  TH1F *KalibriFileParser::hist(const TString &name, const TString &newName) const {
    TH1F *h = 0;
    HistIt it = hists_.find(name);
    if( it == hists_.end() ) {
      std::cerr << "WARNING: No histogram with name '" << name << "'" << std::endl;
    } else {
      // This is weird: a simple TH1F::Clone(newName) to get
      // h produces a crash if the next TFile is opened i.e.
      // if a new object of KalibriFileParser is created...
      // (Even though TH1F::SetDirectory(0) was set above.)
      // Who understands ROOT?!
      TString title = ";";
      title += it->second->GetXaxis()->GetTitle();
      title += ";";
      title += it->second->GetYaxis()->GetTitle();
      h = new TH1F(newName,title,
		   it->second->GetNbinsX(),
		   it->second->GetXaxis()->GetXmin(),
		   it->second->GetXaxis()->GetXmax());
      for(int bin = 1; bin <= h->GetNbinsX(); bin++) {
	h->SetBinContent(bin,it->second->GetBinContent(bin));
	if( it->first != "hRespFit_0" && it->first != "hTruthPDF") {
	  h->SetBinError(bin,it->second->GetBinError(bin));
	} else {
	  h->SetBinError(bin,0.);
	}
      }
      h->SetLineColor(it->second->GetLineColor());
      h->SetLineWidth(it->second->GetLineWidth());
    }

    return h;
  }



  //! This method does the actual file parsing i.e.
  //! initializing the attributes such as \p value_
  //! from the numbers stored in the file. It is
  //! intended to be called only once by constructors.
  // --------------------------------------------
  int KalibriFileParser::parse(const TString &fileName) {
    int ioError = 0;
    if( verbose_ == 2 ) std::cout << "Parsing file '" << fileName << "'" << std::endl;
    
    // Opening file    
    if( verbose_ == 2 ) std::cout << "  Opening file... " << std::flush;
    TFile file(fileName,"READ");
    if( file.IsZombie() ) {
      std::cerr << "  ERROR: Error opening file." << std::endl;
      ioError = -1;
    }
    if( verbose_ == 2 ) std::cout << "ok" << std::endl;


    if( !ioError ) {
      // Read fitted values and statistical uncertainties from file
      if( verbose_ == 2 ) std::cout << "  Getting fitted values... " << std::flush;
      TH1D *h = 0;
      file.GetObject("hAbsoluteParameters",h);
      if( !h ) {
	std::cerr << "  ERROR: 'hAbsoluteParameters' not found." << std::endl;
	ioError = -2;
      } else {
	if( verbose_ == 2 ) std::cout << "ok" << std::endl;
	h->SetDirectory(0);
	for(int i = 0; i < h->GetNbinsX(); i++) {
	  values_.push_back(h->GetBinContent(1+i));
	  statUncert_.push_back(h->GetBinError(1+i));
	  if( verbose_ == 2 ) {
	    std::cout << "  Value " << i << ": " << values_.back() << std::flush;
	    std::cout << " +/- " << statUncert_.back() << std::endl;
	  }
	}
      }

      // Read histograms from file
      if( verbose_ == 2 ) std::cout << "  Getting histograms from file... " << std::flush;
      for(std::map<TString,TH1F*>::iterator it = hists_.begin();
	  it != hists_.end(); it++) {
	TH1F *h = 0;
	TString name = it->first;
	file.GetObject(name,h);
	if( !h ) {
	  std::cerr << "  ERROR: '" << name << "' not found." << std::endl;
	  ioError = -2;
	} else {
	  h->SetDirectory(0);
	  it->second = h;
	}
      }
      if( verbose_ == 2 ) std::cout << "ok" << std::endl;
    }

    file.Close();

    return ioError;
  }



  //! This method intended to be called only once
  //! by constructors after the execution of \p parse().
  // --------------------------------------------
  void KalibriFileParser::setMeanPt() {
    if( verbose_ == 2 ) std::cout << "Setting mean pt values" << std::endl;

    for(HistIt it = hists_.begin(); it != hists_.end(); it++) {
      if( it->first == "hPtGen" ) meanPtGen_ = it->second->GetMean();
      else if( it->first == "hPtDijet" ) meanPtDijet_ = it->second->GetMean();
      else if( it->first == "hTruthPDF" ) meanPdfPtTrue_ = it->second->GetMean();
    }
  }
}
