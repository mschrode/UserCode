// $Id: KalibriFileParser.cc,v 1.10 2010/08/24 09:36:43 mschrode Exp $

#include "KalibriFileParser.h"

#include <cmath>
#include <iostream>

#include "TFile.h"
#include "TH1.h"
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
  KalibriFileParser::KalibriFileParser(const TString &fileName, int verbose, bool readFittedValues)
    : verbose_(verbose), readFittedValues_(readFittedValues) {
    // Histograms to be read from file
    hists_["hPtGen"] = 0;
    hists_["hPtGenJet1"] = 0;
    hists_["hPtDijet"] = 0;
    hists_["hTruthPDF"] = 0;
    hists_["hRespMeas_0"] = 0;
    hists_["hRespFit_0"] = 0;
    hists_["hPtAsym_0"] = 0;
    hists_["hFitPtAsym_0"] = 0;
    hists_["hPtGenAsym_0"] = 0;
    hists_["hPtJet1"] = 0;
    hists_["hPtJet2"] = 0;

    // Parse file
    if( parse(fileName) ) exit(-1);

    // Calculate mean pt values
    meanPtGen_ = 0.;
    meanPtDijet_ = 0.;
    meanPdfPtTrue_ = 0.;
    meanPtGenUncert_ = 0.;	
    meanPtDijetUncert_ = 0.;
    meanPdfPtTrueUncert_ = 0.;
    setMeanPt();
  }


  // --------------------------------------------
  KalibriFileParser::~KalibriFileParser() {
    for(std::map<TString,TH1*>::iterator it = hists_.begin();
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
  TH1 *KalibriFileParser::hist(const TString &name, const TString &newName, bool abs) const {
    TH1 *h = 0;
    HistIt it = hists_.find(name);
    if( it == hists_.end() ) {
      std::cerr << "WARNING: No histogram with name '" << name << "'" << std::endl;
    } else {
      // This is weird: a simple TH1::Clone(newName) to get
      // h produces a crash if the next TFile is opened i.e.
      // if a new object of KalibriFileParser is created...
      // (Even though TH1::SetDirectory(0) was set above.)
      // Who understands ROOT?!
      TString title = ";";
      title += it->second->GetXaxis()->GetTitle();
      title += ";";
      title += it->second->GetYaxis()->GetTitle();
      h = new TH1D(newName,title,
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

      if(abs) h->Scale(it->second->GetEntries()*it->second->GetBinWidth(1));
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
      if( readFittedValues_ ) {
	if( verbose_ == 2 ) std::cout << "  Getting fitted values... " << std::flush;
	TH1 *h = 0;
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
      }

      // Read histograms from file
      if( verbose_ == 2 ) std::cout << "  Getting histograms from file... " << std::flush;
      for(std::map<TString,TH1*>::iterator it = hists_.begin();
	  it != hists_.end(); it++) {
	TH1 *h = 0;
	TString name = it->first;
	file.GetObject(name,h);
	if( !h ) {
	  std::cerr << "  ERROR: '" << name << "' not found." << std::endl;
	  ioError = -2;
	} else {
	  h->SetDirectory(0);
	  h->UseCurrentStyle();
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
      if( it->first == "hPtGen" ) {
	meanPtGen_ = it->second->GetMean();
	meanPtGenUncert_ = it->second->GetMeanError();
      } else if( it->first == "hPtDijet" ) {
	meanPtDijet_ = it->second->GetMean();
	meanPtDijetUncert_ = it->second->GetMeanError();
      } else if( it->first == "hTruthPDF" ) {
	meanPdfPtTrue_ = it->second->GetMean();
	meanPdfPtTrueUncert_ = standardDeviation(it->second);
      }
    }

    if( verbose_ == 2 ) {
      std::cout << "  meanPtGen_      =  " << meanPtGen_ << std::endl;
      std::cout << "  meanPtAve_      =  " << meanPtDijet_ << std::endl;
      std::cout << "  meanPdfPtTrue_  =  " << meanPdfPtTrue_ << std::endl;
    }
  }


  
  double KalibriFileParser::standardDeviation(const TH1 *h) const {
    double sigma = 0.;
    double integral = h->Integral("width");
    if( integral ) {
      double meanX = 0.;
      double meanX2 = 0.;
      for(int bin = 1; bin <= h->GetNbinsX(); ++bin) {
	double x = h->GetBinCenter(bin);
	meanX += x*h->GetBinContent(bin);
	meanX2 += x*x*h->GetBinContent(bin);
      }
      meanX *= h->GetBinWidth(1);
      meanX2 *= h->GetBinWidth(1);
      meanX /= integral;
      meanX2 /= integral;

      sigma = sqrt(meanX2 - meanX*meanX);
    }
    return sigma;
  }
}
