// $Id: KalibriFileParser.cc,v 1.16 2010/12/02 14:32:16 mschrode Exp $

#include "KalibriFileParser.h"

#include <cmath>
#include <iostream>

#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TString.h"

#include "../util/utils.h"

namespace resolutionFit {

  //! Parse the ROOT file of name \p fileName.
  //!
  //! The parameter \p verbose defines the level
  //! of verbosity
  //!  - 0: no output
  //!  - 1: some useful information (default)
  //!  - 2: a lot of information for debugging
  // --------------------------------------------
  KalibriFileParser::KalibriFileParser(const TString &fileName, unsigned int ptBin, int verbose, bool readFittedValues)
    : binId_("_Eta2_Pt"+util::toTString(ptBin)), verbose_(verbose), readFittedValues_(readFittedValues) {

    // Histograms to be read from file
    hists_[("hPtGen"+binId_)] = 0;
    hists_[("hPtGenJet1"+binId_)] = 0;
    hists_[("hPtAveAbs"+binId_)] = 0;
    hists_[("hTruthPDF"+binId_)] = 0;
    hists_[("hRespMeas"+binId_)] = 0;
    hists_[("hRespFit"+binId_)] = 0;
    hists_[("hPtAsym"+binId_)] = 0;
    hists_[("hFitPtAsym"+binId_)] = 0;
    hists_[("hPtGenAsym"+binId_)] = 0;
    hists_[("hPtJet1"+binId_)] = 0;
    hists_[("hPtJet2"+binId_)] = 0;
    hists_[("hPtJet3"+binId_)] = 0;
    hists_[("hPtJet4"+binId_)] = 0;
    hists_[("hPJet3"+binId_)] = 0;
    hists_[("hPJet3Rel"+binId_)] = 0;
    hists_[("hPJet3GenRel"+binId_)] = 0;
    hists_[("hPSJ"+binId_)] = 0;
    hists_[("hPSJRel"+binId_)] = 0;
    hists_[("hPSJGenRel"+binId_)] = 0;
    hists_[("hEta"+binId_)] = 0;
    hists_[("hDeltaPhi12"+binId_)] = 0;
    hists_[("hDeltaPtJet12"+binId_)] = 0;


    // Parse file
    if( parse(fileName,ptBin) ) exit(-1);

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

    TString histName = name;
    histName += binId_;

    TH1 *h = 0;
    HistIt it = hists_.find(histName);
    if( it == hists_.end() ) {
      std::cerr << "WARNING (KalibriFileParser): No histogram with name '" << histName << "'" << std::endl;
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
	if( (it->first).Contains("hRespFit_") || (it->first).Contains("hTruthPDF_") ) {
	  h->SetBinError(bin,0.);
	} else {
	  h->SetBinError(bin,it->second->GetBinError(bin));
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
  int KalibriFileParser::parse(const TString &fileName, unsigned int ptBin) {
    int ioError = 0;
    if( verbose_ == 2 ) std::cout << "(KalibriFileParser) Parsing file '" << fileName << "'" << std::endl;
    
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
	file.GetObject("hAbsoluteParameters"+binId_,h);
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
	      std::cout << "  Value: " << values_.back() << std::flush;
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
	  std::cerr << "  WARNING: '" << name << "' not found." << std::endl;
	  it->second = new TH1D(name,"",1,0,1);
	} else {
	  h->SetDirectory(0);
	  h->UseCurrentStyle();
	  it->second = h;
	  if( verbose_ == 2 ) std::cout << "\n   (h->GetName() == '" << h->GetName() << "')... " << std::endl;
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
      if( it->first == ("hPtGen"+binId_) ) {
	meanPtGen_ = it->second->GetMean();
	meanPtGenUncert_ = it->second->GetMeanError();
      } else if( it->first == ("hPtAveAbs"+binId_) ) {
	meanPtDijet_ = it->second->GetMean();
	meanPtDijetUncert_ = it->second->GetMeanError();
      } else if( it->first == ("hTruthPDF"+binId_) ) {
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
