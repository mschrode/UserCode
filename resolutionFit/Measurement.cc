// $Id: Measurement.cc,v 1.2 2011/02/28 10:53:15 mschrode Exp $

#include "Measurement.h"

#include <cmath>
#include <iostream>

#include "TFile.h"
#include "TH1.h"


namespace resolutionFit {

  unsigned int Measurement::HIST_COUNT = 0;


  // -------------------------------------------------------------------
  Measurement::Measurement(const TString &fileName, const TString &histNameSuffix, double ptSoft, unsigned int verbosity)
    : histNameSuffix_(histNameSuffix), ptSoft_(ptSoft), verbosity_(verbosity) {
    init(fileName,true);
  }

  // -------------------------------------------------------------------
  Measurement::Measurement(const TString &fileName, const TString &histNameSuffix, double ptSoft, bool hasFittedParameters, unsigned int verbosity)
    : histNameSuffix_(histNameSuffix), ptSoft_(ptSoft), verbosity_(verbosity) {
    init(fileName,hasFittedParameters);
  }
  

  // -------------------------------------------------------------------
  void Measurement::init(const TString &fileName, bool hasFittedParameters) {

    // Histograms to be read from file
    hists_[("hPtGen"+histNameSuffix_)] = 0;
    hists_[("hPtGenJet1"+histNameSuffix_)] = 0;
    hists_[("hPtAveAbs"+histNameSuffix_)] = 0;
    hists_[("hTruthPDF"+histNameSuffix_)] = 0;
    hists_[("hPtAsym"+histNameSuffix_)] = 0;
    hists_[("hFitPtAsym"+histNameSuffix_)] = 0;
    hists_[("hRespMeas"+histNameSuffix_)] = 0;
    hists_[("hRespFit"+histNameSuffix_)] = 0;
    hists_[("hPtGenAsym"+histNameSuffix_)] = 0;
    hists_[("hPtJet1"+histNameSuffix_)] = 0;
    hists_[("hPtJet2"+histNameSuffix_)] = 0;
    hists_[("hPtJet3"+histNameSuffix_)] = 0;
    hists_[("hPtJet4"+histNameSuffix_)] = 0;
    hists_[("hPJet3"+histNameSuffix_)] = 0;
    hists_[("hPJet3Rel"+histNameSuffix_)] = 0;
    hists_[("hPJet3GenRel"+histNameSuffix_)] = 0;
    hists_[("hPSJ"+histNameSuffix_)] = 0;
    hists_[("hPSJRel"+histNameSuffix_)] = 0;
    hists_[("hPSJGenRel"+histNameSuffix_)] = 0;
    hists_[("hEta"+histNameSuffix_)] = 0;
    hists_[("hDeltaPhi12"+histNameSuffix_)] = 0;
    hists_[("hDeltaPtJet12"+histNameSuffix_)] = 0;
    
    // Parse file
    if( !parse(fileName,hasFittedParameters) ) exit(-1);

    // Calculate mean pt values
    if( verbosity_ == 2 ) std::cout << "Setting mean pt values" << std::endl;
    setMean("hPtAveAbs"+histNameSuffix_,meanPtAve_,meanPtAveUncert_);
    setMean("hPtGen"+histNameSuffix_,meanPtGen_,meanPtGenUncert_);
    setMean("hTruthPDF"+histNameSuffix_,meanPdfPtTrue_,meanPdfPtTrueUncert_);
    if( verbosity_ == 2 ) {
      std::cout << "  meanPtAve_      =  " << meanPtAve_ << std::endl;
      std::cout << "  meanPtGen_      =  " << meanPtGen_ << std::endl;
      std::cout << "  meanPdfPtTrue_  =  " << meanPdfPtTrue_ << std::endl;
    }
  }



  // -------------------------------------------------------------------
  Measurement::~Measurement() {
    for(HistMapIt it = hists_.begin(); it != hists_.end(); ++it) {
      delete it->second;
    }
    hists_.clear();
  }


  // -------------------------------------------------------------------
  TH1* Measurement::getClone(const TString &name) const {
    TH1* h = 0;
    TString histName = name+histNameSuffix_;
    std::map<TString,TH1*>::const_iterator it = hists_.find(histName);
    if( it != hists_.end() ) {
      TString newName = name+"_UID";
      newName += HIST_COUNT;
      h = static_cast<TH1*>(it->second->Clone(newName));
      ++HIST_COUNT;
    } else {
      std::cerr << "ERROR (Measurement): No histogram with name '" << histName << "'" << std::endl;
      exit(1);
    }

    return h;
  }


  // -------------------------------------------------------------------
  bool Measurement::parse(const TString &fileName, bool hasFittedParameters) {
    bool statusIsGood = true;
    if( verbosity_ == 2 ) std::cout << "(Measurement) Parsing file '" << fileName << "'" << std::endl;
    
    // Opening file    
    if( verbosity_ == 2 ) std::cout << "  Opening file... " << std::endl;
    TFile file(fileName,"READ");
    if( file.IsZombie() ) {
      std::cerr << "  ERROR: Error opening file." << std::endl;
      statusIsGood = false;
    }


    if( statusIsGood ) {
      if( hasFittedParameters ) {
	// Read fitted values and statistical uncertainties from file
	if( verbosity_ == 2 ) std::cout << "  Getting fitted values... " << std::endl;
	TH1 *h = 0;
	file.GetObject("hAbsoluteParameters"+histNameSuffix_,h);
	if( !h ) {
	  std::cerr << "  ERROR: 'hAbsoluteParameters' not found." << std::endl;
	  statusIsGood = false;
	} else {
	  h->SetDirectory(0);
	  for(int i = 0; i < h->GetNbinsX(); i++) {
	    values_.push_back(h->GetBinContent(1+i));
	    statUncert_.push_back(sqrt(2.)*h->GetBinError(1+i)); // Correct for wrong factor 2 in likelihood for LVMINI
	    if( verbosity_ == 2 ) {
	      std::cout << "  Value " << i << ": " << values_.back() << std::flush;
	      std::cout << "  Value: " << values_.back() << std::flush;
	      std::cout << " +/- " << statUncert_.back() << std::endl;
	    }
	  }
	}
      } else {
	// Set default value
	if( verbosity_ == 2 ) std::cout << "  Setting default fitted values... " << std::endl;
	values_.push_back(0.);
	statUncert_.push_back(0.);
      }

      // Read histograms from file
      if( verbosity_ == 2 ) std::cout << "  Getting histograms from file... " << std::endl;
      for(HistMapIt it = hists_.begin(); it != hists_.end(); it++) {
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
	  if( verbosity_ == 2 ) std::cout << "   (h->GetName() == '" << h->GetName() << "')... " << std::endl;
	}
      }
    }

    file.Close();

    return statusIsGood;
  }



  // -------------------------------------------------------------------
  void Measurement::setMean(const TString &name, double &mean, double &uncert) {
    HistMapIt it = hists_.find(name);
    if( it != hists_.end() ) {
      mean = it->second->GetMean();
      uncert = it->second->GetMeanError();
    } else {
      mean = 0.;
      uncert = 10000.;
    }
  }

}
