// $Id: Measurement.cc,v 1.7 2011/08/19 08:33:56 mschrode Exp $

#include "Measurement.h"

#include <cmath>
#include <iostream>

#include "TFile.h"
#include "TH1.h"


namespace resolutionFit {

  unsigned int Measurement::HIST_COUNT = 0;

  // -------------------------------------------------------------------
  Measurement::Measurement(const TString &fileName, unsigned int etaBin, unsigned int ptBin, unsigned int ptSoftBin, double ptMin, double ptMax, double ptSoft, unsigned int verbosity)
    : ptMin_(ptMin), ptMax_(ptMax), ptSoft_(ptSoft), verbosity_(verbosity) {
    init(fileName,etaBin,ptBin,ptSoftBin);
  }
  

  // -------------------------------------------------------------------
  void Measurement::init(const TString &fileName, unsigned int etaBin, unsigned int ptBin, unsigned int ptSoftBin) {

    hnPrefix_ = fileName+":///"+histNamePath(etaBin,ptBin,ptSoftBin);
    hnSuffix_ = histNameSuffix(etaBin,ptBin,ptSoftBin);

    // Histograms to be read from file
    // Within Measurement, the histograms are stored with name
    // "HistName_EtaX_PtX_PtSoftX".
    hists_[("hPtGen"+hnSuffix_)] = 0;
    hists_[("hPtGenJet1"+hnSuffix_)] = 0;
    hists_[("hPtAveAbs"+hnSuffix_)] = 0;
    hists_[("hTruthPDF"+hnSuffix_)] = 0;
    hists_[("hPtAbsAsym"+hnSuffix_)] = 0;
    hists_[("hFitPtAsym"+hnSuffix_)] = 0;
    hists_[("hRespMeas"+hnSuffix_)] = 0;
    hists_[("hRespFit"+hnSuffix_)] = 0;
    hists_[("hPtGenAsym"+hnSuffix_)] = 0;
    hists_[("hPtJet1"+hnSuffix_)] = 0;
    hists_[("hPtJet2"+hnSuffix_)] = 0;
    hists_[("hPtJet3"+hnSuffix_)] = 0;
    hists_[("hPtJet4"+hnSuffix_)] = 0;
    hists_[("hPJet3"+hnSuffix_)] = 0;
    hists_[("hPJet3Rel"+hnSuffix_)] = 0;
    hists_[("hPJet3GenRel"+hnSuffix_)] = 0;
    hists_[("hPSJ"+hnSuffix_)] = 0;
    hists_[("hPSJRel"+hnSuffix_)] = 0;
    hists_[("hPSJGenRel"+hnSuffix_)] = 0;
    hists_[("hEta"+hnSuffix_)] = 0;
    hists_[("hDeltaPhi12"+hnSuffix_)] = 0;
    hists_[("hDeltaPtJet12"+hnSuffix_)] = 0;
    hists_[("hWeights"+hnSuffix_)] = 0;
    hists_[("hNumPU"+hnSuffix_)] = 0;
    hists_[("hNumVtx"+hnSuffix_)] = 0;
    
    // Parse file
    if( !parse(fileName) ) exit(-1);

    // Calculate mean pt values
    if( verbosity_ == 2 ) std::cout << "Setting mean pt values" << std::endl;
    setMean("hPtAveAbs"+hnSuffix_,meanPtAve_,meanPtAveUncert_);
    setMean("hPtGen"+hnSuffix_,meanPtGen_,meanPtGenUncert_);
    setMean("hTruthPDF"+hnSuffix_,meanPdfPtTrue_,meanPdfPtTrueUncert_);
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
    delete hSpectrum_;
  }


  //! \brief Path to histogram within input file
  //!
  //! The histograms are expected to be stored in a subdirectory
  //!  "EtaX/PtX/PtSoftX"
  // -------------------------------------------------------------------
  TString Measurement::histNamePath(unsigned int etaBin, unsigned int ptBin, unsigned int ptSoftBin) const {
    TString path = "Eta";
    path += etaBin;
    path += "/Pt";
    path += ptBin;
    path += "/PtSoft";
    path += ptSoftBin;
    path += "/";

    return path;
  }


  //! \brief Suffix of histogram name
  //!
  //! The suffix has the form "EtaX_PtX_PtSoftX"
  // -------------------------------------------------------------------
  TString Measurement::histNameSuffix(unsigned int etaBin, unsigned int ptBin, unsigned int ptSoftBin) const {
    TString name = "_Eta";
    name += etaBin;
    name += "_Pt";
    name += ptBin;
    name += "_PtSoft";
    name += ptSoftBin;
    
    return name;
  }




  // -------------------------------------------------------------------
  TH1* Measurement::getClone(const TString &name) const {
    TH1* h = 0;
    TString histName = name+hnSuffix_;
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


  //! \brief Parse the input file
  //!
  //! Read and store the histograms and fitted parameter values.
  // -------------------------------------------------------------------
  bool Measurement::parse(const TString &fileName) {
    bool statusIsGood = true;
    if( verbosity_ == 2 ) std::cout << "(Measurement) Parsing file '" << fileName << "'" << std::endl;
    
    // Opening file    
    if( verbosity_ == 2 ) std::cout << "  Opening file... " << std::endl;
    TFile file(fileName,"READ");
    if( file.IsZombie() ) {
      std::cerr << "  ERROR: Error opening file " << fileName << std::endl;
      statusIsGood = false;
    }

    if( statusIsGood ) {
      // Read fitted values and statistical uncertainties from file
      if( verbosity_ == 2 ) std::cout << "  Getting fitted values... " << std::endl;
      TH1 *h = 0;
      file.GetObject(hnPrefix_+"hAbsoluteParameters"+hnSuffix_,h);
      if( !h ) {
	std::cerr << "  ERROR: 'hAbsoluteParameters' not found." << std::endl;
	// Set default value
	if( verbosity_ == 2 ) std::cout << "  Setting default fitted values... " << std::endl;
	values_.push_back(0.);
	statUncert_.push_back(0.);
	//statusIsGood = false;
      } else {
	h->SetDirectory(0);
	for(int i = 0; i < h->GetNbinsX(); i++) {
	  values_.push_back(h->GetBinContent(1+i));
	  statUncert_.push_back(h->GetBinError(1+i));
	  
	  if( verbosity_ == 2 ) {
	    std::cout << "  Value " << i << ": " << values_.back() << std::flush;
	    std::cout << "  Value: " << values_.back() << std::flush;
	    std::cout << " +/- " << statUncert_.back() << std::endl;
	  }
	}
      }

      // Read underlying spectrum from file
      hSpectrum_ = 0;
      file.GetObject(hnPrefix_+"hUnderlyingTruthPDF"+hnSuffix_,hSpectrum_);
      if( !hSpectrum_ ) {
	std::cerr << "  WARNING: '" << hnPrefix_+"hUnderlyingTruthPDF"+hnSuffix_ << "' not found." << std::endl;
	++HIST_COUNT;
	TString name = "Measurement::Dummy::";
	name += HIST_COUNT;
	hSpectrum_ = new TH1D(name,"",1,0,1);
	hSpectrum_->SetDirectory(0);
      } else {
	hSpectrum_->SetDirectory(0);
	hSpectrum_->SetName("hSpectrum::"+hnSuffix_);
	hSpectrum_->UseCurrentStyle();
	hSpectrum_->SetTitle("");
	if( hSpectrum_->Integral() ) hSpectrum_->Scale(1./hSpectrum_->Integral("width"));
	if( verbosity_ == 2 ) std::cout << "   (h->GetName() == '" << hSpectrum_->GetName() << "')... " << std::endl;
      }

      // Read histograms from file
      if( verbosity_ == 2 ) std::cout << "  Getting histograms from file... " << std::endl;
      for(HistMapIt it = hists_.begin(); it != hists_.end(); it++) {
	TH1 *h = 0;
	TString name = it->first;
	file.GetObject(hnPrefix_+name,h);
	if( !h ) {
	  std::cerr << "  WARNING: '" << name << "' not found." << std::endl;
	  it->second = new TH1D(name,"",1,0,1);
	  it->second->SetDirectory(0);
	} else {
	  h->SetDirectory(0);
	  h->UseCurrentStyle();
	  h->SetTitle("");
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


  //! \brief Return pdf of assumed underlying dijet spectrum
  //!
  //! This is the underlying ptGen spectrum assumed in the fit, i.e.
  //! it is not modified to describe migration effects.
  // -------------------------------------------------------------------
  double Measurement::pdfPtTrue(double pt) const {
    return hSpectrum_->Interpolate(pt);
  }
}
