// $Id: PtBin.cc,v 1.11 2010/05/26 21:56:36 mschrode Exp $

#include "PtBin.h"

#include <iostream>

#include "TH1.h"
#include "TH1D.h"

#include "KalibriFileParser.h"

namespace resolutionFit {
  PtBin::PtBin(const Parameters::PtBinParameters *par)
    : par_(par) {

    if( par_->verbosity() == 2 ) {
      std::cout << "\nPtBin::PtBin" << std::endl;
      std::cout << " fileNameStdSel: '" << par_->fileNameStdSel() << "'" << std::endl;
      std::cout << " fileNamesCutVar:" << std::endl;
      for(int i = 0; i < par_->nPt3CutVariations(); i++) {
	std::cout << "  " << i << ": '" << par_->fileNamePt3CutVariations(i) << "'" << std::endl;
      }
    }

    // Perform cut variation and extrapolation
    // of fitted values
    for(int parIdx = 0; parIdx < par_->nFittedPars(); ++parIdx) {
      cutVar_.push_back(new CutVariation(par_,parIdx));
      cutVar_[parIdx]->extrapolate();
      extrapolatedVal_.push_back(cutVar_[parIdx]->extrapolatedValue());
    }

    // Standard selection for reference
    KalibriFileParser *parserStdSel = new KalibriFileParser(par_->fileNameStdSel(),par_->verbosity());

    // Mean pt of this bin
    meanPt_ = parserStdSel->meanPt();
    meanPtUncert_ = parserStdSel->meanPtUncert();

    for(int parIdx = 0; parIdx < par_->nFittedPars(); ++parIdx) { // Loop over paramters
      // Reference sigma for unvaried case
      double refS = parserStdSel->value(parIdx);
      // Sum up systematic uncertainties
      Uncertainty *uncertSyst = new Uncertainty("SystematicUncertainty");
      // Calculate relative deviation after variation
      for(int i = 0; i < par_->nSystUncerts(); ++i) {
	KalibriFileParser *parser = new KalibriFileParser(par_->fileNameSystUncertUp(i),par_->verbosity());
	double dUp = parser->value(parIdx) - refS;
	delete parser;
	if( par_->verbosity() == 2 ) {
	  std::cout << " Syst " << meanPt_ << std::flush;
	  std::cout << ":  " << dUp << std::flush;
	}
	if( par_->isRelParValue(parIdx) ) dUp /= meanPt_;
	uncertSyst->addUncertainty(new Uncertainty(par_->labelSystUncert(i),dUp,0.));
      }  
      // Sum up systematic and statistic uncertainty
      uncert_.push_back(new Uncertainty("TotalUncertainty"));
      uncert_[parIdx]->addUncertainty(new Uncertainty("StatisticUncertainty",
						      cutVar_[parIdx]->extrapolatedUncert()));
      uncert_[parIdx]->addUncertainty(uncertSyst);
    } // End of loop over paramters

    // Store spectrum and response histograms
    if( par_->verbosity() == 2 ) std::cout << "Storing spectrum and resolution histograms... " << std::flush;
    TString name = "hPtGen_PtBin";
    name += ptBinIdx();
    hPtGen_ = parserStdSel->hist("hPtGen",name);
    name = "hPtGen_PtBinJet1";
    name += ptBinIdx();
    hPtGenJet1_ = parserStdSel->hist("hPtGenJet1",name);
    name = "hPdfPtTrue_PtBin";
    name += ptBinIdx();
    hPdfPtTrue_ = parserStdSel->hist("hTruthPDF",name);
    name = "hResGen_PtBin";
    name += ptBinIdx();
    hResGen_ = parserStdSel->hist("hRespMeas_0",name);
    name = "hPdfRes_PtBin";
    name += ptBinIdx();
    hPdfRes_ = parserStdSel->hist("hRespFit_0",name);
    name = "hPtAsym_PtBin";
    name += ptBinIdx();
    hPtAsym_ = parserStdSel->hist("hPtAsym_0",name);
    name = "hPdfPtAsym_PtBin";
    name += ptBinIdx();
    hPdfPtAsym_ = parserStdSel->hist("hFitPtAsym_0",name);

    delete parserStdSel;

    if( par_->hasMCClosure() ) {
      KalibriFileParser *parserMCClosure = new KalibriFileParser(par_->fileNameMCClosure(),par_->verbosity());
      name = "hMCRes_PtBin";
      name += ptBinIdx();
      hMCRes_ = parserMCClosure->hist("hRespMeas_0",name);
      delete parserMCClosure;
    }

    if( par_->verbosity() == 2 ) std::cout << "ok" << std::endl;

    if( par_->verbosity() == 2 ) {
      for(int parIdx = 0; parIdx < par_->nFittedPars(); ++parIdx) {
	std::cout << "Is combined uncertainty: " << std::flush;
	std::cout << ( uncert_[parIdx]->isCombined() ? "yes" : "no" ) << std::endl;
	std::cout << "Syst uncert at " << meanPt_ << " GeV: +" << std::flush;
	std::cout << uncertSystUp(parIdx) << ", -" << uncertSystDown(parIdx) << std::endl;
      }
    }
  }

  PtBin::~PtBin() {
    for(int parIdx = 0; parIdx < par_->nFittedPars(); ++parIdx) {
      delete cutVar_[parIdx];
      delete uncert_[parIdx];
    }
    delete hPtGen_;
    delete hPtGenJet1_;
    delete hPdfPtTrue_;
    delete hResGen_;
    delete hPdfRes_;
    delete hPtAsym_;
    delete hPdfPtAsym_;
  }

  
  TH1 *PtBin::getHist(const TString &name, const TString &newName) const {
    TH1 *h = 0;
    if( name == "hPtGen" ) h = static_cast<TH1D*>(hPtGen_->Clone(newName));
    else if( name == "hPtGenJet1" ) h = static_cast<TH1D*>(hPtGenJet1_->Clone(newName));
    else if( name == "hPdfPtTrue" ) h = static_cast<TH1D*>(hPdfPtTrue_->Clone(newName));
    else if( name == "hResGen" ) h = static_cast<TH1D*>(hResGen_->Clone(newName));
    else if( name == "hPdfRes" ) h = static_cast<TH1D*>(hPdfRes_->Clone(newName));
    else if( name == "hPtAsym" ) h = static_cast<TH1D*>(hPtAsym_->Clone(newName));
    else if( name == "hPdfPtAsym" ) h = static_cast<TH1D*>(hPdfPtAsym_->Clone(newName));
    else if( name == "hMCRes" ) h = static_cast<TH1D*>(hMCRes_->Clone(newName));
    else {
      std::cerr << "ERROR PtBin::getHist: No histogram of name '" << name << "'" << std::endl;
      exit(-1);
    }
      
    return h;
  }
}
