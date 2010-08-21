// $Id: PtBin.cc,v 1.15 2010/08/18 16:18:06 mschrode Exp $

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
      for(int i = 0; i < par_->nPt3Cuts(); i++) {
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

    // Perform cut variation and extrapolation
    // of pt asymmetry
    cutVarAsym_ = new CutVariation(par_,0,false);
    cutVarAsym_->extrapolate();
    extrapolatedAsym_ = cutVarAsym_->extrapolatedValue();

    // Standard selection for reference
    KalibriFileParser *parserStdSel = new KalibriFileParser(par_->fileNameStdSel(),par_->verbosity());
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
	  std::cout << " Syst " << meanPt() << std::flush;
	  std::cout << ":  " << dUp << std::flush;
	}
	if( par_->isRelParValue(parIdx) ) dUp /= meanPt();
	uncertSyst->addUncertainty(new Uncertainty(par_->labelSystUncert(i),dUp,0.));
      }  
      // Sum up systematic and statistic uncertainty
      uncert_.push_back(new Uncertainty("TotalUncertainty"));
      uncert_[parIdx]->addUncertainty(new Uncertainty("StatisticUncertainty",
						      cutVar_[parIdx]->extrapolatedUncert()));
      uncert_[parIdx]->addUncertainty(uncertSyst);
    } // End of loop over paramters
    
    uncertAsym_ = new Uncertainty("StatisticUncertainty",cutVarAsym_->extrapolatedUncert());

    // Store spectrum and response histograms
    if( par_->verbosity() == 2 ) std::cout << "Storing spectrum and resolution histograms... " << std::flush;
    TString name = "hPtGen_PtBin";
    name += ptBinIdx();
    hPtGen_ = parserStdSel->hist("hPtGen",name);
    name = "hPtGenJet1_PtBin";
    name += ptBinIdx();
    hPtGenJet1_ = parserStdSel->hist("hPtGenJet1",name);
    name = "hPdfPtTrue_PtBin";
    name += ptBinIdx();
    hPdfPtTrue_ = parserStdSel->hist("hTruthPDF",name);
    name = "hPtAve_PtBin";
    name += ptBinIdx();
    hPtAve_ = parserStdSel->hist("hPtDijet",name);
    name = "hResGen_PtBin";
    name += ptBinIdx();
    hResGen_ = parserStdSel->hist("hRespMeas_0",name);
    name = "hPdfRes_PtBin";
    name += ptBinIdx();
    hPdfRes_ = parserStdSel->hist("hRespFit_0",name);
    name = "hPtGenAsym_PtBin";
    name += ptBinIdx();
    hPtGenAsym_ = parserStdSel->hist("hPtGenAsym_0",name);

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
	std::cout << "Syst uncert at " << meanPt() << " GeV: +" << std::flush;
	std::cout << uncertSystUp(parIdx) << ", -" << uncertSystDown(parIdx) << std::endl;
      }
    }
  }

  PtBin::~PtBin() {
    for(int parIdx = 0; parIdx < par_->nFittedPars(); ++parIdx) {
      delete cutVar_[parIdx];
      delete uncert_[parIdx];
    }
    delete cutVarAsym_;
    delete uncertAsym_;
    delete hPtGen_;
    delete hPtGenJet1_;
    delete hPdfPtTrue_;
    delete hPtAve_;
    delete hResGen_;
    delete hPdfRes_;
  }

  
  TH1 *PtBin::getHist(const TString &name, const TString &newName) const {
    TH1 *h = 0;
    if( name == "hPtGen" ) h = static_cast<TH1D*>(hPtGen_->Clone(newName));
    else if( name == "hPtGenJet1" ) h = static_cast<TH1D*>(hPtGenJet1_->Clone(newName));
    else if( name == "hPtAve" ) h = static_cast<TH1D*>(hPtAve_->Clone(newName));
    else if( name == "hPdfPtTrue" ) h = static_cast<TH1D*>(hPdfPtTrue_->Clone(newName));
    else if( name == "hResGen" ) h = static_cast<TH1D*>(hResGen_->Clone(newName));
    else if( name == "hPdfRes" ) h = static_cast<TH1D*>(hPdfRes_->Clone(newName));
    else if( name == "hPtGenAsym" ) h = static_cast<TH1D*>(hPtGenAsym_->Clone(newName));
    else if( name == "hMCRes" ) h = static_cast<TH1D*>(hMCRes_->Clone(newName));
    else {
      std::cerr << "ERROR PtBin::getHist: No histogram of name '" << name << "'" << std::endl;
      exit(-1);
    }
      
    return h;
  }
}
