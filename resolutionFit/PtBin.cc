// $Id: PtBin.cc,v 1.9 2010/05/14 09:04:15 mschrode Exp $

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
    cutVar_ = new CutVariation(par_);
    cutVar_->extrapolate();

    // Standard selection for reference
    KalibriFileParser *parserStdSel = new KalibriFileParser(par_->fileNameStdSel(),par_->verbosity());

    // Set up members
    relSigma_ = cutVar_->extrapolatedRelSigma();
    meanPt_ = parserStdSel->meanPt();
    meanPtUncert_ = parserStdSel->meanPtUncert();

    // Sum up systematic uncertainties
    Uncertainty *uncertSyst = new Uncertainty("SystematicUncertainty");
    // Reference sigma for unvaried case
    double refS = parserStdSel->value();
    // Calculate relative deviation after variation
    for(int i = 0; i < par_->nSystUncerts(); ++i) {
      KalibriFileParser *parser = new KalibriFileParser(par_->fileNameSystUncertUp(i),par_->verbosity());
      double dUp = parser->value() - refS;
      delete parser;
      if( par_->verbosity() == 2 ) {
	std::cout << " Syst " << meanPt_ << std::flush;
	std::cout << ":  " << dUp << std::flush;
      }
      dUp /= meanPt_;
      uncertSyst->addUncertainty(new Uncertainty(par_->labelSystUncert(i),dUp,0.));
    }    
    // Sum up systematic and statistic uncertainty
    uncert_ = new Uncertainty("TotalUncertainty");
    uncert_->addUncertainty(new Uncertainty("StatisticUncertainty",cutVar_->extrapolatedUncert()));
    uncert_->addUncertainty(uncertSyst);

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

    if( par_->verbosity() == 2 ) std::cout << "ok" << std::endl;

    delete parserStdSel;

    if( par_->verbosity() == 2 ) {
      std::cout << "Is combined uncertainty: " << std::flush;
      std::cout << ( uncert_->isCombined() ? "yes" : "no" ) << std::endl;
      std::cout << "Syst uncert at " << meanPt_ << " GeV: +" << std::flush;
      std::cout << uncertSystUp() << ", -" << uncertSystDown() << std::endl;
    }
  }

  PtBin::~PtBin() {
    delete cutVar_;
    delete uncert_;
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
    else {
      std::cerr << "ERROR PtBin::getHist: No histogram of name '" << name << "'" << std::endl;
      exit(-1);
    }
      
    return h;
  }
}
