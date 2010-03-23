// $Id: PtBin.cc,v 1.4 2010/03/22 19:45:08 mschrode Exp $

#include "PtBin.h"

#include <iostream>

#include "TH1F.h"

#include "KalibriFileParser.h"

namespace resolutionFit {
  int PtBin::nPtBins = 0;

  PtBin::PtBin(const TString &fileNameStdSel,
	       const std::vector<TString> &fileNamesCutVariation, const std::vector<double> &cutValues,
	       const TString &fileNameMCStatUncert,
	       const std::vector<TString> &fileNamesSystUncertUp, 
	       const std::vector<TString> &fileNamesSystUncertDown,
	       const std::vector<TString> &labelsSystUncertainties, double minPt, double maxPt, int verbose)
    : verbose_(verbose) {

    if( verbose_ == 2 ) {
      std::cout << "\nPtBin::PtBin" << std::endl;
      std::cout << " fileNameStdSel: '" << fileNameStdSel << "'" << std::endl;
      std::cout << " fileNamesCutVar:" << std::endl;
      for(size_t i = 0; i < fileNamesCutVariation.size(); i++) {
	std::cout << "  " << i << ": '" << fileNamesCutVariation[i] << "'" << std::endl;
      }
    }

    // Perform cut variation and extrapolation
    cutVar_ = new CutVariation(fileNamesCutVariation,cutValues,fileNameMCStatUncert,verbose);
    cutVar_->extrapolate();

    // Standard selection for reference
    KalibriFileParser *parserStdSel = new KalibriFileParser(fileNameStdSel,verbose_);

    // Set up members
    relSigma_ = cutVar_->extrapolatedRelSigma();
    minPt_ = minPt;
    maxPt_ = maxPt;
    meanPt_ = parserStdSel->meanPdfPtTrue();

    // Sum up systematic uncertainties
    assert( fileNamesSystUncertDown.size() == fileNamesSystUncertUp.size() );
    Uncertainty *uncertSyst = new Uncertainty("SystematicUncertainty");
    // Reference sigma for unvaried case
    double refS = parserStdSel->value();
    // Calculate relative deviation after variation
    for(size_t i = 0; i < fileNamesSystUncertUp.size(); i++) {
      double dUp = 0.;
      double dDown = 0.;
      // Upward variation
      KalibriFileParser *parser = new KalibriFileParser(fileNamesSystUncertUp[i],verbose_);
      double dVar = parser->value() - refS;
      delete parser;
      if( verbose_ == 2 ) {
	std::cout << " Syst " << meanPt_ << std::flush;
	std::cout << ":  " << dVar << " (up)" << std::flush;
      }
      dVar /= meanPt_;
      if( dVar > 0 ) dUp += dVar;
      else dDown += dVar;

      // Downward variation
      parser = new KalibriFileParser(fileNamesSystUncertDown[i],verbose_);
      dVar = parser->value() - refS;
      delete parser;
      if( verbose_ == 2 ) {
	std::cout << "  " << dDown << " (down)" << std::endl;
      }
      dVar /= meanPt_;
      if( dVar > 0 ) dUp += dVar;
      else dDown += dVar;

      uncertSyst->addUncertainty(new Uncertainty(labelsSystUncertainties[i],dUp,dDown));
    }    
    // Sum up systematic and statistic uncertainty
    uncert_ = new Uncertainty("TotalUncertainty");
    uncert_->addUncertainty(new Uncertainty("StatisticUncertainty",cutVar_->extrapolatedUncert()));
    uncert_->addUncertainty(uncertSyst);

    // Store spectrum and response histograms
    if( verbose_ == 2 ) std::cout << "Storing spectrum and resolution histograms... " << std::flush;
    TString name = "hPtGen_MinPt";
    name += PtBin::nPtBins;
    hPtGen_ = parserStdSel->hist("hPtGen",name);
    name = "hPdfPtTrue_MinPt";
    name += minPt_;
    hPdfPtTrue_ = parserStdSel->hist("hTruthPDF",name);
    name = "hResGen_MinPt";
    name += minPt_;
    hResGen_ = parserStdSel->hist("hRespMeas_0",name);
    name = "hPdfRes_MinPt";
    name += minPt_;
    hPdfRes_ = parserStdSel->hist("hRespFit_0",name);
    if( verbose_ == 2 ) std::cout << "ok" << std::endl;

    delete parserStdSel;

    PtBin::nPtBins++;

    if( verbose_ == 2 ) {
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
    delete hPdfPtTrue_;
    delete hResGen_;
    delete hPdfRes_;
  }

  
  TH1F *PtBin::getHist(const TString &name, const TString &newName) const {
    TH1F *h = 0;
    if( name == "hPtGen" ) h = static_cast<TH1F*>(hPtGen_->Clone(newName));
    else if( name == "hPdfPtTrue" ) h = static_cast<TH1F*>(hPdfPtTrue_->Clone(newName));
    else if( name == "hResGen" ) h = static_cast<TH1F*>(hResGen_->Clone(newName));
    else if( name == "hPdfRes" ) h = static_cast<TH1F*>(hPdfRes_->Clone(newName));
    else {
      std::cerr << "ERROR PtBin::getHist: No histogram of name '" << name << "'" << std::endl;
      exit(-1);
    }
      
    return h;
  }
}
