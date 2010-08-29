// $Id: PtBin.cc,v 1.18 2010/08/28 19:34:19 mschrode Exp $

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
      cutVar_.push_back(new CutVariation(par_,parIdx,"MaxLikeFit"));
      cutVar_[parIdx]->extrapolate();
      extrapolatedVal_.push_back(cutVar_[parIdx]->extrapolatedValue());
    }

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

    // Perform cut variation and extrapolation
    // of pt asymmetry
    cutVarAsym_ = new CutVariation(par_,0,"PtAsym");
    cutVarAsym_->extrapolate();
    extrapolatedAsym_ = cutVarAsym_->extrapolatedValue();
    uncertAsym_ = new Uncertainty("StatisticUncertainty",cutVarAsym_->extrapolatedUncert());

    // Perform cut variation and extrapolation
    // of ptGen asymmetry
    cutVarGenAsym_ = new CutVariation(par_,0,"PtGenAsym");
    cutVarGenAsym_->extrapolate();
    extrapolatedGenAsym_ = cutVarGenAsym_->extrapolatedValue();
    uncertGenAsym_ = new Uncertainty("StatisticUncertainty",cutVarGenAsym_->extrapolatedUncert());


    // Store spectrum and response histograms
    if( par_->verbosity() == 2 ) std::cout << "Storing spectrum and resolution histograms... " << std::flush;
    hPtJet1_ = parserStdSel->hist("hPtJet1","hPtJet1_PtBin"+util::toTString(ptBinIdx()));
    hPtJet2_ = parserStdSel->hist("hPtJet2","hPtJet2_PtBin"+util::toTString(ptBinIdx()));
    hPtJet3_ = parserStdSel->hist("hPtJet3","hPtJet3_PtBin"+util::toTString(ptBinIdx()));
    hPtJet4_ = parserStdSel->hist("hPtJet4","hPtJet4_PtBin"+util::toTString(ptBinIdx()));
    hPtGen_ = parserStdSel->hist("hPtGen","hPtGen_PtBin"+util::toTString(ptBinIdx()));
    hPtGenJet1_ = parserStdSel->hist("hPtGenJet1","hPtGenJet1_PtBin"+util::toTString(ptBinIdx()));
    hPdfPtTrue_ = parserStdSel->hist("hTruthPDF","hPdfPtTrue_PtBin"+util::toTString(ptBinIdx()));
    hPtAve_ = parserStdSel->hist("hPtDijet","hPtAve_PtBin"+util::toTString(ptBinIdx()));
    hPdfRes_ = parserStdSel->hist("hRespFit_0","hPdfRes_PtBin"+util::toTString(ptBinIdx()));
    hPtGenAsym_ = parserStdSel->hist("hPtGenAsym_0","hPtGenAsym_PtBin"+util::toTString(ptBinIdx()));
    hPJet3_ = parserStdSel->hist("hPJet3","hPJet3_PtBin"+util::toTString(ptBinIdx()));
    hPJet3Rel_ = parserStdSel->hist("hPJet3Rel","hPJet3Rel_PtBin"+util::toTString(ptBinIdx()));
    hPJet3GenRel_ = parserStdSel->hist("hPJet3GenRel","hPJet3GenRel_PtBin"+util::toTString(ptBinIdx()));
    hPSJ_ = parserStdSel->hist("hPSJ","hPSJ_PtBin"+util::toTString(ptBinIdx()));
    hPSJRel_ = parserStdSel->hist("hPSJRel","hPSJRel_PtBin"+util::toTString(ptBinIdx()));
    hPSJGenRel_ = parserStdSel->hist("hPSJGenRel","hPSJGenRel_PtBin"+util::toTString(ptBinIdx()));
    hEta_ = parserStdSel->hist("hEta","hEta_PtBin"+util::toTString(ptBinIdx()));
    hDeltaPhi12_ = parserStdSel->hist("hDeltaPhi12","hDeltaPhi12_PtBin"+util::toTString(ptBinIdx()));

    delete parserStdSel;

    if( par_->hasMCClosure() ) {
      KalibriFileParser *parserMCClosure = new KalibriFileParser(par_->fileNameMCClosure(),par_->verbosity());
      //      hMCRes_ = parserMCClosure->hist("hRespMeas_0","hMCRes_PtBin"+util::toTString(ptBinIdx()));
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
    delete cutVarGenAsym_;
    delete uncertGenAsym_;
    delete hPtJet1_;
    delete hPtJet2_;
    delete hPtJet3_;
    delete hPtJet4_;
    delete hPtGen_;
    delete hPtGenJet1_;
    delete hPdfPtTrue_;
    delete hPtAve_;
    delete hPdfRes_;
    delete hPJet3_;
    delete hPJet3Rel_;
    delete hPJet3GenRel_;
    delete hPSJ_;
    delete hPSJRel_;
    delete hPSJGenRel_;
    delete hEta_;
    delete hDeltaPhi12_;
  }


  double PtBin::extrapolatedValue(int parIdx, bool corrected) const {
    double val = extrapolatedVal_.at(parIdx);
    if( parIdx == 0 && corrected ) {
      val = sqrt( val*val - extrapolatedGenAsym_*extrapolatedGenAsym_ );
    }
    return val;
  }


  double PtBin::extrapolatedAsym(bool corrected) const {
    double asym = extrapolatedAsym_;
    if( corrected ) {
      asym = sqrt( asym*asym - extrapolatedGenAsym_*extrapolatedGenAsym_ );
    }
    return asym;
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
    else if( name == "hPtJet1" ) h = static_cast<TH1D*>(hPtJet1_->Clone(newName));
    else if( name == "hPtJet2" ) h = static_cast<TH1D*>(hPtJet2_->Clone(newName));
    else if( name == "hPtJet3" ) h = static_cast<TH1D*>(hPtJet3_->Clone(newName));
    else if( name == "hPtJet4" ) h = static_cast<TH1D*>(hPtJet4_->Clone(newName));
    else if( name == "hPJet3" ) h = static_cast<TH1D*>(hPJet3_->Clone(newName));
    else if( name == "hPJet3Rel" ) h = static_cast<TH1D*>(hPJet3Rel_->Clone(newName));
    else if( name == "hPJet3GenRel" ) h = static_cast<TH1D*>(hPJet3GenRel_->Clone(newName));
    else if( name == "hPSJ" ) h = static_cast<TH1D*>(hPSJ_->Clone(newName));
    else if( name == "hPSJRel" ) h = static_cast<TH1D*>(hPSJRel_->Clone(newName));
    else if( name == "hPSJGenRel" ) h = static_cast<TH1D*>(hPSJGenRel_->Clone(newName));
    else if( name == "hEta" ) h = static_cast<TH1D*>(hEta_->Clone(newName));
    else if( name == "hDeltaPhi12" ) h = static_cast<TH1D*>(hDeltaPhi12_->Clone(newName));

    else {
      std::cerr << "ERROR PtBin::getHist: No histogram of name '" << name << "'" << std::endl;
      exit(-1);
    }
      
    return h;
  }
}
