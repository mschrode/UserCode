// $Id: Parameters.cc,v 1.31 2011/11/21 17:18:05 mschrode Exp $

#include "Parameters.h"

#include <iostream>


namespace resolutionFit {

  // -------------------------------------------------------------------------------------
  Parameters::Parameters(const TString &id, const TString &binAdmCfg, unsigned int verbosity)
    : id_(id), verbosity_(verbosity), binAdm_(new sampleTools::BinningAdmin(binAdmCfg)) {

    jetAlgo_ = JetProperties::AK5;
    jetType_ = JetProperties::Calo;
    outPutMode_ = OutputManager::PSAllInOne;
    nEtaBinsUser_ = 0;
    ptSoftAbsMin_ = 10.;
    wpIdx_ = -1;

    // Print binning info
    binAdm_->printBinning();
    for(unsigned int i = 0; i < nPtSoftBins(); ++i) {
      std::cout << "  PtSoftBin " << i << ": " << ptSoftMin(i) << " - " << ptSoftMax(i) << std::endl;
    }
  }


  // -------------------------------------------------------------------------------------
  Parameters::~Parameters() {
    delete binAdm_;
  }


  // -------------------------------------------------------------------------------------
  void Parameters::setJetProperties(JetProperties::Algo algo, JetProperties::Type type) {
    jetAlgo_ = algo;
    jetType_ = type;

    id_ += "_"+JetProperties::toString(type);
  }


  // -------------------------------------------------------------------------------------
  void Parameters::setOutMode(OutputManager::Mode mode) {
    if( OutputManager::isValidMode(mode) ) {
      outPutMode_ = mode;
    }
  }


  // -------------------------------------------------------------------------------------
  void Parameters::setNEtaBinsUser(unsigned int nBins) {
    if( nBins > binAdm_->nEtaBins() ) {
      std::cerr << "ERROR: Number of maximum eta bins set by user larger than in config" << std::endl;
      exit(1);
    }
    nEtaBinsUser_ = nBins;
  }


  // -------------------------------------------------------------------------------------
  void Parameters::useWPExtrapolation(double cutValue) {
    wpIdx_ = 0;
    if( cutValue > ptSoftMax(nPtSoftBins()-1) ) {
      wpIdx_ = nPtSoftBins()-1;
    } else {
      while( cutValue > ptSoftMax(wpIdx_) ) ++wpIdx_;
    }
  }


  // -------------------------------------------------------------------------------------
  bool JetProperties::isValidType(Type type) {
    bool result = true;
    if( type != Calo && type != JPT && type != PF ) {
      result = false;
      std::cerr << "ERROR in JetProperties::isValidType(): No Type '" << type << "'" << std::endl;
      exit(1);
    }
    
    return result;
  }


  // -------------------------------------------------------------------------------------
  TString JetProperties::toString(Type type) {
    TString str;
    if( isValidType(type) ) {
      if( type == Calo ) str = "Calo";
      else if( type == JPT ) str = "JPT";
      else if( type == PF ) str = "PF";
    }

    return str;
  }
}
