// $Id: $

#include "Parameters.h"

#include <iostream>


namespace resolutionFit {

  // -------------------------------------------------------------------------------------
  Parameters::Parameters(const TString &id, const TString &binAdmCfg, unsigned int verbosity)
    : id_(id), verbosity_(verbosity), binAdm_(new sampleTools::BinningAdmin(binAdmCfg)) {

    jetAlgo_ = JetProperties::AK5;
    jetType_ = JetProperties::Calo;

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
