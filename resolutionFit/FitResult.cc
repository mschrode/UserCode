// $Id: $

#include "FitResult.h"

#include <algorithm>
#include <iostream>


namespace resolutionFit {

  // -------------------------------------------------------------------------------------
  bool FitResult::validType(Type type) {
    if( type != FullMaxLike ) {
      std::cerr << "ERROR in FitResult::validType(): Unknown type '" << type << "'" << std::endl;
      exit(1);
    }
    
    return true;
  }


  // -------------------------------------------------------------------------------------
  FitResult* FitResult::createFitResult(Type type, const std::vector<Measurement*> &meas, unsigned int verbosity) {
    FitResult* fr = 0;
    if( validType(type) ) {
      if( type == FullMaxLike ) {
	fr = new FitResultFullMaxLike(meas,verbosity);
      }
      fr->init();
    }

    return fr;
  }



  // -------------------------------------------------------------------------------------  
  TString FitResult::toString(FitResult::Type type) {
    TString str = "Undefined";
    if( type == FullMaxLike ) str = "FullMaxLike";
    else if( type == SimpleMaxLike ) str = "SimpleMaxLike";
    else if( type == PtAsym ) str = "PtAsym";
    else if( type == PtGenAsym ) str = "PtGenAsym";

    return str;
  }


  // -------------------------------------------------------------------------------------  
  FitResult::FitResult(const Meas &meas, unsigned int verbosity)
    : meas_(meas), verbosity_(verbosity) {

    extrapolation_ = 0;
  }


  // -------------------------------------------------------------------------------------  
  FitResult::~FitResult() {
    if( extrapolation_ ) delete extrapolation_;
  }

  
  void FitResult::init() {
    // Store fitted values and mean pt
    initResult();

    // Perform extrapolation
    if( !extrapolate() ) std::cerr << "WARNING in FitResult: Error in extrapolation" << std::endl;
  }


  bool FitResult::extrapolate() {
    if( verbosity_ == 1 ) std::cout << "Extrapolating" << std::endl;
    Extrapolation extra;

    return extra(meanPt_,ptSoft_,values_,statUncerts_,extrapolation_);
  }




  FitResultFullMaxLike::FitResultFullMaxLike(const Meas meas, unsigned int verbosity)
    : FitResult(meas,verbosity) { }


  bool FitResultFullMaxLike::initResult() {
    values_.clear();
    statUncerts_.clear();
    ptSoft_.clear();

    // Set ptSoft cut values and find smalles ptSoft for meanPt
    double ptSoftSmall = 1000.;
    for(MeasIt it = meas_.begin() ; it != meas_.end(); ++it) {
      ptSoft_.push_back((*it)->ptSoft());
      if( (*it)->ptSoft() < ptSoftSmall ) {
	ptSoftSmall = (*it)->ptSoft();
	
	std::cout << "FitResultFullMaxLike::storeResult(): TODO MeanPt from spectrum convoluted with extrapolated *sigma*!" << std::endl;
	meanPt_ = (*it)->meanPdfPtTrue();
	meanPtUncert_ = (*it)->meanPdfPtTrueUncert();
      }
    }

    // Set fitted values
    for(MeasIt it = meas_.begin(); it != meas_.end(); ++it) {
      values_.push_back((*it)->fittedValue(0)/meanPt_);
      statUncerts_.push_back((*it)->fittedUncert(0)/meanPt_);
    }
  }
}
