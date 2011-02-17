// $Id: FitResult.cc,v 1.1 2011/02/15 18:22:25 mschrode Exp $

#include "FitResult.h"

#include <algorithm>
#include <iostream>


namespace resolutionFit {

  // -------------------------------------------------------------------------------------
  bool FitResult::validType(Type type) {
    if( type != FullMaxLikeRel && type != FullMaxLikeAbs ) {
      std::cerr << "ERROR in FitResult::validType(): Unknown type '" << type << "'" << std::endl;
      exit(1);
    }
    
    return true;
  }


  // -------------------------------------------------------------------------------------
  FitResult* FitResult::createFitResult(Type type, const std::vector<Measurement*> &meas, unsigned int verbosity) {
    FitResult* fr = 0;
    if( validType(type) ) {
      if( type == FullMaxLikeRel ) {
	fr = new FitResultFullMaxLikeRel(meas,verbosity);
      } else if( type == FullMaxLikeAbs ) {
	fr = new FitResultFullMaxLikeAbs(meas,verbosity);
      }
      if( !fr->init() ) {
	std::cerr << "ERROR in FitResult::createFitResult(): initialization failure. Das ist ganz schlecht!" << std::endl;
	exit(1);
      }
    }

    return fr;
  }



  // -------------------------------------------------------------------------------------  
  TString FitResult::toString(FitResult::Type type) {
    TString str = "Undefined";
    if( type == FullMaxLikeRel ) str = "FullMaxLikeRel";
    else if( type == FullMaxLikeAbs ) str = "FullMaxLikeAbs";
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

  
  FitResultFullMaxLikeRel::FitResultFullMaxLikeRel(const Meas meas, unsigned int verbosity)
    : FitResult(meas,verbosity) { }


  bool FitResultFullMaxLikeRel::init() {
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
    
    // Perform extrapolation
    return extrapolate();
  }


  bool FitResultFullMaxLikeRel::extrapolate() {
    if( verbosity_ == 1 ) std::cout << "Extrapolating" << std::endl;
    Extrapolation extra(meanPt_);
    bool result = extra(ptSoft_,values_,statUncerts_,extrapolation_);
    extrapolatedValue_ = extrapolation_->GetParameter(0);
    extrapolatedStatUncert_ = extrapolation_->GetParError(0);

    return result;
  }




  FitResultFullMaxLikeAbs::FitResultFullMaxLikeAbs(const Meas meas, unsigned int verbosity)
    : FitResult(meas,verbosity) { }
  

  bool FitResultFullMaxLikeAbs::init() {
    values_.clear();
    statUncerts_.clear();
    ptSoft_.clear();
    
    for(MeasIt it = meas_.begin() ; it != meas_.end(); ++it) {
      ptSoft_.push_back((*it)->ptSoft());
    }


// 	meanPt_ = (*it)->meanPdfPtTrue();
// 	meanPtUncert_ = (*it)->meanPdfPtTrueUncert();
//       }
//     }

//     // Set fitted values
//     for(MeasIt it = meas_.begin(); it != meas_.end(); ++it) {
//       values_.push_back((*it)->fittedValue(0)/meanPt_);
//       statUncerts_.push_back((*it)->fittedUncert(0)/meanPt_);
//     }
  }
}
