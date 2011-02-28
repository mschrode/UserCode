// $Id: FitResult.cc,v 1.4 2011/02/25 19:50:21 mschrode Exp $

#include "FitResult.h"

#include <algorithm>
#include <iostream>

#include "../util/HistOps.h"


namespace resolutionFit {

  // -------------------------------------------------------------------------------------
  bool FitResult::validType(Type type) {
    if( type != FullMaxLikeRel && type != FullMaxLikeAbs && type != PtAsym && type != PtGenAsym ) {
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
      } else if( type == PtAsym ) {
	fr = new FitResultPtAsym(meas,verbosity);
      } else if( type == PtGenAsym ) {
	fr = new FitResultPtGenAsym(meas,verbosity);
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



  // -------------------------------------------------------------------------------------  
  FitResultFullMaxLikeRel::FitResultFullMaxLikeRel(const Meas meas, unsigned int verbosity)
    : FitResult(meas,verbosity) { }


  // -------------------------------------------------------------------------------------  
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


  // -------------------------------------------------------------------------------------  
  bool FitResultFullMaxLikeRel::extrapolate() {
    if( verbosity_ == 1 ) std::cout << "Extrapolating" << std::endl;
    Extrapolation extra(meanPt_);
    bool result = extra(ptSoft_,values_,statUncerts_,extrapolation_,extrapolatedSystUncert_);
    extrapolatedValue_ = extrapolation_->GetParameter(0);
    extrapolatedStatUncert_ = extrapolation_->GetParError(0);

    return result;
  }




  // -------------------------------------------------------------------------------------  
  FitResultFullMaxLikeAbs::FitResultFullMaxLikeAbs(const Meas meas, unsigned int verbosity)
    : FitResult(meas,verbosity) { }
  

  // -------------------------------------------------------------------------------------  
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




  // -------------------------------------------------------------------------------------  
  FitResultPtAsym::FitResultPtAsym(const Meas meas, unsigned int verbosity)
    : FitResult(meas,verbosity) { }


  // -------------------------------------------------------------------------------------  
  bool FitResultPtAsym::init() {
    values_.clear();
    statUncerts_.clear();
    ptSoft_.clear();

    // Set ptSoft cut values and find smalles ptSoft for meanPt
    double ptSoftSmall = 1000.;
    for(MeasIt it = meas_.begin() ; it != meas_.end(); ++it) {
      ptSoft_.push_back((*it)->ptSoft());
      if( (*it)->ptSoft() < ptSoftSmall ) {
	ptSoftSmall = (*it)->ptSoft();
	meanPt_ = (*it)->meanPtAve();
	meanPtUncert_ = (*it)->meanPtAveUncert();
      }
    }

    // Set fitted values
    for(MeasIt it = meas_.begin(); it != meas_.end(); ++it) {
      double width = 0.;
      double widthErr = 1000.;
      double rms = 0.;
      double rmsErr = 0.;
      TH1* hPtAsym = (*it)->histPtAsym();
      hPtAsym->GetXaxis()->SetRangeUser(-1.,1.);
      if( util::HistOps::fitCoreWidth(hPtAsym,2.,width,widthErr,rms,rmsErr) ) {
	values_.push_back(sqrt(2.)*width);
 	statUncerts_.push_back(sqrt(2.)*widthErr);
      } else {
	values_.push_back(0.);
	statUncerts_.push_back(1000.);
      }
      delete hPtAsym;
    }
    
    // Perform extrapolation
    return extrapolate();
  }


  // -------------------------------------------------------------------------------------  
  bool FitResultPtAsym::extrapolate() {
    if( verbosity_ == 1 ) std::cout << "Extrapolating" << std::endl;
    Extrapolation extra(meanPt_);
    bool result = extra(ptSoft_,values_,statUncerts_,extrapolation_,extrapolatedSystUncert_);
    extrapolatedValue_ = extrapolation_->GetParameter(0);
    extrapolatedStatUncert_ = extrapolation_->GetParError(0);

    return result;
  }




  // -------------------------------------------------------------------------------------  
  FitResultPtGenAsym::FitResultPtGenAsym(const Meas meas, unsigned int verbosity)
    : FitResult(meas,verbosity) { }


  // -------------------------------------------------------------------------------------  
  bool FitResultPtGenAsym::init() {
    values_.clear();
    statUncerts_.clear();
    ptSoft_.clear();

    // Set ptSoft cut values and find smalles ptSoft for meanPt
    double ptSoftSmall = 1000.;
    for(MeasIt it = meas_.begin() ; it != meas_.end(); ++it) {
      ptSoft_.push_back((*it)->ptSoft());
      if( (*it)->ptSoft() < ptSoftSmall ) {
	ptSoftSmall = (*it)->ptSoft();
	meanPt_ = (*it)->meanPtAve();
	meanPtUncert_ = (*it)->meanPtAveUncert();
      }
    }

    // Set fitted values
    bool hasPtGenAsym = true;
    for(MeasIt it = meas_.begin(); it != meas_.end(); ++it) {
      double width = 0.;
      double widthErr = 1000.;
      double rms = 0.;
      double rmsErr = 0.;
      TH1* hPtGenAsym = (*it)->histPtGenAsym();
      hPtGenAsym->GetXaxis()->SetRangeUser(-1.,1.);
      if( hPtGenAsym->Integral(1,hPtGenAsym->GetNbinsX()) > 0. ) {
	if( util::HistOps::fitCoreWidth(hPtGenAsym,2.,width,widthErr,rms,rmsErr) ) {
	  values_.push_back(sqrt(2.)*rms);
	  statUncerts_.push_back(sqrt(2.)*rmsErr);
	} else {
	  values_.push_back(0.);
	  statUncerts_.push_back(1000.);
	}
      } else {
	hasPtGenAsym = false;
	values_.push_back(0.);
	statUncerts_.push_back(1000.);
      }
      delete hPtGenAsym;
    }
    extrapolatedValue_ = 0.;
    extrapolatedStatUncert_ = 0.;
    extrapolatedSystUncert_ = 0.;
    
    // Perform extrapolation
    return hasPtGenAsym ? extrapolate() : true;
  }


  // -------------------------------------------------------------------------------------  
  bool FitResultPtGenAsym::extrapolate() {
    if( verbosity_ == 1 ) std::cout << "Extrapolating" << std::endl;
    Extrapolation extra(meanPt_);
    bool result = extra(ptSoft_,values_,statUncerts_,extrapolation_,extrapolatedSystUncert_);
    extrapolatedValue_ = extrapolation_->GetParameter(0);
    extrapolatedStatUncert_ = extrapolation_->GetParError(0);

    return result;
  }
}

