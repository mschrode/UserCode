// $Id: FitResult.cc,v 1.12 2012/06/07 21:10:55 mschrode Exp $

#include "FitResult.h"

#include <algorithm>
#include <iostream>

#include "TH1.h"
#include "TH1D.h"

#include "../util/HistOps.h"


namespace resolutionFit {

  unsigned int FitResult::HIST_COUNT = 0;
  
  // -------------------------------------------------------------------------------------
  bool FitResult::validType(Type type) {
    if( type != MaxLikeKSoftRel && type != FullMaxLikeRel && type != FullMaxLikeAbs && type != PtAsym && type != PtGenAsym ) {
      std::cerr << "ERROR in FitResult::validType(): Unknown type '" << type << "'" << std::endl;
      exit(1);
    }
    
    return true;
  }


  // -------------------------------------------------------------------------------------
  FitResult* FitResult::createFitResult(Type type, const std::vector<Measurement*> &meas, double minPt3, int wpIdx, unsigned int verbosity) {
    FitResult* fr = 0;
    if( validType(type) ) {
      if( type == MaxLikeKSoftRel ) {
	fr = new FitResultMaxLikeKSoftRel(meas,minPt3,wpIdx,verbosity);
      } else if( type == FullMaxLikeRel ) {
	fr = new FitResultFullMaxLikeRel(meas,minPt3,wpIdx,verbosity);
      } else if( type == FullMaxLikeAbs ) {
	fr = new FitResultFullMaxLikeAbs(meas,minPt3,wpIdx,verbosity);
      } else if( type == PtAsym ) {
	fr = new FitResultPtAsym(meas,minPt3,wpIdx,verbosity);
      } else if( type == PtGenAsym ) {
	fr = new FitResultPtGenAsym(meas,minPt3,wpIdx,verbosity);
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
    if( type == MaxLikeKSoftRel ) str = "MaxLikeKSoftRel";
    else if( type == FullMaxLikeRel ) str = "FullMaxLikeRel";
    else if( type == FullMaxLikeAbs ) str = "FullMaxLikeAbs";
    else if( type == SimpleMaxLike ) str = "SimpleMaxLike";
    else if( type == PtAsym ) str = "PtAsym";
    else if( type == PtGenAsym ) str = "PtGenAsym";

    return str;
  }


  // -------------------------------------------------------------------------------------  
  FitResult::FitResult(const Meas &meas, double minPt3, int wpIdx, unsigned int verbosity)
    : meas_(meas), verbosity_(verbosity), workingPointBin_(3), minPt3_(minPt3), wpIdx_(wpIdx) {

    extrapolation_ = 0;
    kSoftFit_ = 0;
    firstPointInExtrapolation_ = 0;
  }


  // -------------------------------------------------------------------------------------  
  FitResult::~FitResult() {
    if( extrapolation_ ) delete extrapolation_;
    if( kSoftFit_ ) delete kSoftFit_;
  }


  // -------------------------------------------------------------------------------------  
  bool FitResult::extrapolate() {
    if( verbosity_ == 1 ) std::cout << "Extrapolating" << std::endl;
    Extrapolation extra(minPt(),maxPt(),minPt3_,wpIdx_);
    bool result = extra(ptSoft_,values_,statUncerts_,extrapolation_,firstPointInExtrapolation_,extrapolatedValue_,extrapolatedStatUncert_,extrapolatedSystUncert_);

    return result;
  }


//   // -------------------------------------------------------------------------------------  
//   TH1* FitResult::spectrum() const {
//     return meas_.front()->histPdfPtTrue();
//   }


  // -------------------------------------------------------------------------------------  
  double FitResult::kSoftSlope() const {
    double yWP = value(workingPointBin_);
    double yEx = extrapolatedValue_;
    
    return yWP != 0. ? yEx/yWP : 0.;
  }


  // -------------------------------------------------------------------------------------  
  double FitResult::kSoftSlopeStatUncert() const {
    double yWP = value(workingPointBin_);
    double yEx = extrapolatedValue_;
    double yWPErr = statUncert(workingPointBin_);
    double yExErr = extrapolatedStatUncert_;
    double err2 = 0.;
    if( yWPErr > 0. && yExErr > 0. ) {
      err2 = yExErr*yExErr/yWP/yWP;
      err2 += yEx*yEx*yExErr*yExErr/yWP/yWP/yWP/yWP;
    }

    return sqrt(err2);
  }


  // -------------------------------------------------------------------------------------  
  void FitResult::setKSoftFit(const TF1* fit) {
    if( kSoftFit_ ) delete kSoftFit_;
    TString name = extrapolation_->GetName();
    name += "_KSoftFit";
    kSoftFit_ = static_cast<TF1*>(fit->Clone(name));
  }




  // -------------------------------------------------------------------------------------  
  FitResultMaxLikeKSoftRel::FitResultMaxLikeKSoftRel(const Meas meas, double minPt3, int wpIdx, unsigned int verbosity)
    : FitResult(meas,minPt3,wpIdx,verbosity) { }


  // -------------------------------------------------------------------------------------  
  bool FitResultMaxLikeKSoftRel::init() {
    values_.clear();
    statUncerts_.clear();
    ptSoft_.clear();

    // Set ptSoft cut values and find smalles ptSoft for meanPt
    double ptSoftSmall = 1000.;
    for(MeasIt it = meas_.begin() ; it != meas_.end(); ++it) {
      ptSoft_.push_back((*it)->ptSoft());
      if( (*it)->ptSoft() < ptSoftSmall ) {
	ptSoftSmall = (*it)->ptSoft();
	
	std::cout << "FitResultMaxLikeKSoftRel::storeResult(): TODO MeanPt from spectrum convoluted with extrapolated *sigma*!" << std::endl;
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
  void FitResultMaxLikeKSoftRel::setKSoftFit(const TF1* fit) {
    if( kSoftFit_ ) delete kSoftFit_;
    TString name = extrapolation_->GetName();
    name += "_KSoftFit";
    kSoftFit_ = static_cast<TF1*>(fit->Clone(name));
    
    double sigFirstPointInFit = 2.*(extrapolatedSystUncert_ + 0.5*extrapolatedValue_);
    extrapolatedSystUncert_ = 0.5*std::abs(sigFirstPointInFit-value(workingPointBin_));
  }





  // -------------------------------------------------------------------------------------  
  FitResultFullMaxLikeRel::FitResultFullMaxLikeRel(const Meas meas, double minPt3, int wpIdx, unsigned int verbosity)
    : FitResult(meas,minPt3,wpIdx,verbosity) { }


  // -------------------------------------------------------------------------------------  
  bool FitResultFullMaxLikeRel::init() {
    values_.clear();
    statUncerts_.clear();
    ptSoft_.clear();

    for(size_t i = 0; i < meas_.size(); ++i) {
      ptSoft_.push_back(meas_.at(i)->ptSoft());
    }
    meanPt_ = meas_.front()->meanPdfPtTrue();
    meanPtUncert_ = meas_.front()->meanPdfPtTrueUncert();

//     // Hack specific for Spring11 PF sample
//     size_t i = 0;
//     if( meas_.at(0)->fittedValue(0)/meas_.at(1)->fittedValue(0) > 2. ) {
//       std::cerr << "\nWARNING: Resolution a first point large; using second point.\n" << std::endl;
//       i = 1;
//     }
//     meanPt_ = meas_.at(i)->meanPdfPtTrue();
//     meanPtUncert_ = meas_.at(i)->meanPdfPtTrueUncert();


    // Set fitted values
    for(MeasIt it = meas_.begin(); it != meas_.end(); ++it) {
      values_.push_back((*it)->fittedValue(0)/meanPt_);
      statUncerts_.push_back((*it)->fittedUncert(0)/meanPt_);
    }
    
    // Perform extrapolation
    return extrapolate();
  }



  //! \brief Extrapolate absolute sigma and convolute extrapolated resolution
  //! with spectrum to obtain mean pt in each bin
  // -------------------------------------------------------------------------------------  
  FitResultFullMaxLikeAbs::FitResultFullMaxLikeAbs(const Meas meas, double minPt3, int wpIdx, unsigned int verbosity)
    : FitResult(meas,minPt3,wpIdx,verbosity) {
    ++HIST_COUNT;
    TString name = "FitResultFullMaxLikeAbs::Spectrum::";
    name += HIST_COUNT;
    spectrum_ = new TH1D(name,"",5000,0.,3000.);
  }
  

  // -------------------------------------------------------------------------------------  
  FitResultFullMaxLikeAbs::~FitResultFullMaxLikeAbs() {
    delete spectrum_;
  }


  // -------------------------------------------------------------------------------------  
  bool FitResultFullMaxLikeAbs::init() {
    values_.clear();
    statUncerts_.clear();
    ptSoft_.clear();
    
    // Copy ptSoft threholds to local data members
    for(MeasIt it = meas_.begin() ; it != meas_.end(); ++it) {
      ptSoft_.push_back((*it)->ptSoft());
    }

    // Set fitted values, i.e. the absolute resolution, and
    // perform extrapolation
    for(MeasIt it = meas_.begin(); it != meas_.end(); ++it) {
      values_.push_back((*it)->fittedValue(0));
      statUncerts_.push_back((*it)->fittedUncert(0));
    }
    meanPt_ = 1.;
    meanPtUncert_ = 1.;

    if( extrapolate() ) {
//       // Convolute extrapolated absolute resolution with spectrum
//       spectrum_->Reset();
//       for(int tBin = 1; tBin <= spectrum_->GetNbinsX(); ++tBin) {
// 	double ptTrue = spectrum_->GetBinCenter(tBin);
// 	// Convolution with cuts on ptAve
// 	double s = extrapolatedValue_/sqrt(2.);	// This is still the absolute resolution!
// 	double c = sqrt(M_PI/2.)*s*( erf((maxPt()-ptTrue)/s/sqrt(2.)) - erf((minPt()-ptTrue)/s/sqrt(2.)) );
// 	double pdf = c*meas_.front()->pdfPtTrue(ptTrue); // Get truth pdf from bin with lowest ptSoft; this should be closest to real dijet spectrum
	
// 	// Store (un-normalized) truth pdf for
// 	// this value of ptTrue in hash table
// 	spectrum_->SetBinContent(tBin,pdf);
//       } // End of loop over ptTrue
//       // Normalise values of truth pdf
//       if( spectrum_->Integral() ) spectrum_->Scale(1./spectrum_->Integral("width"));

//       // Get mean of pttrue spectrum
//       meanPt_ = spectrum_->GetMean();
//       meanPtUncert_ = meas_.front()->meanPdfPtTrueUncert();

      double meanPtUp = computeMeanPtTrue(std::abs(extrapolatedValue_+extrapolatedStatUncert_));
      double meanPtDn = computeMeanPtTrue(std::abs(extrapolatedValue_-extrapolatedStatUncert_));
      meanPt_ = computeMeanPtTrue(extrapolatedValue_);
      meanPtUncert_ = 0.5*( std::abs(meanPt_-meanPtDn) + std::abs(meanPt_-meanPtUp) );
      
      // Set relative resolution
      values_.clear();
      statUncerts_.clear();
      for(MeasIt it = meas_.begin(); it != meas_.end(); ++it) {
	values_.push_back((*it)->fittedValue(0)/meanPt_);
	statUncerts_.push_back((*it)->fittedUncert(0)/meanPt_);
      }
      if( extrapolation_->GetParameter(1) > 0. ) { // want positive slope
	extrapolatedValue_ /= meanPt_;
	extrapolatedStatUncert_ /= meanPt_;
	extrapolatedSystUncert_ /= meanPt_;
      } else {
	meanPt_ = 1.;
	extrapolatedValue_ = 0.;
	extrapolatedStatUncert_ = 1000.;
	extrapolatedSystUncert_ = 1000.;
      }
      return true;
    } else {
      return false;
    }    
  }

  // Convolute absolute resolution with spectrum
  // to otbain meanPtTrue
  // -------------------------------------------------------------------------------------
  double FitResultFullMaxLikeAbs::computeMeanPtTrue(double sigma) {
    spectrum_->Reset();
    for(int tBin = 1; tBin <= spectrum_->GetNbinsX(); ++tBin) {
      double ptTrue = spectrum_->GetBinCenter(tBin);
      // Convolution with cuts on ptAve
      double s = sigma/sqrt(2.);	// This is still the absolute resolution!
      double c = sqrt(M_PI/2.)*s*( erf((maxPt()-ptTrue)/s/sqrt(2.)) - erf((minPt()-ptTrue)/s/sqrt(2.)) );
      double pdf = c*meas_.front()->pdfPtTrue(ptTrue); // Get truth pdf from bin with lowest ptSoft; this should be closest to real dijet spectrum
   
      // Store (un-normalized) truth pdf for
      // this value of ptTrue in hash table
      spectrum_->SetBinContent(tBin,pdf);
    } // End of loop over ptTrue
      // Normalise values of truth pdf
    if( spectrum_->Integral() ) spectrum_->Scale(1./spectrum_->Integral("width"));

    return spectrum_->GetMean();
  }


  
  //! \brief Get assumed truth pdf for this EtaPt bin
  // -------------------------------------------------------------------------------------  
//   TH1* FitResultFullMaxLikeAbs::spectrum() const {
//     ++HIST_COUNT;
//     TString name = "FitResultFullMaxLikeAbs::spectrum::";
//     name += HIST_COUNT;

//     return static_cast<TH1*>(spectrum_->Clone(name));
//   }



  // -------------------------------------------------------------------------------------  
  FitResultPtAsym::FitResultPtAsym(const Meas meas, double minPt3, int wpIdx, unsigned int verbosity)
    : FitResult(meas,minPt3,wpIdx,verbosity) { }


  // -------------------------------------------------------------------------------------  
  bool FitResultPtAsym::init() {
    values_.clear();
    statUncerts_.clear();
    ptSoft_.clear();

    // Set ptSoft cut values
    for(MeasIt it = meas_.begin() ; it != meas_.end(); ++it) {
      ptSoft_.push_back((*it)->ptSoft());
    }
    meanPt_ = 0.;
    meanPtUncert_ = 1000.;
    for(MeasIt it = meas_.begin() ; it != meas_.end(); ++it) {
      if( (*it)->meanPtAve() > 3. ) {
	meanPt_ = (*it)->meanPtAve();
	meanPtUncert_ = (*it)->meanPtAveUncert();
	break;
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
  FitResultPtGenAsym::FitResultPtGenAsym(const Meas meas, double minPt3, int wpIdx, unsigned int verbosity)
    : FitResult(meas,minPt3,wpIdx,verbosity) { }


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
	meanPt_ = (*it)->meanPtGen();
	meanPtUncert_ = (*it)->meanPtGenUncert();
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


}

