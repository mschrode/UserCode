#include "Sample.h"

#include <iostream>


namespace resolutionFit {

  // -------------------------------------------------------------------
  bool Sample::validType(Type type) {
    if( type != Data && type != MC && type != MCTruth ) {
      std::cerr << "ERROR in Sample::validType(): Unknown type '" << type << "'" << std::endl;
      exit(1);
    }
    
    return true;
  }


  // -------------------------------------------------------------------
  TString Sample::toString(Type type) {
    TString str = "";
    if( validType(type) ) {
      if( type == Data ) str = "Data";
      else if( type == MC ) str = "MC";
      else if( type == MCTruth ) str = "MC Truth";
    }

    return str;
  }


  // -------------------------------------------------------------------
  Sample::Sample(const TString &label, const std::vector<TString> &fileNames, const std::vector<double> &ptSoft, const TString &histNameSuffix, unsigned int verbosity)
    : label_(label), verbosity_(verbosity) {

    if( fileNames.size() != ptSoft.size() ) {
      std::cerr << "ERROR in Sample: Number of specified file names not equal to number of specified ptSoft thresholds" << std::endl;
      exit(1);
    }

    // Read measurements from file
    for(unsigned int i = 0; i < ptSoft.size(); ++i) {
      meas_.push_back(new Measurement(fileNames.at(i),histNameSuffix,ptSoft.at(i),verbosity));
    }
  }


  // -------------------------------------------------------------------
  Sample::~Sample() {
    for(std::vector<Measurement*>::iterator it = meas_.begin(); it != meas_.end(); ++it) {
      delete *it;
    }
    for(std::map<FitResult::Type,FitResult*>::iterator it = fitResult_.begin(); it != fitResult_.end(); ++it) {
      delete it->second;
    }
  }



  // -------------------------------------------------------------------
  bool Sample::addFitResult(FitResult::Type type) {
    bool result = true;
    FitResultMapIt it = fitResult_.find(type);
    if( it != fitResult_.end() ) {
      result = false;
      if( verbosity_ ) std::cerr << "ERROR in Sample: FitResult of type '" << FitResult::toString(type) << "' does already exist" << std::endl;
    } else {
      if( verbosity_ ) std::cout << "Sample::addFitResult(): Adding FitResult '" << FitResult::toString(type) << "'" << std::endl;
      fitResult_[type] = FitResult::createFitResult(type,meas_,verbosity_);
    }

    return result;
  }


  // -------------------------------------------------------------------
  void Sample::ptSoft(std::vector<double> &ptSoftVals) const {
    ptSoftVals.clear();
    for(unsigned int ptSoftBin = 0; ptSoftBin < nPtSoftBins(); ++ptSoftBin) {
      ptSoftVals.push_back(ptSoft(ptSoftBin));
    }
  }
  
  
  // -------------------------------------------------------------------
  void Sample::values(FitResult::Type type, std::vector<double> &val, std::vector<double> &uncert) const {
    val.clear();
    uncert.clear();
    for(unsigned int ptSoftBin = 0; ptSoftBin < nPtSoftBins(); ++ptSoftBin) {
      val.push_back(fittedValue(type,ptSoftBin));
      uncert.push_back(fittedUncert(type,ptSoftBin));
    }
  }


  // -------------------------------------------------------------------
  double Sample::meanPt(FitResult::Type type) const {
    double result = 0.;
    FitResult* fr = 0;
    if( findFitResult(type,fr) ) result = fr->meanPt();
    
    return result;
  }


  // -------------------------------------------------------------------
  double Sample::meanPtStatUncert(FitResult::Type type) const {
    double result = 0.;
    FitResult* fr = 0;
    if( findFitResult(type,fr) ) result = fr->meanPtUncert();
    
    return result;
  }


  // -------------------------------------------------------------------
  double Sample::fittedValue(FitResult::Type type, unsigned int ptSoftBin) const {
    double result = 0.;
    FitResult* fr = 0;
    if( findFitResult(type,fr) ) result = fr->value(ptSoftBin);
    
    return result;
  }


  // -------------------------------------------------------------------
  double Sample::fittedUncert(FitResult::Type type, unsigned int ptSoftBin) const {
    double result = 0.;
    FitResult* fr = 0;
    if( findFitResult(type,fr) ) result = fr->statUncert(ptSoftBin);
    
    return result;
  }


  // -------------------------------------------------------------------
  TF1* Sample::extrapolationFunction(FitResult::Type type, const TString &name) const {
    TF1* result = 0;
    FitResult* fr = 0;
    if( findFitResult(type,fr) ) result = fr->extrapolationFunction(name);
    
    return result;
  }


  // -------------------------------------------------------------------
  double Sample::extrapolatedValue(FitResult::Type type) const {
    double result = 0.;
    FitResult* fr = 0;
    if( findFitResult(type,fr) ) result = fr->extrapolatedValue();
    
    return result;
  }


  // -------------------------------------------------------------------
  double Sample::extrapolatedStatUncert(FitResult::Type type) const {
    double result = 0.;
    FitResult* fr = 0;
    if( findFitResult(type,fr) ) result = fr->extrapolatedStatUncert();
    
    return result;
  }


  // -------------------------------------------------------------------
  double Sample::extrapolatedSystUncert(FitResult::Type type) const {
    double result = 0.;
    FitResult* fr = 0;
    if( findFitResult(type,fr) ) result = fr->extrapolatedSystUncert();
    
    return result;
  }


  // -------------------------------------------------------------------
  bool Sample::findFitResult(FitResult::Type type, FitResult* &fitResult) const {
    bool result = false;
    FitResultMapIt it = fitResult_.find(type);
    if( it != fitResult_.end() ) {
      result = true;
      fitResult = it->second;
    } else {
      std::cerr << "WARNING in Sample::findFitResult(): no fit result type '" << FitResult::toString(type) << "'" << std::endl;
    }

    return result;
  }


  // -------------------------------------------------------------------
  DataSample::DataSample(const TString &label, const std::vector<TString> &fileNames, const std::vector<double> &ptSoft, const TString &histNameSuffix, unsigned int verbosity)
    : Sample(label,fileNames,ptSoft,histNameSuffix,verbosity) {
  }



  // -------------------------------------------------------------------
  MCSample::MCSample(const TString &label, const std::vector<TString> &fileNames, const std::vector<double> &ptSoft, const TString &histNameSuffix, unsigned int verbosity)
    : Sample(label,fileNames,ptSoft,histNameSuffix,verbosity) {
  }
}
